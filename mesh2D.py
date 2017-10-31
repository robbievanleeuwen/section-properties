import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import meshpy.triangle as triangle
import time
import elementDefinitions
import femFunctions

class triMesh:
    '''
    Contains elements within the triangular mesh and computes and
    stores section properties
    '''

    def __init__(self, genMesh, nu=0):
        triElements = [] # list holding all element objects
        pointArray = np.array(genMesh.points) # save points to numpy array
        elementArray = np.array(genMesh.elements) # save elements to numpy array
        # swap mid-node order to retain node ordering consistency
        elementArray[:,[3,4,5]] = elementArray[:,[5,3,4]]

        for tri in elementArray:
            x1 = pointArray[tri[0]][0]
            y1 = pointArray[tri[0]][1]
            x2 = pointArray[tri[1]][0]
            y2 = pointArray[tri[1]][1]
            x3 = pointArray[tri[2]][0]
            y3 = pointArray[tri[2]][1]
            x4 = pointArray[tri[3]][0]
            y4 = pointArray[tri[3]][1]
            x5 = pointArray[tri[4]][0]
            y5 = pointArray[tri[4]][1]
            x6 = pointArray[tri[5]][0]
            y6 = pointArray[tri[5]][1]
            vertices = np.array([[x1, x2, x3, x4, x5, x6], [y1, y2, y3, y4, y5, y6]])
            triElements.append(elementDefinitions.tri6(vertices, tri, nu))

        self.elements = triElements # store element list in mesh2D object

        # store the mesh arrays
        self.pointArray = pointArray
        self.elementArray = elementArray

        self.nu = nu # poissons ratio of material
        self.noNodes = len(genMesh.points) # total number of nodes in mesh

    def computeSimpleProperties(self):
        # initialise variables
        totalArea = 0
        totalQx = 0
        totalQy = 0
        totalIxx_g = 0
        totalIyy_g = 0
        totalIxy_g = 0

        # loop through all elements where summing over all elements is required
        for el in self.elements:
            # calculate total area
            totalArea += el.area()

            # calculate first moments of area about global axis
            totalQx += el.Qx()
            totalQy += el.Qy()

            # calculate second moments of area about global axis
            totalIxx_g += el.ixx()
            totalIyy_g += el.iyy()
            totalIxy_g += el.ixy()

        # ----------------------------------------------------------------------
        # GLOBAL xy AXIS PROPERTIES:
        # ----------------------------------------------------------------------
        # assign total area
        self.area = totalArea

        # assign first moments of area
        self.Qx = totalQx
        self.Qy = totalQy

        # calculate centroids
        self.cx = totalQy / totalArea
        self.cy = totalQx / totalArea

        # assign global axis second moments of area
        self.ixx_g = totalIxx_g
        self.iyy_g = totalIyy_g
        self.ixy_g = totalIxy_g

        # ----------------------------------------------------------------------
        # CENTROIDAL xy AXIS PROPERTIES:
        # ----------------------------------------------------------------------
        # calculate second moments of area about the centroidal xy axis
        self.ixx_c = totalIxx_g - totalQx ** 2 / totalArea
        self.iyy_c = totalIyy_g - totalQy ** 2 / totalArea
        self.ixy_c = totalIxy_g - totalQx * totalQy / totalArea

        # calculate section modulii about the centroidal xy axis
        self.centroidalSectionModulii()

        # calculate radii of gyration about centroidal xy axis
        self.rx_c = (self.ixx_c / totalArea) ** 0.5
        self.ry_c = (self.iyy_c / totalArea) ** 0.5

        # ----------------------------------------------------------------------
        # PRCINCIPAL AXIS PROPERTIES:
        # ----------------------------------------------------------------------
        # calculate prinicpal second moments of area about the centroidal xy axis
        Delta = (((self.ixx_c - self.iyy_c) / 2) ** 2 + self.ixy_c ** 2) ** 0.5
        self.i11_c = (self.ixx_c + self.iyy_c) / 2 + Delta
        self.i22_c = (self.ixx_c + self.iyy_c) / 2 - Delta

        # calculate initial principal axis angle
        self.phi = np.arctan2(self.ixx_c - self.i11_c, self.ixy_c) * 180 / np.pi

        # calculate section modulii about the principal axis
        self.principalSectionModulii()

        # calculate radii of gyration about centroidal principal axis
        self.r1_c = (self.i11_c / totalArea) ** 0.5
        self.r2_c = (self.i22_c / totalArea) ** 0.5

    def computeWarpingProperties(self, A, ixx, iyy, ixy):
        # initialise variables
        shearK = np.zeros((self.noNodes, self.noNodes))
        torsionF = np.zeros(self.noNodes)
        shearFPsi = np.zeros(self.noNodes)
        shearFPhi = np.zeros(self.noNodes)
        shearCentreXInt = 0
        shearCentreYInt = 0
        Q_omega = 0
        i_omega = 0
        i_xomega = 0
        i_yomega = 0

        # loop through all elements where summing over all elements is required
        for el in self.elements:
            # assemble stiffness matrix and load vector for warping function
            indxs = np.ix_(el.nodes, el.nodes)
            shearK[indxs] += el.shearKe()
            torsionF[el.nodes] += el.torsionFe()

        # ----------------------------------------------------------------------
        # TORSION PROPERTIES:
        # ----------------------------------------------------------------------
        # calculate warping constant and torsion constant
        (self.omega, error) = lgMultSolve(shearK, torsionF)
        self.J = ixx + iyy - self.omega.dot(shearK).dot(np.transpose(self.omega))

        # ----------------------------------------------------------------------
        # SHEAR PROPERTIES:
        # ----------------------------------------------------------------------
        # calculate shear functions, shear centre integrals and warping moments
        for el in self.elements:
            shearFPsi[el.nodes] += el.shearFePsi(ixx, ixy)
            shearFPhi[el.nodes] += el.shearFePhi(iyy, ixy)
            shearCentreXInt += el.shearCentreXInt(iyy, ixy)
            shearCentreYInt += el.shearCentreYInt(ixx, ixy)
            Q_omega += el.Q_omega(self.omega[el.nodes])
            i_omega += el.i_omega(self.omega[el.nodes])
            i_xomega += el.i_xomega(self.omega[el.nodes])
            i_yomega += el.i_yomega(self.omega[el.nodes])

        # solve for shear functions
        (self.Psi, error) = lgMultSolve(shearK, shearFPsi)
        (self.Phi, error) = lgMultSolve(shearK, shearFPhi)

        # calculate shear centres (elasticity)
        Delta_s = 2 * (1 + self.nu) * (ixx * iyy - ixy ** 2)
        self.x_se = ((1 / Delta_s) * ((self.nu / 2 * shearCentreXInt) -
            torsionF.dot(self.Phi)))
        self.y_se = ((1 / Delta_s) * ((self.nu / 2 * shearCentreYInt) +
            torsionF.dot(self.Psi)))

        # calculate shear centres (Trefftz's)
        self.x_st = (ixy * i_xomega - iyy * i_yomega) / (ixx * iyy - ixy ** 2)
        self.y_st = (ixx * i_xomega - ixy * i_yomega) / (ixx * iyy - ixy ** 2)

        # calculate warping constant
        self.Gamma = (i_omega - Q_omega ** 2 / A - self.y_se * i_xomega +
            self.x_se * i_yomega)

        # # ----------------------------------------------------------------------
        # # STRESS CALCULATION:
        # # ----------------------------------------------------------------------
        # # calculate torsion stress
        # self.torsionStress()
        # self.shearStress()

    def centroidalSectionModulii(self):
        # determine extreme values of the cartesian co-ordinates
        xmax = self.pointArray[:,0].max()
        xmin = self.pointArray[:,0].min()
        ymax = self.pointArray[:,1].max()
        ymin = self.pointArray[:,1].min()

        # evaluate section modulii
        self.zxx_plus = self.ixx_c / (ymax - self.cy)
        self.zxx_minus = self.ixx_c / (self.cy - ymin)
        self.zyy_plus = self.iyy_c / (xmax - self.cx)
        self.zyy_minus = self.iyy_c / (self.cx - xmin)

    def principalSectionModulii(self):
         # unit vectors in the direction of the principal axes
        u1 = np.array([np.cos(self.phi * np.pi / 180), np.sin(self.phi * np.pi / 180)])
        u2 = np.array([-np.sin(self.phi * np.pi / 180), np.cos(self.phi * np.pi / 180)])

        # initialise min/max distance variables
        self.d1max = self.d1min = self.d2max = self.d2min = 0

        # loop through all co-ordinates to determine extreme values
        # in the principal directions
        for vertex in self.pointArray:
            # vector from point to centroid
            PQ = np.array([self.cx - vertex[0], self.cy - vertex[1]])
            # perpendicular distance from point to 1 and 2 axes
            d1 = np.linalg.norm(np.cross(PQ, u1))
            d2 = np.linalg.norm(np.cross(PQ, u2))

            # check to see if point is above or below axes
            if np.cross(-PQ, u1) < 0: # point is above 1 axis
                self.d1max = max(self.d1max, d1)
            else: # point is below 1 axis
                self.d1min = min(self.d1min, -d1)
            if np.cross(-PQ, u2) < 0: # point is above 2 axis
                self.d2max = max(self.d2max, d2)
            else: # point is below 2 axis
                self.d2min = min(self.d2min, -d2)

        # evaluate principal section modulii
        self.z11_plus = self.i11_c / abs(self.d1max)
        self.z11_minus = self.i11_c / abs(self.d1min)
        self.z22_plus = self.i22_c / abs(self.d2max)
        self.z22_minus = self.i22_c / abs(self.d2min)

    # def torsionStress(self):
    #     self.tau_zx_torsion = np.zeros(self.noNodes)
    #     self.tau_zy_torsion = np.zeros(self.noNodes)
    #     node_count = np.zeros(self.noNodes)
    #
    #     for el in self.elements:
    #         tau_torsion = el.torsionStress(1000, self.omega[el.nodes], self.J)
    #         self.tau_zx_torsion[el.nodes] += tau_torsion[:,0]
    #         self.tau_zy_torsion[el.nodes] += tau_torsion[:,1]
    #         node_count[el.nodes] += 1
    #
    #     self.tau_zx_torsion *= 1 / node_count
    #     self.tau_zy_torsion *= 1 / node_count
    #     self.tau_torsion = (self.tau_zx_torsion ** 2 + self.tau_zy_torsion ** 2) ** 0.5
    #
    # def shearStress(self):
    #     self.tau_zx_shear = np.zeros(self.noNodes)
    #     self.tau_zy_shear = np.zeros(self.noNodes)
    #     node_count = np.zeros(self.noNodes)
    #
    #     for el in self.elements:
    #         tau_shear = el.shearStress(0, 0.001, self.Psi[el.nodes], self.Phi[el.nodes], self.ixx_c, self.iyy_c, self.ixy_c)
    #         self.tau_zx_shear[el.nodes] += tau_shear[:,0]
    #         self.tau_zy_shear[el.nodes] += tau_shear[:,1]
    #         node_count[el.nodes] += 1
    #
    #     self.tau_zx_shear *= 1 / node_count
    #     self.tau_zy_shear *= 1 / node_count
    #     self.tau_shear = (self.tau_zx_shear ** 2 + self.tau_zy_shear ** 2) ** 0.5

    def contourPlot(self, principalAxis = False, z = None, nodes = False, plotTitle = ''):
        plt.figure()
        plt.gca().set_aspect('equal')
        plt.triplot(self.pointArray[:,0], self.pointArray[:,1], self.elementArray[:,0:3], lw=0.5, color='black')
        plt.title(plotTitle)

        if principalAxis:
            d1 = max(abs(self.d2max), abs(self.d2min))
            d2 = max(abs(self.d1max), abs(self.d1min))
            plt.plot([self.cx - d1 * np.cos(self.phi * np.pi / 180), self.cx + d1 * np.cos(self.phi * np.pi / 180)], [self.cy - d1 * np.sin(self.phi * np.pi / 180), self.cy + d1 * np.sin(self.phi * np.pi / 180)])
            plt.plot([self.cx - d2 * np.cos(self.phi * np.pi / 180 + np.pi / 2), self.cx + d2 * np.cos(self.phi * np.pi / 180 + np.pi / 2)], [self.cy - d2 * np.sin(self.phi * np.pi / 180 + np.pi / 2), self.cy + d2 * np.sin(self.phi * np.pi / 180 + np.pi / 2)])

        if z is not None:
            cmap = cm.get_cmap(name = 'jet')
            # v = np.linspace(-10, 10, 15, endpoint=True)
            trictr = plt.tricontourf(self.pointArray[:,0], self.pointArray[:,1], self.elementArray[:,0:3], z, cmap=cmap)
            cbar = plt.colorbar(trictr)

        if nodes:
            plt.plot(self.pointArray[:,0], self.pointArray[:,1], 'ko', markersize = 1)

        plt.show()

    def quiverPlot(self, u, v):
        plt.figure()
        plt.gca().set_aspect('equal')
        plt.triplot(self.pointArray[:,0], self.pointArray[:,1], self.elementArray, lw=0.5, color='black')
        c = np.hypot(u, v)
        cmap = cm.get_cmap(name='jet')
        quiv = plt.quiver(self.pointArray[:,0], self.pointArray[:,1], u, v, c, cmap=cmap)
        cbar = plt.colorbar(quiv)
        plt.show()

    def printSimpleResults(self):
        print "-------------------------"
        print "Global xy Axis Properties"
        print "-------------------------"
        print "Area = {}".format(self.area)
        print "Qx = {}".format(self.Qx)
        print "Qy = {}".format(self.Qy)
        print "cx = {}".format(self.cx)
        print "cy = {}".format(self.cy)
        print "Ixx_g = {}".format(self.ixx_g)
        print "Iyy_g = {}".format(self.iyy_g)
        print "Ixy_g = {}".format(self.ixy_g)
        print ""
        print "-----------------------------"
        print "Centroidal xy Axis Properties"
        print "-----------------------------"
        print "Ixx_c = {}".format(self.ixx_c)
        print "Iyy_c = {}".format(self.iyy_c)
        print "Ixy_c = {}".format(self.ixy_c)
        print "Zxx_plus = {}".format(self.zxx_plus)
        print "Zxx_minus = {}".format(self.zxx_minus)
        print "Zyy_plus = {}".format(self.zyy_plus)
        print "Zyy_minus = {}".format(self.zyy_minus)
        print "rx_c = {}".format(self.rx_c)
        print "ry_c = {}".format(self.ry_c)
        print ""
        print "-----------------------------"
        print "Principal Axis Properties"
        print "-----------------------------"
        print "phi = {}".format(self.phi)
        print "I11_c = {}".format(self.i11_c)
        print "I22_c = {}".format(self.i22_c)
        print "Z11_plus = {}".format(self.z11_plus)
        print "Z11_minus = {}".format(self.z11_minus)
        print "Z22_plus = {}".format(self.z22_plus)
        print "Z22_minus = {}".format(self.z22_minus)
        print "r1_c = {}".format(self.r1_c)
        print "r2_c = {}".format(self.r2_c)
        print ""

    def printWarpingResults(self):
        print "-----------------------------"
        print "Torsional Properties"
        print "-----------------------------"
        print "J = {}".format(self.J)
        print "Iw = {}".format(self.Gamma)
        print ""
        print "-----------------------------"
        print "Shear Properties"
        print "-----------------------------"
        print "x_s,e = {}".format(self.x_se)
        print "y_s,e = {}".format(self.y_se)
        print "x_s,t = {}".format(self.x_st)
        print "y_s,t = {}".format(self.y_st)

def createMesh(points, facets, holes=[], maxArea=[], minAngle=30, meshOrder=2):
    info = triangle.MeshInfo()
    info.set_points(points)
    info.set_holes(holes)
    info.set_facets(facets)

    return triangle.build(info, max_volume = maxArea, min_angle = minAngle, mesh_order = meshOrder)

def shiftGeometry(points, facets, holes, cx, cy):
    # initialise shifted points and holes lists
    shiftedPoints = []
    shiftedHoles =[]

    # shift points by centroid
    for point in points:
        shiftedPoints.append((point[0] - cx, point[1] - cy))

    for hole in holes:
        shiftedHoles.append((hole[0] - cx, hole[1] - cy))

    return (shiftedPoints, shiftedHoles)

def lgMultSolve(K, f):
    Nvec1 = np.ones((K.shape[0], 1))
    Nvec2 = np.ones((1, K.shape[0] + 1))
    Nvec2[:,-1] = 0

    K = np.concatenate((K, Nvec1), axis=1)
    K = np.concatenate((K, Nvec2), axis=0)
    f = np.append(f, 0)

    u = np.linalg.solve(K, f)
    return (u[:-1], u[-1])

def functionTimer(function):
    start_time = time.clock()
    function()
    print("--- %s completed in %s seconds ---" % (function.__name__, time.clock() - start_time))
