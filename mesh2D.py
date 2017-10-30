import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import elementDefinitions
import time
import femFunctions

class triMesh:
    '''
    Contains elements within the triangular mesh and computes and
    stores section properties
    '''

    def __init__(self, genMesh, nu):
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
        self.initialise()

    def initialise(self):
        # initialise variables
        totalArea = totalQx = totalQy = totalIxx_g = totalIyy_g = totalIxy_g = 0
        shearK = np.zeros((self.noNodes, self.noNodes))
        torsionF = shearFPsi = shearFPhi = np.zeros(self.noNodes)

        # loop through all elements where summing over all elements is required
        for el in self.elements:
            # calculate total area
            totalArea += el.area

            # calculate first moments of area about global axis
            totalQx += el.Qx
            totalQy += el.Qy

            # calculate second moments of area about global axis
            totalIxx_g += el.ixx
            totalIyy_g += el.iyy
            totalIxy_g += el.ixy

            # assemble stiffness matrix and load vector for warping constant
            indxs = np.ix_(el.nodes, el.nodes)
            shearK[indxs] += el.shearKe
            torsionF[el.nodes] += el.torsionFe

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

        # ----------------------------------------------------------------------
        # TORSION PROPERTIES:
        # ----------------------------------------------------------------------
        # calculate warping constant and torsion constant
        self.omega = np.linalg.solve(shearK, torsionF)
        self.J = (self.ixx_g + self.iyy_g -
                self.omega.dot(shearK).dot(np.transpose(self.omega)))

        # # ----------------------------------------------------------------------
        # # SHEAR PROPERTIES:
        # # ----------------------------------------------------------------------
        # for el in self.elements:
        #     shearFPsi[el.nodes] += el.shearFePsi(self.ixx_c, self.ixy_c)
        #     shearFPhi[el.nodes] += el.shearFePhi(self.iyy_c, self.ixy_c)
        #
        # self.Psi = np.linalg.solve(shearK, shearFPsi)
        # self.Phi = np.linalg.solve(shearK, shearFPhi)
        #
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

    def contourPlot(self, principalAxis = False, z = False, nodes = False):
        plt.figure()
        plt.gca().set_aspect('equal')
        plt.triplot(self.pointArray[:,0], self.pointArray[:,1], self.elementArray[:,0:3], lw=0.5, color='black')

        if principalAxis:
            d1 = max(abs(self.d2max), abs(self.d2min))
            d2 = max(abs(self.d1max), abs(self.d1min))
            plt.plot([self.cx - d1 * np.cos(self.phi * np.pi / 180), self.cx + d1 * np.cos(self.phi * np.pi / 180)], [self.cy - d1 * np.sin(self.phi * np.pi / 180), self.cy + d1 * np.sin(self.phi * np.pi / 180)])
            plt.plot([self.cx - d2 * np.cos(self.phi * np.pi / 180 + np.pi / 2), self.cx + d2 * np.cos(self.phi * np.pi / 180 + np.pi / 2)], [self.cy - d2 * np.sin(self.phi * np.pi / 180 + np.pi / 2), self.cy + d2 * np.sin(self.phi * np.pi / 180 + np.pi / 2)])

        if z.any():
            cmap = cm.get_cmap(name = 'jet')
            # v = np.linspace(-10, 10, 15, endpoint=True)
            trictr = plt.tricontourf(self.pointArray[:,0], self.pointArray[:,1], self.elementArray[:,0:3], z, cmap=cmap)
            cbar = plt.colorbar(trictr)

        if nodes:
            plt.plot(self.pointArray[:,0], self.pointArray[:,1], 'ko', markersize = 2)

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

def functionTimer(function):
    start_time = time.clock()
    function()
    print("--- %s completed in %s seconds ---" % (function.__name__, time.clock() - start_time))
