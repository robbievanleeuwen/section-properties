import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import elementDefinitions
import femFunctions

class triMesh:
    '''
    Contains elements within the triangular mesh and computes and
    stores section properties
    '''

    def __init__(self, genMesh, nu=0, geometricMesh=None):
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

        # load results of geometricMesh
        if geometricMesh is not None:
            self.geometricMesh = geometricMesh
        else:
            self.geometricMesh = None

    # --------------------------------------------------------------------------
    # SECTION PROPERTY COMPUTATION:
    # --------------------------------------------------------------------------
    def computeGeometricProperties(self):
        # initialise variables
        totalArea = 0
        totalQx = 0
        totalQy = 0
        totalIxx_g = 0
        totalIyy_g = 0
        totalIxy_g = 0

        for el in self.elements:
            # calculate total area
            (elArea, elQx, elQy, elIxx_g, elIyy_g, elIxy_g) = el.geometricProperties()

            totalArea += elArea
            totalQx += elQx
            totalQy += elQy
            totalIxx_g += elIxx_g
            totalIyy_g += elIyy_g
            totalIxy_g += elIxy_g

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
        tol = 1e-9
        if abs(self.ixx_c - self.i11_c) < tol * self.i11_c:
            self.phi = 0
        else:
            self.phi = np.arctan2(self.ixx_c - self.i11_c, self.ixy_c) * 180 / np.pi

        # calculate section modulii about the principal axis
        self.principalSectionModulii()

        # calculate radii of gyration about centroidal principal axis
        self.r1_c = (self.i11_c / totalArea) ** 0.5
        self.r2_c = (self.i22_c / totalArea) ** 0.5

    def computeGlobalPlasticProperties(self, points, facets, holes):
        tol = 1e-6
        ux = np.array([1, 0])
        uy = np.array([0, 1])

        # compute plastic centroids and plastic section modulii
        (x_a_n, topArea, botArea, topCentroidx, botCentroidx) = (plasticCentroidAlgorithm(tol,
            100, uy, self.cx, self.cy, self.xmin - self.cx, self.xmax - self.cx,
            points, facets, holes, self.pointArray, self.elementArray))
        (y_a_n, topArea, botArea, topCentroidy, botCentroidy) = (plasticCentroidAlgorithm(tol,
            100, ux, self.cx, self.cy, self.ymin - self.cy, self.ymax - self.cy,
            points, facets, holes, self.pointArray, self.elementArray))

        self.x_pc = self.cx + x_a_n
        self.y_pc = self.cy + y_a_n
        self.Sxx = self.area / 2 * abs(topCentroidy[1] - botCentroidy[1])
        self.Syy = self.area / 2 * abs(topCentroidx[0] - botCentroidx[0])

    def computePrincipalPlasticProperties(self, points, facets, holes):
        tol = 1e-6
        u1 = np.array([np.cos(self.phi * np.pi / 180), np.sin(self.phi * np.pi / 180)])
        u2 = np.array([-np.sin(self.phi * np.pi / 180), np.cos(self.phi * np.pi / 180)])

        # compute plastic centroids and plastic section modulii
        (x_1_a_n, topArea, botArea, topCentroid1, botCentroid1) = (plasticCentroidAlgorithm(tol,
            100, u1, self.cx, self.cy, self.x_1min, self.x_1max, points, facets,
            holes, self.pointArray, self.elementArray))
        (y_2_a_n, topArea, botArea, topCentroid2, botCentroid2) = (plasticCentroidAlgorithm(tol,
            100, u2, self.cx, self.cy, self.y_2min, self.y_2max, points, facets,
            holes, self.pointArray, self.elementArray))

        (tc1_1, tc1_2) = femFunctions.principalCoordinate(self.phi, topCentroid1[0] - self.cx, topCentroid1[1] - self.cy)
        (bc1_1, bc1_2) = femFunctions.principalCoordinate(self.phi, botCentroid1[0] - self.cx, botCentroid1[1] - self.cy)
        (tc2_1, tc2_2) = femFunctions.principalCoordinate(self.phi, topCentroid2[0] - self.cx, topCentroid2[1] - self.cy)
        (bc2_1, bc2_2) = femFunctions.principalCoordinate(self.phi, botCentroid2[0] - self.cx, botCentroid2[1] - self.cy)

        self.x_1_pc = self.cx + x_1_a_n * u2[0] + y_2_a_n * u1[0]
        self.y_2_pc = self.cy + x_1_a_n * u2[1] + y_2_a_n * u1[1]
        self.S11 = self.area / 2 * abs(tc1_2 - bc1_2)
        self.S22 = self.area / 2 * abs(tc2_1 - bc2_1)

    def computeWarpingProperties(self):
        # load areas and second moments of area
        A = self.geometricMesh.area
        ixx = self.geometricMesh.ixx_c
        iyy = self.geometricMesh.iyy_c
        ixy = self.geometricMesh.ixy_c
        phi = self.geometricMesh.phi

        # calculate stiffness matrix and load vector for warping function
        processText = 'Assembling stiffness matrix and torsion load vector...'
        (shearK, torsionF) = (femFunctions.functionTimer(processText,
            self.assembleTorsionMatrices))

        # invert stiffness matrix
        processText = 'Inverting {} by {} stiffness matrix...'.format(shearK.shape[0], shearK.shape[0])
        invShearK = (femFunctions.functionTimer(processText,
            self.invertStiffnessMatrix, shearK))

        # ----------------------------------------------------------------------
        # TORSION PROPERTIES:
        # ----------------------------------------------------------------------
        # calculate warping constant and torsion constant
        self.omega = invShearK.dot(np.append(torsionF, 0))[:-1]
        self.J = ixx + iyy - self.omega.dot(shearK).dot(np.transpose(self.omega))

        # ----------------------------------------------------------------------
        # SHEAR PROPERTIES:
        # ----------------------------------------------------------------------
        # calculate shear vectors, shear centre integrals and warping moments
        processText = 'Assembling shear vectors and integrals...'
        ((shearFPsi, shearFPhi, shearCentreXInt, shearCentreYInt, Q_omega,
            i_omega, i_xomega, i_yomega)) = (femFunctions.functionTimer(processText,
            self.assembleShearVectors, ixx, iyy, ixy))

        # solve for shear functions
        self.Psi = invShearK.dot(np.append(shearFPsi, 0))[:-1]
        self.Phi = invShearK.dot(np.append(shearFPhi, 0))[:-1]

        # calculate shear centres (elasticity)
        self.Delta_s = 2 * (1 + self.nu) * (ixx * iyy - ixy ** 2)
        self.x_se = ((1 / self.Delta_s) * ((self.nu / 2 * shearCentreXInt) -
            torsionF.dot(self.Phi)))
        self.y_se = ((1 / self.Delta_s) * ((self.nu / 2 * shearCentreYInt) +
            torsionF.dot(self.Psi)))
        (self.x_1_se, self.y_2_se) = femFunctions.principalCoordinate(phi, self.x_se, self.y_se)

        # calculate shear centres (Trefftz's)
        self.x_st = (ixy * i_xomega - iyy * i_yomega) / (ixx * iyy - ixy ** 2)
        self.y_st = (ixx * i_xomega - ixy * i_yomega) / (ixx * iyy - ixy ** 2)

        # calculate shear deformation coefficients
        processText = 'Assembling shear deformation coefficients...'
        (kappa_x, kappa_y, kappa_xy) = (femFunctions.functionTimer(processText,
            self.assembleShearCoefficients, ixx, iyy, ixy))

        # calculate shear areas wrt global axis
        self.A_sx = self.Delta_s ** 2 / kappa_x
        self.A_sy = self.Delta_s ** 2 / kappa_y
        self.A_sxy = self.Delta_s ** 2 / kappa_xy

        # calculate shear areas wrt principal bending axis
        alpha_xx = kappa_x * self.geometricMesh.area / self.Delta_s ** 2
        alpha_yy = kappa_y * self.geometricMesh.area / self.Delta_s ** 2
        alpha_xy = kappa_xy * self.geometricMesh.area / self.Delta_s ** 2
        phi_rad = self.geometricMesh.phi * np.pi / 180
        R = np.array([[np.cos(phi_rad), np.sin(phi_rad)], [-np.sin(phi_rad), np.cos(phi_rad)]])
        rotatedAlpha = R.dot(np.array([[alpha_xx, alpha_xy],[alpha_xy,alpha_yy]])).dot(np.transpose(R))
        self.A_s11 = self.geometricMesh.area / rotatedAlpha[0,0]
        self.A_s22 = self.geometricMesh.area / rotatedAlpha[1,1]

        # calculate warping constant
        self.Gamma = (i_omega - Q_omega ** 2 / A - self.y_se * i_xomega +
            self.x_se * i_yomega)

    # --------------------------------------------------------------------------
    # SECTION PROPERTY METHODS:
    # --------------------------------------------------------------------------
    def invertStiffnessMatrix(self, K):
        Nvec1 = np.ones((K.shape[0], 1))
        Nvec2 = np.ones((1, K.shape[0] + 1))
        Nvec2[:,-1] = 0

        K = np.concatenate((K, Nvec1), axis=1)
        K = np.concatenate((K, Nvec2), axis=0)

        Kinv = np.linalg.inv(K)
        return Kinv

    def assembleTorsionMatrices(self):
        # initialise variables
        shearK = np.zeros((self.noNodes, self.noNodes))
        torsionF = np.zeros(self.noNodes)

        for el in self.elements:
            indxs = np.ix_(el.nodes, el.nodes)
            (elK, elF) = el.torsionProperties()
            shearK[indxs] += elK
            torsionF[el.nodes] += elF

        return (shearK, torsionF)

    def assembleShearVectors(self, ixx, iyy, ixy):
        # initialise variables
        shearFPsi = np.zeros(self.noNodes)
        shearFPhi = np.zeros(self.noNodes)
        shearCentreXInt = 0
        shearCentreYInt = 0
        Q_omega = 0
        i_omega = 0
        i_xomega = 0
        i_yomega = 0

        for el in self.elements:
            ((elShearFPsi, elShearFPhi, elShearCentreXInt, elShearCentreYInt, elQ_omega,
                elI_omega, elI_xomega, elI_yomega)) = (el.shearProperties(ixx, iyy,
                ixy, self.omega[el.nodes]))

            shearFPsi[el.nodes] += elShearFPsi
            shearFPhi[el.nodes] += elShearFPhi
            shearCentreXInt += elShearCentreXInt
            shearCentreYInt += elShearCentreYInt
            Q_omega += elQ_omega
            i_omega += elI_omega
            i_xomega += elI_xomega
            i_yomega += elI_yomega

        return ((shearFPsi, shearFPhi, shearCentreXInt, shearCentreYInt, Q_omega,
            i_omega, i_xomega, i_yomega))

    def assembleShearCoefficients(self, ixx, iyy, ixy):
        # initialise variables
        kappa_x = 0
        kappa_y = 0
        kappa_xy = 0

        for el in self.elements:
            (elKappa_x, elKappa_y, elKappa_xy) = (el.shearCoefficients(ixx,
                iyy, ixy, self.Psi[el.nodes], self.Phi[el.nodes]))

            kappa_x += elKappa_x
            kappa_y += elKappa_y
            kappa_xy += elKappa_xy

        return (kappa_x, kappa_y, kappa_xy)

    def computeAreaSegments(self, u, px, py):
        '''
        Computes the area above and below a line defined by unit vector u and
        point (px,py)
        '''
        # allocate area variables
        topArea = 0
        botArea = 0
        topQx = 0
        topQy = 0
        botQx = 0
        botQy = 0
        topCentroid = 0
        botCentroid = 0

        for el in self.elements:
            # calculate area of element and centroid
            (elArea, Qx, Qy) = el.areaProperties()

            if elArea != 0:
                elCentroid = [Qy / elArea, Qx / elArea]
            else:
                elCentroid = [0, 0]

            # determine location of element and allocate element areas and
            # first moments of area accordingly
            if (femFunctions.pointAboveLine(u, px, py, elCentroid[0], elCentroid[1])):
                topArea += elArea
                topQx += Qx
                topQy += Qy
            else:
                botArea += elArea
                botQx += Qx
                botQy += Qy

        if topArea != 0 and botArea != 0:
            topCentroid = np.array([topQy / topArea, topQx / topArea])
            botCentroid = np.array([botQy / botArea, botQx / botArea])

        return (topArea, botArea, topCentroid, botCentroid)

    def centroidalSectionModulii(self):
        # determine extreme values of the cartesian co-ordinates
        self.xmax = self.pointArray[:,0].max()
        self.xmin = self.pointArray[:,0].min()
        self.ymax = self.pointArray[:,1].max()
        self.ymin = self.pointArray[:,1].min()

        # evaluate section modulii
        self.zxx_plus = self.ixx_c / (self.ymax - self.cy)
        self.zxx_minus = self.ixx_c / (self.cy - self.ymin)
        self.zyy_plus = self.iyy_c / (self.xmax - self.cx)
        self.zyy_minus = self.iyy_c / (self.cx - self.xmin)

    def principalSectionModulii(self):
        # loop through all co-ordinates to determine extreme values
        for (i, vertex) in enumerate(self.pointArray):
            (x_1, y_2) = femFunctions.principalCoordinate(self.phi, vertex[0] - self.cx, vertex[1] - self.cy)

            if i == 0: # initialise min, max variables
                self.x_1max = x_1
                self.x_1min = x_1
                self.y_2max = y_2
                self.y_2min = y_2

            self.x_1max = max(self.x_1max, x_1)
            self.x_1min = min(self.x_1min, x_1)
            self.y_2max = max(self.y_2max, y_2)
            self.y_2min = min(self.y_2min, y_2)

        # evaluate principal section modulii
        self.z11_plus = self.i11_c / abs(self.y_2max)
        self.z11_minus = self.i11_c / abs(self.y_2min)
        self.z22_plus = self.i22_c / abs(self.x_1max)
        self.z22_minus = self.i22_c / abs(self.x_1min)

    # --------------------------------------------------------------------------
    # STRESS CALCULATION:
    # --------------------------------------------------------------------------
    def unitStress(self):
        # calculate stresses due to unit loading
        processText = 'Calculating cross-section stresses...'
        femFunctions.functionTimer(processText, self.calculateStress)

    def calculateStress(self):
        # load geometric propeties
        area = self.geometricMesh.area
        ixx = self.geometricMesh.ixx_c
        iyy = self.geometricMesh.iyy_c
        ixy = self.geometricMesh.ixy_c
        i11 = self.geometricMesh.i11_c
        i22 = self.geometricMesh.i22_c
        phi = self.geometricMesh.phi

        # allocate stress vectors
        self.sigma_zz_axial = np.zeros(self.noNodes)
        self.sigma_zz_bending_xx = np.zeros(self.noNodes)
        self.sigma_zz_bending_yy = np.zeros(self.noNodes)
        self.sigma_zz_bending_11 = np.zeros(self.noNodes)
        self.sigma_zz_bending_22 = np.zeros(self.noNodes)
        self.tau_zx_torsion = np.zeros(self.noNodes)
        self.tau_zy_torsion = np.zeros(self.noNodes)
        self.tau_zx_shear_x = np.zeros(self.noNodes)
        self.tau_zy_shear_x = np.zeros(self.noNodes)
        self.tau_zx_shear_y = np.zeros(self.noNodes)
        self.tau_zy_shear_y = np.zeros(self.noNodes)

        # allocate nodal count vector for nodal averaging
        node_count = np.zeros(self.noNodes)

        for el in self.elements:
            # evaluate stresses at nodes
            ((elSigma_zz_axial, elSigma_zz_bending_xx, elSigma_zz_bending_yy,
                elSigma_zz_bending_11, elSigma_zz_bending_22, elTau_zx_torsion,
                elTau_zy_torsion, elTau_shear_zx_x, elTau_shear_zy_x,
                elTau_shear_zx_y, elTau_shear_zy_y)) = (el.calculateStress(area,
                ixx, iyy, ixy, i11, i22, phi, self.omega[el.nodes], self.J,
                self.Psi[el.nodes], self.Phi[el.nodes], self.Delta_s))

            # add stresses to global vectors
            self.sigma_zz_axial[el.nodes] += elSigma_zz_axial[:,0]
            self.sigma_zz_bending_xx[el.nodes] += elSigma_zz_bending_xx
            self.sigma_zz_bending_yy[el.nodes] += elSigma_zz_bending_yy
            self.sigma_zz_bending_11[el.nodes] += elSigma_zz_bending_11
            self.sigma_zz_bending_22[el.nodes] += elSigma_zz_bending_22
            self.tau_zx_torsion[el.nodes] += elTau_zx_torsion
            self.tau_zy_torsion[el.nodes] += elTau_zy_torsion
            self.tau_zx_shear_x[el.nodes] += elTau_shear_zx_x
            self.tau_zy_shear_x[el.nodes] += elTau_shear_zy_x
            self.tau_zx_shear_y[el.nodes] += elTau_shear_zx_y
            self.tau_zy_shear_y[el.nodes] += elTau_shear_zy_y

            # increment the nodal count vector
            node_count[el.nodes] += 1

        # nodal averaging
        self.sigma_zz_axial *= 1 / node_count
        self.sigma_zz_bending_xx *= 1 / node_count
        self.sigma_zz_bending_yy *= 1 / node_count
        self.sigma_zz_bending_11 *= 1 / node_count
        self.sigma_zz_bending_22 *= 1 / node_count
        self.tau_zx_torsion *= 1 / node_count
        self.tau_zy_torsion *= 1 / node_count
        self.tau_zx_shear_x *= 1 / node_count
        self.tau_zy_shear_x *= 1 / node_count
        self.tau_zx_shear_y *= 1 / node_count
        self.tau_zy_shear_y *= 1 / node_count

    def evaluateSectionStress(self, Nzz, Mxx, Myy, M11, M22, Mzz, Vx, Vy):
        # scale unit stresses by design actions
        self.axialStress = self.sigma_zz_axial * Nzz
        self.bendingStress = (self.sigma_zz_bending_xx * Mxx +
            self.sigma_zz_bending_yy * Myy + self.sigma_zz_bending_11 * M11 +
            self.sigma_zz_bending_xx * M22)
        self.torsionStress_zx = self.tau_zx_torsion * Mzz
        self.torsionStress_zy = self.tau_zy_torsion * Mzz
        self.torsionStress = ((self.torsionStress_zx ** 2 +
            self.torsionStress_zy ** 2) ** 0.5)
        self.shearStress_zx = self.tau_zx_shear_x * Vx + self.tau_zx_shear_y * Vy
        self.shearStress_zy = self.tau_zy_shear_x * Vx + self.tau_zy_shear_y * Vy
        self.shearStress = ((self.shearStress_zx ** 2 +
            self.shearStress_zy ** 2) ** 0.5)

        # compute combined stresses
        self.sigma_zz = self.axialStress + self.bendingStress
        self.tau_zx = self.torsionStress_zx + self.shearStress_zx
        self.tau_zy = self.torsionStress_zy + self.shearStress_zy
        self.tau =  (self.tau_zx ** 2 + self.tau_zy ** 2) ** 0.5
        self.vonMises = ((self.sigma_zz ** 2 + 3 * (self.tau ** 2)) ** 0.5)

    # --------------------------------------------------------------------------
    # POST METHODS:
    # --------------------------------------------------------------------------
    def contourPlot(self, principalAxis=False, z=None, nodes=False, plotTitle='', centroids=False):
        plt.gca().set_aspect('equal')
        # ax = plt.subplot(111)
        plt.triplot(self.pointArray[:,0], self.pointArray[:,1], self.elementArray[:,0:3], lw=0.5, color='black')
        plt.title(plotTitle)
        plt.xlabel('x')
        plt.ylabel('y')

        if principalAxis:
            start_11 = femFunctions.globalCoordinate(self.geometricMesh.phi, self.geometricMesh.x_1min, 0)
            end_11 = femFunctions.globalCoordinate(self.geometricMesh.phi, self.geometricMesh.x_1max, 0)
            start_22 = femFunctions.globalCoordinate(self.geometricMesh.phi, 0, self.geometricMesh.y_2min)
            end_22 = femFunctions.globalCoordinate(self.geometricMesh.phi, 0, self.geometricMesh.y_2max)

            (plt.plot([start_11[0], end_11[0]], [start_11[1], end_11[1]], label='11 axis'))
            (plt.plot([start_22[0], end_22[0]], [start_22[1], end_22[1]], label='22 axis'))

        if centroids:
            plt.scatter(0, 0, facecolors='None', edgecolors='k', marker='o', s=100, label='Elastic Centroid')
            plt.scatter(self.geometricMesh.x_pc - self.geometricMesh.cx, self.geometricMesh.y_pc - self.geometricMesh.cy, c='k', marker='x', s=100, label='Global Plastic Centroid')
            plt.scatter(self.geometricMesh.x_1_pc - self.geometricMesh.cx, self.geometricMesh.y_2_pc - self.geometricMesh.cy, facecolors='None', edgecolors='k', marker='s', s=100, label='Principal Plastic Centroid')
            plt.scatter(self.x_se, self.y_se, c='k', marker='+', s=100, label='Shear Centre')
            plt.legend()

        if z is not None:
            cmap = cm.get_cmap(name = 'jet')
            # v = np.linspace(-10, 10, 15, endpoint=True)
            trictr = plt.tricontourf(self.pointArray[:,0], self.pointArray[:,1], self.elementArray[:,0:3], z, cmap=cmap)
            cbar = plt.colorbar(trictr, label='Stress')

        if nodes:
            plt.plot(self.pointArray[:,0], self.pointArray[:,1], 'ko', markersize = 1)

        plt.grid(True)
        plt.show()

    def quiverPlot(self, u, v, plotTitle=''):
        plt.figure()
        plt.gca().set_aspect('equal')
        plt.triplot(self.pointArray[:,0], self.pointArray[:,1], self.elementArray[:,0:3], lw=0.5, color='black')
        plt.title(plotTitle)
        c = np.hypot(u, v)
        cmap = cm.get_cmap(name='jet')
        quiv = plt.quiver(self.pointArray[:,0], self.pointArray[:,1], u, v, c, cmap=cmap)
        cbar = plt.colorbar(quiv)
        plt.show()

    def printGeometricResults(self):
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

    def printPlasticResults(self):
        print "-----------------------------"
        print "Plastic Properties"
        print "-----------------------------"
        print "x_pc = {}".format(self.x_pc)
        print "y_pc = {}".format(self.y_pc)
        print "Sxx = {}".format(self.Sxx)
        print "Syy = {}".format(self.Syy)
        print "x_1_pc = {}".format(self.x_1_pc)
        print "y_2_pc = {}".format(self.y_2_pc)
        print "S11 = {}".format(self.S11)
        print "S22 = {}".format(self.S22)
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
        print "x_1_s,e = {}".format(self.x_1_se)
        print "y_2_s,e = {}".format(self.y_2_se)
        print "A_s,x = {}".format(self.A_sx)
        print "A_s,y = {}".format(self.A_sy)
        print "A_s,11 = {}".format(self.A_s11)
        print "A_s,22 = {}".format(self.A_s22)
        print ""

def plasticCentroidAlgorithm(tol, maxIt, u, cx, cy, dmin, dmax, points, facets, holes, pointArray, elementArray):
    '''
    Algorithm to find plastic centroid (point at which top area = bot area):
        tol = convergence tolerance
        maxIt = maximum iterations
        u = unit vector in direction of axis
        (cx,cy) = centroid of section
        start = initial guess of the plastic axis location
        (dmin,dmax) = distance from centroid to extreme fibre of section
        points = input points list
        facets = input facets list
        holes = input holes list
        pointArray = np array containing mesh points
        elementArray = np array containing element vertices

        a_n = perpendicular distance from centroid to p.c.
    '''
     # initialise iteration variables
    areaConvergence_n = np.random.rand() * 0.01
    a_n1 = 0
    iterationCount = 0
    u_perp = np.array([u[1], u[0]]) # u vector rotated  90 degrees

    # algorithm
    while ((abs(areaConvergence_n) > tol or iterationCount < 3) and (iterationCount < maxIt)):
        if iterationCount < 3:
            # first two iterations uses a stepping approach
            a_n = a_n1 + areaConvergence_n * (dmax - dmin) / 5 * abs(areaConvergence_n) # compute new trial axis
        else:
            # secant method
            a_n = (a_n2 * areaConvergence_n - a_n1 * areaConvergence_n1) / (areaConvergence_n - areaConvergence_n1)

        # ensure trial axis is within section depth
        if a_n > dmax:
            a_n = 0.95 * dmax
        elif a_n < dmin:
            a_n = 0.95 * dmin

        # determine points (p1,p2) on trial axis
        p1 = np.array([cx + a_n * u_perp[0], cy + a_n * u_perp[1]])
        p2 = np.array([p1[0] + u[0], p1[1] + u[1]])

        # remesh with new trial axis included
        (points_new, facets_new) = (femFunctions.divideMesh(points[:],
            facets[:], pointArray, elementArray, p1[0], p1[1], p2[0], p2[1]))

        newMesh = femFunctions.createMesh(points_new, facets_new, holes, minAngle=None, qualityMeshing=False)

        # create triMesh object with new trial mesh
        meshTrial = triMesh(newMesh)
        # meshTrial.contourPlot(nodes=True)

        # calculate area above and below axis
        (topArea, botArea, topCentroid, botCentroid) = meshTrial.computeAreaSegments(u, p1[0], p1[1])

        # update convergence and solution data
        areaConvergence_n1 = areaConvergence_n
        areaConvergence_n = topArea / botArea - 1 # recalculate convergence
        a_n2 = a_n1
        a_n1 = a_n
        iterationCount += 1 # increment iterations

    return (a_n, topArea, botArea, topCentroid, botCentroid)
