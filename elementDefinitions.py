import femFunctions
import numpy as np

class tri6:
    '''
    Six noded triangle, input a 2 x 6 array containing the co-ordinates of the vertices
    '''

    def __init__(self, vertices, nodes, nu):
        self.xy = vertices # triangle vertex co-ordinates [2 x 6]
        self.nodes = nodes # array of node numbers [1 x 6]
        self.nu = nu # poissons ratio of material
        self.initialise()

    def initialise(self):
        # compute properties of element
        self.area = self._area()
        self.Qx = self._Qx()
        self.Qy = self._Qy()
        self.centroid = self._centroid()
        self.ixx = self._ixx()
        self.iyy = self._iyy()
        self.ixy = self._ixy()
        self.shearKe = self._shearKe()
        self.torsionFe = self._torsionFe()

    def _area(self):
        '''
        Area of element: integral of the determinant of the jacobian over the element
        '''
        # jacobian is constant for a superparametric tri6 element,
        # therefore use 1 point Gaussian integration
        gps = femFunctions.gaussPoints(1) # Gauss point for 1 point Gaussian integration
        # loop through each gauss point to evaluate the integral at each gauss point
        total = 0
        for gp in gps:
            (N, B, j) = femFunctions.shapeFunction(self.xy, gp)
            total += gp[0] * j

        return total

    def _Qx(self):
        '''
        First moment of area about the global x-axis:
            - integral of [N * y * det(J)] over element
        '''
        # N * y * det(J) is quadratic for a tri6 element,
        # therefore use 3 point Gaussian integration
        gps = femFunctions.gaussPoints(3) # Gauss points for 3 point Gaussian integration
        total = 0
        for gp in gps:
            (N, B, j) = femFunctions.shapeFunction(self.xy, gp)
            total += gp[0] * np.dot(N, np.transpose(self.xy[1,:])) * j

        return total

    def _Qy(self):
        '''
        First moment of area about the global y-axis:
            - integral of [N * x * det(J)] over element
        '''
        # N * x * det(J) is quadratic for a tri6 element,
        # therefore use 3 point Gaussian integration
        gps = femFunctions.gaussPoints(3) # Gauss points for 3 point Gaussian integration
        total = 0
        for gp in gps:
            (N, B, j) = femFunctions.shapeFunction(self.xy, gp)
            total += gp[0] * np.dot(N, np.transpose(self.xy[0,:])) * j

        return total

    def _centroid(self):
        ''' Centroid of element '''
        return np.array([self.Qy / self.area, self.Qx / self.area])

    def _ixx(self):
        '''
        Second moment of area about the global x-axis:
            - integral of [(N * y)^2 * det(J)] over element
        '''
        # (N * y)^2 * det(J) is quartic in 2D over a tri6 element,
        # therefore use 6 point Gaussian integration
        gps = femFunctions.gaussPoints(6) # Gauss points for 6 point Gaussian integration
        total = 0
        for gp in gps:
            (N, B, j) = femFunctions.shapeFunction(self.xy, gp)
            total += gp[0] * (np.dot(N, np.transpose(self.xy[1,:]))) ** 2 * j

        return total

    def _iyy(self):
        '''
        Second moment of area about the global y-axis:
            - integral of [(N * x)^2 * det(J)] over element
        '''
        # (N * y)^2 * det(J) is quartic in 2D over a tri6 element,
        # therefore use 6 point Gaussian integration
        gps = femFunctions.gaussPoints(6) # Gauss points for 6 point Gaussian integration
        total = 0
        for gp in gps:
            (N, B, j) = femFunctions.shapeFunction(self.xy, gp)
            total += gp[0] * (np.dot(N, np.transpose(self.xy[0,:]))) ** 2 * j

        return total

    def _ixy(self):
        '''
        Product of inertia about the global xy-axis:
            - integral of [N * y * N * x * det(J)] over element
        '''
        # N * y * N * x * det(J) is quartic in 2D over a tri6 element,
        # therefore use 6 point Gaussian integration
        gps = femFunctions.gaussPoints(6) # Gauss points for 6 point Gaussian integration
        total = 0
        for gp in gps:
            (N, B, j) = femFunctions.shapeFunction(self.xy, gp)
            total += (gp[0] * np.dot(N, np.transpose(self.xy[1,:])) *
                    np.dot(N, np.transpose(self.xy[0,:])) * j)

        return total

    def _shearKe(self):
        '''
        Element stiffness matrix for solving for shear functions and
        torsional warping constant
            - integral of [B' * B * det(J)] over element
        '''
        # B' * B * det(J) is quadratic in 2D over a tri6 element,
        # therefore use 3 point Gaussian integration
        gps = femFunctions.gaussPoints(3) # Gauss points for 3 point Gaussian integration
        total = 0
        for gp in gps:
            (N, B, j) = femFunctions.shapeFunction(self.xy, gp)
            total += gp[0] * np.dot(np.transpose(B), B) * j

        return total

    def _torsionFe(self):
        '''
        Element load vector for solving for the torsional warping constant
            - integral of [B' * [N * y; -N * x] * det(J)] over element
        '''
        # B' * [Ny; -Nx] * det(J) is cubic in 2D over a tri6 element,
        # therefore use 6 point Gaussian integration
        gps = femFunctions.gaussPoints(6) # Gauss points for 6 point Gaussian integration
        total = 0
        for gp in gps:
            (N, B, j) = femFunctions.shapeFunction(self.xy, gp)
            Nx = np.dot(N, np.transpose(self.xy[0,:]))
            Ny = np.dot(N, np.transpose(self.xy[1,:]))

            total += (gp[0] * np.dot(np.transpose(B),
                    np.transpose(np.array([Ny, -Nx]))) * j)

        return total

    # def shearFePsi(self, ixx, ixy):
    #     '''
    #     Element load vector for solving for shear function (Psi) for x-direction
    #         - integral of [nu/2 * B' * [d1; d2] + 2 * (1 + nu) * N' * (Ix * N * x - Ixy * N * y)] * det(J) over element
    #         where:
    #             - d1 = Ix * r - Ixy * q
    #             - d2 = Ixy * r + Ix * q
    #             - r = (N * x) ^ 2 - (N * y)^2
    #             - q = 2 * N * x * N * y
    #     '''
    #     # Function is quadratic in 2D over a tri3 element, therefore use 3 point Gaussian integration
    #     gps = femFunctions.gaussPoints(3) # Gauss points for 3 point Gaussian integration
    #     shearFePsi = np.zeros(3) # allocate shear function load vector
    #
    #     # loop through each gauss point to evaluate the integral at each gauss point
    #     for gp in gps:
    #         (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3')
    #         Nx = np.dot(np.transpose(N), self.xy[:,0])
    #         Ny = np.dot(np.transpose(N), self.xy[:,1])
    #         r = Nx ** 2 - Ny ** 2
    #         q = 2 * Nx * Ny
    #         d1 = ixx * r - ixy * q
    #         d2 = ixy * r + ixx * q
    #
    #         shearFePsi += gp[0] * (self.nu / 2 * dN.dot(np.array([d1, d2])) + 2 * (1 + self.nu) * N * (ixx * Nx - ixy * Ny)) * j
    #
    #     self._shearFePsi = shearFePsi
    #     return shearFePsi
    #
    # def shearFePhi(self, iyy, ixy):
    #     '''
    #     Element load vector for solving for shear function (Phi) for y-direction
    #         - integral of [nu/2 * B' * [h1; h2] + 2 * (1 + nu) * N' * (Iy * N * y - Ixy * N * x)] * det(J) over element
    #         where:
    #             - h1 = -Ixy * r + Iy * q
    #             - h2 = -Iy * r - Ixy * q
    #             - r = (N * x) ^ 2 - (N * y)^2
    #             - q = 2 * N * x * N * y
    #     '''
    #     # Function is quadratic in 2D over a tri3 element, therefore use 3 point Gaussian integration
    #     gps = femFunctions.gaussPoints(3) # Gauss points for 3 point Gaussian integration
    #     shearFePhi = np.zeros(3) # allocate shear load vector
    #
    #     # loop through each gauss point to evaluate the integral at each gauss point
    #     for gp in gps:
    #         (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3')
    #         Nx = np.dot(np.transpose(N), self.xy[:,0])
    #         Ny = np.dot(np.transpose(N), self.xy[:,1])
    #         r = Nx ** 2 - Ny ** 2
    #         q = 2 * Nx * Ny
    #         h1 = -ixy * r + iyy * q
    #         h2 = -iyy * r - ixy * q
    #
    #         shearFePhi += gp[0] * (self.nu / 2 * dN.dot(np.array([h1, h2])) + 2 * (1 + self.nu) * N * (iyy * Ny - ixy * Nx)) * j
    #
    #     self._shearFePhi = shearFePhi
    #     return shearFePhi
    #
    # def torsionStress(self, Mz, omega, J):
    #     '''
    #     Stress due to a unit twisting moment
    #         - tau = Mz / J [B * w - h] at integration points
    #     '''
    #     # B is constant over the 2D tri3 element, therefore evaluate at one integration point only
    #     gp = femFunctions.gaussPoints(1)[0] # Gauss point for 1 point Gaussian integration
    #     (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3')
    #     Nx = np.dot(np.transpose(N), self.xy[:,0])
    #     Ny = np.dot(np.transpose(N), self.xy[:,1])
    #
    #     tau_torsion = Mz / J * (np.transpose(dN).dot(omega) - np.transpose(np.array([Ny, -Nx])))
    #
    #     # extrapolate results to nodes
    #     return femFunctions.extrapolateToNodes(tau_torsion, 'tri3', 1)
    #
    # def shearStress(self, Vx, Vy, Psi, Phi, ixx, iyy, ixy):
    #     '''
    #     Stress due to a shear force Vx and Vy
    #         - tau = Vx / delta * (B * Psi - nu / 2 * [d1; d2]) + Vy / delta * (B * Phi - nu / 2 * [h1; h2]) at integration points
    #     '''
    #     # 3 integration points used to compute stress functions, therefore compute at each integration point
    #     gps = femFunctions.gaussPoints(3) # Gauss points for 3 point Gaussian integration
    #     tau_shear = np.zeros((3,2))
    #
    #     # loop through integration points
    #     for i, gp in enumerate(gps):
    #         (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3')
    #         Nx = np.dot(np.transpose(N), self.xy[:,0])
    #         Ny = np.dot(np.transpose(N), self.xy[:,1])
    #         r = Nx ** 2 - Ny ** 2
    #         q = 2 * Nx * Ny
    #         d1 = ixx * r - ixy * q
    #         d2 = ixy * r + ixx * q
    #         h1 = -ixy * r + iyy * q
    #         h2 = -iyy * r - ixy * q
    #         Delta = 2 * (1 + self.nu) * (ixx * iyy - ixy ** 2)
    #
    #         tau_shear[i,:] = Vx / Delta * (np.transpose(dN).dot(Psi) - self.nu / 2 * np.array([d1, d2])) + Vy / Delta * (np.transpose(dN).dot(Phi) - self.nu / 2 * np.array([h1, h2]))
    #
    #     # extrapolation to nodes results in the nodes taking the values at the gauss points
    #     tau_shear_zx = femFunctions.extrapolateToNodes(tau_shear[:,0], 'tri3', 3)
    #     tau_shear_zy = femFunctions.extrapolateToNodes(tau_shear[:,1], 'tri3', 3)
    #     return np.transpose(np.array([tau_shear_zx, tau_shear_zy]))
