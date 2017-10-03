import femFunctions
import numpy as np

class tri3:
    '''
    Three noded triangle, input a 3 x 2 array containing the co-ordinates of the vertices
    '''

    def __init__(self, vertices, nodes):
        self.xy = vertices # triangle vertex co-ordinates [3 x 2]
        self.nodes = nodes # array of node numbers [1 x 3]
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
        self.torsionKe = self._torsionKe()
        self.torsionFe = self._torsionFe()

    def _area(self):
        '''
        Area of element: integral of the determinant of the jacobian over the element
        '''
        # jacobian is constant for a tri3 element, therefore use 1 point Gaussian integration
        gp = femFunctions.gaussPoints(1)[0] # Gauss point for 1 point Gaussian integration
        (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3') # shape functions etc. for 1 point Gaussian integration
        return gp[0] * j  # Gauss point weight * determinant of jacobian at Gauss point

    def _Qx(self):
        '''
        First moment of area about the global x-axis:
            - integral of [N * y * det(J)] over element
        '''
        # N * y * det(J) is linear for a tri3 element, therefore use 1 point Gaussian integration
        gp = femFunctions.gaussPoints(1)[0] # Gauss point for 1 point Gaussian integration
        (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3') # shape functions etc. for 1 point Gaussian integration
        return gp[0] * np.dot(np.transpose(N), self.xy[:,1]) * j

    def _Qy(self):
        '''
        First moment of area about the global y-axis:
            - integral of [N * x * det(J)] over element
        '''
        # N * x * det(J) is linear for a tri3 element, therefore use 1 point Gaussian integration
        gp = femFunctions.gaussPoints(1)[0] # Gauss point for 1 point Gaussian integration
        (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3') # shape functions etc. for 1 point Gaussian integration
        return gp[0] * np.dot(np.transpose(N), self.xy[:,0]) * j

    def _centroid(self):
        ''' Centroid of element '''
        return np.array([self.Qy / self.area, self.Qx / self.area])

    def _ixx(self):
        '''
        Second moment of area about the global x-axis:
            - integral of [(N * y)^2 * det(J)] over element
        '''
        # (N * y)^2 * det(J) is quadratic in 2D over a tri3 element, therefore use 3 point Gaussian integration
        gps = femFunctions.gaussPoints(3) # Gauss point for 1 point Gaussian integration
        # loop through each gauss point to evaluate the integral at each gauss point
        total = 0
        for gp in gps.T:
            (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3')
            total += gp[0] * (np.dot(np.transpose(N), self.xy[:,1])) ** 2 * j

        return total

    def _iyy(self):
        '''
        Second moment of area about the global y-axis:
            - integral of [(N * x)^2 * det(J)] over element
        '''
        # (N * x)^2 * det(J) is quadratic in 2D over a tri3 element, therefore use 3 point Gaussian integration
        gps = femFunctions.gaussPoints(3) # Gauss point for 1 point Gaussian integration
        # loop through each gauss point to evaluate the integral at each gauss point
        total = 0
        for gp in gps.T:
            (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3')
            total += gp[0] * (np.dot(np.transpose(N), self.xy[:,0])) ** 2 * j

        return total

    def _ixy(self):
        '''
        Product of inertia about the global xy-axis:
            - integral of [N * x * N * y * det(J)] over element
        '''
        # N * x * N * y * det(J) is quadratic in 2D over a tri3 element, therefore use 3 point Gaussian integration
        gps = femFunctions.gaussPoints(3) # Gauss point for 1 point Gaussian integration
        # loop through each gauss point to evaluate the integral at each gauss point
        total = 0
        for gp in gps.T:
            (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3')
            total += gp[0] * np.transpose(N).dot(self.xy[:,0]) * np.transpose(N).dot(self.xy[:,1]) * j

        return total

    def _torsionKe(self):
        '''
        Element stiffness matrix for solving for torsional warping constant
            - integral of [B' * B * det(J)] over element
        '''
        # B' * B * det(J) is constant in 2D over a tri3 element, therefore use 1 point Gaussian integration
        gp = femFunctions.gaussPoints(1)[0] # Gauss point for 1 point Gaussian integration
        (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3')
        return gp[0] * np.dot(dN, np.transpose(dN)) * j

    def _torsionFe(self):
        '''
        Element load vector for solving for torsional warping constant
            - integral of [B' * [N * y; -N * x] * det(J)] over element
        '''
        # B' * [Ny; -Nx] * det(J) is linear in 2D over a tri3 element, therefore use 1 point Gaussian integration
        gp = femFunctions.gaussPoints(1)[0] # Gauss point for 1 point Gaussian integration
        (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3')
        Nx = np.dot(np.transpose(N), self.xy[:,0])
        Ny = np.dot(np.transpose(N), self.xy[:,1])
        return gp[0] * np.dot(dN, np.transpose(np.array([Ny, -Nx]))) * j

    def torsionStress(self, w, J):
        '''
        Stress due to a unit twisting moment
            - tau = Mz / J [B * w - h] at integration points
        '''
        # B is constant over the 2D tri3 element, therefore evaluate at one integration point only
        gp = femFunctions.gaussPoints(1)[0] # Gauss point for 1 point Gaussian integration
        (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3')
        Nx = np.dot(np.transpose(N), self.xy[:,0])
        Ny = np.dot(np.transpose(N), self.xy[:,1])
        tau = 1 / J * (np.transpose(dN).dot(w) - np.transpose(np.array([Ny, -Nx])))
        # extrapolation to nodes results in the nodes taking the values at the gauss points
        return np.array([tau, tau, tau])
