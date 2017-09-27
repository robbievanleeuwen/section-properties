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
        gp1 = femFunctions.gaussPoints(1) # Gauss point for 1 point Gaussian integration
        gp3 = femFunctions.gaussPoints(3) # Gauss points for 3 point Gaussian integration

        self.area = self._area(gp1)
        self.Qx = self._Qx(gp1)
        self.Qy = self._Qy(gp1)
        self.centroid = self._centroid()
        self.ixx = self._ixx(gp3)
        self.iyy = self._iyy(gp3)
        self.ixy = self._ixy(gp3)
        self.torsionKe = self._torsionKe(gp1)
        self.torsionFe = self._torsionFe(gp1)

    def _area(self, gp):
        '''
        Area of element: integral of the determinant of the jacobian over the element
        '''
        # jacobian is constant for a tri3 element, therefore use 1 point Gaussian integration
        (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3') # shape functions etc. for 1 point Gaussian integration
        return gp[0] * j  # Gauss point weight * determinant of jacobian at Gauss point

    def _Qx(self, gp):
        '''
        First moment of area about the global x-axis:
            - integral of [N * y * det(J)] over element
        '''
        # N * y * det(J) is linear for a tri3 element, therefore use 1 point Gaussian integration
        (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3') # shape functions etc. for 1 point Gaussian integration
        return gp[0] * np.dot(np.transpose(N), self.xy[:,1]) * j

    def _Qy(self, gp):
        '''
        First moment of area about the global y-axis:
            - integral of [N * x * det(J)] over element
        '''
        # N * x * det(J) is linear for a tri3 element, therefore use 1 point Gaussian integration
        (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3') # shape functions etc. for 1 point Gaussian integration
        return gp[0] * np.dot(np.transpose(N), self.xy[:,0]) * j

    def _centroid(self):
        ''' Centroid of element '''
        return np.array([self.Qy / self.area, self.Qx / self.area])

    def _ixx(self, gps):
        '''
        Second moment of area about the global x-axis:
            - integral of [(N * y)^2 * det(J)] over element
        '''
        # (N * y)^2 * det(J) is quadratic in 2D over a tri3 element, therefore use 3 point Gaussian integration
        # loop through each gauss point to evaluate the integral at each gauss point
        total = 0
        for gp in gps.T:
            (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3')
            total += gp[0] * (np.dot(np.transpose(N), self.xy[:,1])) ** 2 * j

        return total

    def _iyy(self, gps):
        '''
        Second moment of area about the global y-axis:
            - integral of [(N * x)^2 * det(J)] over element
        '''
        # (N * x)^2 * det(J) is quadratic in 2D over a tri3 element, therefore use 3 point Gaussian integration
        # loop through each gauss point to evaluate the integral at each gauss point
        total = 0
        for gp in gps.T:
            (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3')
            total += gp[0] * (np.dot(np.transpose(N), self.xy[:,0])) ** 2 * j

        return total

    def _ixy(self, gps):
        '''
        Product of inertia about the global xy-axis:
            - integral of [N * x * N * y * det(J)] over element
        '''
        # N * x * N * y * det(J) is quadratic in 2D over a tri3 element, therefore use 3 point Gaussian integration
        # loop through each gauss point to evaluate the integral at each gauss point
        total = 0
        for gp in gps.T:
            (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3')
            total += gp[0] * np.transpose(N).dot(self.xy[:,0]) * np.transpose(N).dot(self.xy[:,1]) * j

        return total

    def _torsionKe(self, gp):
        '''
        Element stiffness matrix for solving for torsional warping constant
            - integral of [B' * B * det(J)] over element
        '''
        # B' * B * det(J) is constant in 2D over a tri3 element, therefore use 1 point Gaussian integration
        (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3')
        return gp[0] * np.dot(dN, np.transpose(dN)) * j

    def _torsionFe(self, gp):
        '''
        Element load vector for solving for torsional warping constant
            - integral of [B' * [N * y; -N * x] * det(J)] over element
        '''
        # B' * [Ny; -Nx] * det(J) is linear in 2D over a tri3 element, therefore use 1 point Gaussian integration
        (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp, 'tri3')
        Nx = np.dot(np.transpose(N), self.xy[:,0])
        Ny = np.dot(np.transpose(N), self.xy[:,1])
        return gp[0] * np.dot(dN, np.transpose(np.array([Ny, -Nx]))) * j
