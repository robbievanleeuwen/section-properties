import numpy as np
import femFunctions

class tri3:
    '''
    Three noded triangle, input a 3 x 2 array containing the co-ordinates of the vertices
    '''

    def __init__(self, vertices):
        self.xy = vertices # triangle vertex co-ordinates [3 x 2]
        self.initialise()

    def initialise(self):
        # compute properties of element
        gp1 = femFunctions.gaussPoints(1) # Gauss point for 1 point Gaussian integration
        gp3 = femFunctions.gaussPoints(3) # Gauss points for 3 point Gaussian integration
        (N1, dN1, j1) = femFunctions.shapeFunction(np.transpose(self.xy), gp1) # shape functions etc. for 1 point Gaussian integration

        self.area = self._area(gp1, j1)
        self.Qx = self._Qx(gp1, N1, j1)
        self.Qy = self._Qy(gp1, N1, j1)
        self.centroid = self._centroid()
        self.ixx = self._ixx(gp3)
        self.iyy = self._iyy(gp3)
        self.ixy = self._ixy(gp3)

    def _area(self, gp, j):
        '''
        Area of element: integral of the determinant of the jacobian over the element
        '''
        # jacobian is constant for a tri3 element, therefore use 1 point Gaussian integration
        return gp[0] * j  # Gauss point weight * determinant of jacobian at Gauss point

    def _Qx(self, gp, N, j):
        '''
        First moment of area about the global x-axis:
            - integral of [N * y * det(J)] over element
        '''
        # N * y * det(J) is linear for a tri3 element, therefore use 1 point Gaussian integration
        return gp[0] * np.dot(np.transpose(N), self.xy[:,1]) * j

    def _Qy(self, gp, N, j):
        '''
        First moment of area about the global y-axis:
            - integral of [N * x * det(J)] over element
        '''
        # N * x * det(J) is linear for a tri3 element, therefore use 1 point Gaussian integration
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
            (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp)
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
            (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp)
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
            (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp)
            total += gp[0] * np.transpose(N).dot(self.xy[:,0]) * np.transpose(N).dot(self.xy[:,1]) * j

        return total

class triMesh:
    '''
    triMesh object contains attribute .elements which an array of tri3 objects
    '''

    def __init__(self,genMesh):
        triElements = []

        for tri in genMesh['triangles']:
            x1 = genMesh['vertices'][tri[0]][0]
            y1 = genMesh['vertices'][tri[0]][1]
            x2 = genMesh['vertices'][tri[1]][0]
            y2 = genMesh['vertices'][tri[1]][1]
            x3 = genMesh['vertices'][tri[2]][0]
            y3 = genMesh['vertices'][tri[2]][1]
            vertex = np.array([[x1,y1], [x2,y2], [x3,y3]])
            triElements.append(tri3(vertex))

        self.elements = triElements
        self.initialise()

    def initialise(self):
        # initialise variables
        totalArea = totalQx = totalQy = totalIxx_g = totalIyy_g = totalIxy_g = 0

        # loop through all elements
        for el in self.elements:
            totalArea += el.area
            totalQx += el.Qx
            totalQy += el.Qy
            totalIxx_g += el.ixx
            totalIyy_g += el.iyy
            totalIxy_g += el.ixy

        self.area = totalArea
        self.Qx = totalQx
        self.Qy = totalQy
        self.cx = totalQy / totalArea
        self.cy = totalQx / totalArea
        self.ixx_g = totalIxx_g
        self.iyy_g = totalIyy_g
        self.ixy_g = totalIxy_g
        self.ixx_c = totalIxx_g - totalQx ** 2 / totalArea
        self.iyy_c = totalIyy_g - totalQy ** 2 / totalArea
        self.ixy_c = totalIxy_g - totalQx * totalQy / totalArea
        self.rx = (self.ixx_c / totalArea) ** 0.5
        self.ry = (self.iyy_c / totalArea) ** 0.5
