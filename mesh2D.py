import numpy as np
import femFunctions

class tri3:
    '''
    Three noded triangle, input a 3 x 2 numpy array containing the co-ordinates of the vertices
    '''

    def __init__(self, vertices):
        self.xy = vertices # triangle vertex co-ordinates [3 x 2]
        self.area = self._area()
        self.Qx = self._Qx()
        self.Qy = self._Qy()
        self.centroid = self._centroid()
        self.ixx = self._ixx()
        self.iyy = self._iyy()

    def _area(self):
        '''
        Area of element: integral of the determinant of the jacobian over the element
        '''
        # jacobian is constant for a tri3 element, therefore use 1 point Gaussian integration
        gp = femFunctions.gaussPoints(1)
        (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp)
        return gp[0] * j  # Gauss point weight * determinant of jacobian at Gauss point

    def _Qx(self):
        '''
        First moment of area about the global x-axis:
            - integral of [N * y * det(J)] over element
        '''
        # N * y * det(J) is linear for a tri3 element, therefore use 1 point Gaussian integration
        gp = femFunctions.gaussPoints(1)
        (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp)
        return gp[0] * np.dot(np.transpose(N), self.xy[:,1]) * j

    def _Qy(self):
        '''
        First moment of area about the global y-axis:
            - integral of [N * x * det(J)] over element
        '''
        # N * x * det(J) is linear for a tri3 element, therefore use 1 point Gaussian integration
        gp = femFunctions.gaussPoints(1)
        (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp)
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
        gps = femFunctions.gaussPoints(3)

        # loop through each gauss point to evaluate the integral at each gauss point
        total = 0
        for gp in gps.T:
            (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp)
            total += gp[0] * (np.dot(np.transpose(N), self.xy[:,1])) ** 2 * j

        return total

    def _iyy(self):
        '''
        Second moment of area about the global y-axis:
            - integral of [(N * x)^2 * det(J)] over element
        '''
        # (N * x)^2 * det(J) is quadratic in 2D over a tri3 element, therefore use 3 point Gaussian integration
        gps = femFunctions.gaussPoints(3)

        # loop through each gauss point to evaluate the integral at each gauss point
        total = 0
        for gp in gps.T:
            (N, dN, j) = femFunctions.shapeFunction(np.transpose(self.xy), gp)
            total += gp[0] * (np.dot(np.transpose(N), self.xy[:,0])) ** 2 * j

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
        self.area = self._area()
        self.Qx = self._Qx()
        self.Qy = self._Qy()
        self.cx = self._cx()
        self.cy = self._cy()
        self.ixx_g = self._ixx_g()
        self.ixx_c = self._ixx_c()
        self.iyy_g = self._iyy_g()
        self.iyy_c = self._iyy_c()
        self.rx = self._rx()
        self.ry = self._ry()

    def _area(self):
        ''' Total area of the cross section'''
        total = 0

        for el in self.elements:
            total += el.area

        return total

    def _Qx(self):
        ''' First moment of area about the global x-axis'''
        total = 0

        for el in self.elements:
            total += el.Qx

        return total

    def _Qy(self):
        ''' First moment of area about the global y-axis'''
        total = 0

        for el in self.elements:
            total += el.Qy

        return total

    def _cx(self):
        ''' Cross-section x-position of centroid'''
        return self.Qy / self.area

    def _cy(self):
        ''' Cross-section y-position of centroid'''
        return self.Qx / self.area

    def _ixx_g(self):
        ''' Second moment of area of about the global x-axis'''
        total = 0

        for el in self.elements:
            total += el.ixx

        return total

    def _iyy_g(self):
        ''' Second moment of area of about the global y-axis'''
        total = 0

        for el in self.elements:
            total += el.iyy

        return total

    def _ixx_c(self):
        ''' Second moment of area of about the centroidal x-axis'''
        return self.ixx_g - self.Qx ** 2 / self.area

    def _iyy_c(self):
        ''' Second moment of area of about the centroidal y-axis'''
        return self.iyy_g - self.Qy ** 2 / self.area

    def _rx(self):
        return (self.ixx_c / self.area) ** 0.5

    def _ry(self):
        return (self.iyy_c / self.area) ** 0.5
