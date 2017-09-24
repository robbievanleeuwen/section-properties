import numpy as np

class tri3:
    '''
    Three noded triangle, input a 3 x 2 numpy array containing the co-ordinates of the vertices
    '''

    def __init__(self, vertices):
        self.xy = vertices
        self.area = self._area()
        self.centroid = self._centroid()
        self.ixx = self._ixx()
        self.iyy = self._iyy()

    def _area(self):
        return 0.5 * abs((self.xy[0,0] - self.xy[2,0]) * (self.xy[1,1] - self.xy[0,1]) - (self.xy[0,0] - self.xy[1,0]) * (self.xy[2,1] - self.xy[0,1]))

    def _centroid(self):
        return np.array([1.0/3 * (np.sum(self.xy[:,0])), 1.0/3 * (np.sum(self.xy[:,1]))])

    def _ixx(self):
        sortedY = self.xy[self.xy[:, 1].argsort()]

        x1 = sortedY[0,0]
        y1 = sortedY[0,1]
        x2 = sortedY[1,0]
        y2 = sortedY[1,1]
        x3 = sortedY[2,0]
        y3 = sortedY[2,1]
        x4 = sortedY[0,0] + (y2 - y1) / (y3 - y1) * (x3 - x1)
        y4 = sortedY[1,1]

        triA = np.array([[x1,y1],[x2,y2],[x4,y4]])
        triB = np.array([[x2,y2],[x3,y3],[x4,y4]])

        triA_h = triA[1,1] - triA[0,1]
        triB_h = triB[1,1] - triB[0,1]
        triA_b = abs(triA[1,0] - triA[2,0])
        triB_b = triA_b

        triA_ixx = triA_b * triA_h ** 3 / 36
        triB_ixx = triB_b * triB_h ** 3 / 36
        triA_cy = triA[1,1] - triA_h / 3
        triB_cy = triB[0,1] + triB_h / 3
        triA_A = 0.5 * abs((triA[0,0] - triA[2,0]) * (triA[1,1] - triA[0,1]) - (triA[0,0] - triA[1,0]) * (triA[2,1] - triA[0,1]))
        triB_A = 0.5 * abs((triB[0,0] - triB[2,0]) * (triB[1,1] - triB[0,1]) - (triB[0,0] - triB[1,0]) * (triB[2,1] - triB[0,1]))

        return triA_ixx + triB_ixx + triA_A * (triA_cy - self.centroid[1]) ** 2 + triB_A * (triB_cy - self.centroid[1]) ** 2

    def _iyy(self):
        sortedX = self.xy[self.xy[:, 0].argsort()]

        x1 = sortedX[0,0]
        y1 = sortedX[0,1]
        x2 = sortedX[1,0]
        y2 = sortedX[1,1]
        x3 = sortedX[2,0]
        y3 = sortedX[2,1]
        x4 = sortedX[1,0]
        y4 = sortedX[0,1] + (x2 - x1) / (x3 - x1) * (y3 - y1)

        triA = np.array([[x1,y1],[x2,y2],[x4,y4]])
        triB = np.array([[x2,y2],[x3,y3],[x4,y4]])

        triA_h = triA[1,0] - triA[0,0]
        triB_h = triB[1,0] - triB[0,0]
        triA_b = abs(triA[1,1] - triA[2,1])
        triB_b = triA_b

        triA_iyy = triA_b * triA_h ** 3 / 36
        triB_iyy = triB_b * triB_h ** 3 / 36
        triA_cx = triA[1,0] - triA_h / 3
        triB_cx = triB[0,0] + triB_h / 3
        triA_A = 0.5 * abs((triA[0,0] - triA[2,0]) * (triA[1,1] - triA[0,1]) - (triA[0,0] - triA[1,0]) * (triA[2,1] - triA[0,1]))
        triB_A = 0.5 * abs((triB[0,0] - triB[2,0]) * (triB[1,1] - triB[0,1]) - (triB[0,0] - triB[1,0]) * (triB[2,1] - triB[0,1]))

        return triA_iyy + triB_iyy + triA_A * (triA_cx - self.centroid[0]) ** 2 + triB_A * (triB_cx - self.centroid[0]) ** 2

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
        self.cx = self._cx()
        self.cy = self._cy()
        self.ixx = self._ixx()
        self.iyy = self._iyy()
        self.rx = self._rx()
        self.ry = self._ry()

    def _area(self):
        totalArea = 0

        for el in self.elements:
            totalArea += el.area

        return totalArea

    def _cx(self):
        num = 0

        for el in self.elements:
            num += el.area * el.centroid[0]

        return num / self.area

    def _cy(self):
        num = 0

        for el in self.elements:
            num += el.area * el.centroid[1]

        return num / self.area

    def _ixx(self):
        num = 0

        for el in self.elements:
            num += el.ixx + el.area * (el.centroid[1] - self.cy) ** 2

        return num

    def _iyy(self):
        num = 0

        for el in self.elements:
            num += el.iyy + el.area * (el.centroid[0] - self.cx) ** 2

        return num

    def _rx(self):
        return (self.ixx / self.area) ** 0.5

    def _ry(self):
        return (self.iyy / self.area) ** 0.5
