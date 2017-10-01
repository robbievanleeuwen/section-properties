import numpy as np
import elementDefinitions

class triMesh:
    '''
    Contains elements within the triangular mesh and computes and stores section properties
    '''

    def __init__(self, genMesh, elementType):
        triElements = [] # list holding all element objects

        if elementType == 'tri3':
            for tri in genMesh['triangles']:
                x1 = genMesh['vertices'][tri[0]][0]
                y1 = genMesh['vertices'][tri[0]][1]
                x2 = genMesh['vertices'][tri[1]][0]
                y2 = genMesh['vertices'][tri[1]][1]
                x3 = genMesh['vertices'][tri[2]][0]
                y3 = genMesh['vertices'][tri[2]][1]
                vertices = np.array([[x1,y1], [x2,y2], [x3,y3]])
                triElements.append(elementDefinitions.tri3(vertices, tri)) # add triangle to element list
        else:
            print 'Element type not programmed'

        self.elements = triElements # store element list in triMesh object
        self.noNodes = len((genMesh)['vertices']) # total number of nodes in mesh
        self.initialise(genMesh)

    def initialise(self, mesh):
        # initialise variables
        totalArea = totalQx = totalQy = totalIxx_g = totalIyy_g = totalIxy_g = 0
        torsionK = np.zeros((self.noNodes, self.noNodes))
        torsionF = np.transpose(np.zeros(self.noNodes))

        # loop through all elements
        for i, el in enumerate(self.elements):
            totalArea += el.area
            totalQx += el.Qx
            totalQy += el.Qy
            totalIxx_g += el.ixx
            totalIyy_g += el.iyy
            totalIxy_g += el.ixy

            # assemble torsion stiffness matrix and load vector
            indxs = np.ix_(el.nodes, el.nodes)
            torsionK[indxs] += el.torsionKe
            torsionF[el.nodes] += el.torsionFe

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
        self.w = np.linalg.solve(torsionK, torsionF)
        self.tf = torsionF
        self.J = self.ixx_g + self.iyy_g - self.w.dot(torsionK).dot(np.transpose(self.w))
