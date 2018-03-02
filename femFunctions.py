'''
This module stores the functions used for generic finite element calculations.
'''

import numpy as np
import meshpy.triangle as triangle


def gaussPoints(n):
    '''
    This function returns an [n x 4] matrix consisting of the weight and eta,
    xi and zeta locations for each Gauss point for a tri6 element.
    '''
    if n == 1:
        # one point gaussian integration
        return np.array([[1, 1.0 / 3, 1.0 / 3, 1.0 / 3]])
    elif n == 3:
        # three point gaussian integration
        return (np.array([[1.0 / 3, 2.0 / 3, 1.0 / 6, 1.0 / 6],
                          [1.0 / 3, 1.0 / 6, 2.0 / 3, 1.0 / 6],
                          [1.0 / 3, 1.0 / 6, 1.0 / 6, 2.0 / 3]]))
    elif n == 6:
        # six point gaussian integration
        g1 = 1.0 / 18 * (8 - np.sqrt(10) + np.sqrt(38 - 44 * np.sqrt(2.0 / 5)))
        g2 = 1.0 / 18 * (8 - np.sqrt(10) - np.sqrt(38 - 44 * np.sqrt(2.0 / 5)))
        w1 = (620 + np.sqrt(213125 - 53320 * np.sqrt(10))) / 3720
        w2 = (620 - np.sqrt(213125 - 53320 * np.sqrt(10))) / 3720
        return (np.array([[w2, 1 - 2 * g2, g2, g2],
                          [w2, g2, 1 - 2 * g2, g2],
                          [w2, g2, g2, 1 - 2 * g2],
                          [w1, g1, g1, 1 - 2 * g1],
                          [w1, 1 - 2 * g1, g1, g1],
                          [w1, g1, 1 - 2 * g1, g1]]))


def shapeFunction(xy, gaussPoint):
    '''
    This function computes shape functions and its derivative, and the
    determinant of the Jacobian matrix for a tri 6 element at a given
    Gauss point.
        INPUT:
        xy          = global co-ordinates of triangle vertics [2 x 6]
        gaussPoint  = isoparametric location of Gauss point [4 x 1]

        OUTPUT:
        N(i)        = value of shape function at given Gauss point [1 x 6]
        B(i,j)      = derivative of shape function i in direction j in global
        co-ordinate system [2 x 6]
        j            = determinant of the Jacobian matrix [1 x 1]
    '''
    B = np.zeros((2, 6))  # allocate B matrix
    # location of isoparametric co-ordinates for each Gauss point
    eta = gaussPoint[1]
    xi = gaussPoint[2]
    zeta = gaussPoint[3]

    # value of the shape functions
    N = (np.array([eta * (2 * eta - 1), xi * (2 * xi - 1),
                   zeta * (2 * zeta - 1), 4 * eta * xi, 4 * xi * zeta,
                   4 * eta * zeta]))
    # derivatives of the shape functions wrt the isoparametric co-ordinates
    B_iso = (np.array([[4 * eta - 1, 0, 0, 4 * xi, 0, 4 * zeta],
                       [0, 4 * xi - 1, 0, 4 * eta, 4 * zeta, 0],
                       [0, 0, 4 * zeta - 1, 0, 4 * xi, 4 * eta]]))

    # form Jacobian matrix
    J_upper = np.array([[1, 1, 1]])
    J_lower = np.dot(xy, np.transpose(B_iso))
    J = np.vstack((J_upper, J_lower))

    # calculate the jacobian -
    try:
        j = 0.5 * np.linalg.det(J)
    except ValueError:
        # handle warning if area is zero during plastic centroid algorithm
        j = 0

    # cacluate the P matrix
    if j != 0:
        P = np.dot(np.linalg.inv(J), np.array([[0, 0], [1, 0], [0, 1]]))
        # calculate the B matrix in terms of cartesian co-ordinates
        B = np.transpose(np.dot(np.transpose(B_iso), P))

    return (N, B, j)


def extrapolateToNodes(w):
    '''
    This function extrapolate results (w) at 6 Gauss points to 6 nodal points
    '''
    H_inv = np.array(
        [[1.87365927351160,	0.138559587411935, 0.138559587411935,
          -0.638559587411936, 0.126340726488397, -0.638559587411935],
         [0.138559587411935, 1.87365927351160, 0.138559587411935,
          -0.638559587411935, -0.638559587411935, 0.126340726488397],
         [0.138559587411935, 0.138559587411935, 1.87365927351160,
          0.126340726488396, -0.638559587411935, -0.638559587411935],
         [0.0749010751157440, 0.0749010751157440, 0.180053080734478,
          1.36051633430762,	-0.345185782636792, -0.345185782636792],
         [0.180053080734478, 0.0749010751157440, 0.0749010751157440,
          -0.345185782636792, 1.36051633430762, -0.345185782636792],
         [0.0749010751157440, 0.180053080734478,  0.0749010751157440,
          -0.345185782636792, -0.345185782636792, 1.36051633430762]])

    return H_inv.dot(w)


def createMesh(points, facets, holes=[], maxArea=[], minAngle=30, meshOrder=2,
               qualityMeshing=True):
    '''
    This function creates a quadratic triangular mesh using the meshpy module,
    which utilises the code 'Triangle', by Jonathan Shewchuk.
    '''
    info = triangle.MeshInfo()
    info.set_points(points)
    info.set_holes(holes)
    info.set_facets(facets)

    return triangle.build(
        info, max_volume=maxArea, min_angle=minAngle, mesh_order=meshOrder,
        quality_meshing=qualityMeshing)
