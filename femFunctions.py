import numpy as np

def gaussPoints(n):
    '''
    Returns an [n x 4] matrix consisting of the weight and eta, xi and zeta
    locations for each Gauss point
    '''
    if n == 1:
        return np.array([[1, 1.0/3, 1.0/3, 1.0/3]])
    elif n == 3:
        return (np.array([[1.0/3, 2.0/3, 1.0/6, 1.0/6],
                          [1.0/3, 1.0/6, 2.0/3, 1.0/6],
                          [1.0/3, 1.0/6, 1.0/6, 2.0/3]]))
    elif n == 6:
        g1 = 1.0 / 18 * (8 - np.sqrt(10) + np.sqrt(38 - 44 * np.sqrt(2.0 / 5)))
        g2 = 1.0 / 18 * (8 - np.sqrt(10) - np.sqrt(38 - 44 * np.sqrt(2.0 / 5)))
        w1 = (620 + np.sqrt(213125 - 53320 * np.sqrt(10))) / 3720
        w2 = (620 - np.sqrt(213125 - 53320 * np.sqrt(10))) / 3720
        return (np.array([[w1, 1 - 2 * g1, g1, g1],
                          [w1, g1, 1 - 2 * g1, g1],
                          [w1, g1, g1, 1 - 2 * g1],
                          [w2, 1 - 2 * g2, g2, g2],
                          [w2, g2, 1 - 2 * g2, g2],
                          [w2, g2, g2, 1 - 2 * g2]]))

def shapeFunction(xy, gaussPoint):
    '''
    Compute shape functions and determinant of the Jacobian matrix for an
    element at a given Gauss point
    INPUT:
        xy          = global co-ordinates of triangle vertics [2 x 6]
        gaussPoint  = isoparametric location of Gauss point [4 x 1]
    OUTPUT:
        N(i)    = value of shape function at given Gauss point [1 x 6]
        B(i,j) = derivative of shape function i in direction j in global
        co-ordinate system [2 x 6]
        j       = determinant of the Jacobian matrix [1 x 1]
    '''
    # location of isoparametric co-ordinates for each Gauss point
    eta  = gaussPoint[1]
    xi = gaussPoint[2]
    zeta = gaussPoint[3]

    # value of the shape functions
    N = (np.array([eta * (2 * eta - 1), xi * (2 * xi - 1),
            zeta * (2 * zeta - 1), 4 * eta * xi, 4 * xi * zeta, 4 * eta * zeta]))
    # derivatives of the shape functions with respect to the isoparametric co-ordinates
    B_iso = (np.array([[4 * eta - 1, 0, 0, 4 * xi, 0, 4 * zeta],
                       [0, 4 * xi - 1, 0, 4 * eta, 4 * zeta, 0],
                       [0, 0, 4 * zeta - 1, 0, 4 * xi, 4 * eta]]))

    # form Jacobian matrix
    J_upper = np.array([[1, 1, 1]])
    J_lower = np.dot(xy, np.transpose(B_iso))
    J = np.vstack((J_upper, J_lower))

    # calculate the jacobian
    j = 0.5 * np.linalg.det(J)

    if j < 0:
        print "Warning: negative Jacobian"

    # cacluate the P matrix
    P = np.dot(np.linalg.inv(J), np.array([[0, 0], [1, 0], [0, 1]]))
    # calculate the B matrix in terms of cartesian co-ordinates
    B = np.transpose(np.dot(np.transpose(B_iso), P))

    return (N, B, j)

# TODO:
def extrapolateToNodes(w, elementType, noGp):
    '''
    Extrapolate reults (w) at Gauss points to nodal points
    '''
    if elementType == 'tri3':
        if noGp == 1:
            return np.array([w, w, w])
        elif noGp == 3:
            H = np.array([[1, 1, -1], [1, -1, 1], [-1, 1, 1]])
            return H.dot(w)
    else:
        print 'Element type not yet programmed'
