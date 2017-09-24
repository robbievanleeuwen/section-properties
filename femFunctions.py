import numpy as np

def gaussPoints(n):
    '''
    Returns an [3 x n] matrix consisting of the weight, x-location and y-location for each Gauss point
    '''
    if n == 1:
        return np.transpose(np.array([0.5, 1.0/3, 1.0/3]))
    elif n == 3:
        return np.array([[1.0/6, 1.0/6, 1.0/6], [0, 0.5, 0.5], [0.5, 0, 0.5]])
    else:
        print "Number of Gauss points not compatible with element type"

def shapeFunction(xy, gaussPoint):
    '''
    Compute shape functions and determinant of the Jacobian matrix for a tri3 element at a given Gauss point
    INPUT:
        xy          = global co-ordinates of triangle vertics [2 x 3]
        gaussPoint  = isoparametric location of Gauss point [3 x 1]
    OUTPUT:
        N(i)    = value of shape function at given Gauss point [3 x 1]
        dN(i,j) = derivative of shape function i in direction j in global co-ordinate system [3 x 2]
        j       = determinant of the Jacobian matrix [1 x 1]
    '''
    # location of isoparametric co-ordinates for each Gauss point
    eta  = gaussPoint[1]
    zeta = gaussPoint[2]

    N = np.transpose(np.array([1.0 - eta - zeta, eta, zeta])) # value of the shape functions
    B = np.array([[-1,-1],[1,0],[0,1]]) # derivatives of the shape functions

    J = np.dot(xy, B) # form Jacobian matrix
    j = np.linalg.det(J) # determinant of Jacobain matrix

    if j < 0:
        print "Warning: negative Jacobian"

    dN = np.dot(B, np.linalg.inv(J)) # transform derivatives to global system

    return (N, dN, j)
