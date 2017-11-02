import numpy as np
import meshpy.triangle as triangle

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
        return (np.array([[w2, 1 - 2 * g2, g2, g2],
                          [w2, g2, 1 - 2 * g2, g2],
                          [w2, g2, g2, 1 - 2 * g2],
                          [w1, g1, g1, 1 - 2 * g1],
                          [w1, 1 - 2 * g1, g1, g1],
                          [w1, g1, 1 - 2 * g1, g1]]))

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

def extrapolateToNodes(w):
    '''
    Extrapolate results (w) at 6 Gauss points to 6 nodal points
    '''
    H_inv = (np.array([[1.87365927351160,	0.138559587411935,	0.138559587411935,	-0.638559587411936,	0.126340726488397,	-0.638559587411935],
                       [0.138559587411935,	1.87365927351160,	0.138559587411935,	-0.638559587411935,	-0.638559587411935,	0.126340726488397],
                       [0.138559587411935,	0.138559587411935,	1.87365927351160,	0.126340726488396,	-0.638559587411935,	-0.638559587411935],
                       [0.0749010751157440,	0.0749010751157440,	0.180053080734478,	1.36051633430762,	-0.345185782636792,	-0.345185782636792],
                       [0.180053080734478,	0.0749010751157440,	0.0749010751157440,	-0.345185782636792,	1.36051633430762,	-0.345185782636792],
                       [0.0749010751157440,	0.180053080734478,	0.0749010751157440,	-0.345185782636792,	-0.345185782636792,	1.36051633430762]]))

    return H_inv.dot(w)

def createMesh(points, facets, holes=[], maxArea=[], minAngle=30, meshOrder=2, qualityMeshing=True):
    info = triangle.MeshInfo()
    info.set_points(points)
    info.set_holes(holes)
    info.set_facets(facets)

    return triangle.build(info, max_volume = maxArea, min_angle = minAngle, mesh_order = meshOrder, quality_meshing= qualityMeshing)

def shiftGeometry(points, holes, cx, cy):
    # initialise shifted points and holes lists
    shiftedPoints = []
    shiftedHoles =[]

    # shift points by centroid
    for point in points:
        shiftedPoints.append((point[0] - cx, point[1] - cy))

    for hole in holes:
        shiftedHoles.append((hole[0] - cx, hole[1] - cy))

    return (shiftedPoints, shiftedHoles)

def divideMesh(points, facets, pointArray, elementArray, x1, y1, x2, y2):
    '''
    Loops through each facet to check for an intersection point with the line defined
    by (x1,y1) and (x2,y2). If so, adds a point to the mesh and then adds facets
    between the added points.
    '''
    # allocate lists for intersection points
    xIntPoints = []
    yIntPoints = []
    facetIntersections = []
    facetIndices = []

    numPoints = len(points) # number of points

    # determine values governed by (x1,y1) & (x2,y2) only
    num_11 = x1 * y2 - y1 * x2

    for (i, line) in enumerate(facets):
        x3 = points[line[0]][0]
        y3 = points[line[0]][1]
        x4 = points[line[1]][0]
        y4 = points[line[1]][1]
        tol = 1e-6

        # calculate denominator
        den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)

        # check to see if there is an intersection
        if den != 0:
            # determine remaining values
            num21 = x3 * y4 - y3 * x4

            # determine intersecting points
            xInt = (num_11 * (x3 - x4) - (x1 - x2) * num21) / den
            yInt = (num_11 * (y3 - y4) - (y1 - y2) * num21) / den

            # check to see if points lie within segment
            if (min(x3,x4) - tol <= xInt <= max(x3,x4) + tol and
                min(y3,y4) - tol <= yInt <= max(y3,y4) + tol):
                xIntPoints.append(xInt)
                yIntPoints.append(yInt)
                facetIndices.append(i)

    # sort intersection points based on x value
    if len(xIntPoints) > 0:
        xIntPoints, yIntPoints, facetIndices = (list(t) for t in zip(*sorted(zip(xIntPoints, yIntPoints, facetIndices))))


    for (i, pt) in enumerate(xIntPoints):
        # add intersection points to point list
        points.append((pt, yIntPoints[i]))
        facetIntersections.append(i)

        # add facets to facet list
        if i != 0:
            # check to see if midpoint of facet lies within an element of mesh
            px = 0.5 * (xIntPoints[i] + xIntPoints[i - 1])
            py = 0.5 * (yIntPoints[i] + yIntPoints[i - 1])
            facetInDomain = False
            for tri in elementArray:
                x1 = pointArray[tri[0]][0]
                y1 = pointArray[tri[0]][1]
                x2 = pointArray[tri[1]][0]
                y2 = pointArray[tri[1]][1]
                x3 = pointArray[tri[2]][0]
                y3 = pointArray[tri[2]][1]

                alpha = (((y2 - y3) * (px - x3) + (x3 - x2) * (py - y3)) /
                    ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)))
                beta = (((y3 - y1) * (px - x3) + (x1 - x3) * (py - y3)) /
                    ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)))
                gamma = 1.0 - alpha - beta

                if alpha >= 0 and beta >= 0 and gamma >= 0:
                    facetInDomain = True

            if facetInDomain:
                facets.append((numPoints + i - 1, numPoints + i))

    newFacets = [] # allocate new facet list

    for (i, facet) in enumerate(facets):
        facetOriginal = True
        for (counter, j) in enumerate(facetIndices):
            if i == j:
                newFacets.append((facet[0], numPoints + facetIntersections[counter]))
                newFacets.append((numPoints + facetIntersections[counter], facet[1]))
                facetOriginal = False

        if facetOriginal:
            newFacets.append((facet[0], facet[1]))

    return (points, newFacets)

def lgMultSolve(K, f):
    Nvec1 = np.ones((K.shape[0], 1))
    Nvec2 = np.ones((1, K.shape[0] + 1))
    Nvec2[:,-1] = 0

    K = np.concatenate((K, Nvec1), axis=1)
    K = np.concatenate((K, Nvec2), axis=0)
    f = np.append(f, 0)

    u = np.linalg.solve(K, f)
    return (u[:-1], u[-1])

def principalCoordinate(u1, u2, cx, cy, x, y):
    '''
    Determines the coordinates of the point (x,y) in the principal axis system
    given unit vectors (u1,u2) defining the prinicpal axis and centroid (cx,cy).
    '''
    # vector from point to centroid
    PQ = np.array([cx - x, cy - y])
    # perpendicular distance from point to 1 and 2 axes
    d1 = np.linalg.norm(np.cross(PQ, u1))
    d2 = np.linalg.norm(np.cross(PQ, u2))

    # check to see if point is first quadrant
    if not (pointAboveLine(u1, cx, cy, x, y)): # point is below 1 axis
        d1 = -d1
    if not (pointAboveLine(u2, cx, cy, x, y)): # point is below 2 axis
        d2 = -d2

    return (d1, d2)

def pointAboveLine(u, px, py, x, y):
    '''
    Determines whether a point (x,y) is a above or below a line defined by
    unit vecotr u and point (px,py).
    '''
    # vector from point to point on line
    PQ = np.array([px - x, py - y])
    return np.cross(PQ, u) > 0

def functionTimer(function):
    start_time = time.clock()
    function()
    print("--- %s completed in %s seconds ---" % (function.__name__, time.clock() - start_time))
