"""
Contains various methods used in the cross-section analysis program.
"""

import time
import numpy as np
import matplotlib.pyplot as plt

from CrossSectionAnalysis import CrossSectionAnalysis
from femUtilities import createMesh


def pcAlgorithm(tol, maxIt, u, dmin, dmax, points, facets, holes,
                controlPoints, nodes, elements, materials, dir):
    """
    Algorithm to find plastic centroid (point at which top area = bot area):
        INPUT:
        tol = convergence tolerance
        maxIt = maximum iterations
        u = unit vector in direction of axis
        (dmin,dmax) = distance from centroid to extreme fibre of section
        points = input points list
        facets = input facets list
        holes = input holes list
        pointArray = np array containing mesh points
        elementArray = np array containing element vertices
        dir = 1 or 2 depending on axis direction (x or y; 11 or 22)

        OUTPUT:
        a_n = perpendicular distance from centroid to p.c.
    """

    # initialise iteration variables
    areaConvergence_n = 0  # area convergence of step n
    areaConvergence_n1 = 0  # area convergence of step n - 1
    areaConvergence_n2 = 0  # area convergence of step n - 2
    a_n = 0  # distance from centroid to pc of step n
    a_n1 = 0  # distance from centroid to pc of step n - 1
    a_n2 = 0  # distance from centroid to pc of step n - 2
    iterationCount = 1

    # determine vectors perpendicular to the current axis
    if (dir == 1):
        u_perp = np.array([u[1], -u[0]])  # u vector rotated  -90 degrees
    elif (dir == 2):
        u_perp = np.array([-u[1], u[0]])  # u vector rotated  90 degrees

    # iterative algorithm
    while ((abs(areaConvergence_n) > tol or iterationCount < 3) and
           (iterationCount < maxIt)):
        if iterationCount < 3:
            # first two steps to setup secant method:
            # random number between -0.5 and 0.5 multiplied by 20% of the depth
            a_n = (np.random.rand() - 0.5) * 0.2 * (dmax - dmin)
        else:
            # secant method
            a_n = (a_n2 * areaConvergence_n1 - a_n1 * areaConvergence_n2) / (
                areaConvergence_n1 - areaConvergence_n2)

        # ensure trial axis is within section depth
        if a_n > dmax:
            a_n = dmax - 0.1 * abs(dmax - dmin)
        elif a_n < dmin:
            a_n = dmin + 0.1 * abs(dmax - dmin)

        # console reporting for debugging purposes
        # print("a_n = {}".format(a_n))
        # print("dmin = {}".format(dmin))
        # print("dmax = {}".format(dmax))

        # determine points (p1,p2) on trial axis
        p1 = [a_n * u_perp[0], a_n * u_perp[1]]
        p2 = [p1[0] + u[0], p1[1] + u[1]]
        # remesh with new trial axis included
        (newPoints, newFacets) = divideMesh(
            points.copy(), facets.copy(), nodes, elements,
            p1[0], p1[1], p2[0], p2[1], abs(dmax - dmin))

        ms = np.ones(len(controlPoints))  # dummy mesh sizes list
        mesh = createMesh(
            newPoints, newFacets, holes, controlPoints, meshSizes=ms,
            minAngle=None, meshOrder=2, qualityMeshing=False,
            volumeConstraints=False, settings=[])

        # create section analysis object with new mesh
        section = CrossSectionAnalysis(mesh, materials, settings=[])
        # section.contourPlot(nodes=True)  # plot for debugging purposes

        # calculate area above and below axis
        (topA, botA, topCen, botCen) = section.computeAreaSegments(
            u, p1[0], p1[1])

        # calculate area convergence
        areaConvergence_n = topA / botA - 1
        # console reporting for debugging purposes
        # print("convergence = {}".format(areaConvergence_n))
        # print("---")

        # update convergence and solution data
        areaConvergence_n2 = areaConvergence_n1
        areaConvergence_n1 = areaConvergence_n
        a_n2 = a_n1
        a_n1 = a_n
        iterationCount += 1

    if (abs(areaConvergence_n) > tol):
        print("WARNING: Plastic centroid algorithm did not converge for the " +
              "axis in the direction x:{0:.3f}; y:{1:.3f}\n".format(
                  u[0], u[1]))

    return (a_n, topA, botA, topCen, botCen)


def divideMesh(points, facets, nodes, elements, x1, y1, x2, y2, d):
    """
    This method loops through each facet to check for an intersection point
    with the line defined by (x1,y1) and (x2,y2). If so, a point is added to
    the mesh at the interesection point. Facets are then added between the new
    points if the new facet is within the mesh domain (and not within a hole).
    """

    # allocate lists for intersection points
    xIntPoints = []  # list of x locations for intersection points
    yIntPoints = []  # list of y locations for intersection points
    facetIndices = []  # list of facet indices that have been intersected
    numPoints = len(points)
    tol = 1e-6 * d  # tolerance for zipping nodes

    # find all intersections between the input line and the facets
    (xIntPoints, yIntPoints, facetIndices) = determineFacetIntersections(
        points, facets, x1, y1, x2, y2, tol)

    # ordered list of point indices that lie along the intersection axis
    facetIntersections = []

    # build new facets along the axis of intersection
    (points, newFacets, facetIntersections) = addNewFacets(
        xIntPoints, yIntPoints, points, facets, nodes, elements, numPoints)

    # reconstruct facet list and subdivide facets at all intersection points
    finalFacets = rebuildFacetList(facets, facetIndices, facetIntersections,
                                   numPoints)
    # plotGeometry(points, newFacets, []) # plot for debugging purposes

    return (points, finalFacets)


def determineFacetIntersections(points, facets, x1, y1, x2, y2, tol):
    """
    This method determines all intersections between the 'facets' and line
    defined by (x1,y1) and (x2,y2). Returned is a sorted list (by x position or
    by y position if intersection line is the y-axis) of intersection locations
    ('xIntPoints' and 'yIntPoints') and a list of facet indices that have been
    intersected.
    """

    # allocate lists for intersection points
    xIntPoints = []  # list of x locations for intersection points
    yIntPoints = []  # list of y locations for intersection points
    facetIndices = []  # list of facet indices that have been intersected

    # calculate intersection points using determinant approach determine
    # values governed by (x1,y1) & (x2,y2) only
    numerator_11 = x1 * y2 - y1 * x2

    # loop through each facet in the mesh to # find all intersections between
    # the input line and the facets
    for (i, line) in enumerate(facets):
        # start and end points of the current facet
        x3 = points[line[0]][0]
        y3 = points[line[0]][1]
        x4 = points[line[1]][0]
        y4 = points[line[1]][1]

        # calculate denominator for determinant approach
        den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)

        # check to see if there is an intersection
        if den != 0:
            # determine remaining values
            numerator_21 = x3 * y4 - y3 * x4

            # determine intersection points
            xInt = (numerator_11 * (x3 - x4) - (x1 - x2) * numerator_21) / den
            yInt = (numerator_11 * (y3 - y4) - (y1 - y2) * numerator_21) / den

            # check to see if points lie within extents of facet
            if (min(x3, x4) - tol <= xInt <= max(x3, x4) + tol and
                    min(y3, y4) - tol <= yInt <= max(y3, y4) + tol):

                # check to see if point already added
                isAdded = False
                for (j, xj) in enumerate(xIntPoints):
                    if ((abs(xj - xInt) < tol) and (
                            abs(yIntPoints[j] - yInt) < tol)):
                        isAdded = True

                # if we are adding a new point
                if isAdded is False:
                    # if the new point is very close to an end-point, take the
                    # value of the end point (zip nodes)
                    if abs(xInt - x3) < tol:
                        xInt = x3
                    elif abs(xInt - x4) < tol:
                        xInt = x4
                    if abs(yInt - y3) < tol:
                        yInt = y3
                    elif abs(yInt - y4) < tol:
                        yInt = y4

                    # add point to intersection list
                    xIntPoints.append(xInt)
                    yIntPoints.append(yInt)

                    # add facet to intersection list
                    facetIndices.append(i)

    # sort intersection lists and facet list based on x or y value
    if len(xIntPoints) > 0:
        # if we are working with the y-axis
        if x1 == x2:
            # sort by y
            (yIntPoints, xIntPoints, facetIndices) = (list(t) for t in zip(
                *sorted(zip(yIntPoints, xIntPoints, facetIndices))))
        else:
            # sort by x
            (xIntPoints, yIntPoints, facetIndices) = (list(t) for t in zip(
                *sorted(zip(xIntPoints, yIntPoints, facetIndices))))

    return (xIntPoints, yIntPoints, facetIndices)


def addNewFacets(xIntPoints, yIntPoints, points, facets, nodes, elements,
                 numPoints):
    """
    This method appends a list of new facets defined by the intersection points
    to the existing list of 'facets' and also appends a list of new point
    indices to the existing list 'points'. New facets are checked to ensure
    that they lie within the mesh domain (i.e. do not lie within a hole).
    Returned is the new list of points and facets and an ordered list of point
    indices that lie along the interseciton axis ('facetIntersections').
    """

    # ordered list of point indices that lie along the intersection axis
    facetIntersections = []

    # loop through all found intersection points to build new facets along axis
    for (i, pt) in enumerate(xIntPoints):
        # add intersection points to geometry point list
        points.append([pt, yIntPoints[i]])

        # add index of new point to the facetIntersections list
        facetIntersections.append(i)

        # add connecting facets to facet list
        if i != 0:  # start by joining new point 2 to new point 1
            # check to see if midpoint of facet lies within an element of mesh,
            # i.e. we are not in a hole
            px = 0.5 * (xIntPoints[i] + xIntPoints[i - 1])
            py = 0.5 * (yIntPoints[i] + yIntPoints[i - 1])

            if (pointWithinElement(px, py, nodes, elements)):
                # add connecting facet along line of intersection
                facets.append([numPoints + i - 1, numPoints + i])

    return (points, facets, facetIntersections)


def rebuildFacetList(facets, facetIndices, facetIntersections, numPoints):
    """
    This method rebuilds the facet list by splitting intersected facets into
    two new facets. Returned is a list of new facets ('newFacets').
    """

    newFacets = []  # allocate new facet list

    # loop through all facets to reconstruct facet list and subdivide facets
    # at all intersection points
    for (i, facet) in enumerate(facets):
        # assume that the facet has already been added (either in the original
        # geometry, or added in the addNewFacets method)
        facetOriginal = True

        # loop through all facets that have been intersected
        for (counter, j) in enumerate(facetIndices):
            # if the current facet is being interesected with a new facet
            if i == j:
                # subdivide this current facet into two facets:
                # facet[0] = start point of original facet
                # facet[1] = end point of original facet
                # facetIntersections[counter] = intersection point index
                newFacets.append([facet[0], numPoints +
                                  facetIntersections[counter]])
                newFacets.append([numPoints + facetIntersections[counter],
                                  facet[1]])
                facetOriginal = False

        # if the current facet has not been intersected
        if facetOriginal:
            newFacets.append([facet[0], facet[1]])

    return newFacets


def pointWithinElement(px, py, nodes, elements):
    """
    This method determines whether the point (px,py) lies within any tri6
    element defined by 'nodes' and 'elements'.
    """

    # loop through all elements in mesh
    for el in elements:
        x1 = nodes[el[0]][0]
        y1 = nodes[el[0]][1]
        x2 = nodes[el[1]][0]
        y2 = nodes[el[1]][1]
        x3 = nodes[el[2]][0]
        y3 = nodes[el[2]][1]

        # compute variables alpha, beta, gamma, to determine if current
        # point lies within the triangular element
        alpha = (((y2 - y3) * (px - x3) + (x3 - x2) * (py - y3)) /
                 ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)))
        beta = (((y3 - y1) * (px - x3) + (x1 - x3) * (py - y3)) /
                ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)))
        gamma = 1.0 - alpha - beta

        # if the point lies within an element
        if alpha >= 0 and beta >= 0 and gamma >= 0:
            return True

    # if the point does not lie within an element
    return False


def principalCoordinate(phi, x, y):
    """
    This method determines the coordinates of the cartesian point (x,y) in
    the principal axis system given an axis rotation angle phi.
    """

    # convert principal axis angle to radians
    phi_rad = phi * np.pi / 180

    # form rotation matrix
    R = np.array([[np.cos(phi_rad), np.sin(phi_rad)],
                  [-np.sin(phi_rad), np.cos(phi_rad)]])

    # calculate rotated x and y coordinates
    x_rotated = R.dot(np.array([x, y]))

    return (x_rotated[0], x_rotated[1])


def globalCoordinate(phi, x_1, y_2):
    '''
    This function determines the coordinates of the point (x_1,y_2) in the
    global axis system given rotation phi.
    '''
    # convert principal axis angle to radians
    phi_rad = phi * np.pi / 180
    # form transposed rotation matrix
    R = (np.array([[np.cos(phi_rad), -np.sin(phi_rad)], [np.sin(phi_rad),
                                                         np.cos(phi_rad)]]))
    # calculate rotated x_1 and y_2 coordinates
    x_rotated = R.dot(np.array([x_1, y_2]))

    return (x_rotated[0], x_rotated[1])


def pointAboveLine(u, px, py, x, y):
    """
    This method determines whether a point (x,y) is a above or below a line
    defined by unit vector u and point (px,py).
    """

    # vector from point to point on line
    PQ = np.array([px - x, py - y])
    return np.cross(PQ, u) > 0


def plotGeometry(points, facets, holes, controlPoints):
    """
    This function plots the geometry defined by points, facets, holes and
    controlPoints.
    """

    fig, ax = plt.subplots()
    plt.ion()  # interactive mode enabled
    plt.show()  # show the plot
    ax.set_aspect("equal")  # set the scale on the x and y axes equal

    # plot the title and axis labels
    ax.set_title("Input Geometry")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    for f in facets:
        # plot the facets
        ax.plot([points[f[0]][0], points[f[1]][0]],
                [points[f[0]][1], points[f[1]][1]],
                'ko-', markersize=2)

    for h in holes:
        # plot the holes
        ax.plot(h[0], h[1], 'rx', markerSize=5)

    for cp in controlPoints:
        # plot the controlPoitns
        ax.plot(cp[0], cp[1], 'bo', markerSize=5)

    ax.grid(True)
    plt.draw()  # render the figure
    plt.pause(0.001)

    return fig


def functionTimer(text, function, *args):
    """
    This method displays the message 'text' and returns the time taken for a
    function, with arguments args, to execute. The value returned by the
    timed function is also returned
    """

    start_time = time.time()

    if text != "":
        print(text)

    result = function(*args)

    if text != "":
        print("---- {0} completed in {1:.6f} seconds ---".format(
            function.__name__, time.time() - start_time))

    return result


class LoadData:
    """
    This class parses the input load data and stores the load values.
    """

    def __init__(self, loads):
        self.containsLoads = False
        try:
            self.Nzz = loads["nzz"]
            self.containsLoads = True
        except KeyError:
            self.Nzz = 0.0

        try:
            self.Vx = loads["vx"]
            self.containsLoads = True
        except KeyError:
            self.Vx = 0.0

        try:
            self.Vy = loads["vy"]
            self.containsLoads = True
        except KeyError:
            self.Vy = 0.0

        try:
            self.Mxx = loads["mxx"]
            self.containsLoads = True
        except KeyError:
            self.Mxx = 0.0

        try:
            self.Myy = loads["myy"]
            self.containsLoads = True
        except KeyError:
            self.Myy = 0.0

        try:
            self.M11 = loads["m11"]
            self.containsLoads = True
        except KeyError:
            self.M11 = 0.0

        try:
            self.M22 = loads["m22"]
            self.containsLoads = True
        except KeyError:
            self.M22 = 0.0

        try:
            self.Mzz = loads["mzz"]
            self.containsLoads = True
        except KeyError:
            self.Mzz = 0.0


class CrossSectionSettings:
    """
    This class contains the settings used for the cross-section analysis.
    """

    def __init__(self, settings):
        # load default settings
        self.checkGeometry = True
        self.checkMesh = True
        self.outputLog = True
        self.outputSettings = True
        self.outputResults = True
        self.plasticAnalysis = True
        self.numberFormat = ".2f"
        self.solverType = "cgs"
        self.tol = 1e-5
        self.plots = []

        # load custom settings
        self.applySettings(settings)

    def applySettings(self, settings):
        # read all valid settings from the dictionary settings
        try:
            testBool = (settings["general"]["check-geometry"].lower() in
                        ["true"])
            self.checkGeometry = testBool
        except KeyError:
            pass

        try:
            testBool = (settings["general"]["check-mesh"].lower() in
                        ["true"])
            self.checkMesh = testBool
        except KeyError:
            pass

        try:
            testBool = (settings["general"]["output-log"].lower() in
                        ["true"])
            self.outputLog = testBool
        except KeyError:
            pass

        try:
            testBool = (settings["general"]["output-settings"].lower() in
                        ["true"])
            self.outputSettings = testBool
        except KeyError:
            pass

        try:
            testBool = (settings["general"]["output-results"].lower() in
                        ["true"])
            self.outputResults = testBool
        except KeyError:
            pass

        try:
            testBool = (settings["general"]["plastic-analysis"].lower() in
                        ["true"])
            self.plasticAnalysis = testBool
        except KeyError:
            pass

        try:
            width = int(settings["number-format"]["width"])
            precision = int(settings["number-format"]["precision"])
            numType = str(settings["number-format"]["type"])
            self.numberFormat = str(width) + "." + str(precision) + numType
        except KeyError:
            pass

        try:
            solverType = settings["solver"]["type"]

            if (solverType.lower() == "cgs"):
                self.solverType = "cgs"
            elif (solverType.lower() == "direct"):
                self.solverType = "direct"
        except KeyError:
            pass

        try:
            self.tol = settings["solver"]["tol"]
        except KeyError:
            pass

        try:
            self.plots = settings["plots"]
        except KeyError:
            pass

    def printSettings(self):
        """
        This method prints the current settings to the console.
        """

        print("\n-----------------------------")
        print("Program Settings")
        print("-----------------------------")
        print("General Settings:")
        print("\tcheck-geometry:\t{}".format(self.checkGeometry))
        print("\tcheck-mesh:\t{}".format(self.checkMesh))
        print("\toutput-log:\t{}".format(self.outputLog))
        print("\toutput-setting:\t{}".format(self.outputSettings))
        print("\toutput-results:\t{}".format(self.outputResults))
        print("Output Settings:")
        print("\tnumber-format:\t{}".format(self.numberFormat))
        print("Solver Settings:")
        print("\ttype:\t\t{}".format(self.solverType))
        print("\ttol:\t\t{}".format(self.tol))
        print("Plot Settings:")
        print("\tplots:\t\t{}\n".format(self.plots))
