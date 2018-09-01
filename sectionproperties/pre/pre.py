import meshpy.triangle as triangle


class Material:
    """a

    aaa

    :param string name: Material name
    :param float elastic_modulus: Material modulus of elasticity
    :param float poissons_ratio: Material Poisson's ratio
    :param float yield_strength: Material yield strength
    :param color: Material color for rendering
    :type color: :class:`matplotlib.colors`

    :cvar string name: Material name
    :cvar float elastic_modulus: Material modulus of elasticity
    :cvar float poissons_ratio: Material Poisson's ratio
    :cvar float yield_strength: Material yield strength
    :cvar color: Material color for rendering
    :vartype color: :class:`matplotlib.colors`
    """

    def __init__(self, name, elastic_modulus, poissons_ratio, yield_strength,
                 color='w'):
        """Inits the Material class"""

        self.name = name
        self.elastic_modulus = elastic_modulus
        self.poissons_ratio = poissons_ratio
        self.yield_strength = yield_strength
        self.color = color


def create_mesh(points, facets, holes, control_points, mesh_sizes):
    """Creates a quadratic triangular mesh using the meshpy module, which
    utilises the code 'Triangle', by Jonathan Shewchuk.

    :param points: List of points *(x, y)* defining the vertices of the
        cross-section
    :type points: list[tuple(float, float)]
    :param facets: List of point index pairs *(p1, p2)* defining the edges of
        the cross-section
    :type points: list[tuple(int, int)]
    :param holes: List of points *(x, y)* defining the locations of holes
        within the cross-section. If there are no holes, provide an empty list
        [].
    :type holes: list[tuple(float, float)]
    :param control_points: A list of points *(x, y)* that define different
        regions of the cross-section. A control point is an arbitrary point
        within a region enclosed by facets.
    :type control_points: list[tuple(float, float)]
    :param mesh_sizes: List of maximum element areas for each region defined by
        a control point
    :type mesh_sizes: list[float]

    :return: Object containing generated mesh data
    :rtype: :class:`meshpy.triangle.MeshInfo`
    """

    mesh = triangle.MeshInfo()  # create mesh info object
    mesh.set_points(points)  # set points
    mesh.set_facets(facets)  # set facets
    mesh.set_holes(holes)  # set holes

    # set regions
    mesh.regions.resize(len(control_points))  # resize regions list
    region_id = 0  # initialise region ID variable

    for (i, cp) in enumerate(control_points):
        mesh.regions[i] = [cp[0], cp[1], region_id, mesh_sizes[i]]
        region_id += 1

    mesh = triangle.build(
        mesh, min_angle=30, mesh_order=2, quality_meshing=True,
        attributes=True, volume_constraints=True)

    return mesh


# def divideMesh(points, facets, nodes, elements, x1, y1, x2, y2, d):
#     """
#     This method loops through each facet to check for an intersection point
#     with the line defined by (x1,y1) and (x2,y2). If so, a point is added to
#     the mesh at the interesection point. Facets are then added between the new
#     points if the new facet is within the mesh domain (and not within a hole).
#     """
#
#     # allocate lists for intersection points
#     xIntPoints = []  # list of x locations for intersection points
#     yIntPoints = []  # list of y locations for intersection points
#     facetIndices = []  # list of facet indices that have been intersected
#     numPoints = len(points)
#     tol = 1e-6 * d  # tolerance for zipping nodes
#
#     # find all intersections between the input line and the facets
#     (xIntPoints, yIntPoints, facetIndices) = determineFacetIntersections(
#         points, facets, x1, y1, x2, y2, tol)
#
#     # ordered list of point indices that lie along the intersection axis
#     facetIntersections = []
#
#     # build new facets along the axis of intersection
#     (points, newFacets, facetIntersections) = addNewFacets(
#         xIntPoints, yIntPoints, points, facets, nodes, elements, numPoints)
#
#     # reconstruct facet list and subdivide facets at all intersection points
#     finalFacets = rebuildFacetList(facets, facetIndices, facetIntersections,
#                                    numPoints)
#     # plotGeometry(points, newFacets, []) # plot for debugging purposes
#
#     return (points, finalFacets)
#
#
# def determineFacetIntersections(points, facets, x1, y1, x2, y2, tol):
#     """
#     This method determines all intersections between the 'facets' and line
#     defined by (x1,y1) and (x2,y2). Returned is a sorted list (by x position or
#     by y position if intersection line is the y-axis) of intersection locations
#     ('xIntPoints' and 'yIntPoints') and a list of facet indices that have been
#     intersected.
#     """
#
#     # allocate lists for intersection points
#     xIntPoints = []  # list of x locations for intersection points
#     yIntPoints = []  # list of y locations for intersection points
#     facetIndices = []  # list of facet indices that have been intersected
#
#     # calculate intersection points using determinant approach determine
#     # values governed by (x1,y1) & (x2,y2) only
#     numerator_11 = x1 * y2 - y1 * x2
#
#     # loop through each facet in the mesh to # find all intersections between
#     # the input line and the facets
#     for (i, line) in enumerate(facets):
#         # start and end points of the current facet
#         x3 = points[line[0]][0]
#         y3 = points[line[0]][1]
#         x4 = points[line[1]][0]
#         y4 = points[line[1]][1]
#
#         # calculate denominator for determinant approach
#         den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
#
#         # check to see if there is an intersection
#         if den != 0:
#             # determine remaining values
#             numerator_21 = x3 * y4 - y3 * x4
#
#             # determine intersection points
#             xInt = (numerator_11 * (x3 - x4) - (x1 - x2) * numerator_21) / den
#             yInt = (numerator_11 * (y3 - y4) - (y1 - y2) * numerator_21) / den
#
#             # check to see if points lie within extents of facet
#             if (min(x3, x4) - tol <= xInt <= max(x3, x4) + tol and
#                     min(y3, y4) - tol <= yInt <= max(y3, y4) + tol):
#
#                 # check to see if point already added
#                 isAdded = False
#                 for (j, xj) in enumerate(xIntPoints):
#                     if ((abs(xj - xInt) < tol) and (
#                             abs(yIntPoints[j] - yInt) < tol)):
#                         isAdded = True
#
#                 # if we are adding a new point
#                 if isAdded is False:
#                     # if the new point is very close to an end-point, take the
#                     # value of the end point (zip nodes)
#                     if abs(xInt - x3) < tol:
#                         xInt = x3
#                     elif abs(xInt - x4) < tol:
#                         xInt = x4
#                     if abs(yInt - y3) < tol:
#                         yInt = y3
#                     elif abs(yInt - y4) < tol:
#                         yInt = y4
#
#                     # add point to intersection list
#                     xIntPoints.append(xInt)
#                     yIntPoints.append(yInt)
#
#                     # add facet to intersection list
#                     facetIndices.append(i)
#
#     # sort intersection lists and facet list based on x or y value
#     if len(xIntPoints) > 0:
#         # if we are working with the y-axis
#         if x1 == x2:
#             # sort by y
#             (yIntPoints, xIntPoints, facetIndices) = (list(t) for t in zip(
#                 *sorted(zip(yIntPoints, xIntPoints, facetIndices))))
#         else:
#             # sort by x
#             (xIntPoints, yIntPoints, facetIndices) = (list(t) for t in zip(
#                 *sorted(zip(xIntPoints, yIntPoints, facetIndices))))
#
#     return (xIntPoints, yIntPoints, facetIndices)
#
#
# def addNewFacets(xIntPoints, yIntPoints, points, facets, nodes, elements,
#                  numPoints):
#     """
#     This method appends a list of new facets defined by the intersection points
#     to the existing list of 'facets' and also appends a list of new point
#     indices to the existing list 'points'. New facets are checked to ensure
#     that they lie within the mesh domain (i.e. do not lie within a hole).
#     Returned is the new list of points and facets and an ordered list of point
#     indices that lie along the interseciton axis ('facetIntersections').
#     """
#
#     # ordered list of point indices that lie along the intersection axis
#     facetIntersections = []
#
#     # loop through all found intersection points to build new facets along axis
#     for (i, pt) in enumerate(xIntPoints):
#         # add intersection points to geometry point list
#         points.append([pt, yIntPoints[i]])
#
#         # add index of new point to the facetIntersections list
#         facetIntersections.append(i)
#
#         # add connecting facets to facet list
#         if i != 0:  # start by joining new point 2 to new point 1
#             # check to see if midpoint of facet lies within an element of mesh,
#             # i.e. we are not in a hole
#             px = 0.5 * (xIntPoints[i] + xIntPoints[i - 1])
#             py = 0.5 * (yIntPoints[i] + yIntPoints[i - 1])
#
#             if (pointWithinElement(px, py, nodes, elements)):
#                 # add connecting facet along line of intersection
#                 facets.append([numPoints + i - 1, numPoints + i])
#
#     return (points, facets, facetIntersections)
#
#
# def rebuildFacetList(facets, facetIndices, facetIntersections, numPoints):
#     """
#     This method rebuilds the facet list by splitting intersected facets into
#     two new facets. Returned is a list of new facets ('newFacets').
#     """
#
#     newFacets = []  # allocate new facet list
#
#     # loop through all facets to reconstruct facet list and subdivide facets
#     # at all intersection points
#     for (i, facet) in enumerate(facets):
#         # assume that the facet has already been added (either in the original
#         # geometry, or added in the addNewFacets method)
#         facetOriginal = True
#
#         # loop through all facets that have been intersected
#         for (counter, j) in enumerate(facetIndices):
#             # if the current facet is being interesected with a new facet
#             if i == j:
#                 # subdivide this current facet into two facets:
#                 # facet[0] = start point of original facet
#                 # facet[1] = end point of original facet
#                 # facetIntersections[counter] = intersection point index
#                 newFacets.append([facet[0], numPoints +
#                                   facetIntersections[counter]])
#                 newFacets.append([numPoints + facetIntersections[counter],
#                                   facet[1]])
#                 facetOriginal = False
#
#         # if the current facet has not been intersected
#         if facetOriginal:
#             newFacets.append([facet[0], facet[1]])
#
#     return newFacets
#
#
# def pointWithinElement(px, py, nodes, elements):
#     """
#     This method determines whether the point (px,py) lies within any tri6
#     element defined by 'nodes' and 'elements'.
#     """
#
#     # loop through all elements in mesh
#     for el in elements:
#         x1 = nodes[el[0]][0]
#         y1 = nodes[el[0]][1]
#         x2 = nodes[el[1]][0]
#         y2 = nodes[el[1]][1]
#         x3 = nodes[el[2]][0]
#         y3 = nodes[el[2]][1]
#
#         # compute variables alpha, beta, gamma, to determine if current
#         # point lies within the triangular element
#         alpha = (((y2 - y3) * (px - x3) + (x3 - x2) * (py - y3)) /
#                  ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)))
#         beta = (((y3 - y1) * (px - x3) + (x1 - x3) * (py - y3)) /
#                 ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)))
#         gamma = 1.0 - alpha - beta
#
#         # if the point lies within an element
#         if alpha >= 0 and beta >= 0 and gamma >= 0:
#             return True
#
#     # if the point does not lie within an element
#     return False
#
#
# class LoadData:
#     """
#     This class parses the input load data and stores the load values.
#     """
#
#     def __init__(self, loads):
#         self.containsLoads = False
#         try:
#             self.Nzz = loads["nzz"]
#             self.containsLoads = True
#         except KeyError:
#             self.Nzz = 0.0
#
#         try:
#             self.Vx = loads["vx"]
#             self.containsLoads = True
#         except KeyError:
#             self.Vx = 0.0
#
#         try:
#             self.Vy = loads["vy"]
#             self.containsLoads = True
#         except KeyError:
#             self.Vy = 0.0
#
#         try:
#             self.Mxx = loads["mxx"]
#             self.containsLoads = True
#         except KeyError:
#             self.Mxx = 0.0
#
#         try:
#             self.Myy = loads["myy"]
#             self.containsLoads = True
#         except KeyError:
#             self.Myy = 0.0
#
#         try:
#             self.M11 = loads["m11"]
#             self.containsLoads = True
#         except KeyError:
#             self.M11 = 0.0
#
#         try:
#             self.M22 = loads["m22"]
#             self.containsLoads = True
#         except KeyError:
#             self.M22 = 0.0
#
#         try:
#             self.Mzz = loads["mzz"]
#             self.containsLoads = True
#         except KeyError:
#             self.Mzz = 0.0
#
#
# class CrossSectionSettings:
#     """
#     This class contains the settings used for the cross-section analysis.
#     """
#
#     def __init__(self, settings):
#         # load default settings
#         self.checkGeometry = True
#         self.checkMesh = True
#         self.outputLog = True
#         self.outputSettings = True
#         self.outputResults = True
#         self.plasticAnalysis = True
#         self.numberFormat = ".2f"
#         self.solverType = "cgs"
#         self.tol = 1e-5
#         self.plots = []
#
#         # load custom settings
#         self.applySettings(settings)
#
#     def applySettings(self, settings):
#         # read all valid settings from the dictionary settings
#         try:
#             testBool = (settings["general"]["check-geometry"].lower() in
#                         ["true"])
#             self.checkGeometry = testBool
#         except KeyError:
#             pass
#
#         try:
#             testBool = (settings["general"]["check-mesh"].lower() in
#                         ["true"])
#             self.checkMesh = testBool
#         except KeyError:
#             pass
#
#         try:
#             testBool = (settings["general"]["output-log"].lower() in
#                         ["true"])
#             self.outputLog = testBool
#         except KeyError:
#             pass
#
#         try:
#             testBool = (settings["general"]["output-settings"].lower() in
#                         ["true"])
#             self.outputSettings = testBool
#         except KeyError:
#             pass
#
#         try:
#             testBool = (settings["general"]["output-results"].lower() in
#                         ["true"])
#             self.outputResults = testBool
#         except KeyError:
#             pass
#
#         try:
#             testBool = (settings["general"]["plastic-analysis"].lower() in
#                         ["true"])
#             self.plasticAnalysis = testBool
#         except KeyError:
#             pass
#
#         try:
#             width = int(settings["number-format"]["width"])
#             precision = int(settings["number-format"]["precision"])
#             numType = str(settings["number-format"]["type"])
#             self.numberFormat = str(width) + "." + str(precision) + numType
#         except KeyError:
#             pass
#
#         try:
#             solverType = settings["solver"]["type"]
#
#             if (solverType.lower() == "cgs"):
#                 self.solverType = "cgs"
#             elif (solverType.lower() == "direct"):
#                 self.solverType = "direct"
#         except KeyError:
#             pass
#
#         try:
#             self.tol = settings["solver"]["tol"]
#         except KeyError:
#             pass
#
#         try:
#             self.plots = settings["plots"]
#         except KeyError:
#             pass
#
#     def printSettings(self):
#         """
#         This method prints the current settings to the console.
#         """
#
#         print("\n-----------------------------")
#         print("Program Settings")
#         print("-----------------------------")
#         print("General Settings:")
#         print("\tcheck-geometry:\t{}".format(self.checkGeometry))
#         print("\tcheck-mesh:\t{}".format(self.checkMesh))
#         print("\toutput-log:\t{}".format(self.outputLog))
#         print("\toutput-setting:\t{}".format(self.outputSettings))
#         print("\toutput-results:\t{}".format(self.outputResults))
#         print("Output Settings:")
#         print("\tnumber-format:\t{}".format(self.numberFormat))
#         print("Solver Settings:")
#         print("\ttype:\t\t{}".format(self.solverType))
#         print("\ttol:\t\t{}".format(self.tol))
#         print("Plot Settings:")
#         print("\tplots:\t\t{}\n".format(self.plots))
