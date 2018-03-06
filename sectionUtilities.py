"""
This module provides a method to construct the geometry for a structural
cross-section comprising of a number of arbitrary base sections.
"""

import numpy as np


def sectionParse(sectionTypes, sectionData, settings):
    """
    Generates the geometry for the structural cross-section to be analysed,
    defined by a number of different sectionTypes, containing various
    sectionData. Note that there must be connectivity between all sections
    (i.e. there cannot be isolated sections) or the meshing and/or
    cross-section analysis will not work.
    """

    # initialise output variables
    points = []
    facets = []
    holes = []
    controlPoints = []

    # initialise pointCount variable
    pointCount = 0

    # loop through each section
    for (i, section) in enumerate(sectionTypes):
        # generate new section depending on section type
        if (section == "custom"):
            # load data from current sectionData
            try:
                pointData = sectionData[i]["points"]
                facetData = sectionData[i]["facets"]
                holeData = sectionData[i]["holes"]
                x = sectionData[i]["x"]
                y = sectionData[i]["y"]
                controlPointData = sectionData[i]["control-point"]
            except KeyError as err:
                handleKeyError(err, section)

            # generate a new section
            newSection = generateCustom(pointData, facetData, holeData, x, y,
                                        controlPointData)

        elif (section == "rectangle"):
            try:
                # load data from current sectionData
                d = sectionData[i]["d"]
                b = sectionData[i]["b"]
                x = sectionData[i]["x"]
                y = sectionData[i]["y"]
            except KeyError as err:
                handleKeyError(err, section)

            # if there is a control-point, load it
            try:
                controlPointData = sectionData[i]["control-point"]
            # if there is no control-point, set it to None
            except KeyError:
                controlPointData = None

            # generate a new section
            newSection = generateRectangle(d, b, x, y, controlPointData)

        elif (section == "circle"):
            # load data from current sectionData
            try:
                d = sectionData[i]["d"]
                n = sectionData[i]["n"]
                x = sectionData[i]["x"]
                y = sectionData[i]["y"]
            except KeyError as err:
                handleKeyError(err, section)

            # if there is a control-point, load it
            try:
                controlPointData = sectionData[i]["control-point"]
            # if there is no control-point, set it to None
            except KeyError:
                controlPointData = None

            # generate a new section
            newSection = generateCircle(d, n, x, y, controlPointData)

        elif (section == "chs"):
            # load data from current sectionData
            try:
                d = sectionData[i]["d"]
                t = sectionData[i]["t"]
                n = sectionData[i]["n"]
                x = sectionData[i]["x"]
                y = sectionData[i]["y"]
            except KeyError as err:
                handleKeyError(err, section)

            # if there is a control-point, load it
            try:
                controlPointData = sectionData[i]["control-point"]
            # if there is no control-point, set it to None
            except KeyError:
                controlPointData = None

            # generate a new section
            newSection = generateCHS(d, t, n, x, y, controlPointData)

        elif (section == "rhs"):
            # load data from current sectionData
            try:
                d = sectionData[i]["d"]
                b = sectionData[i]["b"]
                t = sectionData[i]["t"]
                r_out = sectionData[i]["r_out"]
                n_r = sectionData[i]["n_r"]
                x = sectionData[i]["x"]
                y = sectionData[i]["y"]
            except KeyError as err:
                handleKeyError(err, section)

            # if there is a control-point, load it
            try:
                controlPointData = sectionData[i]["control-point"]
            # if there is no control-point, set it to None
            except KeyError:
                controlPointData = None

            # generate a new section
            newSection = generateRHS(
                d, b, t, r_out, n_r, x, y, controlPointData)

        elif (section == "i-section"):
            # load data from current sectionData
            try:
                d = sectionData[i]["d"]
                b = sectionData[i]["b"]
                tf = sectionData[i]["tf"]
                tw = sectionData[i]["tw"]
                r = sectionData[i]["r"]
                n_r = sectionData[i]["n_r"]
                x = sectionData[i]["x"]
                y = sectionData[i]["y"]
            except KeyError as err:
                handleKeyError(err, section)

            # if there is a control-point, load it
            try:
                controlPointData = sectionData[i]["control-point"]
            # if there is no control-point, set it to None
            except KeyError:
                controlPointData = None

            # generate a new section
            newSection = generateISection(
                d, b, tf, tw, r, n_r, x, y, controlPointData)

        elif (section == "pfc"):
            # load data from current sectionData
            try:
                d = sectionData[i]["d"]
                b = sectionData[i]["b"]
                tf = sectionData[i]["tf"]
                tw = sectionData[i]["tw"]
                r = sectionData[i]["r"]
                n_r = sectionData[i]["n_r"]
                x = sectionData[i]["x"]
                y = sectionData[i]["y"]
            except KeyError as err:
                handleKeyError(err, section)

            # if there is a control-point, load it
            try:
                controlPointData = sectionData[i]["control-point"]
            # if there is no control-point, set it to None
            except KeyError:
                controlPointData = None

            # generate a new section
            newSection = generatePFCSection(
                d, b, tf, tw, r, n_r, x, y, controlPointData)

        elif (section == "tee"):
            # load data from current sectionData
            try:
                d = sectionData[i]["d"]
                b = sectionData[i]["b"]
                tf = sectionData[i]["tf"]
                tw = sectionData[i]["tw"]
                r = sectionData[i]["r"]
                n_r = sectionData[i]["n_r"]
                x = sectionData[i]["x"]
                y = sectionData[i]["y"]
            except KeyError as err:
                handleKeyError(err, section)

            # if there is a control-point, load it
            try:
                controlPointData = sectionData[i]["control-point"]
            # if there is no control-point, set it to None
            except KeyError:
                controlPointData = None

            # generate a new section
            newSection = generateTeeSection(
                d, b, tf, tw, r, n_r, x, y, controlPointData)

        elif (section == "angle"):
            # load data from current sectionData
            try:
                d = sectionData[i]["d"]
                b = sectionData[i]["b"]
                t = sectionData[i]["t"]
                r_root = sectionData[i]["r_root"]
                r_toe = sectionData[i]["r_toe"]
                n_r = sectionData[i]["n_r"]
                x = sectionData[i]["x"]
                y = sectionData[i]["y"]
            except KeyError as err:
                handleKeyError(err, section)

            # if there is a control-point, load it
            try:
                controlPointData = sectionData[i]["control-point"]
            # if there is no control-point, set it to None
            except KeyError:
                controlPointData = None

            # generate a new section
            newSection = generateAngleSection(
                d, b, t, r_root, r_toe, n_r, x, y, controlPointData)

        elif (section == "cee"):
            # load data from current sectionData
            try:
                d = sectionData[i]["d"]
                b = sectionData[i]["b"]
                lip = sectionData[i]["l"]
                t = sectionData[i]["t"]
                r_out = sectionData[i]["r_out"]
                n_r = sectionData[i]["n_r"]
                x = sectionData[i]["x"]
                y = sectionData[i]["y"]
            except KeyError as err:
                handleKeyError(err, section)

            # if there is a control-point, load it
            try:
                controlPointData = sectionData[i]["control-point"]
            # if there is no control-point, set it to None
            except KeyError:
                controlPointData = None

            # generate a new section
            newSection = generateCeeSection(
                d, b, lip, t, r_out, n_r, x, y, controlPointData)

        elif (section == "zed"):
            # load data from current sectionData
            try:
                d = sectionData[i]["d"]
                b1 = sectionData[i]["b1"]
                b2 = sectionData[i]["b2"]
                lip = sectionData[i]["l"]
                t = sectionData[i]["t"]
                r_out = sectionData[i]["r_out"]
                n_r = sectionData[i]["n_r"]
                x = sectionData[i]["x"]
                y = sectionData[i]["y"]
            except KeyError as err:
                handleKeyError(err, section)

            # if there is a control-point, load it
            try:
                controlPointData = sectionData[i]["control-point"]
            # if there is no control-point, set it to None
            except KeyError:
                controlPointData = None

            # generate a new section
            newSection = generateZedSection(
                d, b1, b2, lip, t, r_out, n_r, x, y, controlPointData)

        elif (section == "cruciform"):
            # load data from current sectionData
            try:
                d = sectionData[i]["d"]
                b = sectionData[i]["b"]
                t = sectionData[i]["t"]
                r = sectionData[i]["r"]
                n_r = sectionData[i]["n_r"]
                x = sectionData[i]["x"]
                y = sectionData[i]["y"]
            except KeyError as err:
                handleKeyError(err, section)

            # if there is a control-point, load it
            try:
                controlPointData = sectionData[i]["control-point"]
            # if there is no control-point, set it to None
            except KeyError:
                controlPointData = None

            # generate a new section
            newSection = generateCruciform(
                d, b, t, r, n_r, x, y, controlPointData)

        else:
            print("Error: section type '{}' is not defined.".format(section))
            quit()

        # get points, facets, holes and controlpoint from newSection
        (newPoints, newFacets, newHoles,
         newControlPoint) = newSection.returnSection()

        # loop through the facets in the newSection and append to the list
        for f in newFacets:
            facets.append([f[0] + pointCount, f[1] + pointCount])

        # loop through the points in the newSection and append to the list
        for p in newPoints:
            pointCount += 1
            points.append([p[0], p[1]])

        # loop through the holes in the newSection and append to the list
        for h in newHoles:
            holes.append([h[0], h[1]])

        # append the controlPoint from the newSection
        controlPoints.append([newControlPoint[0], newControlPoint[1]])

    if (settings.outputLog):
        print("-- Loaded {0} points, {1} facets and {2} holes ".format(
            len(points), len(facets), len(holes)) +
            "from {0} sections.".format(len(sectionTypes)))

    return (points, facets, holes, controlPoints)


def handleKeyError(err, section):
    """
    Displays an error message if the correct keys are not provided for the
    current section and quits the program.
    """

    print(
        "Error: Required key {0} not found for section type '{1}'.".format(
            err, section) +
        " Refer to the documentation for the required keys.")
    quit()


class CrossSection:
    """
    Parent class for a cross section generator. The cross section is translated
    by (x,y) before the section is returned in the method returnSection().
    """

    def __init__(self, x, y, controlPoint):
        self.x = x
        self.y = y
        self.controlPoint = controlPoint
        self.points = []
        self.facets = []
        self.holes = []

    def returnSection(self):
        self.translateSection()
        return (self.points, self.facets, self.holes, self.controlPoint)

    def translateSection(self):
        for point in self.points:
            point[0] += self.x
            point[1] += self.y

        for hole in self.holes:
            hole[0] += self.x
            hole[1] += self.y

        self.controlPoint[0] += self.x
        self.controlPoint[1] += self.y


class generateCustom(CrossSection):
    """
    Constructs a cross section from a list of points, facets and holes. The
    user must define the controlPoint (a point lying only within this section).
    """

    def __init__(self, points, facets, holes, x, y, controlPoint):
        super().__init__(x, y, controlPoint)
        self.points = points
        self.facets = facets
        self.holes = holes


class generateRectangle(CrossSection):
    """
    Constructs a rectangular section with depth d and width b. If the user does
    not specify a control point, the point (0,0) is used.
    """

    def __init__(self, d, b, x, y, controlPoint):
        super().__init__(x, y, controlPoint)

        # assign default control point if not specified
        if controlPoint is None:
            self.controlPoint = [0, 0]

        # construct the points and facets
        self.points = [[-0.5 * b, -0.5 * d], [0.5 * b, -0.5 * d],
                       [0.5 * b, 0.5 * d], [-0.5 * b, 0.5 * d]]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 0]]


class generateCircle(CrossSection):
    """
    Constructs a solid circle with diameter d, using n points to construct
    the circle. If the user does not specify a control point, the point
    (0,0) is used.
    """

    def __init__(self, d, n, x, y, controlPoint):
        super().__init__(x, y, controlPoint)

        # assign default control point if not specified
        if controlPoint is None:
            self.controlPoint = [0, 0]

        # loop through each point on the circle
        for i in range(n):
            # determine polar angle
            theta = i * 2 * np.pi * 1.0 / n

            # calculate location of the point
            x = 0.5 * d * np.cos(theta)
            y = 0.5 * d * np.sin(theta)

            # append the current point to the points list
            self.points.append([x, y])

            # if we are not at the last point
            if i != n - 1:
                self.facets.append([i, i + 1])
            # if we are at the last point, complete the circle
            else:
                self.facets.append([i, 0])


class generateCHS(CrossSection):
    """
    Constructs a circular hollow section with diameter d, thickness t, using
    n points to construct the inner and outer circles. If the user does not
    specify a control point, the point (d/2-t/2,0) is used.
    """

    def __init__(self, d, t, n, x, y, controlPoint):
        super().__init__(x, y, controlPoint)

        # assign default control point if not specified
        if controlPoint is None:
            self.controlPoint = [d * 0.5 - t * 0.5, 0]

        # specify a hole in the centre of the CHS
        self.holes = [[0, 0]]

        # loop through each point of the CHS
        for i in range(n):
            # determine polar angle
            theta = i * 2 * np.pi * 1.0 / n

            # calculate location of outer and inner points
            x_outer = 0.5 * d * np.cos(theta)
            y_outer = 0.5 * d * np.sin(theta)
            x_inner = (0.5 * d - t) * np.cos(theta)
            y_inner = (0.5 * d - t) * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_outer, y_outer])
            self.points.append([x_inner, y_inner])

            # if we are not at the last point
            if i != n - 1:
                self.facets.append([i * 2, i * 2 + 2])
                self.facets.append([i * 2 + 1, i * 2 + 3])
            # if we are at the last point, complete the circle
            else:
                self.facets.append([i * 2, 0])
                self.facets.append([i * 2 + 1, 1])


class generateRHS(CrossSection):
    """
    Constructs a rectangular hollow section with depth d, width b, thickness t,
    outer radius r_out, using n_r points to construct the inner and outer
    radii. If the user does not specify a control point, the point
    (b-t/2,d/2) is used.
    """

    def __init__(self, d, b, t, r_out, n_r, x, y, controlPoint):
        super().__init__(x, y, controlPoint)

        # assign default control point if not specified
        if controlPoint is None:
            self.controlPoint = [b - t * 0.5, d * 0.5]

        # specify a hole in the centre of the RHS
        self.holes = [[b * 0.5, d * 0.5]]

        r_in = r_out - t  # calculate internal radius

        # construct the bottom left radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi + i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_outer = r_out + r_out * np.cos(theta)
            y_outer = r_out + r_out * np.sin(theta)
            x_inner = r_out + r_in * np.cos(theta)
            y_inner = r_out + r_in * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_outer, y_outer])
            self.points.append([x_inner, y_inner])

        # construct the bottom right radius
        for i in range(n_r):
            # determine polar angle
            theta = 3.0 / 2 * np.pi + i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_outer = b - r_out + r_out * np.cos(theta)
            y_outer = r_out + r_out * np.sin(theta)
            x_inner = b - r_out + r_in * np.cos(theta)
            y_inner = r_out + r_in * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_outer, y_outer])
            self.points.append([x_inner, y_inner])

        # construct the top right radius
        for i in range(n_r):
            # determine polar angle
            theta = i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_outer = b - r_out + r_out * np.cos(theta)
            y_outer = d - r_out + r_out * np.sin(theta)
            x_inner = b - r_out + r_in * np.cos(theta)
            y_inner = d - r_out + r_in * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_outer, y_outer])
            self.points.append([x_inner, y_inner])

        # construct the top left radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi * 0.5 + i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_outer = r_out + r_out * np.cos(theta)
            y_outer = d - r_out + r_out * np.sin(theta)
            x_inner = r_out + r_in * np.cos(theta)
            y_inner = d - r_out + r_in * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_outer, y_outer])
            self.points.append([x_inner, y_inner])

        # build the facet list
        for i in range(int(len(self.points) / 2)):
            # if we are not at the last point
            if i != int(len(self.points) / 2 - 1):
                self.facets.append([i * 2, i * 2 + 2])
                self.facets.append([i * 2 + 1, i * 2 + 3])
            # if we are at the last point, complete the loop
            else:
                self.facets.append([i * 2, 0])
                self.facets.append([i * 2 + 1, 1])


class generateISection(CrossSection):
    """
    Constructs an I-section with depth d, width b, flange thickness tf, web
    thickness tw, root radius r, using n_r points to construct the root radius.
    If the user does not specify a control point, the point (b/2,d/2) is used.
    """

    def __init__(self, d, b, tf, tw, r, n_r, x, y, controlPoint):
        super().__init__(x, y, controlPoint)

        # assign default control point if not specified
        if controlPoint is None:
            self.controlPoint = [b * 0.5, d * 0.5]

        # add first three points
        self.points.append([0, 0])
        self.points.append([b, 0])
        self.points.append([b, tf])

        # construct the bottom right radius
        for i in range(n_r):
            # determine polar angle
            theta = 3.0 / 2 * np.pi * (1 - i * 1.0 / max(1, n_r - 1) * 1.0 / 3)

            # calculate the locations of the radius points
            x = b * 0.5 + tw * 0.5 + r + r * np.cos(theta)
            y = tf + r + r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # construct the top right radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi * (1 - i * 1.0 / max(1, n_r - 1) * 0.5)

            # calculate the locations of the radius points
            x = b * 0.5 + tw * 0.5 + r + r * np.cos(theta)
            y = d - tf - r + r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # add the next four points
        self.points.append([b, d - tf])
        self.points.append([b, d])
        self.points.append([0, d])
        self.points.append([0, d - tf])

        # construct the top left radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi * 0.5 * (1 - i * 1.0 / max(1, n_r - 1))

            # calculate the locations of the radius points
            x = b * 0.5 - tw * 0.5 - r + r * np.cos(theta)
            y = d - tf - r + r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # construct the bottom left radius
        for i in range(n_r):
            # determine polar angle
            theta = -np.pi * i * 1.0 / max(1, n_r - 1) * 0.5

            # calculate the locations of the radius points
            x = b * 0.5 - tw * 0.5 - r + r * np.cos(theta)
            y = tf + r + r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # add the last point
        self.points.append([0, tf])

        # build the facet list
        for i in range(len(self.points)):
            # if we are not at the last point
            if i != len(self.points) - 1:
                self.facets.append((i, i + 1))
            # if we are at the last point, complete the loop
            else:
                self.facets.append((len(self.points) - 1, 0))


class generatePFCSection(CrossSection):
    """
    Constructs a PFC section with depth d, width b, flange thickness tf, web
    thickness tw, root radius r, using n_r points to construct the root radius.
    If the user does not specify a control point, the point (tw/2,d/2) is used.
    """

    def __init__(self, d, b, tf, tw, r, n_r, x, y, controlPoint):
        super().__init__(x, y, controlPoint)

        # assign default control point if not specified
        if controlPoint is None:
            self.controlPoint = [tw * 0.5, d * 0.5]

        # add first three points
        self.points.append([0, 0])
        self.points.append([b, 0])
        self.points.append([b, tf])

        # construct the bottom right radius
        for i in range(n_r):
            # determine polar angle
            theta = 3.0 / 2 * np.pi * (1 - i * 1.0 / max(1, n_r - 1) * 1.0 / 3)

            # calculate the locations of the radius points
            x = tw + r + r * np.cos(theta)
            y = tf + r + r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # construct the top right radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi * (1 - i * 1.0 / max(1, n_r - 1) * 0.5)

            # calculate the locations of the radius points
            x = tw + r + r * np.cos(theta)
            y = d - tf - r + r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # add last three points
        self.points.append([b, d - tf])
        self.points.append([b, d])
        self.points.append([0, d])

        # build the facet list
        for i in range(len(self.points)):
            # if we are not at the last point
            if i != len(self.points) - 1:
                self.facets.append((i, i + 1))
            # if we are at the last point, complete the loop
            else:
                self.facets.append((len(self.points) - 1, 0))


class generateTeeSection(CrossSection):
    """
    Constructs a Tee section with depth d, width b, flange thickness tf, web
    thickness tw, root radius r, using n_r points to construct the root radius.
    If the user does not specify a control point, the point (b/2,d-tf/2) is
    used.
    """

    def __init__(self, d, b, tf, tw, r, n_r, x, y, controlPoint):
        super().__init__(x, y, controlPoint)

        # assign default control point if not specified
        if controlPoint is None:
            self.controlPoint = [b * 0.5, d - tf * 0.5]

        # add first two points
        self.points.append([b * 0.5 - tw * 0.5, 0])
        self.points.append([b * 0.5 + tw * 0.5, 0])

        # construct the top right radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi * (1 - i * 1.0 / max(1, n_r - 1) * 0.5)

            # calculate the locations of the radius points
            x = b * 0.5 + tw * 0.5 + r + r * np.cos(theta)
            y = d - tf - r + r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # add next four points
        self.points.append([b, d - tf])
        self.points.append([b, d])
        self.points.append([0, d])
        self.points.append([0, d - tf])

        # construct the top left radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi * 0.5 * (1 - i * 1.0 / max(1, n_r - 1))

            # calculate the locations of the radius points
            x = b * 0.5 - tw * 0.5 - r + r * np.cos(theta)
            y = d - tf - r + r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # build the facet list
        for i in range(len(self.points)):
            # if we are not at the last point
            if i != len(self.points) - 1:
                self.facets.append((i, i + 1))
            # if we are at the last point, complete the loop
            else:
                self.facets.append((len(self.points) - 1, 0))


class generateAngleSection(CrossSection):
    """
    Constructs an angle section with depth d, width b, thickness t, root radius
    r_root, toe radius r_toe using n_r points to construct the root radius.
    If the user does not specify a control point, the point (t/2,t/2) is used.
    """

    def __init__(self, d, b, t, r_root, r_toe, n_r, x, y, controlPoint):
        super().__init__(x, y, controlPoint)

        # assign default control point if not specified
        if controlPoint is None:
            self.controlPoint = [t * 0.5, t * 0.5]

        # add first two points
        self.points.append([0, 0])
        self.points.append([b, 0])

        # construct the bottom toe radius
        for i in range(n_r):
            # determine polar angle
            theta = i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate the locations of the radius points
            x = b - r_toe + r_toe * np.cos(theta)
            y = t - r_toe + r_toe * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # construct the root radius
        for i in range(n_r):
            # determine polar angle
            theta = 3.0 / 2 * np.pi * (1 - i * 1.0 / max(1, n_r - 1) * 1.0 / 3)

            # calculate the locations of the radius points
            x = t + r_root + r_root * np.cos(theta)
            y = t + r_root + r_root * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # construct the top toe radius
        for i in range(n_r):
            # determine polar angle
            theta = i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate the locations of the radius points
            x = t - r_toe + r_toe * np.cos(theta)
            y = d - r_toe + r_toe * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # add the next point
        self.points.append([0, d])

        # build the facet list
        for i in range(len(self.points)):
            # if we are not at the last point
            if i != len(self.points) - 1:
                self.facets.append((i, i + 1))
            # if we are at the last point, complete the loop
            else:
                self.facets.append((len(self.points) - 1, 0))


class generateCeeSection(CrossSection):
    """
    Constructs a Cee section with depth d, width b, lip l, thickness t, outer
    radius r_out, using n_r points to construct the root radius. If the user
    does not specify a control point, the point (t/2,d/2) is used.
    """

    def __init__(self, d, b, l, t, r_out, n_r, x, y, controlPoint):
        super().__init__(x, y, controlPoint)

        # assign default control point if not specified
        if controlPoint is None:
            self.controlPoint = [t * 0.5, d * 0.5]

        r_in = r_out - t  # calculate internal radius

        # construct the outer bottom left radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi + i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_outer = r_out + r_out * np.cos(theta)
            y_outer = r_out + r_out * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_outer, y_outer])

        # construct the outer bottom right radius
        for i in range(n_r):
            # determine polar angle
            theta = 3.0 / 2 * np.pi + i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_outer = b - r_out + r_out * np.cos(theta)
            y_outer = r_out + r_out * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_outer, y_outer])

        # add next two points
        self.points.append([b, l])
        self.points.append([b - t, l])

        # construct the inner bottom right radius
        for i in range(n_r):
            # determine polar angle
            theta = 2 * np.pi - i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_inner = b - r_out + r_in * np.cos(theta)
            y_inner = r_out + r_in * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_inner, y_inner])

        # construct the inner bottom left radius
        for i in range(n_r):
            # determine polar angle
            theta = 3.0 / 2 * np.pi - i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_inner = r_out + r_in * np.cos(theta)
            y_inner = r_out + r_in * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_inner, y_inner])

        # construct the inner top left radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi - i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_inner = r_out + r_in * np.cos(theta)
            y_inner = d - r_out + r_in * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_inner, y_inner])

        # construct the inner top right radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi * 0.5 - i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_inner = b - r_out + r_in * np.cos(theta)
            y_inner = d - r_out + r_in * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_inner, y_inner])

        # add next two points
        self.points.append([b - t, d - l])
        self.points.append([b, d - l])

        # construct the outer top right radius
        for i in range(n_r):
            # determine polar angle
            theta = i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_outer = b - r_out + r_out * np.cos(theta)
            y_outer = d - r_out + r_out * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_outer, y_outer])

        # construct the outer top left radius
        for i in range(n_r):
            # determine polar angle
            theta = 0.5 * np.pi + i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_outer = r_out + r_out * np.cos(theta)
            y_outer = d - r_out + r_out * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_outer, y_outer])

        # build the facet list
        for i in range(len(self.points)):
            # if we are not at the last point
            if i != len(self.points) - 1:
                self.facets.append((i, i + 1))
            # if we are at the last point, complete the loop
            else:
                self.facets.append((len(self.points) - 1, 0))


class generateZedSection(CrossSection):
    """
    Constructs a Zed section with depth d, left flange width b1, right flange
    width b2, lip l, thickness t, outer radius r_out, using n_r points to
    construct the root radius. If the user does not specify a control point,
    the point (t/2,d/2) is used.
    """

    def __init__(self, d, b1, b2, l, t, r_out, n_r, x, y, controlPoint):
        super().__init__(x, y, controlPoint)

        # assign default control point if not specified
        if controlPoint is None:
            self.controlPoint = [t * 0.5, d * 0.5]

        r_in = r_out - t  # calculate internal radius

        # construct the outer bottom left radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi + i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_outer = r_out + r_out * np.cos(theta)
            y_outer = r_out + r_out * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_outer, y_outer])

        # construct the outer bottom right radius
        for i in range(n_r):
            # determine polar angle
            theta = 3.0 / 2 * np.pi + i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_outer = b1 - r_out + r_out * np.cos(theta)
            y_outer = r_out + r_out * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_outer, y_outer])

        # add next two points
        self.points.append([b1, l])
        self.points.append([b1 - t, l])

        # construct the inner bottom right radius
        for i in range(n_r):
            # determine polar angle
            theta = 2 * np.pi - i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_inner = b1 - r_out + r_in * np.cos(theta)
            y_inner = r_out + r_in * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_inner, y_inner])

        # construct the inner bottom left radius
        for i in range(n_r):
            # determine polar angle
            theta = 3.0 / 2 * np.pi - i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_inner = r_out + r_in * np.cos(theta)
            y_inner = r_out + r_in * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_inner, y_inner])

        # construct the outer top right radius
        for i in range(n_r):
            # determine polar angle
            theta = i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_outer = t - r_out + r_out * np.cos(theta)
            y_outer = d - r_out + r_out * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_outer, y_outer])

        # construct the outer top left radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi * 0.5 + i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_outer = t - b2 + r_out + r_out * np.cos(theta)
            y_outer = d - r_out + r_out * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_outer, y_outer])

        # add the next two points
        self.points.append([t - b2, d - l])
        self.points.append([t - b2 + t, d - l])

        # construct the inner top left radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi - i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_inner = t - b2 + r_out + r_in * np.cos(theta)
            y_inner = d - r_out + r_in * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_inner, y_inner])

        # construct the inner top right radius
        for i in range(n_r):
            # determine polar angle
            theta = 0.5 * np.pi - i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_inner = t - r_out + r_in * np.cos(theta)
            y_inner = d - r_out + r_in * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_inner, y_inner])

        # build the facet list
        for i in range(len(self.points)):
            # if we are not at the last point
            if i != len(self.points) - 1:
                self.facets.append((i, i + 1))
            # if we are at the last point, complete the loop
            else:
                self.facets.append((len(self.points) - 1, 0))


class generateCruciform(CrossSection):
    """
    Constructs a cruciform section with depth d, width b, thickness t, root
    radius r, using n_r points to construct the root radius. If the user does
    not specify a control point, the point (0,0) is used.
    """

    def __init__(self, d, b, t, r, n_r, x, y, controlPoint):
        super().__init__(x, y, controlPoint)

        # assign default control point if not specified
        if controlPoint is None:
            self.controlPoint = [0, 0]

        # add first two points
        self.points.append([-t * 0.5, -d * 0.5])
        self.points.append([t * 0.5, -d * 0.5])

        # construct the bottom right radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi - i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate the locations of the radius points
            x = 0.5 * t + r + r * np.cos(theta)
            y = -0.5 * t - r + r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # add the next two points
        self.points.append([0.5 * b, -t * 0.5])
        self.points.append([0.5 * b, t * 0.5])

        # construct the top right radius
        for i in range(n_r):
            # determine polar angle
            theta = 1.5 * np.pi - i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate the locations of the radius points
            x = 0.5 * t + r + r * np.cos(theta)
            y = 0.5 * t + r + r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # add the next two points
        self.points.append([t * 0.5, 0.5 * d])
        self.points.append([-t * 0.5, 0.5 * d])

        # construct the top left radius
        for i in range(n_r):
            # determine polar angle
            theta = -i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate the locations of the radius points
            x = -0.5 * t - r + r * np.cos(theta)
            y = 0.5 * t + r + r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # add the next two points
        self.points.append([-0.5 * b, t * 0.5])
        self.points.append([-0.5 * b, -t * 0.5])

        # construct the bottom left radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi * 0.5 - i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate the locations of the radius points
            x = -0.5 * t - r + r * np.cos(theta)
            y = -0.5 * t - r + r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # build the facet list
        for i in range(len(self.points)):
            # if we are not at the last point
            if i != len(self.points) - 1:
                self.facets.append((i, i + 1))
            # if we are at the last point, complete the loop
            else:
                self.facets.append((len(self.points) - 1, 0))
