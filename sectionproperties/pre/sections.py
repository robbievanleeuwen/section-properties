import numpy as np
import matplotlib.pyplot as plt
import sectionproperties.pre.pre as pre
import sectionproperties.post.post as post


# TODO: ensure dimensions are floats

class Geometry:
    """Parent class for a cross-section geometry input.

    Provides an interface for the user to specify the geometry defining a
    cross-section. A method is provided for generating a triangular mesh, for
    translating the cross-section by *(x, y)* and for plotting the geometry.

    :cvar points: List of points *(x, y)* defining the vertices of the
        cross-section
    :vartype points: list[list[float, float]]
    :cvar facets: List of point index pairs *(p1, p2)* defining the edges of
        the cross-section
    :vartype facets: list[list[int, int]]
    :cvar holes: List of points *(x, y)* defining the locations of holes within
        the cross-section. If there are no holes, provide an empty list [].
    :vartype holes: list[list[float, float]]
    :cvar control_points: A list of points *(x, y)* that define different
        regions of the cross-section. A control point is an arbitrary point
        within a region enclosed by facets.
    :vartype control_points: list[list[float, float]]
    :cvar shift: Vector that shifts the cross-section by *(x, y)*
    :vartype shift: list[float, float]
    """

    def __init__(self, control_points, shift):
        """Inits the Geometry class."""

        self.control_points = control_points
        self.shift = shift
        self.points = []
        self.facets = []
        self.holes = []

    def create_mesh(self, mesh_sizes):
        """Creates a quadratic triangular mesh from the Geometry object.

        :param mesh_sizes: A list of maximum element areas corresponding to
            each region within the cross-section geometry.
        :type mesh_size: list[float]

        :return: Object containing generated mesh data
        :rtype: :class:`meshpy.triangle.MeshInfo`

        :raises AssertionError: If the number of mesh sizes does not match the
            number of regions

        The following example creates a circular cross-section with a diameter
        of 50 with 64 points, and generates a mesh with a maximum triangular
        area of 2.5::

            import sectionproperties.pre.sections as sections

            geometry = sections.CircularSection(d=50, n=64)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])

        ..  figure:: ../images/sections/circle_mesh.png
            :align: center
            :scale: 75 %

            Mesh generated from the above geometry.
        """

        str = "Number of mesh_sizes ({0}), ".format(len(mesh_sizes))
        str += "should match the number of regions "
        str += "({0}).".format(len(self.control_points))
        assert(len(mesh_sizes) == len(self.control_points)), str

        return pre.create_mesh(self.points, self.facets, self.holes,
                               self.control_points, mesh_sizes)

    def shift_section(self):
        """Shifts the cross-section parameters by the class variable vector
        *shift*."""

        for point in self.points:
            point[0] += self.shift[0]
            point[1] += self.shift[1]

        for hole in self.holes:
            hole[0] += self.shift[0]
            hole[1] += self.shift[1]

        for cp in self.control_points:
            cp[0] += self.shift[0]
            cp[1] += self.shift[1]

    def rotate_section(self, angle, rot_point=None):
        """Rotates the geometry and specified angle about a point. If the
        rotation point is not provided, rotates the section about the first
        control point in the list of control points of the
        :class:`~sectionproperties.pre.sections.Geometry` object.

        :param float angle: Angle (degrees) by which to rotate the section. A
            positive angle leads to a counter-clockwise rotation.
        :param rot_point: Point *(x, y)* about which to rotate the section
        :type rot_point: list[float, float]

        The following example rotates a 200UB25 section clockwise by 30
        degrees::

            import sectionproperties.pre.sections as sections

            geometry = sections.ISection(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8)
            geometry.rotate_section(angle=-30)
        """

        # convert angle to radians
        rot_phi = angle * np.pi / 180

        def get_r(pt1, pt2):
            """Returns the distance between two points."""

            return ((pt1[0] - pt2[0]) ** 2 + (pt1[1] - pt2[1]) ** 2) ** 0.5

        def get_phi(pt1, pt2):
            """Returns the angle between two points."""

            return np.arctan2(pt1[1] - pt2[1], pt1[0] - pt2[0])

        def rotate_point(pt, rot_point, rot_phi):
            """Rotates a point given a rotation point and rotation angle."""

            r = get_r(pt, rot_point)
            phi = get_phi(pt, rot_point)

            pt[0] = r * np.cos(phi + rot_phi) + rot_point[0]
            pt[1] = r * np.sin(phi + rot_phi) + rot_point[1]

        # use the first control point if no rotation point is specified
        if rot_point is None:
            rot_point = self.control_points[0]

        # rotate all the points
        for point in self.points:
            rotate_point(point, rot_point, rot_phi)

        # rotate all the holes
        for hole in self.holes:
            rotate_point(hole, rot_point, rot_phi)

        # rotate all the control points
        for cp in self.control_points:
            rotate_point(cp, rot_point, rot_phi)

    def mirror_section(self, axis='x', mirror_point=None):
        """Mirrors the geometry about a point on either the x or y-axis. If no
        point is provided, mirrors the geometry about the first control point
        in the list of control points of the
        :class:`~sectionproperties.pre.sections.Geometry` object.

        :param string axis: Axis about which to mirror the geometry, *'x'* or
            *'y'*
        :param mirror_point: Point about which to mirror the geometry *(x, y)*
        :type mirror_point: list[float, float]

        The following example mirrors a 200PFC section about the y-axis and the
        point (0, 0)::

            import sectionproperties.pre.sections as sections

            geometry = sections.PfcSection(d=200, b=75, t_f=12, t_w=6, r=12, n_r=8)
            geometry.mirror_section(axis='y', mirror_point=[0, 0])
        """

        # use the first control point if no mirror point is specified
        if mirror_point is None:
            mirror_point = self.control_points[0]

        # select the axis to mirror
        if axis == 'x':
            i = 1
        elif axis == 'y':
            i = 0
        else:
            raise RuntimeError("Enter a valid axis: 'x' or 'y'")

        # mirror all points
        for point in self.points:
            point[i] = 2 * mirror_point[i] - point[i]

        # mirror all holes
        for hole in self.holes:
            hole[i] = 2 * mirror_point[i] - hole[i]

        # mirror all control points
        for cp in self.control_points:
            cp[i] = 2 * mirror_point[i] - cp[i]

    def add_point(self, point):
        """Adds a point to the geometry and returns the added point id.

        :param point: Location of the point
        :type point: list[float, float]
        :return: Point id
        :rtype: int
        """

        self.points.append(point)
        return len(self.points) - 1

    def add_facet(self, facet):
        """Adds a facet to the geometry and returns the added facet id.

        :param facet: Point indices of the facet
        :type facet: list[float, float]
        :return: Facet id
        :rtype: int
        """

        self.facets.append(facet)
        return len(self.facets) - 1

    def add_hole(self, hole):
        """Adds a hole location to the geometry and returns the added hole id.

        :param hole: Location of the hole
        :type hole: list[float, float]
        :return: Hole id
        :rtype: int
        """

        self.holes.append(hole)
        return len(self.holes) - 1

    def add_control_point(self, control_point):
        """Adds a control point to the geometry and returns the added control
        point id.

        :param hole: Location of the control point
        :type hole: list[float, float]
        :return: Control point id
        :rtype: int
        """

        self.control_points.append(control_point)
        return len(self.control_points) - 1

    def clean_geometry(self, verbose=False):
        """Peforms a full clean on the geometry.

        :param bool verbose: If set to true, information related to the
            geometry cleaning process is printed to the terminal.

        ..  note:: Cleaning the geometry is always recommended when creating a
          merged section which may result in overlapping or intersecting
          facets, or duplicate nodes.
        """

        self = pre.GeometryCleaner(self, verbose).clean_geometry()

    def plot_geometry(self, ax=None, pause=True, labels=False):
        """Plots the geometry defined by the input section. If no axes object
        is supplied a new figure and axis is created.

        :param ax: Axes object on which the mesh is plotted
        :type ax: :class:`matplotlib.axes.Axes`
        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        :param bool labels: If set to true, node and facet labels are displayed

        The following example creates a CHS discretised with 64 points, with a
        diameter of 48 and thickness of 3.2, and plots the geometry::

            import sectionproperties.pre.sections as sections

            geometry = sections.Chs(d=48, t=3.2, n=64)
            geometry.plot_geometry()

        ..  figure:: ../images/sections/chs_geometry.png
            :align: center
            :scale: 75 %

            Geometry generated by the above example.
        """

        # if no axes object is supplied, create and setup the plot
        if ax is None:
            ax_supplied = False
            (fig, ax) = plt.subplots()
            post.setup_plot(ax, pause)
        else:
            ax_supplied = True

        for (i, f) in enumerate(self.facets):
            # plot the points and facets
            if i == 0:
                ax.plot([self.points[f[0]][0], self.points[f[1]][0]],
                        [self.points[f[0]][1], self.points[f[1]][1]],
                        'ko-', markersize=2, label='Points & Facets')
            else:
                ax.plot([self.points[f[0]][0], self.points[f[1]][0]],
                        [self.points[f[0]][1], self.points[f[1]][1]],
                        'ko-', markersize=2)

        for (i, h) in enumerate(self.holes):
            # plot the holes
            if i == 0:
                ax.plot(h[0], h[1], 'rx', markerSize=5, label='Holes')
            else:
                ax.plot(h[0], h[1], 'rx', markerSize=5)

        for (i, cp) in enumerate(self.control_points):
            # plot the control points
            if i == 0:
                ax.plot(cp[0], cp[1], 'bo', markerSize=5,
                        label='Control Points')
            else:
                ax.plot(cp[0], cp[1], 'bo', markerSize=5)

        # display the legend
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        # display the labels
        if labels:
            # plot node labels
            for (i, pt) in enumerate(self.points):
                ax.annotate(str(i), xy=pt, color='r')

            # plot facet labels
            for (i, fct) in enumerate(self.facets):
                pt1 = self.points[fct[0]]
                pt2 = self.points[fct[1]]
                xy = [(pt1[0] + pt2[0]) / 2, (pt1[1] + pt2[1]) / 2]

                ax.annotate(str(i), xy=xy, color='b')

        # if no axes object is supplied, finish the plot
        if not ax_supplied:
            post.finish_plot(ax, pause, title='Cross-Section Geometry')

    def calculate_extents(self):
        """Calculates the minimum and maximum x and y-values amongst the list
        of points.

        :return: Minimum and maximum x and y-values
            *(x_min, x_max, y_min, y_max)*
        :rtype: tuple(float, float, float, float)
        """

        # loop through all points
        for (i, pt) in enumerate(self.points):
            x = pt[0]
            y = pt[1]

            # initialise min, max variables
            if i == 0:
                x_min = x
                x_max = x
                y_min = y
                y_max = y

            # update the mins and maxs where necessary
            x_min = min(x_min, x)
            x_max = max(x_max, x)
            y_min = min(y_min, y)
            y_max = max(y_max, y)

        return (x_min, x_max, y_min, y_max)


class CustomSection(Geometry):
    """Constructs a cross-section from a list of points, facets, holes and a
    user specified control point.

    :param points: List of points *(x, y)* defining the vertices of the
        cross-section
    :type points: list[list[float, float]]
    :param facets: List of point index pairs *(p1, p2)* defining the edges of
        the cross-section
    :type facets: list[list[int, int]]
    :param holes: List of points *(x, y)* defining the locations of holes
        within the cross-section. If there are no holes, provide an empty list
        [].
    :type holes: list[list[float, float]]
    :param control_points: A list of points *(x, y)* that define different
        regions of the cross-section. A control point is an arbitrary point
        within a region enclosed by facets.
    :type control_points: list[list[float, float]]
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a hollow trapezium with a base width of 100,
    top width of 50, height of 50 and a wall thickness of 10. A mesh is
    generated with a maximum triangular area of 2.0::

        import sectionproperties.pre.sections as sections

        points = [[0, 0], [100, 0], [75, 50], [25, 50], [15, 10], [85, 10], [70, 40], [30, 40]]
        facets = [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4]]
        holes = [[50, 25]]
        control_points = [[5, 5]]

        geometry = sections.CustomSection(points, facets, holes, control_points)
        mesh = geometry.create_mesh(mesh_sizes=[2.0])

    ..  figure:: ../images/sections/custom_geometry.png
        :align: center
        :scale: 75 %

        Custom section geometry.

    ..  figure:: ../images/sections/custom_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, points, facets, holes, control_points, shift=[0, 0]):
        """Inits the CustomSection class."""

        super().__init__(control_points, shift)

        self.points = points
        self.facets = facets
        self.holes = holes

        self.shift_section()


class RectangularSection(Geometry):
    """Constructs a rectangular section with the bottom left corner at the
    origin *(0, 0)*, with depth *d* and width *b*.

    :param float d: Depth (y) of the rectangle
    :param float b: Width (x) of the rectangle
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a rectangular cross-section with a depth of
    100 and width of 50, and generates a mesh with a maximum triangular area of
    5::

        import sectionproperties.pre.sections as sections

        geometry = sections.RectangularSection(d=100, b=50)
        mesh = geometry.create_mesh(mesh_sizes=[5])

    ..  figure:: ../images/sections/rectangle_geometry.png
        :align: center
        :scale: 75 %

        Rectangular section geometry.

    ..  figure:: ../images/sections/rectangle_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, d, b, shift=[0, 0]):
        """Inits the RectangularSection class."""

        # assign control point
        control_points = [[0.5 * b, 0.5 * d]]

        super().__init__(control_points, shift)

        # construct the points and facets
        self.points = [[0, 0], [b, 0], [b, d], [0, d]]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 0]]

        self.shift_section()


class CircularSection(Geometry):
    """Constructs a solid circle centered at the origin *(0, 0)* with diameter
    *d* and using *n* points to construct the circle.

    :param float d: Diameter of the circle
    :param int n: Number of points discretising the circle
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a circular cross-section with a diameter of
    50 with 64 points, and generates a mesh with a maximum triangular area of
    2.5::

        import sectionproperties.pre.sections as sections

        geometry = sections.CircularSection(d=50, n=64)
        mesh = geometry.create_mesh(mesh_sizes=[2.5])

    ..  figure:: ../images/sections/circle_geometry.png
        :align: center
        :scale: 75 %

        Circular section geometry.

    ..  figure:: ../images/sections/circle_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, d, n, shift=[0, 0]):
        """Inits the CircularSection class."""

        # assign control point
        control_points = [[0, 0]]

        super().__init__(control_points, shift)

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

        self.shift_section()


class Chs(Geometry):
    """Constructs a circular hollow section centered at the origin *(0, 0)*,
    with diameter *d* and thickness *t*, using *n* points to construct the
    inner and outer circles.

    :param float d: Outer diameter of the CHS
    :param float t: Thickness of the CHS
    :param int n: Number of points discretising the inner and outer circles
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a CHS discretised with 64 points, with a
    diameter of 48 and thickness of 3.2, and generates a mesh with a maximum
    triangular area of 1.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.Chs(d=48, t=3.2, n=64)
        mesh = geometry.create_mesh(mesh_sizes=[1.0])

    ..  figure:: ../images/sections/chs_geometry.png
        :align: center
        :scale: 75 %

        CHS geometry.

    ..  figure:: ../images/sections/chs_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, d, t, n, shift=[0, 0]):
        """Inits the Chs class."""

        # assign control point
        control_points = [[d * 0.5 - t * 0.5, 0]]

        super().__init__(control_points, shift)

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

        self.shift_section()


class Rhs(Geometry):
    """Constructs a rectangular hollow section centered at *(b/2, d/2)*, with
    depth *d*, width *b*, thickness *t* and outer radius *r_out*, using *n_r*
    points to construct the inner and outer radii.

    :param float d: Depth of the RHS
    :param float b: Width of the RHS
    :param float t: Thickness of the RHS
    :param float r_out: Outer radius of the RHS
    :param int n_r: Number of points discretising the inner and outer radii
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates an RHS with a depth of 100, a width of 50, a
    thickness of 6 and an outer radius of 9, using 8 points to discretise the
    inner and outer radii. A mesh is generated with a maximum triangular area
    of 2.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.Rhs(d=100, b=50, t=6, r_out=9, n_r=8)
        mesh = geometry.create_mesh(mesh_sizes=[2.0])

    ..  figure:: ../images/sections/rhs_geometry.png
        :align: center
        :scale: 75 %

        RHS geometry.

    ..  figure:: ../images/sections/rhs_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, d, b, t, r_out, n_r, shift=[0, 0]):
        """Inits the Rhs class."""

        # assign control point
        control_points = [[b - t * 0.5, d * 0.5]]

        super().__init__(control_points, shift)

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

        self.shift_section()


class ISection(Geometry):
    """Constructs an I-section centered at *(b/2, d/2)*, with depth *d*, width
    *b*, flange thickness *t_f*, web thickness *t_w*, and root radius *r*,
    using *n_r* points to construct the root radius.

    :param float d: Depth of the I-section
    :param float b: Width of the I-section
    :param float t_f: Flange thickness of the I-section
    :param float t_w: Web thickness of the I-section
    :param float r: Root radius of the I-section
    :param int n_r: Number of points discretising the root radius
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates an I-section with a depth of 203, a width of
    133, a flange thickness of 7.8, a web thickness of 5.8 and a root radius of
    8.9, using 16 points to discretise the root radius. A mesh is generated
    with a maximum triangular area of 3.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.ISection(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=16)
        mesh = geometry.create_mesh(mesh_sizes=[3.0])

    ..  figure:: ../images/sections/isection_geometry.png
        :align: center
        :scale: 75 %

        I-section geometry.

    ..  figure:: ../images/sections/isection_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, d, b, t_f, t_w, r, n_r, shift=[0, 0]):
        """Inits the ISection class."""

        # assign control point
        control_points = [[b * 0.5, d * 0.5]]

        super().__init__(control_points, shift)

        # add first three points
        self.points.append([0, 0])
        self.points.append([b, 0])
        self.points.append([b, t_f])

        # construct the bottom right radius
        for i in range(n_r):
            # determine polar angle
            theta = 3.0 / 2 * np.pi * (1 - i * 1.0 / max(1, n_r - 1) * 1.0 / 3)

            # calculate the locations of the radius points
            x = b * 0.5 + t_w * 0.5 + r + r * np.cos(theta)
            y = t_f + r + r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # construct the top right radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi * (1 - i * 1.0 / max(1, n_r - 1) * 0.5)

            # calculate the locations of the radius points
            x = b * 0.5 + t_w * 0.5 + r + r * np.cos(theta)
            y = d - t_f - r + r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # add the next four points
        self.points.append([b, d - t_f])
        self.points.append([b, d])
        self.points.append([0, d])
        self.points.append([0, d - t_f])

        # construct the top left radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi * 0.5 * (1 - i * 1.0 / max(1, n_r - 1))

            # calculate the locations of the radius points
            x = b * 0.5 - t_w * 0.5 - r + r * np.cos(theta)
            y = d - t_f - r + r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # construct the bottom left radius
        for i in range(n_r):
            # determine polar angle
            theta = -np.pi * i * 1.0 / max(1, n_r - 1) * 0.5

            # calculate the locations of the radius points
            x = b * 0.5 - t_w * 0.5 - r + r * np.cos(theta)
            y = t_f + r + r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # add the last point
        self.points.append([0, t_f])

        # build the facet list
        for i in range(len(self.points)):
            # if we are not at the last point
            if i != len(self.points) - 1:
                self.facets.append((i, i + 1))
            # if we are at the last point, complete the loop
            else:
                self.facets.append((len(self.points) - 1, 0))

        self.shift_section()


class PfcSection(Geometry):
    """Constructs a PFC section with the bottom left corner at the origin
    *(0, 0)*, with depth *d*, width *b*, flange thickness *t_f*, web  thickness
    *t_w* and root radius *r*, using *n_r* points to construct the root radius.

    :param float d: Depth of the PFC section
    :param float b: Width of the PFC section
    :param float t_f: Flange thickness of the PFC section
    :param float t_w: Web thickness of the PFC section
    :param float r: Root radius of the PFC section
    :param int n_r: Number of points discretising the root radius
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a PFC section with a depth of 250, a width of
    90, a flange thickness of 15, a web thickness of 8 and a root radius of
    12, using 8 points to discretise the root radius. A mesh is generated
    with a maximum triangular area of 5.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.PfcSection(d=250, b=90, t_f=15, t_w=8, r=12, n_r=8)
        mesh = geometry.create_mesh(mesh_sizes=[5.0])

    ..  figure:: ../images/sections/pfc_geometry.png
        :align: center
        :scale: 75 %

        PFC geometry.

    ..  figure:: ../images/sections/pfc_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, d, b, t_f, t_w, r, n_r, shift=[0, 0]):
        """Inits the PfcSection class."""

        # assign control point
        control_points = [[t_w * 0.5, d * 0.5]]

        super().__init__(control_points, shift)

        # add first three points
        self.points.append([0, 0])
        self.points.append([b, 0])
        self.points.append([b, t_f])

        # construct the bottom right radius
        for i in range(n_r):
            # determine polar angle
            theta = 3.0 / 2 * np.pi * (1 - i * 1.0 / max(1, n_r - 1) * 1.0 / 3)

            # calculate the locations of the radius points
            x = t_w + r + r * np.cos(theta)
            y = t_f + r + r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # construct the top right radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi * (1 - i * 1.0 / max(1, n_r - 1) * 0.5)

            # calculate the locations of the radius points
            x = t_w + r + r * np.cos(theta)
            y = d - t_f - r + r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # add last three points
        self.points.append([b, d - t_f])
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

        self.shift_section()


class TeeSection(Geometry):
    """Constructs a Tee section with the top left corner at *(0, d)*, with
    depth *d*, width *b*, flange thickness *t_f*, web thickness *t_w* and root
    radius *r*, using *n_r* points to construct the root radius.

    :param float d: Depth of the Tee section
    :param float b: Width of the Tee section
    :param float t_f: Flange thickness of the Tee section
    :param float t_w: Web thickness of the Tee section
    :param float r: Root radius of the Tee section
    :param int n_r: Number of points discretising the root radius
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a Tee section with a depth of 200, a width of
    100, a flange thickness of 12, a web thickness of 6 and a root radius of
    8, using 8 points to discretise the root radius. A mesh is generated
    with a maximum triangular area of 3.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.TeeSection(d=200, b=100, t_f=12, t_w=6, r=8, n_r=8)
        mesh = geometry.create_mesh(mesh_sizes=[3.0])

    ..  figure:: ../images/sections/tee_geometry.png
        :align: center
        :scale: 75 %

        Tee section geometry.

    ..  figure:: ../images/sections/tee_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, d, b, t_f, t_w, r, n_r, shift=[0, 0]):
        """Inits the TeeSection class."""

        # assign control point
        control_points = [[b * 0.5, d - t_f * 0.5]]

        super().__init__(control_points, shift)

        # add first two points
        self.points.append([b * 0.5 - t_w * 0.5, 0])
        self.points.append([b * 0.5 + t_w * 0.5, 0])

        # construct the top right radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi * (1 - i * 1.0 / max(1, n_r - 1) * 0.5)

            # calculate the locations of the radius points
            x = b * 0.5 + t_w * 0.5 + r + r * np.cos(theta)
            y = d - t_f - r + r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # add next four points
        self.points.append([b, d - t_f])
        self.points.append([b, d])
        self.points.append([0, d])
        self.points.append([0, d - t_f])

        # construct the top left radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi * 0.5 * (1 - i * 1.0 / max(1, n_r - 1))

            # calculate the locations of the radius points
            x = b * 0.5 - t_w * 0.5 - r + r * np.cos(theta)
            y = d - t_f - r + r * np.sin(theta)

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

        self.shift_section()


class AngleSection(Geometry):
    """Constructs an angle section with the bottom left corner at the origin
    *(0, 0)*, with depth *d*, width *b*, thickness *t*, root radius *r_r* and
    toe radius *r_t*, using *n_r* points to construct the radii.

    :param float d: Depth of the angle section
    :param float b: Width of the angle section
    :param float t: Thickness of the angle section
    :param float r_r: Root radius of the angle section
    :param float r_t: Toe radius of the angle section
    :param int n_r: Number of points discretising the radii
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates an angle section with a depth of 150, a width
    of 100, a thickness of 8, a root radius of 12 and a toe radius of 5, using
    16 points to discretise the radii. A mesh is generated with a maximum
    triangular area of 2.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.AngleSection(d=150, b=100, t=8, r_r=12, r_t=5, n_r=16)
        mesh = geometry.create_mesh(mesh_sizes=[2.0])

    ..  figure:: ../images/sections/angle_geometry.png
        :align: center
        :scale: 75 %

        Angle section geometry.

    ..  figure:: ../images/sections/angle_mesh.png
        :align: center
        :scale: 75 %
    """

    def __init__(self, d, b, t, r_r, r_t, n_r, shift=[0, 0]):
        """Inits the AngleSection class."""

        # assign control point
        control_points = [[t * 0.5, t * 0.5]]

        super().__init__(control_points, shift)

        # add first two points
        self.points.append([0, 0])
        self.points.append([b, 0])

        # construct the bottom toe radius
        for i in range(n_r):
            # determine polar angle
            theta = i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate the locations of the radius points
            x = b - r_t + r_t * np.cos(theta)
            y = t - r_t + r_t * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # construct the root radius
        for i in range(n_r):
            # determine polar angle
            theta = 3.0 / 2 * np.pi * (1 - i * 1.0 / max(1, n_r - 1) * 1.0 / 3)

            # calculate the locations of the radius points
            x = t + r_r + r_r * np.cos(theta)
            y = t + r_r + r_r * np.sin(theta)

            # append the current points to the points list
            self.points.append([x, y])

        # construct the top toe radius
        for i in range(n_r):
            # determine polar angle
            theta = i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate the locations of the radius points
            x = t - r_t + r_t * np.cos(theta)
            y = d - r_t + r_t * np.sin(theta)

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

        self.shift_section()


class CeeSection(Geometry):
    """Constructs a Cee section with the bottom left corner at the origin
    *(0, 0)*, with depth *d*, width *b*, lip *l*, thickness *t* and outer
    radius *r_out*, using *n_r* points to construct the radius.

    :param float d: Depth of the Cee section
    :param float b: Width of the Cee section
    :param float l: Lip of the Cee section
    :param float t: Thickness of the Cee section
    :param float r_out: Outer radius of the Cee section
    :param int n_r: Number of points discretising the outer radius
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a Cee section with a depth of 125, a width
    of 50, a lip of 30, a thickness of 1.5 and an outer radius of 6, using 8
    points to discretise the radius. A mesh is generated with a maximum
    triangular area of 0.25::

        import sectionproperties.pre.sections as sections

        geometry = sections.CeeSection(d=125, b=50, l=30, t=1.5, r_out=6, n_r=8)
        mesh = geometry.create_mesh(mesh_sizes=[0.25])

    ..  figure:: ../images/sections/cee_geometry.png
        :align: center
        :scale: 75 %

        Cee section geometry.

    ..  figure:: ../images/sections/cee_mesh.png
        :align: center
        :scale: 75 %
    """

    def __init__(self, d, b, l, t, r_out, n_r, shift=[0, 0]):
        """Inits the CeeSection class."""

        # assign control point
        control_points = [[t * 0.5, d * 0.5]]

        super().__init__(control_points, shift)

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

        self.shift_section()


class ZedSection(Geometry):
    """Constructs a Zed section with the bottom left corner at the origin
    *(0, 0)*, with depth *d*, left flange width *b_l*, right flange width
    *b_r*, lip *l*, thickness *t* and outer radius *r_out*, using *n_r* points
    to construct the radius.

    :param float d: Depth of the Zed section
    :param float b_l: Left flange width of the Zed section
    :param float b_r: Right flange width of the Zed section
    :param float l: Lip of the Zed section
    :param float t: Thickness of the Zed section
    :param float r_out: Outer radius of the Zed section
    :param int n_r: Number of points discretising the outer radius
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a Zed section with a depth of 100, a left
    flange width of 40, a right flange width of 50, a lip of 20, a thickness of
    1.2 and an outer radius of 5, using 8 points to discretise the radius.
    A mesh is generated with a maximum triangular area of 0.15::

        import sectionproperties.pre.sections as sections

        geometry = sections.ZedSection(d=100, b_l=40, b_r=50, l=20, t=1.2, r_out=5, n_r=8)
        mesh = geometry.create_mesh(mesh_sizes=[0.15])

    ..  figure:: ../images/sections/zed_geometry.png
        :align: center
        :scale: 75 %

        Zed section geometry.

    ..  figure:: ../images/sections/zed_mesh.png
        :align: center
        :scale: 75 %
    """

    def __init__(self, d, b_l, b_r, l, t, r_out, n_r, shift=[0, 0]):
        """Inits the ZedSection class."""

        # assign control point
        control_points = [[t * 0.5, d * 0.5]]

        super().__init__(control_points, shift)

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
            x_outer = b_r - r_out + r_out * np.cos(theta)
            y_outer = r_out + r_out * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_outer, y_outer])

        # add next two points
        self.points.append([b_r, l])
        self.points.append([b_r - t, l])

        # construct the inner bottom right radius
        for i in range(n_r):
            # determine polar angle
            theta = 2 * np.pi - i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_inner = b_r - r_out + r_in * np.cos(theta)
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
            x_outer = t - b_l + r_out + r_out * np.cos(theta)
            y_outer = d - r_out + r_out * np.sin(theta)

            # append the current points to the points list
            self.points.append([x_outer, y_outer])

        # add the next two points
        self.points.append([t - b_l, d - l])
        self.points.append([t - b_l + t, d - l])

        # construct the inner top left radius
        for i in range(n_r):
            # determine polar angle
            theta = np.pi - i * 1.0 / max(1, n_r - 1) * np.pi * 0.5

            # calculate location of inner and outer points
            x_inner = t - b_l + r_out + r_in * np.cos(theta)
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

        self.shift_section()


class CruciformSection(Geometry):
    """Constructs a cruciform section centered at the origin *(0, 0)*, with
    depth *d*, width *b*, thickness *t* and root radius *r*, using *n_r* points
    to construct the root radius.

    :param float d: Depth of the cruciform section
    :param float b: Width of the cruciform section
    :param float t: Thickness of the cruciform section
    :param float r: Root radius of the cruciform section
    :param int n_r: Number of points discretising the root radius
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a cruciform section with a depth of 250, a
    width of 175, a thickness of 12 and a root radius of 16, using 16 points to
    discretise the radius. A mesh is generated with a maximum triangular area
    of 5.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.CruciformSection(d=250, b=175, t=12, r=16, n_r=16)
        mesh = geometry.create_mesh(mesh_sizes=[5.0])

    ..  figure:: ../images/sections/cruciform_geometry.png
        :align: center
        :scale: 75 %

        Cruciform section geometry.

    ..  figure:: ../images/sections/cruciform_mesh.png
        :align: center
        :scale: 75 %
    """

    def __init__(self, d, b, t, r, n_r, shift=[0, 0]):
        """Inits the CruciformSection class."""

        # assign control point
        control_points = [[0, 0]]

        super().__init__(control_points, shift)

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

        self.shift_section()


class MergedSection(Geometry):
    """Merges a number of section geometries into one geometry. Note that for
    the meshing algorithm to work, there needs to be connectivity between all
    regions of the provided geometries. Overlapping of geometries is permitted.

    :param sections: A list of geometry objects to merge into one
        :class:`~sectionproperties.pre.sections.Geometry` object
    :type sections: list[:class:`~sectionproperties.pre.sections.Geometry`]

    The following example creates a combined cross-section with a 150x100x6 RHS
    placed on its side on top of a 200UB25.4. A mesh is generated with a
    maximum triangle size of 5.0 for the I-section and 2.5 for the RHS::

        import sectionproperties.pre.sections as sections

        isection = sections.ISection(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8)
        box = sections.Rhs(d=100, b=150, t=6, r_out=15, n_r=8, shift=[-8.5, 203])

        geometry = sections.MergedSection([isection, box])
        geometry.clean_geometry()
        mesh = geometry.create_mesh(mesh_sizes=[5.0, 2.5])

    ..  figure:: ../images/sections/merged_geometry.png
        :align: center
        :scale: 75 %

        Merged section geometry.

    ..  figure:: ../images/sections/merged_mesh.png
        :align: center
        :scale: 75 %
    """

    def __init__(self, sections):
        """Inits the MergedSection class."""

        super().__init__([], [0, 0])

        point_count = 0

        # loop through all sections
        for section in sections:
            # add facets
            for facet in section.facets:
                self.facets.append([facet[0] + point_count,
                                    facet[1] + point_count])

            # add points and count points
            for point in section.points:
                self.points.append([point[0], point[1]])
                point_count += 1

            # add holes
            for hole in section.holes:
                self.holes.append([hole[0], hole[1]])

            # add control points
            for control_point in section.control_points:
                self.control_points.append([control_point[0],
                                            control_point[1]])


# def sectionParse(sectionTypes, sectionData, settings):
#     """
#     Generates the geometry for the structural cross-section to be analysed,
#     defined by a number of different sectionTypes, containing various
#     sectionData. Note that there must be connectivity between all sections
#     (i.e. there cannot be isolated sections) or the meshing and/or
#     cross-section analysis will not work.
#     """
#
#     # initialise output variables
#     points = []
#     facets = []
#     holes = []
#     controlPoints = []
#
#     # initialise pointCount variable
#     pointCount = 0
#
#     # loop through each section
#     for (i, section) in enumerate(sectionTypes):
#         # generate new section depending on section type
#         if (section == "custom"):
#             # load data from current sectionData
#             try:
#                 pointData = sectionData[i]["points"]
#                 facetData = sectionData[i]["facets"]
#                 holeData = sectionData[i]["holes"]
#                 x = sectionData[i]["x"]
#                 y = sectionData[i]["y"]
#                 controlPointData = sectionData[i]["control-point"]
#             except KeyError as err:
#                 handleKeyError(err, section)
#
#             # generate a new section
#             newSection = generateCustom(pointData, facetData, holeData, x, y,
#                                         controlPointData)
#
#         elif (section == "rectangle"):
#             try:
#                 # load data from current sectionData
#                 d = sectionData[i]["d"]
#                 b = sectionData[i]["b"]
#                 x = sectionData[i]["x"]
#                 y = sectionData[i]["y"]
#             except KeyError as err:
#                 handleKeyError(err, section)
#
#             # if there is a control-point, load it
#             try:
#                 controlPointData = sectionData[i]["control-point"]
#             # if there is no control-point, set it to None
#             except KeyError:
#                 controlPointData = None
#
#             # generate a new section
#             newSection = generateRectangle(d, b, x, y, controlPointData)
#
#         elif (section == "circle"):
#             # load data from current sectionData
#             try:
#                 d = sectionData[i]["d"]
#                 n = sectionData[i]["n"]
#                 x = sectionData[i]["x"]
#                 y = sectionData[i]["y"]
#             except KeyError as err:
#                 handleKeyError(err, section)
#
#             # if there is a control-point, load it
#             try:
#                 controlPointData = sectionData[i]["control-point"]
#             # if there is no control-point, set it to None
#             except KeyError:
#                 controlPointData = None
#
#             # generate a new section
#             newSection = generateCircle(d, n, x, y, controlPointData)
#
#         elif (section == "chs"):
#             # load data from current sectionData
#             try:
#                 d = sectionData[i]["d"]
#                 t = sectionData[i]["t"]
#                 n = sectionData[i]["n"]
#                 x = sectionData[i]["x"]
#                 y = sectionData[i]["y"]
#             except KeyError as err:
#                 handleKeyError(err, section)
#
#             # if there is a control-point, load it
#             try:
#                 controlPointData = sectionData[i]["control-point"]
#             # if there is no control-point, set it to None
#             except KeyError:
#                 controlPointData = None
#
#             # generate a new section
#             newSection = generateCHS(d, t, n, x, y, controlPointData)
#
#         elif (section == "rhs"):
#             # load data from current sectionData
#             try:
#                 d = sectionData[i]["d"]
#                 b = sectionData[i]["b"]
#                 t = sectionData[i]["t"]
#                 r_out = sectionData[i]["r_out"]
#                 n_r = sectionData[i]["n_r"]
#                 x = sectionData[i]["x"]
#                 y = sectionData[i]["y"]
#             except KeyError as err:
#                 handleKeyError(err, section)
#
#             # if there is a control-point, load it
#             try:
#                 controlPointData = sectionData[i]["control-point"]
#             # if there is no control-point, set it to None
#             except KeyError:
#                 controlPointData = None
#
#             # generate a new section
#             newSection = generateRHS(
#                 d, b, t, r_out, n_r, x, y, controlPointData)
#
#         elif (section == "i-section"):
#             # load data from current sectionData
#             try:
#                 d = sectionData[i]["d"]
#                 b = sectionData[i]["b"]
#                 tf = sectionData[i]["tf"]
#                 tw = sectionData[i]["tw"]
#                 r = sectionData[i]["r"]
#                 n_r = sectionData[i]["n_r"]
#                 x = sectionData[i]["x"]
#                 y = sectionData[i]["y"]
#             except KeyError as err:
#                 handleKeyError(err, section)
#
#             # if there is a control-point, load it
#             try:
#                 controlPointData = sectionData[i]["control-point"]
#             # if there is no control-point, set it to None
#             except KeyError:
#                 controlPointData = None
#
#             # generate a new section
#             newSection = generateISection(
#                 d, b, tf, tw, r, n_r, x, y, controlPointData)
#
#         elif (section == "pfc"):
#             # load data from current sectionData
#             try:
#                 d = sectionData[i]["d"]
#                 b = sectionData[i]["b"]
#                 tf = sectionData[i]["tf"]
#                 tw = sectionData[i]["tw"]
#                 r = sectionData[i]["r"]
#                 n_r = sectionData[i]["n_r"]
#                 x = sectionData[i]["x"]
#                 y = sectionData[i]["y"]
#             except KeyError as err:
#                 handleKeyError(err, section)
#
#             # if there is a control-point, load it
#             try:
#                 controlPointData = sectionData[i]["control-point"]
#             # if there is no control-point, set it to None
#             except KeyError:
#                 controlPointData = None
#
#             # generate a new section
#             newSection = generatePFCSection(
#                 d, b, tf, tw, r, n_r, x, y, controlPointData)
#
#         elif (section == "tee"):
#             # load data from current sectionData
#             try:
#                 d = sectionData[i]["d"]
#                 b = sectionData[i]["b"]
#                 tf = sectionData[i]["tf"]
#                 tw = sectionData[i]["tw"]
#                 r = sectionData[i]["r"]
#                 n_r = sectionData[i]["n_r"]
#                 x = sectionData[i]["x"]
#                 y = sectionData[i]["y"]
#             except KeyError as err:
#                 handleKeyError(err, section)
#
#             # if there is a control-point, load it
#             try:
#                 controlPointData = sectionData[i]["control-point"]
#             # if there is no control-point, set it to None
#             except KeyError:
#                 controlPointData = None
#
#             # generate a new section
#             newSection = generateTeeSection(
#                 d, b, tf, tw, r, n_r, x, y, controlPointData)
#
#         elif (section == "angle"):
#             # load data from current sectionData
#             try:
#                 d = sectionData[i]["d"]
#                 b = sectionData[i]["b"]
#                 t = sectionData[i]["t"]
#                 r_root = sectionData[i]["r_root"]
#                 r_toe = sectionData[i]["r_toe"]
#                 n_r = sectionData[i]["n_r"]
#                 x = sectionData[i]["x"]
#                 y = sectionData[i]["y"]
#             except KeyError as err:
#                 handleKeyError(err, section)
#
#             # if there is a control-point, load it
#             try:
#                 controlPointData = sectionData[i]["control-point"]
#             # if there is no control-point, set it to None
#             except KeyError:
#                 controlPointData = None
#
#             # generate a new section
#             newSection = generateAngleSection(
#                 d, b, t, r_root, r_toe, n_r, x, y, controlPointData)
#
#         elif (section == "cee"):
#             # load data from current sectionData
#             try:
#                 d = sectionData[i]["d"]
#                 b = sectionData[i]["b"]
#                 lip = sectionData[i]["l"]
#                 t = sectionData[i]["t"]
#                 r_out = sectionData[i]["r_out"]
#                 n_r = sectionData[i]["n_r"]
#                 x = sectionData[i]["x"]
#                 y = sectionData[i]["y"]
#             except KeyError as err:
#                 handleKeyError(err, section)
#
#             # if there is a control-point, load it
#             try:
#                 controlPointData = sectionData[i]["control-point"]
#             # if there is no control-point, set it to None
#             except KeyError:
#                 controlPointData = None
#
#             # generate a new section
#             newSection = generateCeeSection(
#                 d, b, lip, t, r_out, n_r, x, y, controlPointData)
#
#         elif (section == "zed"):
#             # load data from current sectionData
#             try:
#                 d = sectionData[i]["d"]
#                 b1 = sectionData[i]["b1"]
#                 b2 = sectionData[i]["b2"]
#                 lip = sectionData[i]["l"]
#                 t = sectionData[i]["t"]
#                 r_out = sectionData[i]["r_out"]
#                 n_r = sectionData[i]["n_r"]
#                 x = sectionData[i]["x"]
#                 y = sectionData[i]["y"]
#             except KeyError as err:
#                 handleKeyError(err, section)
#
#             # if there is a control-point, load it
#             try:
#                 controlPointData = sectionData[i]["control-point"]
#             # if there is no control-point, set it to None
#             except KeyError:
#                 controlPointData = None
#
#             # generate a new section
#             newSection = generateZedSection(
#                 d, b1, b2, lip, t, r_out, n_r, x, y, controlPointData)
#
#         elif (section == "cruciform"):
#             # load data from current sectionData
#             try:
#                 d = sectionData[i]["d"]
#                 b = sectionData[i]["b"]
#                 t = sectionData[i]["t"]
#                 r = sectionData[i]["r"]
#                 n_r = sectionData[i]["n_r"]
#                 x = sectionData[i]["x"]
#                 y = sectionData[i]["y"]
#             except KeyError as err:
#                 handleKeyError(err, section)
#
#             # if there is a control-point, load it
#             try:
#                 controlPointData = sectionData[i]["control-point"]
#             # if there is no control-point, set it to None
#             except KeyError:
#                 controlPointData = None
#
#             # generate a new section
#             newSection = generateCruciform(
#                 d, b, t, r, n_r, x, y, controlPointData)
#
#         else:
#             print("Error: section type '{}' is not defined.".format(section))
#             quit()
#
#         # get points, facets, holes and controlpoint from newSection
#         (newPoints, newFacets, newHoles,
#          newControlPoint) = newSection.returnSection()
#
#         # loop through the facets in the newSection and append to the list
#         for f in newFacets:
#             facets.append([f[0] + pointCount, f[1] + pointCount])
#
#         # loop through the points in the newSection and append to the list
#         for p in newPoints:
#             pointCount += 1
#             points.append([p[0], p[1]])
#
#         # loop through the holes in the newSection and append to the list
#         for h in newHoles:
#             holes.append([h[0], h[1]])
#
#         # append the controlPoint from the newSection
#         controlPoints.append([newControlPoint[0], newControlPoint[1]])
#
#     if (settings.outputLog):
#         print("-- Loaded {0} points, {1} facets and {2} holes ".format(
#             len(points), len(facets), len(holes)) +
#             "from {0} sections.".format(len(sectionTypes)))
#
#     return (points, facets, holes, controlPoints)
#
#
# def handleKeyError(err, section):
#     """
#     Displays an error message if the correct keys are not provided for the
#     current section and quits the program.
#     """
#
#     print(
#         "Error: Required key {0} not found for section type '{1}'.".format(
#             err, section) +
#         " Refer to the documentation for the required keys.")
#     quit()
