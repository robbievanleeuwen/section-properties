"""Basic sections."""

import numpy as np

import sectionproperties.post.post as post
import sectionproperties.pre.pre as pre


class Geometry(pre.GeometryCleanerMixin):
    """Parent class for a cross-section geometry input.

    Provides an interface for the user to specify the geometry defining a cross-section. A method
    is provided for generating a triangular mesh, for translating the cross-section by *(x, y)* and
    for plotting the geometry.

    Attributes
    ----------
    control_points : list[list[float, float]]
        A list of points *(x, y)* that define different regions of the cross-section. A control
        point is an arbitrary point within a region enclosed by facets.
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*
    points : list[list[float, float]]
        List of points *(x, y)* defining the vertices of the cross-section
    facets : list[list[int, int]]
        List of point index pairs *(p1, p2)* defining the edges of the cross-section
    holes : list[list[float, float]]
        List of points *(x, y)* defining the locations of holes within the cross-section. If there
        are no holes, provide an empty list [].
    perimeter : list[int]
        List of facet indices defining the perimeter of the cross-section
    """

    def __init__(self, control_points, shift):
        """Initialises the Geometry class."""
        self.control_points = control_points
        self.shift = shift
        self.points = []
        self.facets = []
        self.holes = []
        self.perimeter = []

    def create_mesh(self, mesh_sizes):
        """Creates a quadratic triangular mesh from the Geometry object.

        Parameters
        ----------
        mesh_sizes : list[float]
            A list of maximum element areas corresponding to each region within the cross-section
            geometry.

        Returns
        -------
        :class:`meshpy.triangle.MeshInfo`
            Object containing generated mesh data

        Raises
        ------
        AssertionError
            If the number of mesh sizes does not match the number of regions

        Examples
        --------
        The following example creates a circular cross-section with a diameter of 50 with 64
        points, and generates a mesh with a maximum triangular area of 2.5

        .. plot::
           :context: reset

           import sectionproperties.pre.sections as sections
           from sectionproperties.analysis.cross_section import CrossSection
           geometry = sections.CircularSection(d=50, n=64)
           geometry.plot_geometry()
           mesh = geometry.create_mesh(mesh_sizes=[2.5])
           section = CrossSection(geometry, mesh)
           section.plot_mesh()
        """
        msg = "Number of mesh_sizes ({0}), should match the number of regions ({1})".format(
            len(mesh_sizes), len(self.control_points)
        )
        assert len(mesh_sizes) == len(self.control_points), msg

        return pre.create_mesh(
            self.points, self.facets, self.holes, self.control_points, mesh_sizes
        )

    def shift_section(self):
        """Shifts the cross-section parameters by the class variable vector *shift*."""
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
        """Rotates the geometry and specified angle about a point. If the rotation point is not
        provided, rotates the section about the first control point in the list of control points
        of the :class:`~sectionproperties.pre.sections.Geometry` object.

        Parameters
        ----------
        angle : float
            Angle (degrees) by which to rotate the section. A positive angle leads to a counter-
            clockwise rotation.
        rot_point : list[float, float]
            Point *(x, y)* about which to rotate the section

        Examples
        --------
        The following example rotates a 200UB25 section clockwise by 30 degrees

        .. plot::
           :context: reset

           import sectionproperties.pre.sections as sections
           geometry = sections.ISection(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8)
           geometry.rotate_section(angle=-30)
           geometry.plot_geometry()
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
        """Mirrors the geometry about a point on either the x or y-axis. If no point is provided,
        mirrors the geometry about the first control point in the list of control points of the
        :class:`~sectionproperties.pre.sections.Geometry` object.

        Parameters
        ----------
        axis : str
            Axis about which to mirror the geometry, *'x'* or *'y'*
        mirror_point : list[float, float]
            Point about which to mirror the geometry *(x, y)*

        Examples
        --------
        The following example mirrors a 200PFC section about the y-axis and the point (0, 0)

        .. plot::
           :context: reset

           import sectionproperties.pre.sections as sections
           geometry = sections.PfcSection(d=200, b=75, t_f=12, t_w=6, r=12, n_r=8)
           geometry.mirror_section(axis='y', mirror_point=[0, 0])
           geometry.plot_geometry()
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

        Parameters
        ----------
        point : list[float, float]
            Location of the point

        Returns
        -------
        int
            Point id
        """
        self.points.append(point)
        return len(self.points) - 1

    def add_facet(self, facet):
        """Adds a facet to the geometry and returns the added facet id.

        Parameters
        ----------
        facet : list[float, float]
            Point indices of the facet

        Returns
        -------
        int
            Facet id
        """
        self.facets.append(facet)
        return len(self.facets) - 1

    def add_hole(self, hole):
        """Adds a hole location to the geometry and returns the added hole id.

        Parameters
        ----------
        hole : list[float, float]
            Location of the hole

        Returns
        -------
        int
            Hole id
        """
        self.holes.append(hole)
        return len(self.holes) - 1

    def add_control_point(self, control_point):
        """Adds a control point to the geometry and returns the added control
        point id.

        Parameters
        ----------
        hole : list[float, float]
            Location of the control point

        Returns
        -------
        int
            Control point id
        """
        self.control_points.append(control_point)
        return len(self.control_points) - 1

    def clean_geometry(self, verbose=False):
        """Performs a full clean on the geometry.

        Parameters
        ----------
        verbose : bool
            If set to true, information related to the geometry cleaning process is printed to the
            terminal.

        Notes
        -----
        Cleaning the geometry is always recommended when creating a merged section, which may result
        in overlapping or intersecting facets, or duplicate nodes.
        """
        self.zip_points(verbose=verbose)
        self.remove_zero_length_facets(verbose)
        self.remove_duplicate_facets(verbose)
        self.remove_overlapping_facets(verbose)
        self.remove_unused_points(verbose)
        self.intersect_facets(verbose)

    def plot_geometry(
        self, labels=False, perimeter=False, title='Cross-Section Geometry', **kwargs
    ):
        r"""Plots the geometry defined by the input section. If no axes object is supplied a new
        figure and axis is created.

        Parameters
        ----------
        labels : bool
            If set to true, node and facet labels are displayed
        perimeter : bool
            If set to true, boldens the perimeter of the cross-section
        title : str
            Plot title
        \**kwargs
            Passed to :func:`~sectionproperties.post.post.plotting_context`

        Returns
        -------
        (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes.Axes`)
            Matplotlib figure and axes objects (fig, ax)

        Examples
        --------
        The following example creates a CHS discretised with 64 points, with a diameter of 48 and
        thickness of 3.2, and plots the geometry

        .. plot::
           :context: reset

           import sectionproperties.pre.sections as sections
           geometry = sections.Chs(d=48, t=3.2, n=64)
           geometry.plot_geometry()
        """
        with post.plotting_context(title=title, **kwargs) as (fig, ax):
            for (i, f) in enumerate(self.facets):
                if perimeter:
                    if i in self.perimeter:
                        linewidth = 3
                    else:
                        linewidth = 1.5
                else:
                    linewidth = 1.5

                # plot the points and facets
                if i == 0:
                    ax.plot(
                        [self.points[f[0]][0], self.points[f[1]][0]],
                        [self.points[f[0]][1], self.points[f[1]][1]],
                        'ko-',
                        markersize=2,
                        linewidth=linewidth,
                        label='Points & Facets',
                    )
                else:
                    ax.plot(
                        [self.points[f[0]][0], self.points[f[1]][0]],
                        [self.points[f[0]][1], self.points[f[1]][1]],
                        'ko-',
                        markersize=2,
                        linewidth=linewidth,
                    )

            for (i, h) in enumerate(self.holes):
                # plot the holes
                if i == 0:
                    ax.plot(h[0], h[1], 'rx', markersize=5, label='Holes')
                else:
                    ax.plot(h[0], h[1], 'rx', markersize=5)

            for (i, cp) in enumerate(self.control_points):
                # plot the control points
                if i == 0:
                    ax.plot(cp[0], cp[1], 'bo', markersize=5, label='Control Points')
                else:
                    ax.plot(cp[0], cp[1], 'bo', markersize=5)

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

        return fig, ax

    def calculate_extents(self):
        """Calculates the minimum and maximum x and y-values amongst the list of points.

        Returns
        -------
        tuple(float, float, float, float)
            Minimum and maximum x and y-values *(x_min, x_max, y_min, y_max)*
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

            # update the mins and maxes where necessary
            x_min = min(x_min, x)
            x_max = max(x_max, x)
            y_min = min(y_min, y)
            y_max = max(y_max, y)

        return (x_min, x_max, y_min, y_max)

    def draw_radius(self, pt, r, theta, n, anti=True):
        """Adds a quarter radius of points to the points list - centered at point *pt*, with radius
        *r*, starting at angle *theta*, with *n* points. If r = 0, adds pt only.

        Parameters
        ----------
        pt : list[float, float]
            Centre of radius *(x,y)*
        r : float
            Radius
        theta : float
            Initial angle
        n : int
            Number of points
        anti : bool
            Anticlockwise rotation?
        """
        if r == 0:
            self.points.append(pt)
            return

        if anti:
            mult = 1
        else:
            mult = -1

        # calculate radius of points
        for i in range(n):
            # determine angle
            t = theta + mult * i * 1.0 / max(1, n - 1) * np.pi * 0.5

            x = pt[0] + r * np.cos(t)
            y = pt[1] + r * np.sin(t)
            self.points.append([x, y])

    def calculate_facet_length(self, facet):
        """Calculates the length of the facet.

        Parameters
        ----------
        facet : list[int, int]
            Point index pair *(p1, p2)* defining a facet

        Returns
        -------
        float
            Facet length
        """
        # get facet points
        p1 = self.points[facet[0]]
        p2 = self.points[facet[1]]

        # calculate distance between two points
        return np.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2)

    def calculate_perimeter(self):
        """Calculates the perimeter of the cross-section by summing the length of all facets in the
        ``perimeter`` class variable.

        Returns
        -------
        float
            Cross-section perimeter, returns 0 if there is no perimeter defined
        """
        # check to see if there are any facets in the perimeter variable
        if len(self.perimeter) == 0:
            return 0

        # initialise perimeter variable
        perimeter = 0

        # loop through all the facets along the perimeter
        for facet_idx in self.perimeter:
            perimeter += self.calculate_facet_length(self.facets[facet_idx])

        return perimeter


class CustomSection(Geometry):
    """Constructs a cross-section from a list of points, facets, holes and a user specified control
    point.

    Parameters
    ----------
    points : list[list[float, float]]
        List of points *(x, y)* defining the vertices of the cross-section
    facets : list[list[int, int]]
        List of point index pairs *(p1, p2)* defining the edges of the cross-section
    holes : list[list[float, float]]
        List of points *(x, y)* defining the locations of holes within the cross-section. If there
        are no holes, provide an empty list [].
    control_points : list[list[float, float]]
        A list of points *(x, y)* that define different regions of the cross-section. A control
        point is an arbitrary point within a region enclosed by facets.
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*
    perimeter : list[int]
        List of facet indices defining the perimeter of the cross-section

    Examples
    --------
    The following example creates a hollow trapezium with a base width of 100, top width of 50,
    height of 50 and a wall thickness of 10. A mesh is generated with a maximum triangular area of
    2.0

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       points = [[0, 0], [100, 0], [75, 50], [25, 50], [15, 10], [85, 10], [70, 40], [30, 40]]
       facets = [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4]]
       holes = [[50, 25]]
       control_points = [[5, 5]]
       perimeter = [0, 1, 2, 3]
       geometry = sections.CustomSection(points, facets, holes, control_points, perimeter=perimeter)
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[2.0])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, points, facets, holes, control_points, shift=(0, 0), perimeter=()):
        """Initialises the CustomSection class."""
        super().__init__(control_points, shift)

        self.points = points
        self.facets = facets
        self.holes = holes
        self.perimeter = list(perimeter)

        self.shift_section()


class RectangularSection(Geometry):
    """Constructs a rectangular section with the bottom left corner at the origin *(0, 0)*, with
    depth *d* and width *b*.

    Parameters
    ----------
    d : float
        Depth (y) of the rectangle
    b : float
        Width (x) of the rectangle
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*

    Examples
    --------
    The following example creates a rectangular cross-section with a depth of 100 and width of 50,
    and generates a mesh with a maximum triangular area of 5

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       geometry = sections.RectangularSection(d=100, b=50)
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[5])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, d, b, shift=(0, 0)):
        """Initialises the RectangularSection class."""
        # assign control point
        control_points = [[0.5 * b, 0.5 * d]]

        super().__init__(control_points, shift)

        # construct the points and facets
        self.points = [[0, 0], [b, 0], [b, d], [0, d]]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 0]]
        self.perimeter = list(range(len(self.facets)))

        self.shift_section()


class CircularSection(Geometry):
    """Constructs a solid circle centered at the origin *(0, 0)* with diameter *d* and using *n*
    points to construct the circle.

    Parameters
    ----------
    d : float
        Diameter of the circle
    n : int
        Number of points discretising the circle
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*

    Examples
    --------
    The following example creates a circular cross-section with a diameter of 50 with 64 points,
    and generates a mesh with a maximum triangular area of 2.5

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       geometry = sections.CircularSection(d=50, n=64)
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[2.5])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, d, n, shift=(0, 0)):
        """Initialises the CircularSection class."""
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

        self.perimeter = list(range(len(self.facets)))

        self.shift_section()


class Chs(Geometry):
    """Constructs a circular hollow section centered at the origin *(0, 0)*, with diameter *d* and
    thickness *t*, using *n* points to construct the inner and outer circles.

    Parameters
    ----------
    d : float
        Outer diameter of the CHS
    t : float
        Thickness of the CHS
    n : int
        Number of points discretising the inner and outer circles
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*

    Examples
    --------
    The following example creates a CHS discretised with 64 points, with a diameter of 48 and
    thickness of 3.2, and generates a mesh with a maximum triangular area of 1.0

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       geometry = sections.Chs(d=48, t=3.2, n=64)
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[1.0])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, d, t, n, shift=(0, 0)):
        """Initialises the Chs class."""
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

        self.perimeter = list(range(0, len(self.facets), 2))

        self.shift_section()


class EllipticalSection(Geometry):
    """Constructs a solid ellipse centered at the origin *(0, 0)* with vertical diameter *d_y* and
    horizontal diameter *d_x*, using *n* points to construct the ellipse.

    Parameters
    ----------
    d_y : float
        Diameter of the ellipse in the y-dimension
    d_x : float
        Diameter of the ellipse in the x-dimension
    n : int
        Number of points discretising the ellipse
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*

    Examples
    --------
    The following example creates an elliptical cross-section with a vertical diameter of 25 and
    horizontal diameter of 50, with 40 points, and generates a mesh with a maximum triangular area
    of 1.0

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       geometry = sections.EllipticalSection(d_y=25, d_x=50, n=40)
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[1.0])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, d_y, d_x, n, shift=(0, 0)):
        """Initialises the EllipticalSection class."""
        # assign control point centered at zero
        control_points = [[0, 0]]

        super().__init__(control_points, shift)

        # loop through each point on the ellipse
        for i in range(n):
            # determine polar angle
            theta = i * 2 * np.pi * 1.0 / n

            # calculate location of the point
            x = 0.5 * d_x * np.cos(theta)
            y = 0.5 * d_y * np.sin(theta)

            # append the current point to the points list
            self.points.append([x, y])

            # if we are not at the last point
            if i != n - 1:
                self.facets.append([i, i + 1])
            # if we are at the last point, complete the ellipse
            else:
                self.facets.append([i, 0])

        self.perimeter = list(range(len(self.facets)))

        self.shift_section()


class Ehs(Geometry):
    """Constructs an elliptical hollow section centered at the origin *(0, 0)*, with outer vertical
    diameter *d_y*, outer horizontal diameter *d_x*, and thickness *t*, using *n* points to
    construct the inner and outer ellipses.

    Parameters
    ----------
    d_y : float
        Diameter of the ellipse in the y-dimension
    d_x : float
        Diameter of the ellipse in the x-dimension
    t : float
        Thickness of the EHS
    n : int
        Number of points discretising the inner and outer ellipses
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*

    Examples
    --------
    The following example creates a EHS discretised with 30 points, with a outer vertical diameter
    of 25, outer horizontal diameter of 50, and thickness of 2.0, and generates a mesh with a
    maximum triangular area of 0.5

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       geometry = sections.Ehs(d_y=25, d_x=50, t=2.0, n=64)
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[0.5])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, d_y, d_x, t, n, shift=(0, 0)):
        """Initialises the Ehs class."""
        # assign control point
        control_points = [[(d_x * 0.5) - (t * 0.5), 0]]

        super().__init__(control_points, shift)

        # specify a hole in the centre of the EHS
        self.holes = [[0, 0]]

        # loop through each point of the EHS
        for i in range(n):
            # determine polar angle
            theta = i * 2 * np.pi * 1.0 / n

            # calculate location of outer and inner points
            x_outer = 0.5 * d_x * np.cos(theta)
            y_outer = 0.5 * d_y * np.sin(theta)
            x_inner = ((0.5 * d_x) - t) * np.cos(theta)
            y_inner = ((0.5 * d_y) - t) * np.sin(theta)

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

        self.perimeter = list(range(0, len(self.facets), 2))

        self.shift_section()


class Rhs(Geometry):
    """Constructs a rectangular hollow section centered at *(b/2, d/2)*, with depth *d*, width *b*,
    thickness *t* and outer radius *r_out*, using *n_r* points to construct the inner and outer
    radii. If the outer radius is less than the thickness of the RHS, the inner radius is set to
    zero.

    Parameters
    ----------
    d : float
        Depth of the RHS
    b : float
        Width of the RHS
    t : float
        Thickness of the RHS
    r_out : float
        Outer radius of the RHS
    n_r : int
        Number of points discretising the inner and outer radii
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*

    Examples
    --------
    The following example creates an RHS with a depth of 100, a width of 50, a thickness of 6 and
    an outer radius of 9, using 8 points to discretise the inner and outer radii. A mesh is
    generated with a maximum triangular area of 2.0

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       geometry = sections.Rhs(d=100, b=50, t=6, r_out=9, n_r=8)
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[2.0])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, d, b, t, r_out, n_r, shift=(0, 0)):
        """Initialises the Rhs class."""
        # assign control point
        control_points = [[b - t * 0.5, d * 0.5]]

        super().__init__(control_points, shift)

        # specify a hole in the centre of the RHS
        self.holes = [[b * 0.5, d * 0.5]]

        # calculate internal radius
        r_in = max(r_out - t, 0)

        # construct the outer radius points
        self.draw_radius([r_out, r_out], r_out, np.pi, n_r)
        self.draw_radius([b - r_out, r_out], r_out, 1.5 * np.pi, n_r)
        self.draw_radius([b - r_out, d - r_out], r_out, 0, n_r)
        self.draw_radius([r_out, d - r_out], r_out, 0.5 * np.pi, n_r)

        # construct the outer radius facet list
        n_outer = len(self.points)
        for i in range(n_outer):
            # if we are not at the last point
            if i != n_outer - 1:
                self.facets.append([i, i + 1])
            # if we are at the last point, complete the loop
            else:
                self.facets.append([i, 0])

        # construct the inner radius points
        self.draw_radius([t + r_in, t + r_in], r_in, np.pi, n_r)
        self.draw_radius([b - t - r_in, t + r_in], r_in, 1.5 * np.pi, n_r)
        self.draw_radius([b - t - r_in, d - t - r_in], r_in, 0, n_r)
        self.draw_radius([t + r_in, d - t - r_in], r_in, 0.5 * np.pi, n_r)

        # construct the inner radius facet list
        n_inner = len(self.points) - n_outer
        for i in range(n_inner):
            # if we are not at the last point
            if i != n_inner - 1:
                self.facets.append([i + n_outer, i + n_outer + 1])
            # if we are at the last point, complete the loop
            else:
                self.facets.append([i + n_outer, n_outer])

        self.perimeter = list(range(int(len(self.facets) / 2)))

        self.shift_section()


class ISection(Geometry):
    """Constructs an I-section centered at *(b/2, d/2)*, with depth *d*, width *b*, flange
    thickness *t_f*, web thickness *t_w*, and root radius *r*, using *n_r* points to construct the
    root radius.

    Parameters
    ----------
    d : float
        Depth of the I-section
    b : float
        Width of the I-section
    t_f : float
        Flange thickness of the I-section
    t_w : float
        Web thickness of the I-section
    r : float
        Root radius of the I-section
    n_r : int
        Number of points discretising the root radius
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*

    Examples
    --------
    The following example creates an I-section with a depth of 203, a width of 133, a flange
    thickness of 7.8, a web thickness of 5.8 and a root radius of 8.9, using 16 points to
    discretise the root radius. A mesh is generated with a maximum triangular area of 3.0

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       geometry = sections.ISection(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=16)
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[3.0])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, d, b, t_f, t_w, r, n_r, shift=(0, 0)):
        """Initialises the ISection class."""
        # assign control point
        control_points = [[b * 0.5, d * 0.5]]

        super().__init__(control_points, shift)

        # add first three points
        self.points.append([0, 0])
        self.points.append([b, 0])
        self.points.append([b, t_f])

        # construct the bottom right radius
        pt = [b * 0.5 + t_w * 0.5 + r, t_f + r]
        self.draw_radius(pt, r, 1.5 * np.pi, n_r, False)

        # construct the top right radius
        pt = [b * 0.5 + t_w * 0.5 + r, d - t_f - r]
        self.draw_radius(pt, r, np.pi, n_r, False)

        # add the next four points
        self.points.append([b, d - t_f])
        self.points.append([b, d])
        self.points.append([0, d])
        self.points.append([0, d - t_f])

        # construct the top left radius
        pt = [b * 0.5 - t_w * 0.5 - r, d - t_f - r]
        self.draw_radius(pt, r, 0.5 * np.pi, n_r, False)

        # construct the bottom left radius
        pt = [b * 0.5 - t_w * 0.5 - r, t_f + r]
        self.draw_radius(pt, r, 0, n_r, False)

        # add the last point
        self.points.append([0, t_f])

        # build the facet list
        for i in range(len(self.points)):
            # if we are not at the last point
            if i != len(self.points) - 1:
                self.facets.append([i, i + 1])
            # if we are at the last point, complete the loop
            else:
                self.facets.append([len(self.points) - 1, 0])

        self.perimeter = list(range(len(self.facets)))

        self.shift_section()


class MonoISection(Geometry):
    """Constructs a monosymmetric I-section centered at *(max(b_t, b_b)/2, d/2)*, with depth *d*,
    top flange width *b_t*, bottom flange width *b_b*, top flange thickness *t_ft*, top flange
    thickness *t_fb*, web thickness *t_w*, and root radius *r*, using *n_r* points to construct the
    root radius.

    Parameters
    ----------
    d : float
        Depth of the I-section
    b_t : float
        Top flange width
    b_b : float
        Bottom flange width
    t_ft : float
        Top flange thickness of the I-section
    t_fb : float
        Bottom flange thickness of the I-section
    t_w : float
        Web thickness of the I-section
    r : float
        Root radius of the I-section
    n_r : int
        Number of points discretising the root radius
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*

    Examples
    --------
    The following example creates a monosymmetric I-section with a depth of 200, a top flange width
    of 50, a top flange thickness of 12, a bottom flange width of 130, a bottom flange thickness of
    8, a web thickness of 6 and a root radius of 8, using 16 points to discretise the root radius.
    A mesh is generated with a maximum triangular area of 3.0

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       geometry = sections.MonoISection(d=200, b_t=50, b_b=130, t_ft=12, t_fb=8, t_w=6, r=8, n_r=16)
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[3.0])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, d, b_t, b_b, t_fb, t_ft, t_w, r, n_r, shift=(0, 0)):
        """Initialises the ISection class."""
        # assign control point
        control_points = [[max(b_t, b_b) * 0.5, d * 0.5]]

        super().__init__(control_points, shift)

        # calculate central axis
        x_central = max(b_t, b_b) * 0.5

        # add first three points
        self.points.append([x_central - b_b * 0.5, 0])
        self.points.append([x_central + b_b * 0.5, 0])
        self.points.append([x_central + b_b * 0.5, t_fb])

        # construct the bottom right radius
        pt = [x_central + t_w * 0.5 + r, t_fb + r]
        self.draw_radius(pt, r, 1.5 * np.pi, n_r, False)

        # construct the top right radius
        pt = [x_central + t_w * 0.5 + r, d - t_ft - r]
        self.draw_radius(pt, r, np.pi, n_r, False)

        # add the next four points
        self.points.append([x_central + b_t * 0.5, d - t_ft])
        self.points.append([x_central + b_t * 0.5, d])
        self.points.append([x_central - b_t * 0.5, d])
        self.points.append([x_central - b_t * 0.5, d - t_ft])

        # construct the top left radius
        pt = [x_central - t_w * 0.5 - r, d - t_ft - r]
        self.draw_radius(pt, r, 0.5 * np.pi, n_r, False)

        # construct the bottom left radius
        pt = [x_central - t_w * 0.5 - r, t_fb + r]
        self.draw_radius(pt, r, 0, n_r, False)

        # add the last point
        self.points.append([x_central - b_b * 0.5, t_fb])

        # build the facet list
        for i in range(len(self.points)):
            # if we are not at the last point
            if i != len(self.points) - 1:
                self.facets.append([i, i + 1])
            # if we are at the last point, complete the loop
            else:
                self.facets.append([len(self.points) - 1, 0])

        self.perimeter = list(range(len(self.facets)))

        self.shift_section()


class TaperedFlangeISection(Geometry):
    """Constructs a Tapered Flange I-section centered at *(b/2, d/2)*, with depth *d*, width *b*,
    mid-flange thickness *t_f*, web thickness *t_w*, root radius *r_r*, flange radius *r_f* and
    flange angle *alpha*, using *n_r* points to construct the radii.

    Parameters
    ----------
    d : float
        Depth of the Tapered Flange I-section
    b : float
        Width of the Tapered Flange I-section
    t_f : float
        Mid-flange thickness of the Tapered Flange I-section (measured at the point equidistant from
        the face of the web to the edge of the flange)
    t_w : float
        Web thickness of the Tapered Flange I-section
    r_r : float
        Root radius of the Tapered Flange I-section
    r_f : float
        Flange radius of the Tapered Flange I-section
    alpha : float
        Flange angle of the Tapered Flange I-section (degrees)
    n_r : int
        Number of points discretising the radii
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*

    Examples
    --------
    The following example creates a Tapered Flange I-section with a depth of 588, a width of 191, a
    mid-flange thickness of 27.2, a web thickness of 15.2, a root radius of 17.8, a flange radius
    of 8.9 and a flange angle of 8°, using 16 points to discretise the radii. A mesh is generated
    with a maximum triangular area of 20.0

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       geometry = sections.TaperedFlangeISection(
           d=588, b=191, t_f=27.2, t_w=15.2, r_r=17.8, r_f=8.9, alpha=8, n_r=16
       )
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[20.0])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, d, b, t_f, t_w, r_r, r_f, alpha, n_r, shift=(0, 0)):
        """Initialises the ISection class."""
        # assign control point
        control_points = [[b * 0.5, d * 0.5]]

        super().__init__(control_points, shift)

        # calculate alpha in radians
        alpha_rad = np.pi * alpha / 180

        # calculate the height of the flange toe and dimensions of the straight
        x1 = b * 0.25 - t_w * 0.25 - r_f * (1 - np.sin(alpha_rad))
        y1 = x1 * np.tan(alpha_rad)
        x2 = b * 0.25 - t_w * 0.25 - r_r * (1 - np.sin(alpha_rad))
        y2 = x2 * np.tan(alpha_rad)
        y_t = t_f - y1 - r_f * np.cos(alpha_rad)

        # add first two points
        self.points.append([0, 0])
        self.points.append([b, 0])

        # construct the bottom right flange toe radius
        if r_f == 0:
            self.points.append([b, y_t])
        else:
            for i in range(n_r):
                # determine polar angle
                theta = i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)

                # calculate the locations of the radius points
                x = b - r_f + r_f * np.cos(theta)
                y = y_t + r_f * np.sin(theta)

                # append the current points to the points list
                self.points.append([x, y])

        # construct the bottom right root radius
        if r_r == 0:
            self.points.append([b * 0.5 + t_w * 0.5, t_f + y2])
        else:
            for i in range(n_r):
                # determine polar angle
                theta = (3.0 / 2 * np.pi - alpha_rad) - (
                    i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)
                )

                # calculate the locations of the radius points
                x = b * 0.5 + t_w * 0.5 + r_r + r_r * np.cos(theta)
                y = t_f + y2 + r_r * np.cos(alpha_rad) + r_r * np.sin(theta)

                # append the current points to the points list
                self.points.append([x, y])

        # construct the top right root radius
        if r_r == 0:
            self.points.append([b * 0.5 + t_w * 0.5, d - t_f - y2])
        else:
            for i in range(n_r):
                # determine polar angle
                theta = np.pi - i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)

                # calculate the locations of the radius points
                x = b * 0.5 + t_w * 0.5 + r_r + r_r * np.cos(theta)
                y = d - t_f - y2 - r_r * np.cos(alpha_rad) + r_r * np.sin(theta)

                # append the current points to the points list
                self.points.append([x, y])

        # construct the top right flange toe radius
        if r_f == 0:
            self.points.append([b, d - y_t])
        else:
            for i in range(n_r):
                # determine polar angle
                theta = (3.0 * np.pi / 2 + alpha_rad) + i * 1.0 / max(1, n_r - 1) * (
                    np.pi * 0.5 - alpha_rad
                )

                # calculate the locations of the radius points
                x = b - r_f + r_f * np.cos(theta)
                y = d - y_t + r_f * np.sin(theta)

                # append the current points to the points list
                self.points.append([x, y])

        # add the next two points
        self.points.append([b, d])
        self.points.append([0, d])

        # construct the top left flange toe radius
        if r_f == 0:
            self.points.append([0, d - y_t])
        else:
            for i in range(n_r):
                # determine polar angle
                theta = np.pi + (i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad))

                # calculate the locations of the radius points
                x = r_f + r_f * np.cos(theta)
                y = d - y_t + r_f * np.sin(theta)

                # append the current points to the points list
                self.points.append([x, y])

        # construct the top left root radius
        if r_r == 0:
            self.points.append([b * 0.5 - t_w * 0.5, d - t_f - y2])
        else:
            for i in range(n_r):
                # determine polar angle
                theta = (np.pi * 0.5 - alpha_rad) - (
                    i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)
                )

                # calculate the locations of the radius points
                x = b * 0.5 - t_w * 0.5 - r_r + r_r * np.cos(theta)
                y = d - t_f - y2 - r_r * np.cos(alpha_rad) + r_r * np.sin(theta)

                # append the current points to the points list
                self.points.append([x, y])

        # construct the bottom left root radius
        if r_r == 0:
            self.points.append([b * 0.5 - t_w * 0.5, t_f + y2])
        else:
            for i in range(n_r):
                # determine polar angle
                theta = -i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)

                # calculate the locations of the radius points
                x = b * 0.5 - t_w * 0.5 - r_r + r_r * np.cos(theta)
                y = t_f + y2 + r_r * np.cos(alpha_rad) + r_r * np.sin(theta)

                # append the current points to the points list
                self.points.append([x, y])

        # construct the bottom left flange toe radius
        if r_f == 0:
            self.points.append([0, y_t])
        else:
            for i in range(n_r):
                # determine polar angle
                theta = (np.pi * 0.5 + alpha_rad) + (
                    i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)
                )

                # calculate the locations of the radius points
                x = r_f + r_f * np.cos(theta)
                y = y_t + r_f * np.sin(theta)

                # append the current points to the points list
                self.points.append([x, y])

        # build the facet list
        for i in range(len(self.points)):
            # if we are not at the last point
            if i != len(self.points) - 1:
                self.facets.append([i, i + 1])
            # if we are at the last point, complete the loop
            else:
                self.facets.append([len(self.points) - 1, 0])

        self.perimeter = list(range(len(self.facets)))

        self.shift_section()


class PfcSection(Geometry):
    """Constructs a PFC section with the bottom left corner at the origin *(0, 0)*, with depth *d*,
    width *b*, flange thickness *t_f*, web  thickness *t_w* and root radius *r*, using *n_r* points
    to construct the root radius.

    Parameters
    ----------
    d : float
        Depth of the PFC section
    b : float
        Width of the PFC section
    t_f : float
        Flange thickness of the PFC section
    t_w : float
        Web thickness of the PFC section
    r : float
        Root radius of the PFC section
    n_r : int
        Number of points discretising the root radius
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*

    Examples
    --------
    The following example creates a PFC section with a depth of 250, a width of 90, a flange
    thickness of 15, a web thickness of 8 and a root radius of 12, using 8 points to discretise the
    root radius. A mesh is generated with a maximum triangular area of 5.0

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       geometry = sections.PfcSection(d=250, b=90, t_f=15, t_w=8, r=12, n_r=8)
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[5.0])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, d, b, t_f, t_w, r, n_r, shift=(0, 0)):
        """Initialises the PfcSection class."""
        # assign control point
        control_points = [[t_w * 0.5, d * 0.5]]

        super().__init__(control_points, shift)

        # add first three points
        self.points.append([0, 0])
        self.points.append([b, 0])
        self.points.append([b, t_f])

        # construct the bottom right radius
        pt = [t_w + r, t_f + r]
        self.draw_radius(pt, r, 1.5 * np.pi, n_r, False)

        # construct the top right radius
        pt = [t_w + r, d - t_f - r]
        self.draw_radius(pt, r, np.pi, n_r, False)

        # add last three points
        self.points.append([b, d - t_f])
        self.points.append([b, d])
        self.points.append([0, d])

        # build the facet list
        for i in range(len(self.points)):
            # if we are not at the last point
            if i != len(self.points) - 1:
                self.facets.append([i, i + 1])
            # if we are at the last point, complete the loop
            else:
                self.facets.append([len(self.points) - 1, 0])

        self.perimeter = list(range(len(self.facets)))

        self.shift_section()


class TaperedFlangeChannel(Geometry):
    """Constructs a Tapered Flange Channel section with the bottom left corner at the origin
    *(0, 0)*, with depth *d*, width *b*, mid-flange thickness *t_f*, web thickness *t_w*, root
    radius *r_r*, flange radius *r_f* and flange angle *alpha*, using *n_r* points to construct the
    radii.

    Parameters
    ----------
    d : float
        Depth of the Tapered Flange Channel section
    b : float
        Width of the Tapered Flange Channel section
    t_f : float
        Mid-flange thickness of the Tapered Flange Channel section (measured at the point
        equidistant from the face of the web to the edge of the flange)
    t_w : float
        Web thickness of the Tapered Flange Channel section
    r_r : float
        Root radius of the Tapered Flange Channel section
    r_f : float
        Flange radius of the Tapered Flange Channel section
    alpha : float
        Flange angle of the Tapered Flange Channel section (degrees)
    n_r : int
        Number of points discretising the radii
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*

    Examples
    --------
    The following example creates a Tapered Flange Channel section with a depth of 10, a width of
    3.5, a mid-flange thickness of 0.575, a web thickness of 0.475, a root radius of 0.575, a
    flange radius of 0.4 and a flange angle of 8°, using 16 points to discretise the radii. A mesh
    is generated with a maximum triangular area of 0.02

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       geometry = sections.TaperedFlangeChannel(
           d=10, b=3.5, t_f=0.575, t_w=0.475, r_r=0.575, r_f=0.4, alpha=8, n_r=16
       )
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[0.02])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, d, b, t_f, t_w, r_r, r_f, alpha, n_r, shift=(0, 0)):
        """Initialises the ISection class."""
        # assign control point
        control_points = [[t_w * 0.5, d * 0.5]]

        super().__init__(control_points, shift)

        # calculate alpha in radians
        alpha_rad = np.pi * alpha / 180

        # calculate the height of the flange toe and dimensions of the straight
        x1 = b * 0.5 - t_w * 0.5 - r_f * (1 - np.sin(alpha_rad))
        y1 = x1 * np.tan(alpha_rad)
        x2 = b * 0.5 - t_w * 0.5 - r_r * (1 - np.sin(alpha_rad))
        y2 = x2 * np.tan(alpha_rad)
        y_t = t_f - y1 - r_f * np.cos(alpha_rad)

        # add first two points
        self.points.append([0, 0])
        self.points.append([b, 0])

        # construct the bottom right flange toe radius
        if r_f == 0:
            self.points.append([b, y_t])
        else:
            for i in range(n_r):
                # determine polar angle
                theta = i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)

                # calculate the locations of the radius points
                x = b - r_f + r_f * np.cos(theta)
                y = y_t + r_f * np.sin(theta)

                # append the current points to the points list
                self.points.append([x, y])

        # construct the bottom right root radius
        if r_r == 0:
            self.points.append([t_w, t_f + y2])
        else:
            for i in range(n_r):
                # determine polar angle
                theta = (3.0 / 2 * np.pi - alpha_rad) - (
                    i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)
                )

                # calculate the locations of the radius points
                x = t_w + r_r + r_r * np.cos(theta)
                y = t_f + y2 + r_r * np.cos(alpha_rad) + r_r * np.sin(theta)

                # append the current points to the points list
                self.points.append([x, y])

        # construct the top right root radius
        if r_r == 0:
            self.points.append([t_w, d - t_f - y2])
        else:
            for i in range(n_r):
                # determine polar angle
                theta = np.pi - i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)

                # calculate the locations of the radius points
                x = t_w + r_r + r_r * np.cos(theta)
                y = d - t_f - y2 - r_r * np.cos(alpha_rad) + r_r * np.sin(theta)

                # append the current points to the points list
                self.points.append([x, y])

        # construct the top right flange toe radius
        if r_f == 0:
            self.points.append([b, d - y_t])
        else:
            for i in range(n_r):
                # determine polar angle
                theta = (3.0 * np.pi / 2 + alpha_rad) + (
                    i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)
                )

                # calculate the locations of the radius points
                x = b - r_f + r_f * np.cos(theta)
                y = d - y_t + r_f * np.sin(theta)

                # append the current points to the points list
                self.points.append([x, y])

        # add the final two points
        self.points.append([b, d])
        self.points.append([0, d])

        # build the facet list
        for i in range(len(self.points)):
            # if we are not at the last point
            if i != len(self.points) - 1:
                self.facets.append([i, i + 1])
            # if we are at the last point, complete the loop
            else:
                self.facets.append([len(self.points) - 1, 0])

        self.perimeter = list(range(len(self.facets)))

        self.shift_section()


class TeeSection(Geometry):
    """Constructs a Tee section with the top left corner at *(0, d)*, with depth *d*, width *b*,
    flange thickness *t_f*, web thickness *t_w* and root radius *r*, using *n_r* points to
    construct the root radius.

    Parameters
    ----------
    d : float
        Depth of the Tee section
    b : float
        Width of the Tee section
    t_f : float
        Flange thickness of the Tee section
    t_w : float
        Web thickness of the Tee section
    r : float
        Root radius of the Tee section
    n_r : int
        Number of points discretising the root radius
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*

    Examples
    --------
    The following example creates a Tee section with a depth of 200, a width of 100, a flange
    thickness of 12, a web thickness of 6 and a root radius of 8, using 8 points to discretise the
    root radius. A mesh is generated with a maximum triangular area of 3.0

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       geometry = sections.TeeSection(d=200, b=100, t_f=12, t_w=6, r=8, n_r=8)
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[3.0])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, d, b, t_f, t_w, r, n_r, shift=(0, 0)):
        """Initialises the TeeSection class."""
        # assign control point
        control_points = [[b * 0.5, d - t_f * 0.5]]

        super().__init__(control_points, shift)

        # add first two points
        self.points.append([b * 0.5 - t_w * 0.5, 0])
        self.points.append([b * 0.5 + t_w * 0.5, 0])

        # construct the top right radius
        pt = [b * 0.5 + t_w * 0.5 + r, d - t_f - r]
        self.draw_radius(pt, r, np.pi, n_r, False)

        # add next four points
        self.points.append([b, d - t_f])
        self.points.append([b, d])
        self.points.append([0, d])
        self.points.append([0, d - t_f])

        # construct the top left radius
        pt = [b * 0.5 - t_w * 0.5 - r, d - t_f - r]
        self.draw_radius(pt, r, 0.5 * np.pi, n_r, False)

        # build the facet list
        for i in range(len(self.points)):
            # if we are not at the last point
            if i != len(self.points) - 1:
                self.facets.append([i, i + 1])
            # if we are at the last point, complete the loop
            else:
                self.facets.append([len(self.points) - 1, 0])

        self.perimeter = list(range(len(self.facets)))

        self.shift_section()


class AngleSection(Geometry):
    """Constructs an angle section with the bottom left corner at the origin *(0, 0)*, with depth
    *d*, width *b*, thickness *t*, root radius *r_r* and toe radius *r_t*, using *n_r* points to
    construct the radii.

    Parameters
    ----------
    d : float
        Depth of the angle section
    b : float
        Width of the angle section
    t : float
        Thickness of the angle section
    r_r : float
        Root radius of the angle section
    r_t : float
        Toe radius of the angle section
    n_r : int
        Number of points discretising the radii
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*

    Examples
    --------
    The following example creates an angle section with a depth of 150, a width of 100, a thickness
    of 8, a root radius of 12 and a toe radius of 5, using 16 points to discretise the radii. A
    mesh is generated with a maximum triangular area of 2.0

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       geometry = sections.AngleSection(d=150, b=100, t=8, r_r=12, r_t=5, n_r=16)
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[2.0])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, d, b, t, r_r, r_t, n_r, shift=(0, 0)):
        """Initialises the AngleSection class."""
        # assign control point
        control_points = [[t * 0.5, t * 0.5]]

        super().__init__(control_points, shift)

        # add first two points
        self.points.append([0, 0])
        self.points.append([b, 0])

        # construct the bottom toe radius
        pt = [b - r_t, t - r_t]
        self.draw_radius(pt, r_t, 0, n_r)

        # construct the root radius
        pt = [t + r_r, t + r_r]
        self.draw_radius(pt, r_r, 1.5 * np.pi, n_r, False)

        # construct the top toe radius
        pt = [t - r_t, d - r_t]
        self.draw_radius(pt, r_t, 0, n_r)

        # add the next point
        self.points.append([0, d])

        # build the facet list
        for i in range(len(self.points)):
            # if we are not at the last point
            if i != len(self.points) - 1:
                self.facets.append([i, i + 1])
            # if we are at the last point, complete the loop
            else:
                self.facets.append([len(self.points) - 1, 0])

        self.perimeter = list(range(len(self.facets)))

        self.shift_section()


class CeeSection(Geometry):
    """Constructs a Cee section with the bottom left corner at the origin *(0, 0)*, with depth *d*,
    width *b*, lip *l*, thickness *t* and outer radius *r_out*, using *n_r* points to construct the
    radius. If the outer radius is less than the thickness of the Cee Section, the inner radius is
    set to zero.

    Parameters
    ----------
    d : float
        Depth of the Cee section
    b : float
        Width of the Cee section
    l : float
        Lip of the Cee section
    t : float
        Thickness of the Cee section
    r_out : float
        Outer radius of the Cee section
    n_r : int
        Number of points discretising the outer radius
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*

    Raises
    ------
    Exception
        Lip length must be greater than the outer radius

    Examples
    --------
    The following example creates a Cee section with a depth of 125, a width of 50, a lip of 30, a
    thickness of 1.5 and an outer radius of 6, using 8 points to discretise the radius. A mesh is
    generated with a maximum triangular area of 0.25

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       geometry = sections.CeeSection(d=125, b=50, l=30, t=1.5, r_out=6, n_r=8)
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[0.25])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, d, b, l, t, r_out, n_r, shift=(0, 0)):
        """Initialises the CeeSection class."""
        # ensure the lip length is greater than the outer radius
        if l < r_out:
            raise Exception('Lip length must be greater than the outer radius')

        # assign control point
        control_points = [[t * 0.5, d * 0.5]]

        super().__init__(control_points, shift)

        # calculate internal radius
        r_in = max(r_out - t, 0)

        # construct the outer bottom left radius
        self.draw_radius([r_out, r_out], r_out, np.pi, n_r)

        # construct the outer bottom right radius
        self.draw_radius([b - r_out, r_out], r_out, 1.5 * np.pi, n_r)

        if r_out != l:
            # add next two points
            self.points.append([b, l])
            self.points.append([b - t, l])

        # construct the inner bottom right radius
        self.draw_radius([b - t - r_in, t + r_in], r_in, 0, n_r, False)

        # construct the inner bottom left radius
        self.draw_radius([t + r_in, t + r_in], r_in, 1.5 * np.pi, n_r, False)

        # construct the inner top left radius
        self.draw_radius([t + r_in, d - t - r_in], r_in, np.pi, n_r, False)

        # construct the inner top right radius
        self.draw_radius([b - t - r_in, d - t - r_in], r_in, 0.5 * np.pi, n_r, False)

        if r_out != l:
            # add next two points
            self.points.append([b - t, d - l])
            self.points.append([b, d - l])

        # construct the outer top right radius
        self.draw_radius([b - r_out, d - r_out], r_out, 0, n_r)

        # construct the outer top left radius
        self.draw_radius([r_out, d - r_out], r_out, 0.5 * np.pi, n_r)

        # build the facet list
        for i in range(len(self.points)):
            # if we are not at the last point
            if i != len(self.points) - 1:
                self.facets.append([i, i + 1])
            # if we are at the last point, complete the loop
            else:
                self.facets.append([len(self.points) - 1, 0])

        self.perimeter = list(range(len(self.facets)))

        self.shift_section()


class ZedSection(Geometry):
    """Constructs a Zed section with the bottom left corner at the origin *(0, 0)*, with depth *d*,
    left flange width *b_l*, right flange width *b_r*, lip *l*, thickness *t* and outer radius
    *r_out*, using *n_r* points to construct the radius. If the outer radius is less than the
    thickness of the Zed Section, the inner radius is set to zero.

    Parameters
    ----------
    d : float
        Depth of the Zed section
    b_l : float
        Left flange width of the Zed section
    b_r : float
        Right flange width of the Zed section
    l : float
        Lip of the Zed section
    t : float
        Thickness of the Zed section
    r_out : float
        Outer radius of the Zed section
    n_r : int
        Number of points discretising the outer radius
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*

    Raises
    ------
    Exception
        Lip length must be greater than the outer radius

    Examples
    --------
    The following example creates a Zed section with a depth of 100, a left flange width of 40, a
    right flange width of 50, a lip of 20, a thickness of 1.2 and an outer radius of 5, using 8
    points to discretise the radius. A mesh is generated with a maximum triangular area of 0.15

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       geometry = sections.ZedSection(d=100, b_l=40, b_r=50, l=20, t=1.2, r_out=5, n_r=8)
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[0.15])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, d, b_l, b_r, l, t, r_out, n_r, shift=(0, 0)):
        """Initialises the ZedSection class."""
        # ensure the lip length is greater than the outer radius
        if l < r_out:
            raise Exception('Lip length must be greater than the outer radius')

        # assign control point
        control_points = [[t * 0.5, d * 0.5]]

        super().__init__(control_points, shift)

        # calculate internal radius
        r_in = max(r_out - t, 0)

        # construct the outer bottom left radius
        self.draw_radius([r_out, r_out], r_out, np.pi, n_r)

        # construct the outer bottom right radius
        self.draw_radius([b_r - r_out, r_out], r_out, 1.5 * np.pi, n_r)

        if r_out != l:
            # add next two points
            self.points.append([b_r, l])
            self.points.append([b_r - t, l])

        # construct the inner bottom right radius
        self.draw_radius([b_r - t - r_in, t + r_in], r_in, 0, n_r, False)

        # construct the inner bottom left radius
        self.draw_radius([t + r_in, t + r_in], r_in, 1.5 * np.pi, n_r, False)

        # construct the outer top right radius
        self.draw_radius([t - r_out, d - r_out], r_out, 0, n_r)

        # construct the outer top left radius
        self.draw_radius([t - b_l + r_out, d - r_out], r_out, 0.5 * np.pi, n_r)

        if r_out != l:
            # add the next two points
            self.points.append([t - b_l, d - l])
            self.points.append([t - b_l + t, d - l])

        # construct the inner top left radius
        self.draw_radius([2 * t - b_l + r_in, d - t - r_in], r_in, np.pi, n_r, False)

        # construct the inner top right radius
        self.draw_radius([-r_in, d - t - r_in], r_in, 0.5 * np.pi, n_r, False)

        # build the facet list
        for i in range(len(self.points)):
            # if we are not at the last point
            if i != len(self.points) - 1:
                self.facets.append([i, i + 1])
            # if we are at the last point, complete the loop
            else:
                self.facets.append([len(self.points) - 1, 0])

        self.perimeter = list(range(len(self.facets)))

        self.shift_section()


class CruciformSection(Geometry):
    """Constructs a cruciform section centered at the origin *(0, 0)*, with depth *d*, width *b*,
    thickness *t* and root radius *r*, using *n_r* points to construct the root radius.

    Parameters
    ----------
    d : float
        Depth of the cruciform section
    b : float
        Width of the cruciform section
    t : float
        Thickness of the cruciform section
    r : float
        Root radius of the cruciform section
    n_r : int
        Number of points discretising the root radius
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*

    Examples
    --------
    The following example creates a cruciform section with a depth of 250, a width of 175, a
    thickness of 12 and a root radius of 16, using 16 points to discretise the radius. A mesh is
    generated with a maximum triangular area of 5.0

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       geometry = sections.CruciformSection(d=250, b=175, t=12, r=16, n_r=16)
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[5.0])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, d, b, t, r, n_r, shift=(0, 0)):
        """Initialises the CruciformSection class."""
        # assign control point
        control_points = [[0, 0]]

        super().__init__(control_points, shift)

        # add first two points
        self.points.append([-t * 0.5, -d * 0.5])
        self.points.append([t * 0.5, -d * 0.5])

        # construct the bottom right radius
        pt = [0.5 * t + r, -0.5 * t - r]
        self.draw_radius(pt, r, np.pi, n_r, False)

        # add the next two points
        self.points.append([0.5 * b, -t * 0.5])
        self.points.append([0.5 * b, t * 0.5])

        # construct the top right radius
        pt = [0.5 * t + r, 0.5 * t + r]
        self.draw_radius(pt, r, 1.5 * np.pi, n_r, False)

        # add the next two points
        self.points.append([t * 0.5, 0.5 * d])
        self.points.append([-t * 0.5, 0.5 * d])

        # construct the top left radius
        pt = [-0.5 * t - r, 0.5 * t + r]
        self.draw_radius(pt, r, 0, n_r, False)

        # add the next two points
        self.points.append([-0.5 * b, t * 0.5])
        self.points.append([-0.5 * b, -t * 0.5])

        # construct the bottom left radius
        pt = [-0.5 * t - r, -0.5 * t - r]
        self.draw_radius(pt, r, 0.5 * np.pi, n_r, False)

        # build the facet list
        for i in range(len(self.points)):
            # if we are not at the last point
            if i != len(self.points) - 1:
                self.facets.append([i, i + 1])
            # if we are at the last point, complete the loop
            else:
                self.facets.append([len(self.points) - 1, 0])

        self.perimeter = list(range(len(self.facets)))

        self.shift_section()


class PolygonSection(Geometry):
    """Constructs a regular hollow polygon section centered at *(0, 0)*, with a pitch circle
    diameter of bounding polygon *d*, thickness *t*, number of sides *n_sides* and an optional
    inner radius *r_in*, using *n_r* points to construct the inner and outer radii (if radii is
    specified).

    Parameters
    ----------
    d : float
        Pitch circle diameter of the outer bounding polygon (i.e. diameter of circle that passes
        through all vertices of the outer polygon)
    t : float
        Thickness of the polygon section wall
    r_in : float
        Inner radius of the polygon corners. By default, if not specified, a polygon with no corner
        radii is generated.
    n_r : int
        Number of points discretising the inner and outer radii, ignored if no inner radii is
        specified
    rot
        Initial counterclockwise rotation in degrees. By default bottom face is aligned with x axis.
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*

    Raises
    ------
    Exception
        Number of sides in polygon must be greater than or equal to 3

    Examples
    --------
    The following example creates an Octagonal section (8 sides) with a diameter of 200, a
    thickness of 6 and an inner radius of 20, using 12 points to discretise the inner and outer
    radii. A mesh is generated with a maximum triangular area of 5

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       geometry = sections.PolygonSection(d=200, t=6, n_sides=8, r_in=20, n_r=12)
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[5])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, d, t, n_sides, r_in=0, n_r=1, rot=0, shift=(0, 0)):
        """Initialises the PolygonSection class."""
        if n_sides < 3:
            msg = 'n_sides required to be greater than 3 for PolygonSection class'
            raise Exception(msg)

        # initial rotation
        rot = rot * np.pi / 180  # radians

        # determine triangular segment angle
        alpha = 2 * np.pi / n_sides  # radians

        # determine distance from origin to point perpendicular on face of side
        a_out = d / 2 * np.cos(alpha / 2)
        a_in = a_out - t

        # determine side length for outer & inner faces neglecting radii
        side_length_out = d * np.sin(alpha / 2)
        side_length_in = a_in / a_out * side_length_out

        # check limit on internal radii, if exceeded then radii merge to circle
        if r_in > a_in:
            r_in = a_in
            circle = True
        else:
            circle = False

        # calculate external radius, if r_in is zero, r_out also is zero
        if r_in == 0:
            r_out = 0
            n_r = 1
        else:
            r_out = r_in + t

        # equivalent side length of half the corner radii triangular segment
        c_out = r_out * (side_length_out / 2) / a_out
        c_in = r_in * (side_length_in / 2) / a_in

        # determine straight side length between corner radii (if present)
        side_length_straight_out = side_length_out - (2 * c_out)
        side_length_straight_in = side_length_in - (2 * c_in)

        # assign control point central on bottom side length & rotate to initial rotation specified
        control_points = [self.rotate([0, -a_out + t / 2], rot)]

        super().__init__(control_points, shift)

        # temp list for repeating geometry
        base_points = []

        # specify a hole in the centre of the Polygon section
        self.holes = [[0, 0]]

        # start at bottom face, constructing one corner radii, then rotate by initial rotation +
        # alpha and repeat for n_side number of times to form full section perimeter

        # construct the first radius (bottom right)
        for i in range(n_r):
            # determine polar angle
            theta = 1 / 2 * np.pi + i * 1.0 / max(1, n_r - 1) * alpha

            # calculate location of inner and outer points
            x_outer = side_length_straight_out / 2 - r_out * np.cos(theta)
            y_outer = -a_out + r_out - r_out * np.sin(theta)
            x_inner = side_length_straight_in / 2 - r_in * np.cos(theta)
            y_inner = -a_in + r_in - r_in * np.sin(theta)

            # append the current temporary points to the temporary points list
            base_points.append([x_outer, y_outer])
            base_points.append([x_inner, y_inner])

        # if radii merged to circle with an outer diameter of a_out then skip last point as causes
        # overlapping end points which causes meshing issues if geometry is not cleaned by user
        if circle:
            base_points = base_points[0:-2]

        # iterate and add subsequent corner radii one point at a time for each side
        for i in range(n_sides):
            for point in base_points:
                point_new = self.rotate(point, alpha * i + rot)
                self.points.append(point_new)

        # build the facet list
        num_points = int(len(self.points) / 2)
        for i in range(num_points):
            # if we are not at the last point
            if i != num_points - 1:
                self.facets.append([i * 2, i * 2 + 2])
                self.facets.append([i * 2 + 1, i * 2 + 3])
            # if we are at the last point, complete the loop
            else:
                self.facets.append([i * 2, 0])
                self.facets.append([i * 2 + 1, 1])

        self.perimeter = list(range(0, len(self.facets), 2))

        self.shift_section()

    @staticmethod
    def rotate(point, angle):
        """Rotate a point counterclockwise by a given angle around origin [0, 0]

        Parameters
        ----------
        point : list
            Point coordinates to be rotated
        angle : float
            Angle to rotate point coordinates

        Returns
        -------
        list[float, float]
            Coordinates of rotated point
        """
        pt_x, pt_y = point

        c = np.cos(angle)
        s = np.sin(angle)

        new_x = c * pt_x - s * pt_y
        new_y = s * pt_x + c * pt_y

        return [new_x, new_y]


class BoxGirderSection(Geometry):
    """Constructs a Box Girder section centered at at *(max(b_t, b_b)/2, d/2)*, with depth *d*, top
    width *b_t*, bottom width *b_b*, top flange thickness *t_ft*, bottom flange thickness *t_fb*
    and web thickness *t_w*.

    Parameters
    ----------
    d : float
        Depth of the Box Girder section
    b_t : float
        Top width of the Box Girder section
    b_b : float
        Bottom width of the Box Girder section
    t_ft : float
        Top flange thickness of the Box Girder section
    t_fb : float
        Bottom flange thickness of the Box Girder section
    t_w : float
        Web thickness of the Box Girder section
    shift : list[float, float]
        Vector that shifts the cross-section by *(x, y)*

    Examples
    --------
    The following example creates a Box Girder section with a depth of 1200, a top width of 1200, a
    bottom width of 400, a top flange thickness of 16, a bottom flange thickness of 12 and a web
    thickness of 8. A mesh is generated with a maximum triangular area of 5.0

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       geometry = sections.BoxGirderSection(d=1200, b_t=1200, b_b=400, t_ft=100, t_fb=80, t_w=50)
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[200.0])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, d, b_t, b_b, t_ft, t_fb, t_w, shift=(0, 0)):
        """Initialises the BoxGirderSection class."""
        # assign control point
        control_points = [[max(b_t, b_b) * 0.5, t_fb * 0.5]]

        super().__init__(control_points, shift)

        # calculate central axis
        x_c = max(b_t, b_b) * 0.5

        # specify a hole in the centre of the Box Girder
        self.holes = [[x_c, d * 0.5]]

        # determine side wall angle
        if b_t < b_b:
            phi_b = np.arctan2(d, 0.5 * (b_b - b_t))
            phi_t = np.pi - phi_b
        else:
            phi_t = np.arctan2(d, 0.5 * (b_t - b_b))
            phi_b = np.pi - phi_t

        # determine inner wall x-offsets
        x_bot = t_fb / np.tan(np.pi - phi_b)
        x_top = t_ft / np.tan(np.pi - phi_t)
        web_x = abs(t_w / np.sin(np.pi - phi_b))

        # add outer points
        self.points.append([x_c - 0.5 * b_b, 0])
        self.points.append([x_c + 0.5 * b_b, 0])
        self.points.append([x_c + 0.5 * b_t, d])
        self.points.append([x_c - 0.5 * b_t, d])

        # add inner points
        self.points.append([x_c - 0.5 * b_b - x_bot + web_x, t_fb])
        self.points.append([x_c + 0.5 * b_b + x_bot - web_x, t_fb])
        self.points.append([x_c + 0.5 * b_t + x_top - web_x, d - t_ft])
        self.points.append([x_c - 0.5 * b_t - x_top + web_x, d - t_ft])

        # build facet list
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4]]
        self.perimeter = [0, 1, 2, 3]

        self.shift_section()


class MergedSection(Geometry):
    """Merges a number of section geometries into one geometry. Note that for the meshing algorithm
    to work, there needs to be connectivity between all regions of the provided geometries.
    Overlapping of geometries is permitted.

    Parameters
    ----------
    sections : list[:class:`~sectionproperties.pre.sections.Geometry`]
        A list of geometry objects to merge into one
        :class:`~sectionproperties.pre.sections.Geometry` object

    Examples
    --------
    The following example creates a combined cross-section with a 150x100x6 RHS placed on its side
    on top of a 200UB25.4. A mesh is generated with a maximum triangle size of 5.0 for the
    I-section and 2.5 for the RHS

    .. plot::
       :context: reset

       import sectionproperties.pre.sections as sections
       from sectionproperties.analysis.cross_section import CrossSection
       isection = sections.ISection(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8)
       box = sections.Rhs(d=100, b=150, t=6, r_out=15, n_r=8, shift=[-8.5, 203])
       geometry = sections.MergedSection([isection, box])
       geometry.clean_geometry()
       geometry.plot_geometry()
       mesh = geometry.create_mesh(mesh_sizes=[5.0, 2.5])
       section = CrossSection(geometry, mesh)
       section.plot_mesh()
    """

    def __init__(self, sections):
        """Initialises the MergedSection class."""
        super().__init__([], [0, 0])

        point_count = 0

        # loop through all sections
        for section in sections:
            # add facets
            for facet in section.facets:
                self.facets.append([facet[0] + point_count, facet[1] + point_count])

            # add points and count points
            for point in section.points:
                self.points.append([point[0], point[1]])
                point_count += 1

            # add holes
            for hole in section.holes:
                self.holes.append([hole[0], hole[1]])

            # add control points
            for control_point in section.control_points:
                self.control_points.append([control_point[0], control_point[1]])
