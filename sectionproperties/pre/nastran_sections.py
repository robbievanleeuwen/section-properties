from sectionproperties.pre.sections import *


#TODO Update documentation.

class BoxSection(Geometry):
    """
    Constructs a rectangular box section with the center at the
    origin *(0, 0)*, with four parameters defining dimensions.
    See Nastran documentation for definition of parameters.

    :param float DIM1: Width (x) of box
    :param float DIM2: Depth (y) of box
    :param float DIM3: Thickness of box in y direction
    :param float DIM4: Thickness of box in x direction
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    Added by JohnDN90.

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

    def __init__(self, DIM1, DIM2, DIM3, DIM4, shift=[0, 0]):
        """Inits the BoxSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0

        # assign control point
        control_points = [[0, DIM2 / 2 - DIM3 / 2]]

        super().__init__(control_points, shift)

        # specify a hole in the centre of the Box
        self.holes = [[0, 0]]

        # construct the points and facets
        self.points = [[-DIM1 / 2., -DIM2 / 2.], [DIM1 / 2., -DIM2 / 2.], [DIM1 / 2., DIM2 / 2.], [-DIM1 / 2., DIM2 / 2.],
                       [-DIM1 / 2. + DIM4, -DIM2 / 2. + DIM3], [DIM1 / 2. - DIM4, -DIM2 / 2. + DIM3],
                       [DIM1 / 2. - DIM4, DIM2 / 2. - DIM3], [-DIM1 / 2. + DIM4, DIM2 / 2. - DIM3]]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4]]

        self.shift_section()


class ISection(Geometry):
    """
    Constructs a I section with the bottom flange's middle center at
    the origin *(0, 0)*, with six parameters defining dimensions.
    See Nastran documentation for definition of parameters.

    :param float DIM1: Depth(y) of the I-section
    :param float DIM2: Width (x) of bottom flange
    :param float DIM3: Width (x) of top flange
    :param float DIM4: Thickness of web
    :param float DIM5: Thickness of bottom web
    :param float DIM6: Thickness of top web
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    Added by JohnDN90.

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

    def __init__(self, DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, shift=[0, 0]):
        """Inits the BoxSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        DIM5 *= 1.0
        DIM6 *= 1.0

        # assign control point
        control_points = [[0.5*DIM2, 0.5*DIM5]]

        shift = [-0.5*DIM2+shift[0], -0.5*DIM1+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        db = (DIM2 - DIM4) / 2
        dt = (DIM3 - DIM4) / 2
        self.points = [[0.,0.], [DIM2,0.], [DIM2, DIM5], [db+DIM4, DIM5], [db + DIM4, DIM1-DIM6], [db+DIM4+dt, DIM1-DIM6],
                       [db+DIM4+dt, DIM1], [db-dt, DIM1], [db-dt, DIM1-DIM6], [db, DIM1-DIM6], [db, DIM5], [0, DIM5]]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4,5], [5,6], [6,7], [7,8], [8,9], [9,10], [10,11], [11,0]]

        self.shift_section()


class TSection(Geometry):
    """
    Constructs a T section with the top flange's middle center at
    the origin *(0, 0)*, with four parameters defining dimensions.
    See Nastran documentation for more details.

    :param float DIM1: Width (x) of top flange
    :param float DIM2: Depth (y) of the T-section.
    :param float DIM3: Thickness of top flange.
    :param float DIM4: Thickness of web.
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    Added by JohnDN90.

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

    def __init__(self, DIM1, DIM2, DIM3, DIM4, shift=[0, 0]):
        """Inits the BoxSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0

        d = DIM2
        b = DIM1
        t_f = DIM3
        t_w = DIM4
        r = 0
        n_r = 1
        shift = [-DIM1/2.0+shift[0], -(DIM2-DIM3/2.0)+shift[1]]

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

        self.shift_section()


class LSection(Geometry):
    """
    Constructs a L section with the intersection's center at
    the origin *(0, 0)*, with four parameters defining dimensions.
    See Nastran documentation for more details.

    :param float DIM1: Width (x) of the L-section.
    :param float DIM2: Depth (y) of the L-section.
    :param float DIM3: Thickness of flange (horizontal portion).
    :param float DIM4: Thickness of web (vertical portion).
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    Added by JohnDN90.

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

    def __init__(self, DIM1, DIM2, DIM3, DIM4, shift=[0, 0]):
        """Inits the BoxSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0

        # assign control point
        control_points = [[0.5*DIM1, 0.5*DIM3]]

        shift = [-0.5*DIM4+shift[0], -DIM3/3.0+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        self.points = [[0,0], [DIM1, 0], [DIM1, DIM3], [DIM4, DIM3], [DIM4, DIM2], [0, DIM2]]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4,5], [5,0]]

        self.shift_section()


class ChanSection(Geometry):
    """
    Constructs a CHAN (C-Channel) section with the web's middle
    center at the origin *(0, 0)*, with six parameters defining
    dimensions. See Nastran documentation for more details.

    :param float DIM1: Width (x) of the Chan section.
    :param float DIM2: Depth (y) of the Chan section.
    :param float DIM3: Thickness of web (vertical portion).
    :param float DIM4: Thickness of flanges (top/bottom portion).
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    Added by JohnDN90.

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

    def __init__(self, DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, shift=[0, 0]):
        """Inits the BoxSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        DIM5 *= 1.0
        DIM6 *= 1.0

        # assign control point
        control_points = [[0.5*DIM1, 0.5*DIM4]]

        shift = [-0.5*DIM3+shift[0], -0.5*DIM2+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        self.points = [[0.,0.], [DIM1, 0.], [DIM1, DIM4], [DIM3, DIM4], [DIM3, DIM2-DIM4], [DIM1, DIM2-DIM4], [DIM1, DIM2], [0., DIM2]]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4,5], [5,6], [6,7], [7,0]]

        self.shift_section()


# class Section(Geometry):
#     """
#     Constructs a I section with the bottom web's middle center at
#     the origin *(0, 0)*, with six parameters defining dimensions.
#     See Nastran documentation for more details.
#
#     :param float DIM1:
#     :param float DIM2:
#     :param float DIM3:
#     :param float DIM4:
#     :param float DIM5:
#     :param float DIM6:
#     :param shift: Vector that shifts the cross-section by *(x, y)*
#     :type shift: list[float, float]
#
#     Added by JohnDN90.
#
#     The following example creates a rectangular cross-section with a depth of
#     100 and width of 50, and generates a mesh with a maximum triangular area of
#     5::
#
#         import sectionproperties.pre.sections as sections
#
#         geometry = sections.RectangularSection(d=100, b=50)
#         mesh = geometry.create_mesh(mesh_sizes=[5])
#
#     ..  figure:: ../images/sections/rectangle_geometry.png
#         :align: center
#         :scale: 75 %
#
#         Rectangular section geometry.
#
#     ..  figure:: ../images/sections/rectangle_mesh.png
#         :align: center
#         :scale: 75 %
#
#         Mesh generated from the above geometry.
#     """
#
#     def __init__(self, DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, shift=[0, 0]):
#         """Inits the BoxSection class."""
#
#         # force dimensions to be floating point values
#         DIM1 *= 1.0
#         DIM2 *= 1.0
#         DIM3 *= 1.0
#         DIM4 *= 1.0
#         DIM5 *= 1.0
#         DIM6 *= 1.0
#
#         # assign control point
#         control_points =
#
#         shift = [+shift[0], +shift[1]]
#         super().__init__(control_points, shift)
#
#         # construct the points and facets
#         self.points =
#         self.facets =
#
#         self.shift_section()