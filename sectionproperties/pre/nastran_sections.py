from sectionproperties.pre.sections import *


#TODO Update documentation.

class BarSection(Geometry):
    """
    Constructs a rectangular section with the center at the
    origin *(0, 0)*, with two parameters defining dimensions.
    See Nastran documentation for definition of parameters.

    :param float DIM1: Width (x) of bar
    :param float DIM2: Depth (y) of bar
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

    def __init__(self, DIM1, DIM2, shift=[0, 0]):
        """Inits the BoxSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0

        # assign control point
        control_points = [[0., 0.]]

        shift = [-0.5*DIM1+shift[0], -0.5*DIM2+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        self.points = [[-DIM1 / 2., -DIM2 / 2.], [DIM1 / 2., -DIM2 / 2.], [DIM1 / 2., DIM2 / 2.], [-DIM1 / 2., DIM2 / 2.]]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 0]]

        self.shift_section()


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
        control_points = [[0., DIM2 / 2.0 - DIM3 / 2.0]]

        super().__init__(control_points, shift)

        # specify a hole in the centre of the Box
        self.holes = [[0., 0.]]

        # construct the points and facets
        self.points = [[-DIM1 / 2., -DIM2 / 2.], [DIM1 / 2., -DIM2 / 2.], [DIM1 / 2., DIM2 / 2.], [-DIM1 / 2., DIM2 / 2.],
                       [-DIM1 / 2. + DIM4, -DIM2 / 2. + DIM3], [DIM1 / 2. - DIM4, -DIM2 / 2. + DIM3],
                       [DIM1 / 2. - DIM4, DIM2 / 2. - DIM3], [-DIM1 / 2. + DIM4, DIM2 / 2. - DIM3]]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4]]

        self.shift_section()


class Box1Section(Geometry):
    """
    Constructs a Box1 section with the bottom web's middle center at
    the origin *(0, 0)*, with six parameters defining dimensions.
    See Nastran documentation for more details.

    :param float DIM1:
    :param float DIM2:
    :param float DIM3:
    :param float DIM4:
    :param float DIM5:
    :param float DIM6:
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

        shift = [-0.5*DIM1+shift[0], -0.5*DIM2+shift[1]]
        super().__init__(control_points, shift)

        # specify a hole in the centre of the Box
        self.holes = [[DIM6 + 0.5*(DIM1-DIM5), DIM4+0.5*(DIM2-DIM3)]]

        # construct the points and facets
        self.points = [ [0.,0.], [DIM1,0.], [DIM1, DIM2], [0., DIM2],
                        [DIM6, DIM4], [DIM1-DIM5, DIM4], [DIM1-DIM5, DIM2-DIM3], [DIM6, DIM2-DIM3]]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4]]

        self.shift_section()


class HexaSection(Geometry):
    """
    Constructs a Hexa section with the bottom web's middle center at
    the origin *(0, 0)*, with six parameters defining dimensions.
    See Nastran documentation for more details.

    :param float DIM1:
    :param float DIM2:
    :param float DIM3:
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

    def __init__(self, DIM1, DIM2, DIM3, shift=[0, 0]):
        """Inits the BoxSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0

        # assign control point
        control_points = [[0.5*DIM2, 0.5*DIM3]]

        shift = [-0.5*DIM2+shift[0], -0.5*DIM3+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        self.points = [[DIM1, 0.], [DIM2-DIM1, 0.], [DIM2, 0.5*DIM3], [DIM2-DIM1, DIM3], [DIM1, DIM3], [0., 0.5*DIM3]]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]

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


class I1Section(Geometry):
    """
    Constructs a I1 section with the right flange's middle center at
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
        """Inits the T1section class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0

        shift = [-0.5*(DIM1+DIM2)+shift[0], -0.5*DIM4+shift[1]]

        # assign control point
        control_points = [[0.5*(DIM1+DIM2), 0.5*DIM4]]

        super().__init__(control_points, shift)

        # construct the points and facets
        t = (DIM4 - DIM3) / 2
        self.points = [[0,0], [DIM1+DIM2, 0], [DIM1+DIM2, t], [0.5*DIM1+DIM2, t], [0.5*DIM1+DIM2, t+DIM3],
                       [DIM1+DIM2, t+DIM3], [DIM1+DIM2, DIM4], [0, DIM4], [0, t+DIM3], [0.5*DIM1, t+DIM3], [0.5*DIM1, t], [0, t]]
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


class T1Section(Geometry):
    """
    Constructs a T1 section with the right flange's middle center at
    the origin *(0, 0)*, with four parameters defining dimensions.
    See Nastran documentation for more details.

    :param float DIM1:
    :param float DIM2:
    :param float DIM3:
    :param float DIM4:
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
        """Inits the T1section class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0

        shift = [-0.5*DIM3+shift[0], -0.5*DIM1+shift[1]]

        # assign control point
        control_points = [[0.5*DIM3, 0.5*DIM1]]

        super().__init__(control_points, shift)

        # construct the points and facets
        d1 = (DIM1 - DIM4) / 2.0
        self.points = [[0, 0], [DIM3, 0], [DIM3, DIM1], [0, DIM1], [0, d1 + DIM4], [-DIM2, d1 + DIM4], [-DIM2, d1], [0, d1]]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 0]]

        self.shift_section()


class T2Section(Geometry):
    """
    Constructs a T2 section with the bottom web's middle center at
    the origin *(0, 0)*, with six parameters defining dimensions.
    See Nastran documentation for more details.

    :param float DIM1:
    :param float DIM2:
    :param float DIM3:
    :param float DIM4:
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

        shift = [-0.5*DIM1+shift[0], -0.5*DIM3+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        d1 = 0.5*(DIM1 - DIM4)
        self.points = [[0.,0.], [DIM1, 0.], [DIM1, DIM3], [DIM1-d1, DIM3],
                       [DIM1-d1, DIM2], [d1, DIM2], [d1, DIM3], [0, DIM3]]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4,5], [5,6], [6,7], [7,0]]

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

    def __init__(self, DIM1, DIM2, DIM3, DIM4, shift=[0, 0]):
        """Inits the BoxSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0

        # assign control point
        control_points = [[0.5*DIM1, 0.5*DIM4]]

        shift = [-0.5*DIM3+shift[0], -0.5*DIM2+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        self.points = [[0.,0.], [DIM1, 0.], [DIM1, DIM4], [DIM3, DIM4], [DIM3, DIM2-DIM4], [DIM1, DIM2-DIM4], [DIM1, DIM2], [0., DIM2]]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4,5], [5,6], [6,7], [7,0]]

        self.shift_section()


class Chan1Section(Geometry):
    """
    Constructs a CHAN1 (C-Channel) section with the web's middle
    center at the origin *(0, 0)*, with six parameters defining
    dimensions. See Nastran documentation for more details.

    :param float DIM1:
    :param float DIM2:
    :param float DIM3:
    :param float DIM4:
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
        """Inits the Chan1Section class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0

        # assign control point
        control_points = [[0.5*DIM1, 0.5*DIM4]]

        shift = [-0.5*DIM2+shift[0], -0.5*DIM4+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        tf = 0.5 * (DIM4 - DIM3)
        self.points = [[0,0], [DIM1+DIM2, 0], [DIM1+DIM2, tf], [DIM2, tf],
                       [DIM2, tf+DIM3], [DIM2+DIM1, tf+DIM3], [DIM2+DIM1, DIM4], [0, DIM4]]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4,5], [5,6], [6,7], [7,0]]

        self.shift_section()


class Chan2Section(Geometry):
    """
    Constructs a Chan2 section with the bottom web's middle center at
    the origin *(0, 0)*, with six parameters defining dimensions.
    See Nastran documentation for more details.

    :param float DIM1:
    :param float DIM2:
    :param float DIM3:
    :param float DIM4:
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
        control_points = [[0.5*DIM1, 0.5*DIM4]]

        shift = [-0.5*DIM4+shift[0], -0.5*DIM2+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        self.points = [[0.,0.], [DIM4,0.], [DIM4,DIM3], [DIM4-DIM1,DIM3],
                       [DIM4-DIM1,DIM2], [DIM1,DIM2], [DIM1,DIM3], [0.,DIM3]]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4,5], [5,6], [6,7], [7,0]]

        self.shift_section()


class RodSection(Geometry):
    """
    Constructs a circular rod section with the center at
    the origin *(0, 0)*, with one parameter defining dimensions.
    See Nastran documentation for more details.

    :param float DIM1: Radius of the circular rod section
    :param int n: Number of points discretising the circle
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

    def __init__(self, DIM1, n, shift=[0, 0]):
        """Inits the BoxSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0

        # assign control point
        control_points = [[0, 0]]

        super().__init__(control_points, shift)

        # loop through each point on the circle
        d = 2.0*DIM1
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


class TubeSection(Geometry):
    """
    Constructs a circular tube section with the center at
    the origin *(0, 0)*, with one parameter defining dimensions.
    See Nastran documentation for more details.

    :param float DIM1: Outer radius of the circular tube section
    :param float DIM2: Inner radius of the circular tube section
    :param int n: Number of points discretising the circle
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

    def __init__(self, DIM1, DIM2, n, shift=[0, 0]):
        """Inits the TubeSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0

        d = 2.0*DIM1
        t = DIM1-DIM2

        # assign control point
        control_points = [[d * 0.5 - t * 0.5, 0]]

        super().__init__(control_points, shift)

        # specify a hole in the centre of the CHS
        self.holes = [[0., 0.]]

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


class Hat1Section(Geometry):
    """
    Constructs a Hat1 section with the bottom plate's bottom center at
    the origin *(0, 0)*, with five parameters defining dimensions.
    See Nastran documentation for definition of parameters.

    :param float DIM1: Width(x) of the Hat1-section
    :param float DIM2: Depth (y) of the Hat1-section
    :param float DIM3: Width (x) of hat's top flange.
    :param float DIM4: Thickness of hat stiffener.
    :param float DIM5: Thicknesss of bottom plate.
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

    def __init__(self, DIM1, DIM2, DIM3, DIM4, DIM5, shift=[0, 0]):
        """Inits the Hat1Section class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        DIM5 *= 1.0

        shift = [-0.5*DIM1+shift[0], shift[1]]

        # create bottom rectangular plate
        bottom_plate = RectangularSection(d=DIM5, b=DIM1, shift=shift)

        # create the hat stiffener
        d1 = DIM2 - DIM5
        d2 = DIM4
        d4 = (DIM1 - DIM3) / 2

        # specify a hole in the combined plate and hat structure
        holes = [[DIM1 / 2, DIM2 / 2]]

        # assign control point
        control_points = [[d4 / 2, DIM5 + DIM4 / 2]]

        super().__init__(control_points, shift)

        # construct the points and facets
        points = [[0, DIM5 + 0], [d4 + d2, DIM5 + 0], [d4 + d2, DIM5 + d1 - d2], [d4 + DIM3 - d2, DIM5 + d1 - d2],
                  [d4 + DIM3 - d2, DIM5 + 0], [2 * d4 + DIM3, DIM5 + 0],
                  [2 * d4 + DIM3, DIM5 + d2], [d4 + DIM3, DIM5 + d2], [d4 + DIM3, DIM5 + d1], [d4, DIM5 + d1],
                  [d4, DIM5 + d2], [0, DIM5 + d2]]
        facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10], [10, 11], [11, 0]]

        hat = CustomSection(points, facets, holes, control_points, shift=shift)

        # Create a list of the sections to merge
        section_list = [bottom_plate, hat]

        # Merge the three sections into one geometry
        geometry = MergedSection(section_list)

        # Clean the geometry and print information to the terminal
        geometry.clean_geometry(verbose=False)

        self.control_points = geometry.control_points
        self.shift = geometry.shift
        self.points = geometry.points
        self.facets = geometry.facets
        self.holes = geometry.holes


    def create_mesh(self, mesh_sizes):
        """Creates a quadratic triangular mesh from the Geometry object.
        This is overloaded here to allow specifying only one mesh_size
        which is used for both regions in the Hat1 section.

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
        mesh_sizes *= 2
        str = "Number of mesh_sizes ({0}), ".format(len(mesh_sizes))
        str += "should match the number of regions "
        str += "({0}).".format(len(self.control_points))
        assert(len(mesh_sizes) == len(self.control_points)), str

        return pre.create_mesh(self.points, self.facets, self.holes,
                               self.control_points, mesh_sizes)


class CruciformSection(Geometry):
    """
    Constructs a cruciform/cross section with the intersection's middle
    center at the origin *(0, 0)*, with four parameters defining
    dimensions. See Nastran documentation for more details.

    :param float DIM1:
    :param float DIM2:
    :param float DIM3:
    :param float DIM4:
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
        control_points = [[0.5*DIM1+0.5*DIM2, 0.5*DIM3]]

        shift = [-(0.5*DIM1+0.5*DIM2)+shift[0], -(0.5*DIM3)+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        d = 0.5*(DIM3 - DIM4)
        self.points = [[0.5*DIM1, 0], [0.5*DIM1+DIM2, 0], [0.5*DIM1+DIM2, d], [DIM1+DIM2, d], [DIM1+DIM2, d+DIM4], [0.5*DIM1+DIM2, d+DIM4],
                       [0.5*DIM1+DIM2, DIM3], [0.5*DIM1, DIM3], [0.5*DIM1, d+DIM4], [0, d+DIM4], [0, d], [0.5*DIM1, d]]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4,5], [5,6], [6,7], [7,8], [8,9], [9,10], [10,11], [11,0]]

        self.shift_section()


class HSection(Geometry):
    """
    Constructs a H section with the bottom web's middle center at
    the origin *(0, 0)*, with six parameters defining dimensions.
    See Nastran documentation for more details.

    :param float DIM1:
    :param float DIM2:
    :param float DIM3:
    :param float DIM4:
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
        """Inits the HSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0

        d1 = 0.5 * (DIM3 - DIM4)
        d2 = 0.5 * DIM2

        # assign control point
        control_points = [[0.5*d2, 0.5*DIM3]]

        shift = [-0.5*(DIM2+DIM1)+shift[0], -0.5*DIM3+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        self.points = [[0,0], [d2,0], [d2,d1], [d2+DIM1, d1], [d2+DIM1,0], [DIM1+DIM2, 0],
                       [DIM1+DIM2, DIM3], [DIM1+DIM2-d2, DIM3], [DIM1+DIM2-d2, d1+DIM4], [d2,d1+DIM4], [d2,DIM3], [0,DIM3]]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4,5], [5,6], [6,7], [7,8], [8,9], [9,10], [10,11], [11,0]]

        self.shift_section()


class ZSection(Geometry):
    """
    Constructs a Z section with the bottom web's middle center at
    the origin *(0, 0)*, with six parameters defining dimensions.
    See Nastran documentation for more details.

    :param float DIM1:
    :param float DIM2:
    :param float DIM3:
    :param float DIM4:
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
        control_points = [[DIM1+0.5*DIM2, 0.5*DIM4]]

        shift = [-0.5*(DIM1+DIM2)+shift[0], -0.5*DIM4+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        t = 0.5*(DIM4 - DIM3)
        self.points = [[DIM1,      0.],    [2.*DIM1+DIM2, 0.],    [2.*DIM1+DIM2, t],       [DIM1+DIM2,     t],
                       [DIM1+DIM2, DIM4], [0.,           DIM4], [0.,           DIM4-t],  [DIM1,        DIM4-t]]
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