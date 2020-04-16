from sectionproperties.pre.sections import (
    Geometry, RectangularSection, CustomSection, MergedSection
)
from sectionproperties.pre.pre import create_mesh
import numpy as np


class BARSection(Geometry):
    """Constructs a BAR section with the center at the origin *(0, 0)*, with two parameters
    defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ [5]_ for definition of
    parameters. Added by JohnDN90.

    :param float DIM1: Width (x) of bar
    :param float DIM2: Depth (y) of bar
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a BAR cross-section with a depth of 1.5 and width of 2.0, and
    generates a mesh with a maximum triangular area of 0.001::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.BARSection(DIM1=2.0, DIM2=1.5)
        mesh = geometry.create_mesh(mesh_sizes=[0.001])

    ..  figure:: ../images/sections/bar_geometry.png
        :align: center
        :scale: 75 %

        BAR section geometry.

    ..  figure:: ../images/sections/bar_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, shift=[0, 0]):
        """Inits the BARSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2

        # assign control point
        control_points = [[0., 0.]]

        # shift = [-0.5*DIM1+shift[0], -0.5*DIM2+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        self.points = [
            [-0.5*DIM1, -0.5*DIM2], [0.5*DIM1, -0.5*DIM2],
            [0.5*DIM1, 0.5*DIM2], [-0.5*DIM1, 0.5*DIM2]
        ]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 0]]

        self.shift_section()

    def getStressPoints(self, shift=(0., 0.)):
        """Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (0.5*self.DIM1-shift[0], 0.5*self.DIM2-shift[1])
        D = (0.5*self.DIM1-shift[0], -0.5*self.DIM2-shift[1])
        E = (-0.5*self.DIM1-shift[0], -0.5*self.DIM2-shift[1])
        F = (-0.5*self.DIM1-shift[0], 0.5*self.DIM2-shift[1])

        return C, D, E, F


class BOXSection(Geometry):
    """ Constructs a BOX section with the center at the origin *(0, 0)*, with four parameters
    defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ [5]_ for definition of
    parameters. Added by JohnDN90.

    :param float DIM1: Width (x) of box
    :param float DIM2: Depth (y) of box
    :param float DIM3: Thickness of box in y direction
    :param float DIM4: Thickness of box in x direction
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a BOX cross-section with a depth of 3.0 and width of 4.0, and
    generates a mesh with a maximum triangular area of 0.001::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.BOXSection(DIM1=4.0, DIM2=3.0, DIM3=0.375, DIM4=0.5)
        mesh = geometry.create_mesh(mesh_sizes=[0.001])

    ..  figure:: ../images/sections/box_geometry.png
        :align: center
        :scale: 75 %

        BOX section geometry.

    ..  figure:: ../images/sections/box_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, DIM3, DIM4, shift=[0, 0]):
        """Inits the BOXSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3
        self.DIM4 = DIM4

        # Ensure dimensions are physically relevant
        np.testing.assert_(2.0*DIM4 < DIM1, "Invalid geometry specified.")
        np.testing.assert_(2.0*DIM3 < DIM2, "Invalid geometry specified.")

        # assign control point
        control_points = [[0., 0.5*DIM2 - 0.5*DIM3]]

        super().__init__(control_points, shift)

        # specify a hole in the centre of the Box
        self.holes = [[0., 0.]]

        # construct the points and facets
        self.points = [
            [-0.5*DIM1, -0.5*DIM2], [0.5*DIM1, -0.5*DIM2], [0.5*DIM1, 0.5*DIM2],
            [-0.5*DIM1, 0.5*DIM2], [-0.5*DIM1 + DIM4, -0.5*DIM2 + DIM3],
            [0.5*DIM1 - DIM4, -0.5*DIM2 + DIM3], [0.5*DIM1 - DIM4, 0.5*DIM2 - DIM3],
            [-0.5*DIM1 + DIM4, 0.5*DIM2 - DIM3]
        ]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4]]

        self.shift_section()

    def getStressPoints(self, shift=(0., 0.)):
        """Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (0.5*self.DIM1-shift[0], 0.5*self.DIM2-shift[1])
        D = (0.5*self.DIM1-shift[0], -0.5*self.DIM2-shift[1])
        E = (-0.5*self.DIM1-shift[0], -0.5*self.DIM2-shift[1])
        F = (-0.5*self.DIM1-shift[0], 0.5*self.DIM2-shift[1])

        return C, D, E, F


class BOX1Section(Geometry):
    """Constructs a BOX1 section with the center at the origin *(0, 0)*, with six parameters
    defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for more details. Added by
    JohnDN90.

    :param float DIM1: Width (x) of box
    :param float DIM2: Depth (y) of box
    :param float DIM3: Thickness of top wall
    :param float DIM4: Thickness of bottom wall
    :param float DIM5: Thickness of left wall
    :param float DIM6: Thickness of right wall
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a BOX1 cross-section with a depth of 3.0 and width of 4.0, and
    generates a mesh with a maximum triangular area of 0.007::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.BOX1Section(
            DIM1=4.0, DIM2=3.0, DIM3=0.375, DIM4=0.5, DIM5=0.25, DIM6=0.75
        )
        mesh = geometry.create_mesh(mesh_sizes=[0.007])

    ..  figure:: ../images/sections/box1_geometry.png
        :align: center
        :scale: 75 %

        BOX1 section geometry.

    ..  figure:: ../images/sections/box1_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, shift=[0, 0]):
        """Inits the Box1Section class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        DIM5 *= 1.0
        DIM6 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3
        self.DIM4 = DIM4
        self.DIM5 = DIM5
        self.DIM6 = DIM6

        # Ensure dimensions are physically relevant
        np.testing.assert_(DIM5+DIM6 < DIM1, "Invalid geometry specified.")
        np.testing.assert_(DIM3+DIM4 < DIM2, "Invalid geometry specified.")

        # assign control point
        control_points = [[0.5*DIM1, 0.5*DIM4]]

        shift = [-0.5*DIM1+shift[0], -0.5*DIM2+shift[1]]
        super().__init__(control_points, shift)

        # specify a hole in the centre of the Box
        self.holes = [[DIM6 + 0.5*(DIM1-DIM5), DIM4+0.5*(DIM2-DIM3)]]

        # construct the points and facets
        self.points = [
            [0., 0.], [DIM1, 0.], [DIM1, DIM2], [0., DIM2], [DIM6, DIM4], [DIM1-DIM5, DIM4],
            [DIM1-DIM5, DIM2-DIM3], [DIM6, DIM2-DIM3]
        ]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4]]

        self.shift_section()

    def getStressPoints(self, shift=(0., 0.)):
        """ Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (0.5*self.DIM1-shift[0], 0.5*self.DIM2-shift[1])
        D = (0.5*self.DIM1-shift[0], -0.5*self.DIM2-shift[1])
        E = (-0.5*self.DIM1-shift[0], -0.5*self.DIM2-shift[1])
        F = (-0.5*self.DIM1-shift[0], 0.5*self.DIM2-shift[1])

        return C, D, E, F


class CHANSection(Geometry):
    """ Constructs a CHAN (C-Channel) section with the web's middle center at the origin *(0, 0)*,
    with four parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for
    more details. Added by JohnDN90.

    :param float DIM1: Width (x) of the CHAN-section
    :param float DIM2: Depth (y) of the CHAN-section
    :param float DIM3: Thickness of web (vertical portion)
    :param float DIM4: Thickness of flanges (top/bottom portion)
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a CHAN cross-section with a depth of 4.0 and width of 2.0, and
    generates a mesh with a maximum triangular area of 0.008::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.CHANSection(DIM1=2.0, DIM2=4.0, DIM3=0.25, DIM4=0.5)
        mesh = geometry.create_mesh(mesh_sizes=[0.008])

    ..  figure:: ../images/sections/chan_geometry.png
        :align: center
        :scale: 75 %

        CHAN section geometry.

    ..  figure:: ../images/sections/chan_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, DIM3, DIM4, shift=[0, 0]):
        """Inits the CHANSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3
        self.DIM4 = DIM4

        # Ensure dimensions are physically relevant
        np.testing.assert_(2.0*DIM4 < DIM2, "Invalid geometry specified.")
        np.testing.assert_(DIM3 < DIM1, "Invalid geometry specified.")

        # assign control point
        control_points = [[0.5*DIM1, 0.5*DIM4]]

        shift = [-0.5*DIM3+shift[0], -0.5*DIM2+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        self.points = [
            [0., 0.], [DIM1, 0.], [DIM1, DIM4], [DIM3, DIM4], [DIM3, DIM2-DIM4], [DIM1, DIM2-DIM4],
            [DIM1, DIM2], [0., DIM2]
        ]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 0]]

        self.shift_section()

    def getStressPoints(self, shift=(0., 0.)):
        """ Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (self.DIM1-0.5*self.DIM3-shift[0], 0.5*self.DIM2-shift[1])
        D = (self.DIM1-0.5*self.DIM3-shift[0], -0.5*self.DIM2-shift[1])
        E = (-0.5*self.DIM3-shift[0], -0.5*self.DIM2-shift[1])
        F = (-0.5*self.DIM3-shift[0], 0.5*self.DIM2-shift[1])

        return C, D, E, F


class CHAN1Section(Geometry):
    """ Constructs a CHAN1 (C-Channel) section with the web's middle center at the origin *(0, 0)*,
    with four parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for
    more details. Added by JohnDN90.

    :param float DIM1: Width (x) of channels
    :param float DIM2: Thicknesss (x) of web
    :param float DIM3: Spacing between channels (length of web)
    :param float DIM4: Depth (y) of CHAN1-section
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a CHAN1 cross-section with a depth of 4.0 and width of 1.75, and
    generates a mesh with a maximum triangular area of 0.01::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.CHAN1Section(DIM1=0.75, DIM2=1.0, DIM3=3.5, DIM4=4.0)
        mesh = geometry.create_mesh(mesh_sizes=[0.01])

    ..  figure:: ../images/sections/chan1_geometry.png
        :align: center
        :scale: 75 %

        CHAN1 section geometry.

    ..  figure:: ../images/sections/chan1_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, DIM3, DIM4, shift=[0, 0]):
        """Inits the CHAN1Section class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3
        self.DIM4 = DIM4

        # Ensure dimensions are physically relevant
        np.testing.assert_(DIM4 > DIM3, "Invalid geometry specified.")

        # assign control point
        control_points = [[0.5*DIM1, 0.5*DIM4]]

        shift = [-0.5*DIM2+shift[0], -0.5*DIM4+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        tf = 0.5 * (DIM4 - DIM3)
        self.points = [
            [0, 0], [DIM1+DIM2, 0], [DIM1+DIM2, tf], [DIM2, tf], [DIM2, tf+DIM3],
            [DIM2+DIM1, tf+DIM3], [DIM2+DIM1, DIM4], [0, DIM4]
        ]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 0]]

        self.shift_section()

    def getStressPoints(self, shift=(0., 0.)):
        """ Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (0.5*self.DIM2+self.DIM1-shift[0], 0.5*self.DIM4-shift[1])
        D = (0.5*self.DIM2+self.DIM1-shift[0], -0.5*self.DIM4-shift[1])
        E = (-0.5*self.DIM2-shift[0], -0.5*self.DIM4-shift[1])
        F = (-0.5*self.DIM2-shift[0], 0.5*self.DIM4-shift[1])

        return C, D, E, F


class CHAN2Section(Geometry):
    """ Constructs a CHAN2 (C-Channel) section with the bottom web's middle center at the origin
    *(0, 0)*, with four parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_
    [4]_ for more details. Added by JohnDN90.

    :param float DIM1: Thickness of channels
    :param float DIM2: Thickness of web
    :param float DIM3: Depth (y) of CHAN2-section
    :param float DIM4: Width (x) of CHAN2-section
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a CHAN2 cross-section with a depth of 2.0 and width of 4.0, and
    generates a mesh with a maximum triangular area of 0.01::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.CHAN2Section(DIM1=0.375, DIM2=0.5, DIM3=2.0, DIM4=4.0)
        mesh = geometry.create_mesh(mesh_sizes=[0.01])

    ..  figure:: ../images/sections/chan2_geometry.png
        :align: center
        :scale: 75 %

        CHAN2 section geometry.

    ..  figure:: ../images/sections/chan2_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, DIM3, DIM4, shift=[0, 0]):
        """Inits the CHAN2Section class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3
        self.DIM4 = DIM4

        # Ensure dimensions are physically relevant
        np.testing.assert_(DIM4 > 2.0*DIM1, "Invalid geometry specified.")
        np.testing.assert_(DIM3 > DIM2, "Invalid geometry specified.")

        # assign control point
        control_points = [[0.5*DIM4, 0.5*DIM2]]

        shift = [-0.5*DIM4+shift[0], -0.5*DIM2+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        self.points = [
            [0., 0.], [DIM4, 0.], [DIM4, DIM3], [DIM4-DIM1, DIM3], [DIM4-DIM1, DIM2], [DIM1, DIM2],
            [DIM1, DIM3], [0., DIM3]
        ]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 0]]

        self.shift_section()

    def getStressPoints(self, shift=(0., 0.)):
        """ Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (0.5*self.DIM4-shift[0], self.DIM3-0.5*self.DIM2-shift[1])
        D = (0.5*self.DIM4-shift[0], -0.5*self.DIM2-shift[1])
        E = (-0.5*self.DIM4-shift[0], -0.5*self.DIM2-shift[1])
        F = (-0.5*self.DIM4-shift[0], self.DIM3-0.5*self.DIM2-shift[1])

        return C, D, E, F


class CROSSSection(Geometry):
    """ Constructs Nastran's cruciform/cross section with the intersection's middle center at the
    origin *(0, 0)*, with four parameters defining dimensions. See Nastran documentation [1]_ [2]_
    [3]_ [4]_ for more details. Added by JohnDN90.

    :param float DIM1: Twice the width of horizontal member protruding from the vertical center
        member
    :param float DIM2: Thickness of the vertical member
    :param float DIM3: Depth (y) of the CROSS-section
    :param float DIM4: Thickness of the horizontal members
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a rectangular cross-section with a depth of 3.0 and width of
    1.875, and generates a mesh with a maximum triangular area of 0.008::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.CROSSSection(DIM1=1.5, DIM2=0.375, DIM3=3.0, DIM4=0.25)
        mesh = geometry.create_mesh(mesh_sizes=[0.008])

    ..  figure:: ../images/sections/cross_geometry.png
        :align: center
        :scale: 75 %

        Cruciform/cross section geometry.

    ..  figure:: ../images/sections/cross_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, DIM3, DIM4, shift=[0, 0]):
        """Inits the CROSSSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3
        self.DIM4 = DIM4

        # Ensure dimensions are physically relevant
        np.testing.assert_(DIM4 < DIM3, "Invalid geometry specified.")

        # assign control point
        control_points = [[0.5*DIM1+0.5*DIM2, 0.5*DIM3]]

        shift = [-(0.5*DIM1+0.5*DIM2)+shift[0], -(0.5*DIM3)+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        d = 0.5*(DIM3 - DIM4)
        self.points = [
            [0.5*DIM1, 0], [0.5*DIM1+DIM2, 0], [0.5*DIM1+DIM2, d], [DIM1+DIM2, d],
            [DIM1+DIM2, d+DIM4], [0.5*DIM1+DIM2, d+DIM4], [0.5*DIM1+DIM2, DIM3], [0.5*DIM1, DIM3],
            [0.5*DIM1, d+DIM4], [0, d+DIM4], [0, d], [0.5*DIM1, d]
        ]
        self.facets = [
            [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
            [10, 11], [11, 0]
        ]

        self.shift_section()

    def getStressPoints(self, shift=(0., 0.)):
        """Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (-shift[0], 0.5*self.DIM3-shift[1])
        D = (0.5*(self.DIM1+self.DIM2)-shift[0], -shift[1])
        E = (-shift[0], -0.5*self.DIM3-shift[1])
        F = (-0.5*(self.DIM1+self.DIM2)-shift[0], -shift[1])

        return C, D, E, F


class FCROSSSection(Geometry):
    """ Constructs a flanged cruciform/cross section with the intersection's middle center at the
    origin *(0, 0)*, with eight parameters defining dimensions. Added by JohnDN90.

    :param float DIM1: Depth (y) of flanged cruciform
    :param float DIM2: Width (x) of flanged cruciform
    :param float DIM3: Thickness of vertical web
    :param float DIM4: Thickness of horizontal web
    :param float DIM5: Length of flange attached to vertical web
    :param float DIM6: Thickness of flange attached to vertical web
    :param float DIM7: Length of flange attached to horizontal web
    :param float DIM8: Thickness of flange attached to horizontal web
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example demonstrates the creation of a flanged cross section::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.FCROSSSection(
            DIM1=9.0, DIM2=6.0, DIM3=0.75, DIM4=0.625, DIM5=2.1, DIM6=0.375, DIM7=4.5, DIM8=0.564
        )
        mesh = geometry.create_mesh(mesh_sizes=[0.03])

    ..  figure:: ../images/sections/fcross_geometry.png
        :align: center
        :scale: 75 %

        Flanged Cruciform/cross section geometry.

    ..  figure:: ../images/sections/fcross_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, shift=[0, 0]):
        """Inits the FCROSSSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        DIM5 *= 1.0
        DIM6 *= 1.0
        DIM7 *= 1.0
        DIM8 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3
        self.DIM4 = DIM4
        self.DIM5 = DIM5
        self.DIM6 = DIM6
        self.DIM7 = DIM7
        self.DIM7 = DIM8

        # Ensure dimensions are physically relevant
        # TODO: Finish dimension checks.
        np.testing.assert_(DIM5 > DIM3, "Invalid geometry specified.")
        np.testing.assert_(DIM7 > DIM4, "Invalid geometry specified.")
        np.testing.assert_(DIM7 < DIM1, "Invalid geometry specified.")
        np.testing.assert_(DIM5 < DIM2, "Invalid geometry specified.")
        np.testing.assert_(DIM8 < (0.5*DIM2-0.5*DIM3), "Invalid geometry specified.")
        np.testing.assert_(DIM6 < (0.5*DIM1-0.5*DIM4), "Invalid geometry specified.")

        # assign control point
        control_points = [[0.0, 0.0]]

        shift = [shift[0], shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        self.points = [
            [0.5*DIM3, -0.5*DIM4], [0.5*DIM2-DIM8, -0.5*DIM4], [0.5*DIM2-DIM8, -0.5*DIM7],
            [0.5*DIM2, -0.5*DIM7], [0.5*DIM2, 0.5*DIM7], [0.5*DIM2-DIM8, 0.5*DIM7],
            [0.5*DIM2-DIM8, 0.5*DIM4], [0.5*DIM3, 0.5*DIM4], [0.5*DIM3, 0.5*DIM1-DIM6],
            [0.5*DIM5, 0.5*DIM1-DIM6], [0.5*DIM5, 0.5*DIM1], [-0.5*DIM5, 0.5*DIM1],
            [-0.5*DIM5, 0.5*DIM1-DIM6],  [-0.5*DIM3, 0.5*DIM1-DIM6], [-0.5*DIM3, 0.5*DIM4],
            [-0.5*DIM2+DIM8, 0.5*DIM4],  [-0.5*DIM2+DIM8, 0.5*DIM7], [-0.5*DIM2, 0.5*DIM7],
            [-0.5*DIM2, -0.5*DIM7], [-0.5*DIM2+DIM8, -0.5*DIM7], [-0.5*DIM2+DIM8, -0.5*DIM4],
            [-0.5*DIM3, -0.5*DIM4], [-0.5*DIM3, -0.5*DIM1+DIM6], [-0.5*DIM5, -0.5*DIM1+DIM6],
            [-0.5*DIM5, -0.5*DIM1], [0.5*DIM5, -0.5*DIM1], [0.5*DIM5, -0.5*DIM1+DIM6],
            [0.5*DIM3, -0.5*DIM1+DIM6]
        ]
        self.facets = [
            [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
            [10, 11], [11, 12], [12, 13], [13, 14], [14, 15], [15, 16], [16, 17], [17, 18],
            [18, 19], [19, 20], [20, 21], [21, 22], [22, 23], [23, 24], [24, 25], [25, 26],
            [26, 27], [27, 0]
        ]

        self.shift_section()

    def getStressPoints(self, shift=(0., 0.)):
        """Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: list[float, float]
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (-shift[0], 0.5*self.DIM1-shift[1])
        D = (0.5*self.DIM2-shift[0], -shift[1])
        E = (-shift[0], -0.5*self.DIM1-shift[1])
        F = (-0.5*self.DIM2-shift[0], -shift[1])

        return C, D, E, F


class DBOXSection(Geometry):
    """ Constructs a DBOX section with the center at the origin *(0, 0)*, with ten parameters
    defining dimensions. See MSC Nastran documentation [1]_ for more details. Added by JohnDN90.

    :param float DIM1: Width (x) of the DBOX-section
    :param float DIM2: Depth (y) of the DBOX-section
    :param float DIM3: Width (x) of left-side box
    :param float DIM4: Thickness of left wall
    :param float DIM5: Thickness of center wall
    :param float DIM6: Thickness of right wall
    :param float DIM7: Thickness of top left wall
    :param float DIM8: Thickness of bottom left wall
    :param float DIM9: Thickness of top right wall
    :param float DIM10: Thickness of bottom right wall
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a DBOX cross-section with a depth of 3.0 and width of 8.0, and
    generates a mesh with a maximum triangular area of 0.01::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.DBOXSection(
            DIM1=8.0, DIM2=3.0, DIM3=3.0, DIM4=0.5, DIM5=0.625, DIM6=0.75, DIM7=0.375, DIM8=0.25,
            DIM9=0.5, DIM10=0.375
        )
        mesh = geometry.create_mesh(mesh_sizes=[0.01])

    ..  figure:: ../images/sections/dbox_geometry.png
        :align: center
        :scale: 75 %

        DBOX section geometry.

    ..  figure:: ../images/sections/dbox_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, DIM7, DIM8, DIM9, DIM10, shift=[0, 0]):
        """Inits the DBOXSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        DIM5 *= 1.0
        DIM6 *= 1.0
        DIM7 *= 1.0
        DIM8 *= 1.0
        DIM9 *= 1.0
        DIM10 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3
        self.DIM4 = DIM4
        self.DIM5 = DIM5
        self.DIM6 = DIM6
        self.DIM7 = DIM7
        self.DIM8 = DIM8
        self.DIM9 = DIM9
        self.DIM10 = DIM10

        # Ensure dimensions are physically relevant
        np.testing.assert_((DIM4+DIM5+DIM6) < DIM1, "Invalid geometry specified.")
        np.testing.assert_((DIM4+0.5*DIM5) < DIM3, "Invalid geometry specified.")
        np.testing.assert_((DIM7+DIM8) < DIM2, "Invalid geometry specified.")
        np.testing.assert_((DIM9+DIM10) < DIM2, "Invalid geometry specified.")

        # assign control point
        control_points = [[0.5*DIM3, 0.5*DIM8]]

        shift = [-0.5*DIM1+shift[0], -0.5*DIM2+shift[1]]
        super().__init__(control_points, shift)

        # specify a hole in the centre of the Box
        d2 = 0.5*(DIM1 - DIM6 - DIM3 - 0.5*DIM5)
        self.holes = [
            [DIM4 + 0.5*(DIM3 - DIM4 - 0.5*DIM5), DIM8 + 0.5*(DIM2 - DIM8 - DIM7)],
            [DIM3 + 0.5*DIM5 + d2, DIM10 + 0.5*(DIM2 - DIM10 - DIM9)]
        ]

        # construct the points and facets
        self.points = [
            [0., 0.], [DIM1, 0.], [DIM1, DIM2], [0., DIM2], [DIM4, DIM8], [DIM3-DIM5/2., DIM8],
            [DIM3-DIM5/2., DIM2-DIM7], [DIM4, DIM2-DIM7], [DIM3+DIM5/2., DIM10],
            [DIM1-DIM6, DIM10], [DIM1-DIM6, DIM2-DIM9], [DIM3+DIM5/2., DIM2-DIM9]
        ]
        self.facets = [
            [0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4], [8, 9], [9, 10],
            [10, 11], [11, 8]
        ]

        self.shift_section()

    def getStressPoints(self, shift=(0., 0.)):
        """
        Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (0.5 * self.DIM1 - shift[0], 0.5 * self.DIM2 - shift[1])
        D = (0.5 * self.DIM1 - shift[0], -0.5 * self.DIM2 - shift[1])
        E = (-0.5 * self.DIM1 - shift[0], -0.5 * self.DIM2 - shift[1])
        F = (-0.5 * self.DIM1 - shift[0], 0.5 * self.DIM2 - shift[1])

        return C, D, E, F


class GBOXSection(Geometry):
    """ Constructs a GBOX section with the center at the origin *(0, 0)*, with six parameters
    defining dimensions. See ASTROS documentation [5]_ for more details. Added by JohnDN90.

    :param float DIM1: Width (x) of the GBOX-section
    :param float DIM2: Depth (y) of the GBOX-section
    :param float DIM3: Thickness of top flange
    :param float DIM4: Thickness of bottom flange
    :param float DIM5: Thickness of webs
    :param float DIM6: Spacing between webs
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a GBOX cross-section with a depth of 2.5 and width of 6.0, and
    generates a mesh with a maximum triangular area of 0.01::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.GBOXSection(
            DIM1=6.0, DIM2=2.5, DIM3=0.375, DIM4=0.25, DIM5=0.625, DIM6=1.0
        )
        mesh = geometry.create_mesh(mesh_sizes=[0.01])

    ..  figure:: ../images/sections/gbox_geometry.png
        :align: center
        :scale: 75 %

        GBOX section geometry.

    ..  figure:: ../images/sections/gbox_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, shift=[0, 0]):
        """Inits the GBOXSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        DIM5 *= 1.0
        DIM6 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3
        self.DIM4 = DIM4
        self.DIM5 = DIM5
        self.DIM6 = DIM6

        # Ensure dimensions are physically relevant
        np.testing.assert_((DIM3+DIM4) < DIM2, "Invalid geometry specified.")
        np.testing.assert_((2.0*DIM5+DIM6) < DIM1, "Invalid geometry specified.")

        # assign control point
        control_points = [[0.5*DIM1, 0.5*DIM4]]

        shift = [-(0.5*DIM1)+shift[0], -(DIM4 + 0.5*(DIM2-DIM3-DIM4))+shift[1]]
        super().__init__(control_points, shift)

        # specify a hole in the centre of the GBOX
        self.holes = [[0.5*DIM1, 0.5*DIM2]]

        # construct the points and facets
        d = 0.5*(DIM1 - DIM6 - 2.0 * DIM5)
        self.points = [
            [0., 0.], [DIM1, 0.], [DIM1, DIM4], [d + 2. * DIM5 + DIM6, DIM4],
            [d + 2. * DIM5 + DIM6, DIM2 - DIM3], [DIM1, DIM2 - DIM3], [DIM1, DIM2], [0., DIM2],
            [0., DIM2 - DIM3], [d, DIM2 - DIM3], [d, DIM4], [0., DIM4], [d + DIM5, DIM4],
            [d + DIM5 + DIM6, DIM4], [d + DIM5 + DIM6, DIM2 - DIM3], [d + DIM5, DIM2 - DIM3]
        ]
        self.facets = [
            [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
            [10, 11], [11, 0], [12, 13], [13, 14], [14, 15], [15, 12]
        ]

        self.shift_section()

    def getStressPoints(self, shift=(0., 0.)):
        """Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (0.5*self.DIM1-shift[0], 0.5*self.DIM2-shift[1])
        D = (0.5*self.DIM1-shift[0], -0.5*self.DIM2-shift[1])
        E = (-0.5*self.DIM1-shift[0], -0.5*self.DIM2-shift[1])
        F = (-0.5*self.DIM1-shift[0], 0.5*self.DIM2-shift[1])

        return C, D, E, F


class HSection(Geometry):
    """Constructs a H section with the middle web's middle center at the origin *(0, 0)*, with four
    parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for more details.
    Added by JohnDN90.

    :param float DIM1: Spacing between vertical flanges (length of web)
    :param float DIM2: Twice the thickness of the vertical flanges
    :param float DIM3: Depth (y) of the H-section
    :param float DIM4: Thickness of the middle web
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a H cross-section with a depth of 3.5 and width of 2.75, and
    generates a mesh with a maximum triangular area of 0.005::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.HSection(DIM1=2.0, DIM2=0.75, DIM3=3.5, DIM4=0.25)
        mesh = geometry.create_mesh(mesh_sizes=[0.005])

    ..  figure:: ../images/sections/h_geometry.png
        :align: center
        :scale: 75 %

        H section geometry.

    ..  figure:: ../images/sections/h_mesh.png
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
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3
        self.DIM4 = DIM4

        # Ensure dimensions are physically relevant
        np.testing.assert_(DIM4 < DIM3, "Invalid geometry specified.")

        d1 = 0.5 * (DIM3 - DIM4)
        d2 = 0.5 * DIM2

        # assign control point
        control_points = [[0.5*d2, 0.5*DIM3]]

        shift = [-0.5*(DIM2+DIM1)+shift[0], -0.5*DIM3+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        self.points = [
            [0, 0], [d2, 0], [d2, d1], [d2+DIM1, d1], [d2+DIM1, 0], [DIM1+DIM2, 0],
            [DIM1+DIM2, DIM3], [DIM1+DIM2-d2, DIM3], [DIM1+DIM2-d2, d1+DIM4], [d2, d1+DIM4],
            [d2, DIM3], [0, DIM3]
        ]
        self.facets = [
            [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
            [10, 11], [11, 0]
        ]

        self.shift_section()

    def getStressPoints(self, shift=(0., 0.)):
        """Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (0.5*(self.DIM1+self.DIM2)-shift[0], 0.5*self.DIM3-shift[1])
        D = (0.5*(self.DIM1+self.DIM2)-shift[0], -0.5*self.DIM3-shift[1])
        E = (-0.5*(self.DIM1+self.DIM2)-shift[0], -0.5*self.DIM3-shift[1])
        F = (-0.5*(self.DIM1+self.DIM2)-shift[0], 0.5*self.DIM3-shift[1])

        return C, D, E, F


class HATSection(Geometry):
    """Constructs a Hat section with the top most section's middle center at the origin *(0, 0)*,
    with four parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for
    more details. Note that HAT in ASTROS is actually HAT1 in this code. Added by JohnDN90.

    :param float DIM1: Depth (y) of HAT-section
    :param float DIM2: Thickness of HAT-section
    :param float DIM3: Width (x) of top most section
    :param float DIM4: Width (x) of bottom sections
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a HAT cross-section with a depth of 1.25 and width of 2.5, and
    generates a mesh with a maximum triangular area of 0.001::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.HATSection(DIM1=1.25, DIM2=0.25, DIM3=1.5, DIM4=0.5)
        mesh = geometry.create_mesh(mesh_sizes=[0.001])

    ..  figure:: ../images/sections/hat_geometry.png
        :align: center
        :scale: 75 %

        HAT section geometry.

    ..  figure:: ../images/sections/hat_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, DIM3, DIM4, shift=[0, 0]):
        """Inits the HATSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3
        self.DIM4 = DIM4

        # Ensure dimensions are physically relevant
        np.testing.assert_(2.0*DIM2 < DIM1, "Invalid geometry specified.")

        # assign control point
        control_points = [[0.5*DIM4, 0.5*DIM2]]

        shift = [-DIM4-0.5*DIM3+shift[0], -DIM1+0.5*DIM2+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        self.points = [
            [0., 0.], [DIM4+DIM2, 0.], [DIM4+DIM2, DIM1-DIM2], [DIM4+DIM3-DIM2, DIM1-DIM2],
            [DIM4+DIM3-DIM2, 0.], [2*DIM4+DIM3, 0.], [2.*DIM4+DIM3, DIM2], [DIM4+DIM3, DIM2],
            [DIM4+DIM3, DIM1], [DIM4, DIM1], [DIM4, DIM2], [0., DIM2]
        ]
        self.facets = [
            [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
            [10, 11], [11, 0]
        ]

        self.shift_section()

    def getStressPoints(self, shift=(0., 0.)):
        """Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the origin by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (0.5*self.DIM3 - shift[0], 0.5*self.DIM2 - shift[1])
        D = (0.5*self.DIM3 + self.DIM4 - shift[0], -self.DIM1 + self.DIM2 - shift[1])
        E = (-0.5*self.DIM3 - self.DIM4 - shift[0], -self.DIM1 + self.DIM2 - shift[1])
        F = (-0.5*self.DIM3 - shift[0], 0.5*self.DIM2 - shift[1])

        return C, D, E, F


class HAT1Section(Geometry):
    """ Constructs a HAT1 section with the bottom plate's bottom center at the origin *(0, 0)*,
    with five parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [5]_ for
    definition of parameters. Note that in ASTROS, HAT1 is called HAT. Added by JohnDN90.

    :param float DIM1: Width(x) of the HAT1-section
    :param float DIM2: Depth (y) of the HAT1-section
    :param float DIM3: Width (x) of hat's top flange
    :param float DIM4: Thickness of hat stiffener
    :param float DIM5: Thicknesss of bottom plate
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a HAT1 cross-section with a depth of 2.0 and width of 4.0, and
    generates a mesh with a maximum triangular area of 0.005::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.HAT1Section(DIM1=4.0, DIM2=2.0, DIM3=1.5, DIM4=0.1875, DIM5=0.375)
        mesh = geometry.create_mesh(mesh_sizes=[0.005])

    ..  figure:: ../images/sections/hat1_geometry.png
        :align: center
        :scale: 75 %

        HAT1 section geometry.

    ..  figure:: ../images/sections/hat1_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, DIM3, DIM4, DIM5, shift=[0, 0]):
        """Inits the HAT1Section class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        DIM5 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3
        self.DIM4 = DIM4
        self.DIM5 = DIM5

        # Ensure dimensions are physically relevant
        np.testing.assert_((2.0*DIM4+DIM5) < DIM2, "Invalid geometry specified.")
        np.testing.assert_(DIM3 < DIM1, "Invalid geometry specified.")

        shift = [-0.5*DIM1+shift[0], shift[1]]

        # create bottom rectangular plate
        bottom_plate = RectangularSection(d=DIM5, b=DIM1, shift=shift)

        # create the hat stiffener
        d1 = DIM2 - DIM5
        d2 = DIM4
        d4 = 0.5*(DIM1 - DIM3)

        # specify a hole in the combined plate and hat structure
        holes = [[0.5*DIM1, 0.5*DIM2]]

        # assign control point
        control_points = [[0.5*d4, DIM5 + 0.5*DIM4]]

        super().__init__(control_points, shift)

        # construct the points and facets
        points = [
            [0., DIM5 + 0.], [d4 + d2, DIM5 + 0.], [d4 + d2, DIM5 + d1 - d2],
            [d4 + DIM3 - d2, DIM5 + d1 - d2], [d4 + DIM3 - d2, DIM5 + 0.],
            [2. * d4 + DIM3, DIM5 + 0.], [2. * d4 + DIM3, DIM5 + d2], [d4 + DIM3, DIM5 + d2],
            [d4 + DIM3, DIM5 + d1], [d4, DIM5 + d1], [d4, DIM5 + d2], [0, DIM5 + d2]
        ]
        facets = [
            [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
            [10, 11], [11, 0]
        ]

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
        """Creates a quadratic triangular mesh from the Geometry object. This is overloaded here to
        allow specifying only one mesh_size which is used for both regions in the Hat1 section.

        :param mesh_sizes: A list of maximum element areas corresponding to each region within the
            cross-section geometry.
        :type mesh_size: list[float]

        :return: Object containing generated mesh data
        :rtype: :class:`meshpy.triangle.MeshInfo`

        :raises AssertionError: If the number of mesh sizes does not match the number of regions
        """

        mesh_sizes *= 2

        str = "Number of mesh_sizes ({0}), should match the number of regions ({1})".format(
            len(mesh_sizes), len(self.control_points)
        )
        assert(len(mesh_sizes) == len(self.control_points)), str

        return create_mesh(self.points, self.facets, self.holes, self.control_points, mesh_sizes)

    def getStressPoints(self, shift=(0., 0.)):
        """ Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the origin by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (-0.5*self.DIM1 - shift[0], -shift[1])
        D = (0.5*self.DIM1 - shift[0],  -shift[1])
        E = (-0.5*self.DIM3 - shift[0], self.DIM2 - shift[1])
        F = (0.5*self.DIM3 - shift[0], self.DIM2 - shift[1])

        return C, D, E, F


class HEXASection(Geometry):
    """ Constructs a HEXA (hexagon) section with the center at the origin *(0, 0)*, with three
    parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for more details.
    Added by JohnDN90.

    :param float DIM1: Spacing between bottom right point and right most point
    :param float DIM2: Width (x) of hexagon
    :param float DIM3: Depth (y) of hexagon
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a rectangular cross-section with a depth of 1.5 and width of 2.0,
    and generates a mesh with a maximum triangular area of 0.005::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.HEXASection(DIM1=0.5, DIM2=2.0, DIM3=1.5)
        mesh = geometry.create_mesh(mesh_sizes=[0.005])

    ..  figure:: ../images/sections/hexa_geometry.png
        :align: center
        :scale: 75 %

        HEXA section geometry.

    ..  figure:: ../images/sections/hexa_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, DIM3, shift=[0, 0]):
        """Inits the HEXASection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3

        # Ensure dimensions are physically relevant
        np.testing.assert_(DIM2 > DIM1, "Invalid geometry specified.")

        # assign control point
        control_points = [[0.5*DIM2, 0.5*DIM3]]

        shift = [-0.5*DIM2+shift[0], -0.5*DIM3+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        self.points = [
            [DIM1, 0.], [DIM2-DIM1, 0.], [DIM2, 0.5*DIM3], [DIM2-DIM1, DIM3], [DIM1, DIM3],
            [0., 0.5*DIM3]
        ]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]

        self.shift_section()

    def getStressPoints(self, shift=(0., 0.)):
        """Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (-shift[0], 0.5*self.DIM3-shift[1])
        D = (-shift[0], -0.5*self.DIM3-shift[1])
        E = (0.5*self.DIM2-shift[0], -shift[1])
        F = (-0.5*self.DIM2-shift[0], -shift[1])

        return C, D, E, F


class NISection(Geometry):
    """Constructs Nastran's I section with the bottom flange's middle center at the origin
    *(0, 0)*, with six parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_
    [4]_ for definition of parameters. Added by JohnDN90.

    :param float DIM1: Depth(y) of the I-section
    :param float DIM2: Width (x) of bottom flange
    :param float DIM3: Width (x) of top flange
    :param float DIM4: Thickness of web
    :param float DIM5: Thickness of bottom web
    :param float DIM6: Thickness of top web
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a Nastran I cross-section with a depth of 5.0, and generates a
    mesh with a maximum triangular area of 0.008::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.NISection(
            DIM1=5.0, DIM2=2.0, DIM3=3.0, DIM4=0.25, DIM5=0.375, DIM6=0.5
        )
        mesh = geometry.create_mesh(mesh_sizes=[0.008])

    ..  figure:: ../images/sections/ni_geometry.png
        :align: center
        :scale: 75 %

        Nastran's I section geometry.

    ..  figure:: ../images/sections/ni_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, DIM3, DIM4, DIM5, DIM6, shift=[0, 0]):
        """Inits the NISection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        DIM5 *= 1.0
        DIM6 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3
        self.DIM4 = DIM4
        self.DIM5 = DIM5
        self.DIM6 = DIM6

        # Ensure dimensions are physically relevant
        np.testing.assert_((DIM5 + DIM6) < DIM1, "Invalid geometry specified.")
        np.testing.assert_(DIM4 < DIM3, "Invalid geometry specified.")
        np.testing.assert_(DIM4 < DIM2, "Invalid geometry specified.")

        # assign control point
        control_points = [[0.5*DIM2, 0.5*DIM5]]

        shift = [-0.5*DIM2+shift[0], -0.5*DIM1+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        db = 0.5*(DIM2 - DIM4)
        dt = 0.5*(DIM3 - DIM4)
        self.points = [
            [0., 0.], [DIM2, 0.], [DIM2, DIM5], [db+DIM4, DIM5], [db + DIM4, DIM1-DIM6],
            [db+DIM4+dt, DIM1-DIM6], [db+DIM4+dt, DIM1], [db-dt, DIM1], [db-dt, DIM1-DIM6],
            [db, DIM1-DIM6], [db, DIM5], [0, DIM5]
        ]
        self.facets = [
            [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
            [10, 11], [11, 0]
        ]

        self.shift_section()

    def getStressPoints(self, shift=(0., 0.)):
        """Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (0.5*self.DIM3-shift[0], 0.5*self.DIM1-shift[1])
        D = (0.5*self.DIM3-shift[0], -0.5*self.DIM1-shift[1])
        E = (-0.5*self.DIM3-shift[0], -0.5*self.DIM1-shift[1])
        F = (-0.5*self.DIM3-shift[0], 0.5*self.DIM1-shift[1])

        return C, D, E, F


class I1Section(Geometry):
    """Constructs a I1 section with the web's middle center at the origin *(0, 0)*, with four
    parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for more details.
    Added by JohnDN90.

    :param float DIM1: Twice distance from web end to flange end
    :param float DIM2: Thickness of web
    :param float DIM3: Length of web (spacing between flanges)
    :param float DIM4: Depth (y) of the I1-section
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a I1 cross-section with a depth of
    5.0 and width of 1.75, and generates a mesh with a maximum triangular area of
    0.02::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.I1Section(DIM1=1.0, DIM2=0.75, DIM3=4.0, DIM4=5.0)
        mesh = geometry.create_mesh(mesh_sizes=[0.02])

    ..  figure:: ../images/sections/i1_geometry.png
        :align: center
        :scale: 75 %

        I1 section geometry.

    ..  figure:: ../images/sections/i1_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, DIM3, DIM4, shift=[0, 0]):
        """Inits the I1section class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3
        self.DIM4 = DIM4

        # Ensure dimensions are physically relevant
        np.testing.assert_(DIM4 > DIM3, "Invalid geometry specified.")

        shift = [-0.5*(DIM1+DIM2)+shift[0], -0.5*DIM4+shift[1]]

        # assign control point
        control_points = [[0.5*(DIM1+DIM2), 0.5*DIM4]]

        super().__init__(control_points, shift)

        # construct the points and facets
        t = 0.5*(DIM4 - DIM3)
        self.points = [
            [0., 0.], [DIM1+DIM2, 0.], [DIM1+DIM2, t], [0.5*DIM1+DIM2, t], [0.5*DIM1+DIM2, t+DIM3],
            [DIM1+DIM2, t+DIM3], [DIM1+DIM2, DIM4], [0., DIM4], [0., t+DIM3], [0.5*DIM1, t+DIM3],
            [0.5*DIM1, t], [0., t]
        ]
        self.facets = [
            [0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 9], [9, 10],
            [10, 11], [11, 0]
        ]

        self.shift_section()

    def getStressPoints(self, shift=(0., 0.)):
        """Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (0.5*(self.DIM1+self.DIM2)-shift[0], 0.5*self.DIM4-shift[1])
        D = (0.5*(self.DIM1+self.DIM2)-shift[0], -0.5*self.DIM4-shift[1])
        E = (-0.5*(self.DIM1+self.DIM2)-shift[0], -0.5*self.DIM4-shift[1])
        F = (-0.5*(self.DIM1+self.DIM2)-shift[0], 0.5*self.DIM4-shift[1])

        return C, D, E, F


class LSection(Geometry):
    """Constructs a L section with the intersection's center at the origin *(0, 0)*, with four
    parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ for more details.
    Added by JohnDN90.

    :param float DIM1: Width (x) of the L-section
    :param float DIM2: Depth (y) of the L-section
    :param float DIM3: Thickness of flange (horizontal portion)
    :param float DIM4: Thickness of web (vertical portion)
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a L cross-section with a depth of 6.0 and width of 3.0, and
    generates a mesh with a maximum triangular area of 0.01::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.LSection(DIM1=3.0, DIM2=6.0, DIM3=0.375, DIM4=0.625)
        mesh = geometry.create_mesh(mesh_sizes=[0.01])

    ..  figure:: ../images/sections/l_geometry.png
        :align: center
        :scale: 75 %

        L section geometry.

    ..  figure:: ../images/sections/l_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, DIM3, DIM4, shift=[0, 0]):
        """Inits the LSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3
        self.DIM4 = DIM4

        # Ensure dimensions are physically relevant
        np.testing.assert_(DIM4 < DIM1, "Invalid geometry specified.")
        np.testing.assert_(DIM3 < DIM2, "Invalid geometry specified.")

        # assign control point
        control_points = [[0.5*DIM1, 0.5*DIM3]]

        shift = [-0.5*DIM4+shift[0], -0.5*DIM3+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        self.points = [[0, 0], [DIM1, 0], [DIM1, DIM3], [DIM4, DIM3], [DIM4, DIM2], [0, DIM2]]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]

        self.shift_section()

    def getStressPoints(self, shift=(0., 0.)):
        """Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (0.5*self.DIM4-shift[0], self.DIM2-0.5*self.DIM3-shift[1])
        D = (self.DIM1-0.5*self.DIM4-shift[0], -0.5*self.DIM3-shift[1])
        E = (-0.5*self.DIM4-shift[0], -0.5*self.DIM3-shift[1])
        F = (-0.5*self.DIM4-shift[0], self.DIM2-0.5*self.DIM3-shift[1])

        return C, D, E, F


class RODSection(Geometry):
    """Constructs a circular rod section with the center at the origin *(0, 0)*, with one parameter
    defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for more details. Added by
    JohnDN90.

    :param float DIM1: Radius of the circular rod section
    :param int n: Number of points discretising the circle
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a circular rod with a radius of 3.0 and 50 points discretizing
    the boundary, and generates a mesh with a maximum triangular area of 0.01::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.RODSection(DIM1=3.0, n=50)
        mesh = geometry.create_mesh(mesh_sizes=[0.01])

    ..  figure:: ../images/sections/rod_geometry.png
        :align: center
        :scale: 75 %

        Rod section geometry.

    ..  figure:: ../images/sections/rod_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, n, shift=[0, 0]):
        """Inits the RODSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        self.DIM1 = DIM1

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

    def getStressPoints(self, shift=(0., 0.)):
        """Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param float DIM1: Radius of the circular rod section
        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (-shift[0], self.DIM1-shift[1])
        D = (self.DIM1-shift[0], -shift[1])
        E = (-shift[0], -self.DIM1-shift[1])
        F = (-self.DIM1-shift[0], -shift[1])

        return C, D, E, F


class TSection(Geometry):
    """Constructs a T section with the top flange's middle center at the origin *(0, 0)*, with four
    parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ [5]_ for more
    details. Added by JohnDN90.

    :param float DIM1: Width (x) of top flange
    :param float DIM2: Depth (y) of the T-section
    :param float DIM3: Thickness of top flange
    :param float DIM4: Thickness of web
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a T cross-section with a depth of 4.0 and width of 3.0, and
    generates a mesh with a maximum triangular area of 0.001::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.TSection(DIM1=3.0, DIM2=4.0, DIM3=0.375, DIM4=0.25)
        mesh = geometry.create_mesh(mesh_sizes=[0.001])

    ..  figure:: ../images/sections/t_geometry.png
        :align: center
        :scale: 75 %

        T section geometry.

    ..  figure:: ../images/sections/t_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, DIM3, DIM4, shift=[0, 0]):
        """Inits the TSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3
        self.DIM4 = DIM4

        # Ensure dimensions are physically relevant
        np.testing.assert_(DIM4 < DIM1, "Invalid geometry specified.")
        np.testing.assert_(DIM3 < DIM2, "Invalid geometry specified.")

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

    def getStressPoints(self, shift=(0., 0.)):
        """
        Returns the coordinates of the stress evaluation points relative to the origin
        of the cross-section. The shift parameter can be used to make the coordinates
        relative to the centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (-shift[0], 0.5*self.DIM3-shift[1])
        D = (0.5*self.DIM1-shift[0], 0.5*self.DIM3-shift[1])
        E = (-shift[0], 0.5*self.DIM3-self.DIM2-shift[1])
        F = (-0.5*self.DIM1-shift[0], 0.5*self.DIM3-shift[1])

        return C, D, E, F


class T1Section(Geometry):
    """Constructs a T1 section with the right flange's middle center at the origin *(0, 0)*, with
    four parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for more
    details. Added by JohnDN90.

    :param float DIM1: Depth (y) of T1-section
    :param float DIM2: Length (x) of web
    :param float DIM3: Thickness of right flange
    :param float DIM4: Thickness of web
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a T1 cross-section with a depth of 3.0 and width of 3.875, and
    generates a mesh with a maximum triangular area of 0.001::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.T1Section(DIM1=3.0, DIM2=3.5, DIM3=0.375, DIM4=0.25)
        mesh = geometry.create_mesh(mesh_sizes=[0.001])

    ..  figure:: ../images/sections/t1_geometry.png
        :align: center
        :scale: 75 %

        T1 section geometry.

    ..  figure:: ../images/sections/t1_mesh.png
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
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3
        self.DIM4 = DIM4

        # Ensure dimensions are physically relevant
        np.testing.assert_(DIM4 < DIM1, "Invalid geometry specified.")

        shift = [-0.5*DIM3+shift[0], -0.5*DIM1+shift[1]]

        # assign control point
        control_points = [[0.5*DIM3, 0.5*DIM1]]

        super().__init__(control_points, shift)

        # construct the points and facets
        d1 = (DIM1 - DIM4) / 2.0
        self.points = [
            [0, 0], [DIM3, 0], [DIM3, DIM1], [0, DIM1], [0, d1 + DIM4], [-DIM2, d1 + DIM4],
            [-DIM2, d1], [0, d1]
        ]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 0]]

        self.shift_section()

    def getStressPoints(self, shift=(0., 0.)):
        """Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (0.5*self.DIM3-shift[0], -shift[1])
        D = (0.5*self.DIM3-shift[0], -0.5*self.DIM1-shift[1])
        E = (-0.5*self.DIM3-self.DIM2-shift[0], -shift[1])
        F = (0.5*self.DIM3-shift[0], 0.5*self.DIM1-shift[1])

        return C, D, E, F


class T2Section(Geometry):
    """Constructs a T2 section with the bottom flange's middle center at the origin *(0, 0)*, with
    four parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for more
    details. Added by JohnDN90.

    :param float DIM1: Width (x) of T2-section
    :param float DIM2: Depth (y) of T2-section
    :param float DIM3: Thickness of bottom flange
    :param float DIM4: Thickness of web
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a T2 cross-section with a depth of 4.0 and width of 3.0, and
    generates a mesh with a maximum triangular area of 0.005::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.T2Section(DIM1=3.0, DIM2=4.0, DIM3=0.375, DIM4=0.5)
        mesh = geometry.create_mesh(mesh_sizes=[0.005])

    ..  figure:: ../images/sections/t2_geometry.png
        :align: center
        :scale: 75 %

        T2 section geometry.

    ..  figure:: ../images/sections/t2_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, DIM3, DIM4, shift=[0, 0]):
        """Inits the T2Section class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3
        self.DIM4 = DIM4

        # Ensure dimensions are physically relevant
        np.testing.assert_(DIM4 < DIM1, "Invalid geometry specified.")
        np.testing.assert_(DIM3 < DIM2, "Invalid geometry specified.")

        # assign control point
        control_points = [[0.5*DIM1, 0.5*DIM3]]

        shift = [-0.5*DIM1+shift[0], -0.5*DIM3+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        d1 = 0.5*(DIM1 - DIM4)
        self.points = [
            [0., 0.], [DIM1, 0.], [DIM1, DIM3], [DIM1-d1, DIM3], [DIM1-d1, DIM2], [d1, DIM2],
            [d1, DIM3], [0, DIM3]
        ]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 0]]

        self.shift_section()

    def getStressPoints(self, shift=(0., 0.)):
        """Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (0.5*self.DIM4-shift[0], self.DIM2-0.5*self.DIM3-shift[1])
        D = (0.5*self.DIM1-shift[0], -0.5*self.DIM3-shift[1])
        E = (-0.5*self.DIM1-shift[0], -0.5*self.DIM3-shift[1])
        F = (-0.5*self.DIM4-shift[0], self.DIM2-0.5*self.DIM3-shift[1])

        return C, D, E, F


class TUBESection(Geometry):
    """Constructs a circular tube section with the center at the origin *(0, 0)*, with two
    parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for more
    details. Added by JohnDN90.

    :param float DIM1: Outer radius of the circular tube section
    :param float DIM2: Inner radius of the circular tube section
    :param int n: Number of points discretising the circle
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a circular tube cross-section with an outer radius of 3.0 and an
    inner radius of 2.5, and generates a mesh with 37 points discretizing the boundaries and a
    maximum triangular area of 0.01::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.TUBESection(DIM1=3.0, DIM2=2.5, n=37)
        mesh = geometry.create_mesh(mesh_sizes=[0.01])

    ..  figure:: ../images/sections/tube_geometry.png
        :align: center
        :scale: 75 %

        TUBE section geometry.

    ..  figure:: ../images/sections/tube_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, n, shift=[0, 0]):
        """Inits the TUBESection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2

        # Ensure dimensions are physically relevant
        np.testing.assert_(DIM2 < DIM1, "Invalid geometry specified.")

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

    def getStressPoints(self, shift=(0., 0.)):
        """Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (-shift[0], self.DIM1-shift[1])
        D = (self.DIM1-shift[0], -shift[1])
        E = (-shift[0], -self.DIM1-shift[1])
        F = (-self.DIM1-shift[0], -shift[1])

        return C, D, E, F


class TUBE2Section(Geometry):
    """Constructs a circular TUBE2 section with the center at the origin *(0, 0)*, with two
    parameters defining dimensions. See MSC Nastran documentation [1]_ for more details. Added by
    JohnDN90.

    :param float DIM1: Outer radius of the circular tube section
    :param float DIM2: Thickness of wall
    :param int n: Number of points discretising the circle
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a ciruclar TUBE2 cross-section with an outer radius of 3.0 and a
    wall thickness of 0.5, and generates a mesh with 37 point discritizing the boundary and a
    maximum triangular area of 0.01::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.TUBE2Section(DIM1=3.0, DIM2=0.5, n=37)
        mesh = geometry.create_mesh(mesh_sizes=[0.01])

    ..  figure:: ../images/sections/tube2_geometry.png
        :align: center
        :scale: 75 %

        TUBE2 section geometry.

    ..  figure:: ../images/sections/tube2_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, n, shift=[0, 0]):
        """Inits the TUBE2Section class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2

        # Ensure dimensions are physically relevant
        np.testing.assert_(DIM2 < DIM1, "Invalid geometry specified.")

        d = 2.0*DIM1
        t = DIM2

        # assign control point
        control_points = [[d * 0.5 - t * 0.5, 0]]

        super().__init__(control_points, shift)

        # specify a hole in the centre of the section
        self.holes = [[0., 0.]]

        # loop through each point of the section
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

    def getStressPoints(self, shift=(0., 0.)):
        """Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (-shift[0], self.DIM1-shift[1])
        D = (self.DIM1-shift[0], -shift[1])
        E = (-shift[0], -self.DIM1-shift[1])
        F = (-self.DIM1-shift[0], -shift[1])

        return C, D, E, F


class ZSection(Geometry):
    """Constructs a Z section with the web's middle center at the origin *(0, 0)*, with four
    parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for more details.
    Added by JohnDN90.

    :param float DIM1: Width (x) of horizontal members
    :param float DIM2: Thickness of web
    :param float DIM3: Spacing between horizontal members (length of web)
    :param float DIM4: Depth (y) of Z-section
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a rectangular cross-section with a depth of 4.0 and width of
    2.75, and generates a mesh with a maximum triangular area of 0.005::

        import sectionproperties.pre.nastran_sections as nsections

        geometry = nsections.ZSection(DIM1=1.125, DIM2=0.5, DIM3=3.5, DIM4=4.0)
        mesh = geometry.create_mesh(mesh_sizes=[0.005])

    ..  figure:: ../images/sections/z_geometry.png
        :align: center
        :scale: 75 %

        Z section geometry.

    ..  figure:: ../images/sections/z_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """

    def __init__(self, DIM1, DIM2, DIM3, DIM4, shift=[0, 0]):
        """Inits the ZSection class."""

        # force dimensions to be floating point values
        DIM1 *= 1.0
        DIM2 *= 1.0
        DIM3 *= 1.0
        DIM4 *= 1.0
        self.DIM1 = DIM1
        self.DIM2 = DIM2
        self.DIM3 = DIM3
        self.DIM4 = DIM4

        # Ensure dimensions are physically relevant
        np.testing.assert_(DIM4 > DIM3, "Invalid geometry specified.")

        # assign control point
        control_points = [[DIM1+0.5*DIM2, 0.5*DIM4]]

        shift = [-0.5*(DIM1+DIM2)+shift[0], -0.5*DIM4+shift[1]]
        super().__init__(control_points, shift)

        # construct the points and facets
        t = 0.5*(DIM4 - DIM3)
        self.points = [
            [DIM1, 0.], [2.*DIM1+DIM2, 0.], [2.*DIM1+DIM2, t], [DIM1+DIM2, t], [DIM1+DIM2, DIM4],
            [0., DIM4], [0., DIM4-t],  [DIM1, DIM4-t]
        ]
        self.facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 0]]

        self.shift_section()

    def getStressPoints(self, shift=(0., 0.)):
        """Returns the coordinates of the stress evaluation points relative to the origin of the
        cross-section. The shift parameter can be used to make the coordinates relative to the
        centroid or the shear center.

        :param shift: Vector that shifts the cross-section by *(x, y)*
        :type shift: tuple(float, float)
        :returns: Stress evaluation points relative to shifted origin - C, D, E, F
        """

        C = (0.5*self.DIM2-shift[0], 0.5*self.DIM4-shift[1])
        D = (0.5*self.DIM2+self.DIM1-shift[0], -0.5*self.DIM4-shift[1])
        E = (-0.5*self.DIM2-shift[0], -0.5*self.DIM4-shift[1])
        F = (-0.5*self.DIM2-self.DIM1-shift[0], 0.5*self.DIM4-shift[1])

        return C, D, E, F
