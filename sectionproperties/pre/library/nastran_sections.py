import numpy as np
from shapely.geometry import Polygon
import sectionproperties.pre.geometry as geometry
import sectionproperties.pre.pre as pre
from sectionproperties.pre.library.utils import draw_radius


def nastran_bar(
    DIM1: float, DIM2: float, material: pre.Material = pre.DEFAULT_MATERIAL
) -> geometry.Geometry:
    """Constructs a BAR section with the center at the origin *(0, 0)*, with two parameters
    defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ [5]_ for definition of
    parameters. Added by JohnDN90.

    :param float DIM1: Width (x) of bar
    :param float DIM2: Depth (y) of bar
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a BAR cross-section with a depth of 1.5 and width of 2.0, and
    generates a mesh with a maximum triangular area of 0.001::

        from sectionproperties.pre.library.nastran_sections import nastran_bar

        geom = nastran_bar(DIM1=2.0, DIM2=1.5)
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
    points = [
        [-0.5 * DIM1, -0.5 * DIM2],
        [0.5 * DIM1, -0.5 * DIM2],
        [0.5 * DIM1, 0.5 * DIM2],
        [-0.5 * DIM1, 0.5 * DIM2],
    ]
    geom = geometry.Geometry(Polygon(points), material)
    C = (0.5 * DIM1, 0.5 * DIM2)
    D = (0.5 * DIM1, -0.5 * DIM2)
    E = (-0.5 * DIM1, -0.5 * DIM2)
    F = (-0.5 * DIM1, 0.5 * DIM2)
    geom.recovery_points = [C, D, E, F]
    return geom


def nastran_box(
    DIM1: float,
    DIM2: float,
    DIM3: float,
    DIM4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a BOX section with the center at the origin *(0, 0)*, with four parameters
    defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ [5]_ for definition of
    parameters. Added by JohnDN90.

    :param float DIM1: Width (x) of box
    :param float DIM2: Depth (y) of box
    :param float DIM3: Thickness of box in y direction
    :param float DIM4: Thickness of box in x direction
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a BOX cross-section with a depth of 3.0 and width of 4.0, and
    generates a mesh with a maximum triangular area of 0.001::

        from sectionproperties.pre.library.nastran_sections import nastran_box

        geom = nastran_box(DIM1=4.0, DIM2=3.0, DIM3=0.375, DIM4=0.5)
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
    # Ensure dimensions are physically relevant
    np.testing.assert_(2.0 * DIM4 < DIM1, "Invalid geometry specified.")
    np.testing.assert_(2.0 * DIM3 < DIM2, "Invalid geometry specified.")
    points_outer = [
        [-0.5 * DIM1, -0.5 * DIM2],
        [0.5 * DIM1, -0.5 * DIM2],
        [0.5 * DIM1, 0.5 * DIM2],
        [-0.5 * DIM1, 0.5 * DIM2],
    ]
    points_inner = [
        [-0.5 * DIM1 + DIM4, -0.5 * DIM2 + DIM3],
        [0.5 * DIM1 - DIM4, -0.5 * DIM2 + DIM3],
        [0.5 * DIM1 - DIM4, 0.5 * DIM2 - DIM3],
        [-0.5 * DIM1 + DIM4, 0.5 * DIM2 - DIM3],
    ]

    inner_box = Polygon(points_inner)
    outer_box = Polygon(points_outer)

    C = (0.5 * DIM1, 0.5 * DIM2)
    D = (0.5 * DIM1, -0.5 * DIM2)
    E = (-0.5 * DIM1, -0.5 * DIM2)
    F = (-0.5 * DIM1, 0.5 * DIM2)

    geom = geometry.Geometry(outer_box - inner_box, material)
    geom.recovery_points = [C, D, E, F]

    return geom


def nastran_box1(
    DIM1: float,
    DIM2: float,
    DIM3: float,
    DIM4: float,
    DIM5: float,
    DIM6: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a BOX1 section with the center at the origin *(0, 0)*, with six parameters
    defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for more details. Added by
    JohnDN90.

    :param float DIM1: Width (x) of box
    :param float DIM2: Depth (y) of box
    :param float DIM3: Thickness of top wall
    :param float DIM4: Thickness of bottom wall
    :param float DIM5: Thickness of left wall
    :param float DIM6: Thickness of right wall
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a BOX1 cross-section with a depth of 3.0 and width of 4.0, and
    generates a mesh with a maximum triangular area of 0.007::

        from sectionproperties.pre.library.nastran_sections import nastran_box1

        geom = nastran_box1(
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
    # Ensure dimensions are physically relevant
    np.testing.assert_(DIM5 + DIM6 < DIM1, "Invalid geometry specified.")
    np.testing.assert_(DIM3 + DIM4 < DIM2, "Invalid geometry specified.")

    exterior_points = [
        [0.0, 0.0],
        [DIM1, 0.0],
        [DIM1, DIM2],
        [0.0, DIM2],
    ]
    interior_points = [
        [DIM6, DIM4],
        [DIM1 - DIM5, DIM4],
        [DIM1 - DIM5, DIM2 - DIM3],
        [DIM6, DIM2 - DIM3],
    ]
    geom = geometry.Geometry(
        Polygon(exterior_points) - Polygon(interior_points), material
    )

    C = (0.5 * DIM1, 0.5 * DIM2)
    D = (0.5 * DIM1, -0.5 * DIM2)
    E = (-0.5 * DIM1, -0.5 * DIM2)
    F = (-0.5 * DIM1, 0.5 * DIM2)
    geom.recovery_points = [C, D, E, F]
    return geom


def nastran_chan(
    DIM1: float,
    DIM2: float,
    DIM3: float,
    DIM4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a CHAN (C-Channel) section with the web's middle center at the origin *(0, 0)*,
    with four parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for
    more details. Added by JohnDN90.

    :param float DIM1: Width (x) of the CHAN-section
    :param float DIM2: Depth (y) of the CHAN-section
    :param float DIM3: Thickness of web (vertical portion)
    :param float DIM4: Thickness of flanges (top/bottom portion)
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a CHAN cross-section with a depth of 4.0 and width of 2.0, and
    generates a mesh with a maximum triangular area of 0.008::

        from sectionproperties.pre.library.nastran_sections import nastran_chan

        geom = nastran_chan(DIM1=2.0, DIM2=4.0, DIM3=0.25, DIM4=0.5)
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
    # Ensure dimensions are physically relevant
    np.testing.assert_(2.0 * DIM4 < DIM2, "Invalid geometry specified.")
    np.testing.assert_(DIM3 < DIM1, "Invalid geometry specified.")

    # construct the points
    points = [
        [0.0, 0.0],
        [DIM1, 0.0],
        [DIM1, DIM4],
        [DIM3, DIM4],
        [DIM3, DIM2 - DIM4],
        [DIM1, DIM2 - DIM4],
        [DIM1, DIM2],
        [0.0, DIM2],
    ]

    geom = geometry.Geometry(Polygon(points), material)

    C = (0.5 * DIM1, 0.5 * DIM2)
    D = (0.5 * DIM1, -0.5 * DIM2)
    E = (-0.5 * DIM1, -0.5 * DIM2)
    F = (-0.5 * DIM1, 0.5 * DIM2)

    geom.recovery_points = [C, D, E, F]
    return geom


def nastran_chan1(
    DIM1: float,
    DIM2: float,
    DIM3: float,
    DIM4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a CHAN1 (C-Channel) section with the web's middle center at the origin *(0, 0)*,
    with four parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for
    more details. Added by JohnDN90.

    :param float DIM1: Width (x) of channels
    :param float DIM2: Thickness (x) of web
    :param float DIM3: Spacing between channels (length of web)
    :param float DIM4: Depth (y) of CHAN1-section
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a CHAN1 cross-section with a depth of 4.0 and width of 1.75, and
    generates a mesh with a maximum triangular area of 0.01::

        from sectionproperties.pre.library.nastran_sections import nastran_chan1

        geom = nastran_chan1(DIM1=0.75, DIM2=1.0, DIM3=3.5, DIM4=4.0)
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
    # Ensure dimensions are physically relevant
    np.testing.assert_(DIM4 > DIM3, "Invalid geometry specified.")

    # construct the points and facets
    tf = 0.5 * (DIM4 - DIM3)
    points = [
        [0, 0],
        [DIM1 + DIM2, 0],
        [DIM1 + DIM2, tf],
        [DIM2, tf],
        [DIM2, tf + DIM3],
        [DIM2 + DIM1, tf + DIM3],
        [DIM2 + DIM1, DIM4],
        [0, DIM4],
    ]
    geom = geometry.Geometry(Polygon(points), material)
    C = (0.5 * DIM2 + DIM1, 0.5 * DIM4)
    D = (0.5 * DIM2 + DIM1, -0.5 * DIM4)
    E = (-0.5 * DIM2, -0.5 * DIM4)
    F = (-0.5 * DIM2, 0.5 * DIM4)
    geom.recovery_points = [C, D, E, F]

    return geom


def nastran_chan2(
    DIM1: float,
    DIM2: float,
    DIM3: float,
    DIM4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a CHAN2 (C-Channel) section with the bottom web's middle center at the origin
    *(0, 0)*, with four parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_
    [4]_ for more details. Added by JohnDN90.

    :param float DIM1: Thickness of channels
    :param float DIM2: Thickness of web
    :param float DIM3: Depth (y) of CHAN2-section
    :param float DIM4: Width (x) of CHAN2-section
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a CHAN2 cross-section with a depth of 2.0 and width of 4.0, and
    generates a mesh with a maximum triangular area of 0.01::

        from sectionproperties.pre.library.nastran_sections import nastran_chan2

        geom = nastran_chan2(DIM1=0.375, DIM2=0.5, DIM3=2.0, DIM4=4.0)
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
    # Ensure dimensions are physically relevant
    np.testing.assert_(DIM4 > 2.0 * DIM1, "Invalid geometry specified.")
    np.testing.assert_(DIM3 > DIM2, "Invalid geometry specified.")

    # construct the points and facets
    points = [
        [0.0, 0.0],
        [DIM4, 0.0],
        [DIM4, DIM3],
        [DIM4 - DIM1, DIM3],
        [DIM4 - DIM1, DIM2],
        [DIM1, DIM2],
        [DIM1, DIM3],
        [0.0, DIM3],
    ]
    geom = geometry.Geometry(Polygon(points), material)
    C = (0.5 * DIM4, DIM3 - 0.5 * DIM2)
    D = (0.5 * DIM4, -0.5 * DIM2)
    E = (-0.5 * DIM4, -0.5 * DIM2)
    F = (-0.5 * DIM4, DIM3 - 0.5 * DIM2)
    geom.recovery_points = [C, D, E, F]

    return geom


def nastran_cross(
    DIM1: float,
    DIM2: float,
    DIM3: float,
    DIM4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs Nastran's cruciform/cross section with the intersection's middle center at the
    origin *(0, 0)*, with four parameters defining dimensions. See Nastran documentation [1]_ [2]_
    [3]_ [4]_ for more details. Added by JohnDN90.

    :param float DIM1: Twice the width of horizontal member protruding from the vertical center
        member
    :param float DIM2: Thickness of the vertical member
    :param float DIM3: Depth (y) of the CROSS-section
    :param float DIM4: Thickness of the horizontal members
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a rectangular cross-section with a depth of 3.0 and width of
    1.875, and generates a mesh with a maximum triangular area of 0.008::

        from sectionproperties.pre.library.nastran_sections import nastran_cross

        geom = nastran_cross(DIM1=1.5, DIM2=0.375, DIM3=3.0, DIM4=0.25)
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
    # Ensure dimensions are physically relevant
    np.testing.assert_(DIM4 < DIM3, "Invalid geometry specified.")

    # construct the points and facets
    d = 0.5 * (DIM3 - DIM4)
    points = [
        [0.5 * DIM1, 0],
        [0.5 * DIM1 + DIM2, 0],
        [0.5 * DIM1 + DIM2, d],
        [DIM1 + DIM2, d],
        [DIM1 + DIM2, d + DIM4],
        [0.5 * DIM1 + DIM2, d + DIM4],
        [0.5 * DIM1 + DIM2, DIM3],
        [0.5 * DIM1, DIM3],
        [0.5 * DIM1, d + DIM4],
        [0, d + DIM4],
        [0, d],
        [0.5 * DIM1, d],
    ]
    geom = geometry.Geometry(Polygon(points), material)
    C = (0, 0.5 * DIM3)
    D = (0.5 * (DIM1 + DIM2), 0)
    E = (0, -0.5 * DIM3)
    F = (-0.5 * (DIM1 + DIM2), 0)
    geom.recovery_points = [C, D, E, F]
    return geom


def nastran_fcross(
    DIM1: float,
    DIM2: float,
    DIM3: float,
    DIM4: float,
    DIM5: float,
    DIM6: float,
    DIM7: float,
    DIM8: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a flanged cruciform/cross section with the intersection's middle center at the
    origin *(0, 0)*, with eight parameters defining dimensions. Added by JohnDN90.

    :param float DIM1: Depth (y) of flanged cruciform
    :param float DIM2: Width (x) of flanged cruciform
    :param float DIM3: Thickness of vertical web
    :param float DIM4: Thickness of horizontal web
    :param float DIM5: Length of flange attached to vertical web
    :param float DIM6: Thickness of flange attached to vertical web
    :param float DIM7: Length of flange attached to horizontal web
    :param float DIM8: Thickness of flange attached to horizontal web
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example demonstrates the creation of a flanged cross section::

        from sectionproperties.pre.library.nastran_sections import nastran_fcross

        geom = nastran_fcross(
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
    # Ensure dimensions are physically relevant
    # TODO: Finish dimension checks.
    np.testing.assert_(DIM5 > DIM3, "Invalid geometry specified.")
    np.testing.assert_(DIM7 > DIM4, "Invalid geometry specified.")
    np.testing.assert_(DIM7 < DIM1, "Invalid geometry specified.")
    np.testing.assert_(DIM5 < DIM2, "Invalid geometry specified.")
    np.testing.assert_(DIM8 < (0.5 * DIM2 - 0.5 * DIM3), "Invalid geometry specified.")
    np.testing.assert_(DIM6 < (0.5 * DIM1 - 0.5 * DIM4), "Invalid geometry specified.")

    # construct the points and facets
    points = [
        [0.5 * DIM3, -0.5 * DIM4],
        [0.5 * DIM2 - DIM8, -0.5 * DIM4],
        [0.5 * DIM2 - DIM8, -0.5 * DIM7],
        [0.5 * DIM2, -0.5 * DIM7],
        [0.5 * DIM2, 0.5 * DIM7],
        [0.5 * DIM2 - DIM8, 0.5 * DIM7],
        [0.5 * DIM2 - DIM8, 0.5 * DIM4],
        [0.5 * DIM3, 0.5 * DIM4],
        [0.5 * DIM3, 0.5 * DIM1 - DIM6],
        [0.5 * DIM5, 0.5 * DIM1 - DIM6],
        [0.5 * DIM5, 0.5 * DIM1],
        [-0.5 * DIM5, 0.5 * DIM1],
        [-0.5 * DIM5, 0.5 * DIM1 - DIM6],
        [-0.5 * DIM3, 0.5 * DIM1 - DIM6],
        [-0.5 * DIM3, 0.5 * DIM4],
        [-0.5 * DIM2 + DIM8, 0.5 * DIM4],
        [-0.5 * DIM2 + DIM8, 0.5 * DIM7],
        [-0.5 * DIM2, 0.5 * DIM7],
        [-0.5 * DIM2, -0.5 * DIM7],
        [-0.5 * DIM2 + DIM8, -0.5 * DIM7],
        [-0.5 * DIM2 + DIM8, -0.5 * DIM4],
        [-0.5 * DIM3, -0.5 * DIM4],
        [-0.5 * DIM3, -0.5 * DIM1 + DIM6],
        [-0.5 * DIM5, -0.5 * DIM1 + DIM6],
        [-0.5 * DIM5, -0.5 * DIM1],
        [0.5 * DIM5, -0.5 * DIM1],
        [0.5 * DIM5, -0.5 * DIM1 + DIM6],
        [0.5 * DIM3, -0.5 * DIM1 + DIM6],
    ]
    geom = geometry.Geometry(Polygon(points), material)

    C = (0, 0.5 * DIM1)
    D = (0.5 * DIM2, 0)
    E = (0, -0.5 * DIM1)
    F = (-0.5 * DIM2, 0)
    geom.recovery_points = [C, D, E, F]
    return geom


def nastran_dbox(
    DIM1: float,
    DIM2: float,
    DIM3: float,
    DIM4: float,
    DIM5: float,
    DIM6: float,
    DIM7: float,
    DIM8: float,
    DIM9: float,
    DIM10: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a DBOX section with the center at the origin *(0, 0)*, with ten parameters
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
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a DBOX cross-section with a depth of 3.0 and width of 8.0, and
    generates a mesh with a maximum triangular area of 0.01::

        from sectionproperties.pre.library.nastran_sections import nastran_dbox

        geom = nastran_dbox(
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
    # Ensure dimensions are physically relevant
    np.testing.assert_((DIM4 + DIM5 + DIM6) < DIM1, "Invalid geometry specified.")
    np.testing.assert_((DIM4 + 0.5 * DIM5) < DIM3, "Invalid geometry specified.")
    np.testing.assert_((DIM7 + DIM8) < DIM2, "Invalid geometry specified.")
    np.testing.assert_((DIM9 + DIM10) < DIM2, "Invalid geometry specified.")

    # construct the points and facets
    exterior_points = [
        [0.0, 0.0],
        [DIM1, 0.0],
        [DIM1, DIM2],
        [0.0, DIM2],
    ]
    interior_points_1 = [
        [DIM4, DIM8],
        [DIM3 - DIM5 / 2.0, DIM8],
        [DIM3 - DIM5 / 2.0, DIM2 - DIM7],
        [DIM4, DIM2 - DIM7],
    ]
    interior_points_2 = [
        [DIM3 + DIM5 / 2.0, DIM10],
        [DIM1 - DIM6, DIM10],
        [DIM1 - DIM6, DIM2 - DIM9],
        [DIM3 + DIM5 / 2.0, DIM2 - DIM9],
    ]
    geom = geometry.Geometry(
        Polygon(exterior_points)
        - Polygon(interior_points_1)
        - Polygon(interior_points_2)
    )
    C = (0.5 * DIM1, 0.5 * DIM2)
    D = (0.5 * DIM1, -0.5 * DIM2)
    E = (-0.5 * DIM1, -0.5 * DIM2)
    F = (-0.5 * DIM1, 0.5 * DIM2)
    geom.recovery_points = [C, D, E, F]
    return geom


def nastran_gbox(
    DIM1: float,
    DIM2: float,
    DIM3: float,
    DIM4: float,
    DIM5: float,
    DIM6: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a GBOX section with the center at the origin *(0, 0)*, with six parameters
    defining dimensions. See ASTROS documentation [5]_ for more details. Added by JohnDN90.

    :param float DIM1: Width (x) of the GBOX-section
    :param float DIM2: Depth (y) of the GBOX-section
    :param float DIM3: Thickness of top flange
    :param float DIM4: Thickness of bottom flange
    :param float DIM5: Thickness of webs
    :param float DIM6: Spacing between webs
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a GBOX cross-section with a depth of 2.5 and width of 6.0, and
    generates a mesh with a maximum triangular area of 0.01::

        from sectionproperties.pre.library.nastran_sections import nastran_gbox

        geom = nastran_gbox(
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
    # Ensure dimensions are physically relevant
    np.testing.assert_((DIM3 + DIM4) < DIM2, "Invalid geometry specified.")
    np.testing.assert_((2.0 * DIM5 + DIM6) < DIM1, "Invalid geometry specified.")

    # construct the points and facets
    d = 0.5 * (DIM1 - DIM6 - 2.0 * DIM5)
    exterior_points = [
        [0.0, 0.0],
        [DIM1, 0.0],
        [DIM1, DIM4],
        [d + 2.0 * DIM5 + DIM6, DIM4],
        [d + 2.0 * DIM5 + DIM6, DIM2 - DIM3],
        [DIM1, DIM2 - DIM3],
        [DIM1, DIM2],
        [0.0, DIM2],
        [0.0, DIM2 - DIM3],
        [d, DIM2 - DIM3],
        [d, DIM4],
        [0.0, DIM4],
    ]
    interior_points = [
        [d + DIM5, DIM4],
        [d + DIM5 + DIM6, DIM4],
        [d + DIM5 + DIM6, DIM2 - DIM3],
        [d + DIM5, DIM2 - DIM3],
    ]
    geom = geometry.Geometry(
        Polygon(exterior_points) - Polygon(interior_points), material
    )

    C = (0.5 * DIM1, 0.5 * DIM2)
    D = (0.5 * DIM1, -0.5 * DIM2)
    E = (-0.5 * DIM1, -0.5 * DIM2)
    F = (-0.5 * DIM1, 0.5 * DIM2)
    geom.recovery_points = [C, D, E, F]

    return geom


def nastran_h(
    DIM1: float,
    DIM2: float,
    DIM3: float,
    DIM4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a H section with the middle web's middle center at the origin *(0, 0)*, with four
    parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for more details.
    Added by JohnDN90.

    :param float DIM1: Spacing between vertical flanges (length of web)
    :param float DIM2: Twice the thickness of the vertical flanges
    :param float DIM3: Depth (y) of the H-section
    :param float DIM4: Thickness of the middle web
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a H cross-section with a depth of 3.5 and width of 2.75, and
    generates a mesh with a maximum triangular area of 0.005::

        from sectionproperties.pre.library.nastran_sections import nastran_h

        geom = nastran_h(DIM1=2.0, DIM2=0.75, DIM3=3.5, DIM4=0.25)
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
    # Ensure dimensions are physically relevant
    np.testing.assert_(DIM4 < DIM3, "Invalid geometry specified.")

    d1 = 0.5 * (DIM3 - DIM4)
    d2 = 0.5 * DIM2

    # construct the points and facets
    points = [
        [0, 0],
        [d2, 0],
        [d2, d1],
        [d2 + DIM1, d1],
        [d2 + DIM1, 0],
        [DIM1 + DIM2, 0],
        [DIM1 + DIM2, DIM3],
        [DIM1 + DIM2 - d2, DIM3],
        [DIM1 + DIM2 - d2, d1 + DIM4],
        [d2, d1 + DIM4],
        [d2, DIM3],
        [0, DIM3],
    ]
    geom = geometry.Geometry(Polygon(points), material)
    C = (0.5 * (DIM1 + DIM2), 0.5 * DIM3)
    D = (0.5 * (DIM1 + DIM2), -0.5 * DIM3)
    E = (-0.5 * (DIM1 + DIM2), -0.5 * DIM3)
    F = (-0.5 * (DIM1 + DIM2), 0.5 * DIM3)
    geom.recovery_points = [C, D, E, F]
    return geom


def nastran_hat(
    DIM1: float,
    DIM2: float,
    DIM3: float,
    DIM4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a Hat section with the top most section's middle center at the origin *(0, 0)*,
    with four parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for
    more details. Note that HAT in ASTROS is actually HAT1 in this code. Added by JohnDN90.

    :param float DIM1: Depth (y) of HAT-section
    :param float DIM2: Thickness of HAT-section
    :param float DIM3: Width (x) of top most section
    :param float DIM4: Width (x) of bottom sections
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a HAT cross-section with a depth of 1.25 and width of 2.5, and
    generates a mesh with a maximum triangular area of 0.001::

        from sectionproperties.pre.library.nastran_sections import nastran_hat

        geom = nastran_hat(DIM1=1.25, DIM2=0.25, DIM3=1.5, DIM4=0.5)
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
    # Ensure dimensions are physically relevant
    np.testing.assert_(2.0 * DIM2 < DIM1, "Invalid geometry specified.")

    # construct the points and facets
    points = [
        [0.0, 0.0],
        [DIM4 + DIM2, 0.0],
        [DIM4 + DIM2, DIM1 - DIM2],
        [DIM4 + DIM3 - DIM2, DIM1 - DIM2],
        [DIM4 + DIM3 - DIM2, 0.0],
        [2 * DIM4 + DIM3, 0.0],
        [2.0 * DIM4 + DIM3, DIM2],
        [DIM4 + DIM3, DIM2],
        [DIM4 + DIM3, DIM1],
        [DIM4, DIM1],
        [DIM4, DIM2],
        [0.0, DIM2],
    ]

    geom = geometry.Geometry(Polygon(points), material)
    C = (0.5 * DIM3, 0.5 * DIM2)
    D = (0.5 * DIM3 + DIM4, -DIM1 + DIM2)
    E = (-0.5 * DIM3 - DIM4, -DIM1 + DIM2)
    F = (-0.5 * DIM3, 0.5 * DIM2)
    geom.recovery_points = [C, D, E, F]
    # geometry.compile_geometry()
    return geom


def nastran_hat1(
    DIM1: float,
    DIM2: float,
    DIM3: float,
    DIM4: float,
    DIM5: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a HAT1 section with the bottom plate's bottom center at the origin *(0, 0)*,
    with five parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [5]_ for
    definition of parameters. Note that in ASTROS, HAT1 is called HAT. Added by JohnDN90.

    :param float DIM1: Width(x) of the HAT1-section
    :param float DIM2: Depth (y) of the HAT1-section
    :param float DIM3: Width (x) of hat's top flange
    :param float DIM4: Thickness of hat stiffener
    :param float DIM5: Thicknesss of bottom plate
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a HAT1 cross-section with a depth of 2.0 and width of 4.0, and
    generates a mesh with a maximum triangular area of 0.005::

        from sectionproperties.pre.library.nastran_sections import nastran_hat1

        geom = nastran_hat1(DIM1=4.0, DIM2=2.0, DIM3=1.5, DIM4=0.1875, DIM5=0.375)
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
    # Ensure dimensions are physically relevant
    np.testing.assert_((2.0 * DIM4 + DIM5) < DIM2, "Invalid geometry specified.")
    np.testing.assert_(DIM3 < DIM1, "Invalid geometry specified.")

    # create bottom rectangular plate
    bottom_plate = nastran_bar(DIM1=DIM1, DIM2=DIM5).shift_section(y_offset=DIM5 / 2)

    # create the hat stiffener
    d1 = DIM2 - DIM5
    d2 = DIM4
    d3 = DIM3
    d4 = 0.5 * (DIM1 - DIM3)

    hat = nastran_hat(DIM1=d1, DIM2=d2, DIM3=d3, DIM4=d4)
    # Merge the two sections into one geometry
    geom = (
        hat.align_center(bottom_plate).align_to(bottom_plate, on="top") + bottom_plate
    )

    C = (-0.5 * DIM1, 0)
    D = (0.5 * DIM1, 0)
    E = (-0.5 * DIM3, DIM2)
    F = (0.5 * DIM3, DIM2)

    geom.recovery_points = [C, D, E, F]
    # geometry.compile_geometry()

    return geom


def nastran_hexa(
    DIM1: float, DIM2: float, DIM3: float, material: pre.Material = pre.DEFAULT_MATERIAL
) -> geometry.Geometry:
    """Constructs a HEXA (hexagon) section with the center at the origin *(0, 0)*, with three
    parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for more details.
    Added by JohnDN90.

    :param float DIM1: Spacing between bottom right point and right most point
    :param float DIM2: Width (x) of hexagon
    :param float DIM3: Depth (y) of hexagon
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a rectangular cross-section with a depth of 1.5 and width of 2.0,
    and generates a mesh with a maximum triangular area of 0.005::

        from sectionproperties.pre.library.nastran_sections import nastran_hexa

        geom = nastran_hexa(DIM1=0.5, DIM2=2.0, DIM3=1.5)
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

    # Ensure dimensions are physically relevant
    np.testing.assert_(DIM2 > DIM1, "Invalid geometry specified.")

    # construct the points and facets
    points = [
        [DIM1, 0.0],
        [DIM2 - DIM1, 0.0],
        [DIM2, 0.5 * DIM3],
        [DIM2 - DIM1, DIM3],
        [DIM1, DIM3],
        [0.0, 0.5 * DIM3],
    ]

    geom = geometry.Geometry(Polygon(points), material)
    C = (0, 0.5 * DIM3)
    D = (0, -0.5 * DIM3)
    E = 0.5 * DIM2
    F = -0.5 * DIM2
    geom.recovery_points = [C, D, E, F]
    return geom


def nastran_i(
    DIM1: float,
    DIM2: float,
    DIM3: float,
    DIM4: float,
    DIM5: float,
    DIM6: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs Nastran's I section with the bottom flange's middle center at the origin
    *(0, 0)*, with six parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_
    [4]_ for definition of parameters. Added by JohnDN90.

    :param float DIM1: Depth(y) of the I Section
    :param float DIM2: Width (x) of bottom flange
    :param float DIM3: Width (x) of top flange
    :param float DIM4: Thickness of web
    :param float DIM5: Thickness of bottom web
    :param float DIM6: Thickness of top web
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a Nastran I cross-section with a depth of 5.0, and generates a
    mesh with a maximum triangular area of 0.008::

        from sectionproperties.pre.library.nastran_sections import nastran_i

        geom = nastran_i(
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
    # Ensure dimensions are physically relevant
    np.testing.assert_((DIM5 + DIM6) < DIM1, "Invalid geometry specified.")
    np.testing.assert_(DIM4 < DIM3, "Invalid geometry specified.")
    np.testing.assert_(DIM4 < DIM2, "Invalid geometry specified.")

    # construct the points and facets
    db = 0.5 * (DIM2 - DIM4)
    dt = 0.5 * (DIM3 - DIM4)
    points = [
        [0.0, 0.0],
        [DIM2, 0.0],
        [DIM2, DIM5],
        [db + DIM4, DIM5],
        [db + DIM4, DIM1 - DIM6],
        [db + DIM4 + dt, DIM1 - DIM6],
        [db + DIM4 + dt, DIM1],
        [db - dt, DIM1],
        [db - dt, DIM1 - DIM6],
        [db, DIM1 - DIM6],
        [db, DIM5],
        [0, DIM5],
    ]
    geom = geometry.Geometry(Polygon(points), material)

    C = (0.5 * DIM3, 0.5 * DIM1)
    D = (0.5 * DIM3, -0.5 * DIM1)
    E = (-0.5 * DIM3, -0.5 * DIM1)
    F = (-0.5 * DIM3, 0.5 * DIM1)

    geom.recovery_points = [C, D, E, F]
    return geom


def nastran_i1(
    DIM1: float,
    DIM2: float,
    DIM3: float,
    DIM4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a I1 section with the web's middle center at the origin *(0, 0)*, with four
    parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for more details.
    Added by JohnDN90.

    :param float DIM1: Twice distance from web end to flange end
    :param float DIM2: Thickness of web
    :param float DIM3: Length of web (spacing between flanges)
    :param float DIM4: Depth (y) of the I1-section
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a I1 cross-section with a depth of
    5.0 and width of 1.75, and generates a mesh with a maximum triangular area of
    0.02::

        from sectionproperties.pre.library.nastran_sections import nastran_i1

        geom = nastran_i1(DIM1=1.0, DIM2=0.75, DIM3=4.0, DIM4=5.0)
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
    # Ensure dimensions are physically relevant
    np.testing.assert_(DIM4 > DIM3, "Invalid geometry specified.")

    # construct the points and facets
    t = 0.5 * (DIM4 - DIM3)
    points = [
        [0.0, 0.0],
        [DIM1 + DIM2, 0.0],
        [DIM1 + DIM2, t],
        [0.5 * DIM1 + DIM2, t],
        [0.5 * DIM1 + DIM2, t + DIM3],
        [DIM1 + DIM2, t + DIM3],
        [DIM1 + DIM2, DIM4],
        [0.0, DIM4],
        [0.0, t + DIM3],
        [0.5 * DIM1, t + DIM3],
        [0.5 * DIM1, t],
        [0.0, t],
    ]
    geom = geometry.Geometry(Polygon(points), material)

    C = (0.5 * (DIM1 + DIM2), 0.5 * DIM4)
    D = (0.5 * (DIM1 + DIM2), -0.5 * DIM4)
    E = (-0.5 * (DIM1 + DIM2), -0.5 * DIM4)
    F = (-0.5 * (DIM1 + DIM2), 0.5 * DIM4)

    geom.recovery_points = [C, D, E, F]
    return geom


def nastran_l(
    DIM1: float,
    DIM2: float,
    DIM3: float,
    DIM4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a L section with the intersection's center at the origin *(0, 0)*, with four
    parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ for more details.
    Added by JohnDN90.

    :param float DIM1: Width (x) of the L-section
    :param float DIM2: Depth (y) of the L-section
    :param float DIM3: Thickness of flange (horizontal portion)
    :param float DIM4: Thickness of web (vertical portion)
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a L cross-section with a depth of 6.0 and width of 3.0, and
    generates a mesh with a maximum triangular area of 0.01::

        from sectionproperties.pre.library.nastran_sections import nastran_l

        geom = nastran_l(DIM1=3.0, DIM2=6.0, DIM3=0.375, DIM4=0.625)
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
    # Ensure dimensions are physically relevant
    np.testing.assert_(DIM4 < DIM1, "Invalid geometry specified.")
    np.testing.assert_(DIM3 < DIM2, "Invalid geometry specified.")

    # construct the points and facets
    points = [[0, 0], [DIM1, 0], [DIM1, DIM3], [DIM4, DIM3], [DIM4, DIM2], [0, DIM2]]

    geom = geometry.Geometry(Polygon(points), material)

    C = (0.5 * DIM4, DIM2 - 0.5 * DIM3)
    D = (DIM1 - 0.5 * DIM4, -0.5 * DIM3)
    E = (-0.5 * DIM4, -0.5 * DIM3)
    F = (-0.5 * DIM4, DIM2 - 0.5 * DIM3)

    geom.recovery_points = [C, D, E, F]
    return geom


def nastran_rod(
    DIM1: float, n: int, material: pre.Material = pre.DEFAULT_MATERIAL
) -> geometry.Geometry:
    """Constructs a circular rod section with the center at the origin *(0, 0)*, with one parameter
    defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for more details. Added by
    JohnDN90.

    :param float DIM1: Radius of the circular rod section
    :param int n: Number of points discretising the circle
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a circular rod with a radius of 3.0 and 50 points discretising
    the boundary, and generates a mesh with a maximum triangular area of 0.01::

        from sectionproperties.pre.library.nastran_sections import nastran_rod

        geom = nastran_rod(DIM1=3.0, n=50)
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

    # loop through each point on the circle
    d = 2.0 * DIM1
    points = []
    for i in range(n):
        # determine polar angle
        theta = i * 2 * np.pi * 1.0 / n

        # calculate location of the point
        x = 0.5 * d * np.cos(theta)
        y = 0.5 * d * np.sin(theta)

        # append the current point to the points list
        points.append([x, y])

    geom = geometry.Geometry(Polygon(points), material)
    C = (0, DIM1)
    D = (DIM1, 0)
    E = (0, -DIM1)
    F = (-DIM1, 0)

    geom.recovery_points = [C, D, E, F]

    return geom


def nastran_tee(
    DIM1: float,
    DIM2: float,
    DIM3: float,
    DIM4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a T section with the top flange's middle center at the origin *(0, 0)*, with four
    parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ [5]_ for more
    details. Added by JohnDN90.

    :param float DIM1: Width (x) of top flange
    :param float DIM2: Depth (y) of the T-section
    :param float DIM3: Thickness of top flange
    :param float DIM4: Thickness of web
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a T cross-section with a depth of 4.0 and width of 3.0, and
    generates a mesh with a maximum triangular area of 0.001::

        from sectionproperties.pre.library.nastran_sections import nastran_tee

        geom = nastran_tee(DIM1=3.0, DIM2=4.0, DIM3=0.375, DIM4=0.25)
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
    # Ensure dimensions are physically relevant
    np.testing.assert_(DIM4 < DIM1, "Invalid geometry specified.")
    np.testing.assert_(DIM3 < DIM2, "Invalid geometry specified.")

    d = DIM2
    b = DIM1
    t_f = DIM3
    t_w = DIM4
    r = 0
    n_r = 1

    points = []

    # add first two points
    points.append([b * 0.5 - t_w * 0.5, 0])
    points.append([b * 0.5 + t_w * 0.5, 0])

    # construct the top right radius
    pt = [b * 0.5 + t_w * 0.5 + r, d - t_f - r]
    points += draw_radius(pt, r, np.pi, n_r, False)

    # add next four points
    points.append([b, d - t_f])
    points.append([b, d])
    points.append([0, d])
    points.append([0, d - t_f])

    # construct the top left radius
    pt = [b * 0.5 - t_w * 0.5 - r, d - t_f - r]
    points += draw_radius(pt, r, 0.5 * np.pi, n_r, False)

    geom = geometry.Geometry(Polygon(points), material)

    C = (0, 0.5 * DIM3)
    D = (0.5 * DIM1, 0.5 * DIM3)
    E = (0, 0.5 * DIM3 - DIM2)
    F = (-0.5 * DIM1, 0.5 * DIM3)
    geom.recovery_points = [C, D, E, F]

    return geom


def nastran_tee1(
    DIM1: float,
    DIM2: float,
    DIM3: float,
    DIM4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a T1 section with the right flange's middle center at the origin *(0, 0)*, with
    four parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for more
    details. Added by JohnDN90.

    :param float DIM1: Depth (y) of T1-section
    :param float DIM2: Length (x) of web
    :param float DIM3: Thickness of right flange
    :param float DIM4: Thickness of web
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a T1 cross-section with a depth of 3.0 and width of 3.875, and
    generates a mesh with a maximum triangular area of 0.001::

        from sectionproperties.pre.library.nastran_sections import nastran_tee1

        geom = nastran_tee1(DIM1=3.0, DIM2=3.5, DIM3=0.375, DIM4=0.25)
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
    # Ensure dimensions are physically relevant
    np.testing.assert_(DIM4 < DIM1, "Invalid geometry specified.")

    # construct the points and facets
    d1 = (DIM1 - DIM4) / 2.0
    points = [
        [0, 0],
        [DIM3, 0],
        [DIM3, DIM1],
        [0, DIM1],
        [0, d1 + DIM4],
        [-DIM2, d1 + DIM4],
        [-DIM2, d1],
        [0, d1],
    ]

    geom = geometry.Geometry(Polygon(points), material)

    C = (0.5 * DIM3, 0)
    D = (0.5 * DIM3, -0.5 * DIM1)
    E = (-0.5 * DIM3 - DIM2, 0)
    F = (0.5 * DIM3, 0.5 * DIM1)
    geom.recovery_points = [C, D, E, F]

    return geom


def nastran_tee2(
    DIM1: float,
    DIM2: float,
    DIM3: float,
    DIM4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a T2 section with the bottom flange's middle center at the origin *(0, 0)*, with
    four parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for more
    details. Added by JohnDN90.

    :param float DIM1: Width (x) of T2-section
    :param float DIM2: Depth (y) of T2-section
    :param float DIM3: Thickness of bottom flange
    :param float DIM4: Thickness of web
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a T2 cross-section with a depth of 4.0 and width of 3.0, and
    generates a mesh with a maximum triangular area of 0.005::

        from sectionproperties.pre.library.nastran_sections import nastran_tee2

        geom = nastran_tee2(DIM1=3.0, DIM2=4.0, DIM3=0.375, DIM4=0.5)
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

    # Ensure dimensions are physically relevant
    np.testing.assert_(DIM4 < DIM1, "Invalid geometry specified.")
    np.testing.assert_(DIM3 < DIM2, "Invalid geometry specified.")

    # construct the points and facets
    d1 = 0.5 * (DIM1 - DIM4)
    points = [
        [0.0, 0.0],
        [DIM1, 0.0],
        [DIM1, DIM3],
        [DIM1 - d1, DIM3],
        [DIM1 - d1, DIM2],
        [d1, DIM2],
        [d1, DIM3],
        [0, DIM3],
    ]

    geom = geometry.Geometry(Polygon(points), material)
    C = (0.5 * DIM4, DIM2 - 0.5 * DIM3)
    D = (0.5 * DIM1, -0.5 * DIM3)
    E = (-0.5 * DIM1, -0.5 * DIM3)
    F = (-0.5 * DIM4, DIM2 - 0.5 * DIM3)
    geom.recovery_points = [C, D, E, F]
    return geom


def nastran_tube(
    DIM1: float, DIM2: float, n: int, material: pre.Material = pre.DEFAULT_MATERIAL
) -> geometry.Geometry:
    """Constructs a circular tube section with the center at the origin *(0, 0)*, with two
    parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for more
    details. Added by JohnDN90.

    :param float DIM1: Outer radius of the circular tube section
    :param float DIM2: Inner radius of the circular tube section
    :param int n: Number of points discretising the circle
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a circular tube cross-section with an outer radius of 3.0 and an
    inner radius of 2.5, and generates a mesh with 37 points discretising the boundaries and a
    maximum triangular area of 0.01::

        from sectionproperties.pre.library.nastran_sections import nastran_tube

        geom = nastran_tube(DIM1=3.0, DIM2=2.5, n=37)
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
    # Ensure dimensions are physically relevant
    np.testing.assert_(DIM2 < DIM1, "Invalid geometry specified.")

    d = 2.0 * DIM1
    t = DIM1 - DIM2
    points_inner = []
    points_outer = []
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
        points_outer.append([x_outer, y_outer])
        points_inner.append([x_inner, y_inner])

    exterior = geometry.Geometry(Polygon(points_outer), material)
    interior = geometry.Geometry(Polygon(points_inner), material)

    geom = exterior - interior

    C = (0, DIM1)
    D = (DIM1, 0)
    E = (0, -DIM1)
    F = (-DIM1, 0)
    geom.recovery_points = [C, D, E, F]
    return geom


def nastran_tube2(
    DIM1: float, DIM2: float, n: float, material: pre.Material = pre.DEFAULT_MATERIAL
) -> geometry.Geometry:
    """Constructs a circular TUBE2 section with the center at the origin *(0, 0)*, with two
    parameters defining dimensions. See MSC Nastran documentation [1]_ for more details. Added by
    JohnDN90.

    :param float DIM1: Outer radius of the circular tube section
    :param float DIM2: Thickness of wall
    :param int n: Number of points discretising the circle
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a circular TUBE2 cross-section with an outer radius of 3.0 and a
    wall thickness of 0.5, and generates a mesh with 37 point discretising the boundary and a
    maximum triangular area of 0.01::

        from sectionproperties.pre.library.nastran_sections import nastran_tube2

        geom = nastran_tube2(DIM1=3.0, DIM2=0.5, n=37)
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
    # Ensure dimensions are physically relevant
    np.testing.assert_(DIM2 < DIM1, "Invalid geometry specified.")

    d = 2.0 * DIM1
    t = DIM2

    points_inner = []
    points_outer = []

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
        points_outer.append([x_outer, y_outer])
        points_inner.append([x_inner, y_inner])

    exterior = geometry.Geometry(Polygon(points_outer), material)
    interior = geometry.Geometry(Polygon(points_inner), material)
    geom = exterior - interior

    C = (0, DIM1)
    D = (DIM1, 0)
    E = (0, -DIM1)
    F = (-DIM1, 0)
    geom.recovery_points = [C, D, E, F]
    return geom


def nastran_zed(
    DIM1: float,
    DIM2: float,
    DIM3: float,
    DIM4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a Z section with the web's middle center at the origin *(0, 0)*, with four
    parameters defining dimensions. See Nastran documentation [1]_ [2]_ [3]_ [4]_ for more details.
    Added by JohnDN90.

    :param float DIM1: Width (x) of horizontal members
    :param float DIM2: Thickness of web
    :param float DIM3: Spacing between horizontal members (length of web)
    :param float DIM4: Depth (y) of Z-section
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a rectangular cross-section with a depth of 4.0 and width of
    2.75, and generates a mesh with a maximum triangular area of 0.005::

        from sectionproperties.pre.library.nastran_sections import nastran_zed

        geom = nastran_zed(DIM1=1.125, DIM2=0.5, DIM3=3.5, DIM4=4.0)
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
    # Ensure dimensions are physically relevant
    np.testing.assert_(DIM4 > DIM3, "Invalid geometry specified.")

    # construct the points and facets
    t = 0.5 * (DIM4 - DIM3)
    points = [
        [DIM1, 0.0],
        [2.0 * DIM1 + DIM2, 0.0],
        [2.0 * DIM1 + DIM2, t],
        [DIM1 + DIM2, t],
        [DIM1 + DIM2, DIM4],
        [0.0, DIM4],
        [0.0, DIM4 - t],
        [DIM1, DIM4 - t],
    ]
    geom = geometry.Geometry(Polygon(points), material)

    C = (0.5 * DIM2, 0.5 * DIM4)
    D = (0.5 * DIM2 + DIM1, -0.5 * DIM4)
    E = (-0.5 * DIM2, -0.5 * DIM4)
    F = (-0.5 * DIM2 - DIM1, 0.5 * DIM4)

    geom.recovery_points = [C, D, E, F]

    return geom
