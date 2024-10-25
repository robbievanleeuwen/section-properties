"""NASTRAN sections library."""

from __future__ import annotations

import numpy as np
from shapely import Polygon

import sectionproperties.pre.geometry as geometry
import sectionproperties.pre.library.utils as sp_utils
import sectionproperties.pre.pre as pre


def nastran_bar(
    dim_1: float,
    dim_2: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a BAR section.

    Constructs a BAR section with the center at the origin ``(0, 0)``, with two
    parameters defining dimensions.

    Args:
        dim_1: Width (x) of bar
        dim_2: Depth (y) of bar
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        BAR section geometry

    Example:
        The following example creates a BAR cross-section with a depth of 1.5 and width
        of 2.0:

        .. plot::
            :include-source: True
            :caption: Bar section geometry

            from sectionproperties.pre.library import nastran_bar

            nastran_bar(dim_1=2.0, dim_2=1.5).plot_geometry()
    """
    points = [
        (-0.5 * dim_1, -0.5 * dim_2),
        (0.5 * dim_1, -0.5 * dim_2),
        (0.5 * dim_1, 0.5 * dim_2),
        (-0.5 * dim_1, 0.5 * dim_2),
    ]

    geom = geometry.Geometry(geom=Polygon(points), material=material)

    c = (0.5 * dim_1, 0.5 * dim_2)
    d = (0.5 * dim_1, -0.5 * dim_2)
    e = (-0.5 * dim_1, -0.5 * dim_2)
    f = (-0.5 * dim_1, 0.5 * dim_2)
    geom.recovery_points = [c, d, e, f]

    return geom


def nastran_box(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    dim_4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a BOX section.

    Constructs a BOX section with the center at the origin ``(0, 0)``, with four
    parameters defining dimensions.

    Args:
        dim_1: Width (x) of box
        dim_2: Depth (y) of box
        dim_3: Thickness of box in y direction
        dim_4: Thickness of box in x direction
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        RuntimeError: If the geometry generation fails

    Returns:
        BOX section geometry

    Example:
        The following example creates a BOX cross-section with a depth of 3.0 and width
        of 4.0:

        .. plot::
            :include-source: True
            :caption: BOX section geometry

            from sectionproperties.pre.library import nastran_box

            nastran_box(dim_1=4.0, dim_2=3.0, dim_3=0.375, dim_4=0.5).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not 2.0 * dim_4 < dim_1 or not 2.0 * dim_3 < dim_2:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    points_outer = [
        (-0.5 * dim_1, -0.5 * dim_2),
        (0.5 * dim_1, -0.5 * dim_2),
        (0.5 * dim_1, 0.5 * dim_2),
        (-0.5 * dim_1, 0.5 * dim_2),
    ]
    points_inner = [
        (-0.5 * dim_1 + dim_4, -0.5 * dim_2 + dim_3),
        (0.5 * dim_1 - dim_4, -0.5 * dim_2 + dim_3),
        (0.5 * dim_1 - dim_4, 0.5 * dim_2 - dim_3),
        (-0.5 * dim_1 + dim_4, 0.5 * dim_2 - dim_3),
    ]

    inner_box = Polygon(points_inner)
    outer_box = Polygon(points_outer)
    poly_sub = outer_box - inner_box

    c = (0.5 * dim_1, 0.5 * dim_2)
    d = (0.5 * dim_1, -0.5 * dim_2)
    e = (-0.5 * dim_1, -0.5 * dim_2)
    f = (-0.5 * dim_1, 0.5 * dim_2)

    if isinstance(poly_sub, Polygon):
        geom = geometry.Geometry(geom=poly_sub, material=material)
        geom.recovery_points = [c, d, e, f]

        return geom

    msg = "Geometry generation failed."
    raise RuntimeError(msg)


def nastran_box1(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    dim_4: float,
    dim_5: float,
    dim_6: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a BOX1 section.

    Constructs a BOX1 section with the center at the origin ``(0, 0)``, with six
    parameters defining dimensions.

    Args:
        dim_1: Width (x) of box
        dim_2: Depth (y) of box
        dim_3: Thickness of top wall
        dim_4: Thickness of bottom wall
        dim_5: Thickness of left wall
        dim_6: Thickness of right wall
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        RuntimeError: If the geometry generation fails

    Returns:
        BOX1 section geometry

    Example:
        The following example creates a BOX1 cross-section with a depth of 3.0 and width
        of 4.0:

        .. plot::
            :include-source: True
            :caption: BOX1 section geometry

            from sectionproperties.pre.library import nastran_box1

            nastran_box1(
                dim_1=4.0, dim_2=3.0, dim_3=0.375, dim_4=0.5, dim_5=0.25, dim_6=0.75
            ).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not dim_5 + dim_6 < dim_1 or not dim_3 + dim_4 < dim_2:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    exterior_points = [
        (0.0, 0.0),
        (dim_1, 0.0),
        (dim_1, dim_2),
        (0.0, dim_2),
    ]
    interior_points = [
        (dim_6, dim_4),
        (dim_1 - dim_5, dim_4),
        (dim_1 - dim_5, dim_2 - dim_3),
        (dim_6, dim_2 - dim_3),
    ]
    poly_sub = Polygon(exterior_points) - Polygon(interior_points)

    c = (0.5 * dim_1, 0.5 * dim_2)
    d = (0.5 * dim_1, -0.5 * dim_2)
    e = (-0.5 * dim_1, -0.5 * dim_2)
    f = (-0.5 * dim_1, 0.5 * dim_2)

    if isinstance(poly_sub, Polygon):
        geom = geometry.Geometry(geom=poly_sub, material=material)
        geom.recovery_points = [c, d, e, f]

        return geom

    msg = "Geometry generation failed."
    raise RuntimeError(msg)


def nastran_chan(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    dim_4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a CHAN section.

    Constructs a CHAN (C-Channel) section with the web's middle center at the origin
    ``(0, 0)``, with four parameters defining dimensions.

    Args:
        dim_1: Width (x) of the CHAN-section
        dim_2: Depth (y) of the CHAN-section
        dim_3: Thickness of web (vertical portion)
        dim_4: Thickness of flanges (top/bottom portion)
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        CHAN section geometry

    Example:
        The following example creates a CHAN cross-section with a depth of 4.0 and width
        of 2.0:

        .. plot::
            :include-source: True
            :caption: CHAN section geometry

            from sectionproperties.pre.library import nastran_chan

            nastran_chan(dim_1=2.0, dim_2=4.0, dim_3=0.25, dim_4=0.5).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not 2.0 * dim_4 < dim_2 or not dim_3 < dim_1:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    # construct the points
    points = [
        (0.0, 0.0),
        (dim_1, 0.0),
        (dim_1, dim_4),
        (dim_3, dim_4),
        (dim_3, dim_2 - dim_4),
        (dim_1, dim_2 - dim_4),
        (dim_1, dim_2),
        (0.0, dim_2),
    ]

    geom = geometry.Geometry(geom=Polygon(points), material=material)

    c = (0.5 * dim_1, 0.5 * dim_2)
    d = (0.5 * dim_1, -0.5 * dim_2)
    e = (-0.5 * dim_1, -0.5 * dim_2)
    f = (-0.5 * dim_1, 0.5 * dim_2)

    geom.recovery_points = [c, d, e, f]

    return geom


def nastran_chan1(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    dim_4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a CHAN1 section.

    Constructs a CHAN1 (C-Channel) section with the web's middle center at the origin
    ``(0, 0)``, with four parameters defining dimensions.

    Args:
        dim_1: Width (x) of channels
        dim_2: Thickness (x) of web
        dim_3: Spacing between channels (length of web)
        dim_4: Depth (y) of CHAN1-section
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        CHAN1 section geometry

    Example:
        The following example creates a CHAN1 cross-section with a depth of 4.0 and
        width of 1.75:

        .. plot::
            :include-source: True
            :caption: CHAN1 section geometry

            from sectionproperties.pre.library import nastran_chan1

            nastran_chan1(dim_1=0.75, dim_2=1.0, dim_3=3.5, dim_4=4.0).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not dim_4 > dim_3:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    # construct the points and facets
    tf = 0.5 * (dim_4 - dim_3)
    points = [
        (0, 0),
        (dim_1 + dim_2, 0),
        (dim_1 + dim_2, tf),
        (dim_2, tf),
        (dim_2, tf + dim_3),
        (dim_2 + dim_1, tf + dim_3),
        (dim_2 + dim_1, dim_4),
        (0, dim_4),
    ]

    geom = geometry.Geometry(geom=Polygon(points), material=material)

    c = (0.5 * dim_2 + dim_1, 0.5 * dim_4)
    d = (0.5 * dim_2 + dim_1, -0.5 * dim_4)
    e = (-0.5 * dim_2, -0.5 * dim_4)
    f = (-0.5 * dim_2, 0.5 * dim_4)
    geom.recovery_points = [c, d, e, f]

    return geom


def nastran_chan2(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    dim_4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a CHAN2 section.

    Constructs a CHAN2 (C-Channel) section with the bottom web's middle center at the
    origin ``(0, 0)``, with four parameters defining dimensions.

    Args:
        dim_1: Thickness of channels
        dim_2: Thickness of web
        dim_3: Depth (y) of CHAN2-section
        dim_4: Width (x) of CHAN2-section
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        CHAN2 section geometry

    Example:
        The following example creates a CHAN2 cross-section with a depth of 2.0 and
        width of 4.0:

        .. plot::
            :include-source: True
            :caption: CHAN2 section geometry

            from sectionproperties.pre.library import nastran_chan2

            nastran_chan2(dim_1=0.375, dim_2=0.5, dim_3=2.0, dim_4=4.0).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not dim_4 > 2.0 * dim_1 or not dim_3 > dim_2:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    # construct the points and facets
    points = [
        (0.0, 0.0),
        (dim_4, 0.0),
        (dim_4, dim_3),
        (dim_4 - dim_1, dim_3),
        (dim_4 - dim_1, dim_2),
        (dim_1, dim_2),
        (dim_1, dim_3),
        (0.0, dim_3),
    ]

    geom = geometry.Geometry(geom=Polygon(points), material=material)

    c = (0.5 * dim_4, dim_3 - 0.5 * dim_2)
    d = (0.5 * dim_4, -0.5 * dim_2)
    e = (-0.5 * dim_4, -0.5 * dim_2)
    f = (-0.5 * dim_4, dim_3 - 0.5 * dim_2)
    geom.recovery_points = [c, d, e, f]

    return geom


def nastran_cross(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    dim_4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs Nastran's cruciform section.

    Constructs Nastran's cruciform cross-section with the intersection's middle center
    at the origin ``(0, 0)``, with four parameters defining dimensions.

    Args:
        dim_1: Twice the width of horizontal member protruding from the vertical center
            member
        dim_2: Thickness of the vertical member
        dim_3: Depth (y) of the CROSS-section
        dim_4: Thickness of the horizontal members
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Cruciform section geometry

    Example:
        The following example creates a cruciform cross-section with a depth of 3.0:

        .. plot::
            :include-source: True
            :caption: Cruciform section geometry

            from sectionproperties.pre.library import nastran_cross

            nastran_cross(dim_1=1.5, dim_2=0.375, dim_3=3.0, dim_4=0.25).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not dim_4 < dim_3:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    # construct the points and facets
    d = 0.5 * (dim_3 - dim_4)
    points = [
        (0.5 * dim_1, 0.0),
        (0.5 * dim_1 + dim_2, 0.0),
        (0.5 * dim_1 + dim_2, d),
        (dim_1 + dim_2, d),
        (dim_1 + dim_2, d + dim_4),
        (0.5 * dim_1 + dim_2, d + dim_4),
        (0.5 * dim_1 + dim_2, dim_3),
        (0.5 * dim_1, dim_3),
        (0.5 * dim_1, d + dim_4),
        (0.0, d + dim_4),
        (0.0, d),
        (0.5 * dim_1, d),
    ]

    geom = geometry.Geometry(geom=Polygon(points), material=material)

    c = (0.0, 0.5 * dim_3)
    d_r = (0.5 * (dim_1 + dim_2), 0.0)
    e = (0.0, -0.5 * dim_3)
    f = (-0.5 * (dim_1 + dim_2), 0.0)
    geom.recovery_points = [c, d_r, e, f]

    return geom


def nastran_fcross(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    dim_4: float,
    dim_5: float,
    dim_6: float,
    dim_7: float,
    dim_8: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a flanged cruciform section.

    Constructs a flanged cruciform cross-section with the intersection's middle
    center at the origin ``(0, 0)``, with eight parameters defining dimensions

    Args:
        dim_1: Depth (y) of flanged cruciform
        dim_2: Width (x) of flanged cruciform
        dim_3: Thickness of vertical web
        dim_4: Thickness of horizontal web
        dim_5: Length of flange attached to vertical web
        dim_6: Thickness of flange attached to vertical web
        dim_7: Length of flange attached to horizontal web
        dim_8: Thickness of flange attached to horizontal web
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Flanged cruciform section geometry

    Example:
        The following example demonstrates the creation of a flanged cross section:

        .. plot::
            :include-source: True
            :caption: Flanged cruciform section geometry

            from sectionproperties.pre.library import nastran_fcross

            nastran_fcross(
                dim_1=9.0, dim_2=6.0, dim_3=0.75, dim_4=0.625, dim_5=2.1, dim_6=0.375,
                dim_7=4.5, dim_8=0.564
            ).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if (
        not dim_5 > dim_3
        or not dim_7 > dim_4
        or not dim_7 < dim_1
        or not dim_5 < dim_2
        or not dim_8 < (0.5 * dim_2 - 0.5 * dim_3)
        or not dim_6 < (0.5 * dim_1 - 0.5 * dim_4)
    ):
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    # construct the points and facets
    points = [
        (0.5 * dim_3, -0.5 * dim_4),
        (0.5 * dim_2 - dim_8, -0.5 * dim_4),
        (0.5 * dim_2 - dim_8, -0.5 * dim_7),
        (0.5 * dim_2, -0.5 * dim_7),
        (0.5 * dim_2, 0.5 * dim_7),
        (0.5 * dim_2 - dim_8, 0.5 * dim_7),
        (0.5 * dim_2 - dim_8, 0.5 * dim_4),
        (0.5 * dim_3, 0.5 * dim_4),
        (0.5 * dim_3, 0.5 * dim_1 - dim_6),
        (0.5 * dim_5, 0.5 * dim_1 - dim_6),
        (0.5 * dim_5, 0.5 * dim_1),
        (-0.5 * dim_5, 0.5 * dim_1),
        (-0.5 * dim_5, 0.5 * dim_1 - dim_6),
        (-0.5 * dim_3, 0.5 * dim_1 - dim_6),
        (-0.5 * dim_3, 0.5 * dim_4),
        (-0.5 * dim_2 + dim_8, 0.5 * dim_4),
        (-0.5 * dim_2 + dim_8, 0.5 * dim_7),
        (-0.5 * dim_2, 0.5 * dim_7),
        (-0.5 * dim_2, -0.5 * dim_7),
        (-0.5 * dim_2 + dim_8, -0.5 * dim_7),
        (-0.5 * dim_2 + dim_8, -0.5 * dim_4),
        (-0.5 * dim_3, -0.5 * dim_4),
        (-0.5 * dim_3, -0.5 * dim_1 + dim_6),
        (-0.5 * dim_5, -0.5 * dim_1 + dim_6),
        (-0.5 * dim_5, -0.5 * dim_1),
        (0.5 * dim_5, -0.5 * dim_1),
        (0.5 * dim_5, -0.5 * dim_1 + dim_6),
        (0.5 * dim_3, -0.5 * dim_1 + dim_6),
    ]

    geom = geometry.Geometry(geom=Polygon(points), material=material)

    c = (0.0, 0.5 * dim_1)
    d = (0.5 * dim_2, 0.0)
    e = (0.0, -0.5 * dim_1)
    f = (-0.5 * dim_2, 0.0)
    geom.recovery_points = [c, d, e, f]

    return geom


def nastran_dbox(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    dim_4: float,
    dim_5: float,
    dim_6: float,
    dim_7: float,
    dim_8: float,
    dim_9: float,
    dim_10: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a DBOX section.

    Constructs a DBOX section with the center at the origin ``(0, 0)``, with ten
    parameters defining dimensions.

    Args:
        dim_1: Width (x) of the DBOX-section
        dim_2: Depth (y) of the DBOX-section
        dim_3: Width (x) of left-side box
        dim_4: Thickness of left wall
        dim_5: Thickness of center wall
        dim_6: Thickness of right wall
        dim_7: Thickness of top left wall
        dim_8: Thickness of bottom left wall
        dim_9: Thickness of top right wall
        dim_10: Thickness of bottom right wall
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        RuntimeError: If the geometry generation fails

    Returns:
        DBOX section geometry

    Example:
        The following example creates a DBOX cross-section with a depth of 3.0 and width
        of 8.0:

        .. plot::
            :include-source: True
            :caption: DBOX section geometry

            from sectionproperties.pre.library import nastran_dbox

            nastran_dbox(
                dim_1=8.0, dim_2=3.0, dim_3=3.0, dim_4=0.5, dim_5=0.625, dim_6=0.75,
                dim_7=0.375, dim_8=0.25, dim_9=0.5, dim_10=0.375
            ).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if (
        not (dim_4 + dim_5 + dim_6) < dim_1
        or not (dim_4 + 0.5 * dim_5) < dim_3
        or not (dim_7 + dim_8) < dim_2
        or not (dim_9 + dim_10) < dim_2
    ):
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    # construct the points and facets
    exterior_points = [
        (0.0, 0.0),
        (dim_1, 0.0),
        (dim_1, dim_2),
        (0.0, dim_2),
    ]
    interior_points_1 = [
        (dim_4, dim_8),
        (dim_3 - dim_5 / 2.0, dim_8),
        (dim_3 - dim_5 / 2.0, dim_2 - dim_7),
        (dim_4, dim_2 - dim_7),
    ]
    interior_points_2 = [
        (dim_3 + dim_5 / 2.0, dim_10),
        (dim_1 - dim_6, dim_10),
        (dim_1 - dim_6, dim_2 - dim_9),
        (dim_3 + dim_5 / 2.0, dim_2 - dim_9),
    ]
    poly_sub = (
        Polygon(exterior_points)
        - Polygon(interior_points_1)
        - Polygon(interior_points_2)
    )

    c = (0.5 * dim_1, 0.5 * dim_2)
    d = (0.5 * dim_1, -0.5 * dim_2)
    e = (-0.5 * dim_1, -0.5 * dim_2)
    f = (-0.5 * dim_1, 0.5 * dim_2)

    if isinstance(poly_sub, Polygon):
        geom = geometry.Geometry(geom=poly_sub, material=material)
        geom.recovery_points = [c, d, e, f]

        return geom

    msg = "Geometry generation failed."
    raise RuntimeError(msg)


def nastran_gbox(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    dim_4: float,
    dim_5: float,
    dim_6: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a GBOX section.

    Constructs a GBOX section with the center at the origin ``(0, 0)``, with six
    parameters defining dimensions.

    Args:
        dim_1: Width (x) of the GBOX-section
        dim_2: Depth (y) of the GBOX-section
        dim_3: Thickness of top flange
        dim_4: Thickness of bottom flange
        dim_5: Thickness of webs
        dim_6: Spacing between webs
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        RuntimeError: If the geometry generation fails

    Returns:
        GBOX section geometry

    Example:
        The following example creates a GBOX cross-section with a depth of 2.5 and width
        of 6.0:

        .. plot::
            :include-source: True
            :caption: GBOX section geometry

            from sectionproperties.pre.library import nastran_gbox

            nastran_gbox(
                dim_1=6.0, dim_2=2.5, dim_3=0.375, dim_4=0.25, dim_5=0.625, dim_6=1.0
            ).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not (dim_3 + dim_4) < dim_2 or not (2.0 * dim_5 + dim_6) < dim_1:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    # construct the points and facets
    d = 0.5 * (dim_1 - dim_6 - 2.0 * dim_5)
    exterior_points = [
        (0.0, 0.0),
        (dim_1, 0.0),
        (dim_1, dim_4),
        (d + 2.0 * dim_5 + dim_6, dim_4),
        (d + 2.0 * dim_5 + dim_6, dim_2 - dim_3),
        (dim_1, dim_2 - dim_3),
        (dim_1, dim_2),
        (0.0, dim_2),
        (0.0, dim_2 - dim_3),
        (d, dim_2 - dim_3),
        (d, dim_4),
        (0.0, dim_4),
    ]
    interior_points = [
        (d + dim_5, dim_4),
        (d + dim_5 + dim_6, dim_4),
        (d + dim_5 + dim_6, dim_2 - dim_3),
        (d + dim_5, dim_2 - dim_3),
    ]
    poly_sub = Polygon(exterior_points) - Polygon(interior_points)

    c = (0.5 * dim_1, 0.5 * dim_2)
    d_r = (0.5 * dim_1, -0.5 * dim_2)
    e = (-0.5 * dim_1, -0.5 * dim_2)
    f = (-0.5 * dim_1, 0.5 * dim_2)

    if isinstance(poly_sub, Polygon):
        geom = geometry.Geometry(geom=poly_sub, material=material)
        geom.recovery_points = [c, d_r, e, f]

        return geom

    msg = "Geometry generation failed."
    raise RuntimeError(msg)


def nastran_h(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    dim_4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a H section.

    Constructs a H section with the middle web's middle center at the origin ``(0, 0)``,
    with four parameters defining dimensions.

    Args:
        dim_1: Spacing between vertical flanges (length of web)
        dim_2: Twice the thickness of the vertical flanges
        dim_3: Depth (y) of the H-section
        dim_4: Thickness of the middle web
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        H section geometry

    Example:
        The following example creates a H cross-section with a depth of 3.5 and width of
        2.75:

        .. plot::
            :include-source: True
            :caption: H section geometry

            from sectionproperties.pre.library import nastran_h

            nastran_h(
                dim_1=2.0, dim_2=0.75, dim_3=3.5, dim_4=0.2
            ).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not dim_4 < dim_3:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    d1 = 0.5 * (dim_3 - dim_4)
    d2 = 0.5 * dim_2

    # construct the points and facets
    points = [
        (0, 0),
        (d2, 0),
        (d2, d1),
        (d2 + dim_1, d1),
        (d2 + dim_1, 0),
        (dim_1 + dim_2, 0),
        (dim_1 + dim_2, dim_3),
        (dim_1 + dim_2 - d2, dim_3),
        (dim_1 + dim_2 - d2, d1 + dim_4),
        (d2, d1 + dim_4),
        (d2, dim_3),
        (0, dim_3),
    ]

    geom = geometry.Geometry(geom=Polygon(points), material=material)

    c = (0.5 * (dim_1 + dim_2), 0.5 * dim_3)
    d = (0.5 * (dim_1 + dim_2), -0.5 * dim_3)
    e = (-0.5 * (dim_1 + dim_2), -0.5 * dim_3)
    f = (-0.5 * (dim_1 + dim_2), 0.5 * dim_3)
    geom.recovery_points = [c, d, e, f]

    return geom


def nastran_hat(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    dim_4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a HAT section.

    Constructs a HAT section with the top most section's middle center at the origin
    ``(0, 0)``, with four parameters defining dimensions.

    Args:
        dim_1: Depth (y) of HAT-section
        dim_2: Thickness of HAT-section
        dim_3: Width (x) of top most section
        dim_4: Width (x) of bottom sections
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        HAT section geometry

    Example:
        The following example creates a HAT cross-section with a depth of 1.25 and width
        of 2.5:

        .. plot::
            :include-source: True
            :caption: HAT section geometry

            from sectionproperties.pre.library import nastran_hat

            nastran_hat(dim_1=1.25, dim_2=0.25, dim_3=1.5, dim_4=0.5).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not 2.0 * dim_2 < dim_1:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    # construct the points and facets
    points = [
        (0.0, 0.0),
        (dim_4 + dim_2, 0.0),
        (dim_4 + dim_2, dim_1 - dim_2),
        (dim_4 + dim_3 - dim_2, dim_1 - dim_2),
        (dim_4 + dim_3 - dim_2, 0.0),
        (2 * dim_4 + dim_3, 0.0),
        (2.0 * dim_4 + dim_3, dim_2),
        (dim_4 + dim_3, dim_2),
        (dim_4 + dim_3, dim_1),
        (dim_4, dim_1),
        (dim_4, dim_2),
        (0.0, dim_2),
    ]

    geom = geometry.Geometry(geom=Polygon(points), material=material)

    c = (0.5 * dim_3, 0.5 * dim_2)
    d = (0.5 * dim_3 + dim_4, -dim_1 + dim_2)
    e = (-0.5 * dim_3 - dim_4, -dim_1 + dim_2)
    f = (-0.5 * dim_3, 0.5 * dim_2)
    geom.recovery_points = [c, d, e, f]

    return geom


def nastran_hat1(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    dim_4: float,
    dim_5: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a HAT1 section.

    Constructs a HAT1 section with the bottom plate's bottom center at the origin
    ``(0, 0)``, with five parameters defining dimensions.

    Args:
        dim_1: Width(x) of the HAT1-section
        dim_2: Depth (y) of the HAT1-section
        dim_3: Width (x) of hat's top flange
        dim_4: Thickness of hat stiffener
        dim_5: Thicknesss of bottom plate
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        HAT1 section geometry

    Example:
        The following example creates a HAT1 cross-section with a depth of 2.0 and width
        of 4.0:

        .. plot::
            :include-source: True
            :caption: HAT1 section geometry

            from sectionproperties.pre.library import nastran_hat1

            nastran_hat1(
                dim_1=4.0, dim_2=2.0, dim_3=1.5, dim_4=0.1875, dim_5=0.375
            ).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not (2.0 * dim_4 + dim_5) < dim_2 or not dim_3 < dim_1:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    # create bottom rectangular plate
    bottom_plate = nastran_bar(
        dim_1=dim_1,
        dim_2=dim_5,
        material=material,
    ).shift_section(y_offset=dim_5 / 2)

    # create the hat stiffener
    d1 = dim_2 - dim_5
    d2 = dim_4
    d3 = dim_3
    d4 = 0.5 * (dim_1 - dim_3)

    hat = nastran_hat(dim_1=d1, dim_2=d2, dim_3=d3, dim_4=d4, material=material)
    # Merge the two sections into one geometry
    geom = (
        hat.align_center(align_to=bottom_plate).align_to(other=bottom_plate, on="top")
        + bottom_plate
    )

    c = (-0.5 * dim_1, 0.0)
    d = (0.5 * dim_1, 0.0)
    e = (-0.5 * dim_3, dim_2)
    f = (0.5 * dim_3, dim_2)

    geom.recovery_points = [c, d, e, f]

    return geom


def nastran_hexa(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a HEXA section.

    Constructs a HEXA (hexagon) section with the center at the origin ``(0, 0)``, with
    three parameters defining dimensions.

    Args:
        dim_1: Spacing between bottom right point and right most point
        dim_2: Width (x) of hexagon
        dim_3: Depth (y) of hexagon
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        HEXA section geometry

    Example:
        The following example creates a HEXA cross-section with a depth of 1.5 and width
        of 2.0:

        .. plot::
            :include-source: True
            :caption: HEXA section geometry

            from sectionproperties.pre.library import nastran_hexa

            nastran_hexa(dim_1=0.5, dim_2=2.0, dim_3=1.5).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not dim_2 > dim_1:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    # construct the points and facets
    points = [
        (dim_1, 0.0),
        (dim_2 - dim_1, 0.0),
        (dim_2, 0.5 * dim_3),
        (dim_2 - dim_1, dim_3),
        (dim_1, dim_3),
        (0.0, 0.5 * dim_3),
    ]

    geom = geometry.Geometry(geom=Polygon(points), material=material)

    c = 0.0, 0.5 * dim_3
    d = 0.0, -0.5 * dim_3
    e = 0.5 * dim_2, 0.0
    f = -0.5 * dim_2, 0.0
    geom.recovery_points = [c, d, e, f]

    return geom


def nastran_i(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    dim_4: float,
    dim_5: float,
    dim_6: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs Nastran's I section.

    Constructs Nastran's I section with the bottom flange's middle center at the origin
    ``(0, 0)``, with six parameters defining dimensions.

    Args:
        dim_1: Depth(y) of the I Section
        dim_2: Width (x) of bottom flange
        dim_3: Width (x) of top flange
        dim_4: Thickness of web
        dim_5: Thickness of bottom web
        dim_6: Thickness of top web
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        I section geometry

    Example:
        The following example creates an I cross-section with a depth of 5.0:

        .. plot::
            :include-source: True
            :caption: I section geometry

            from sectionproperties.pre.library import nastran_i

            nastran_i(
                dim_1=5.0, dim_2=2.0, dim_3=3.0, dim_4=0.25, dim_5=0.375, dim_6=0.5
            ).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not (dim_5 + dim_6) < dim_1 or not dim_4 < dim_3 or not dim_4 < dim_2:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    # construct the points and facets
    db = 0.5 * (dim_2 - dim_4)
    dt = 0.5 * (dim_3 - dim_4)
    points = [
        (0.0, 0.0),
        (dim_2, 0.0),
        (dim_2, dim_5),
        (db + dim_4, dim_5),
        (db + dim_4, dim_1 - dim_6),
        (db + dim_4 + dt, dim_1 - dim_6),
        (db + dim_4 + dt, dim_1),
        (db - dt, dim_1),
        (db - dt, dim_1 - dim_6),
        (db, dim_1 - dim_6),
        (db, dim_5),
        (0, dim_5),
    ]
    geom = geometry.Geometry(geom=Polygon(points), material=material)

    c = (0.5 * dim_3, 0.5 * dim_1)
    d = (0.5 * dim_3, -0.5 * dim_1)
    e = (-0.5 * dim_3, -0.5 * dim_1)
    f = (-0.5 * dim_3, 0.5 * dim_1)
    geom.recovery_points = [c, d, e, f]

    return geom


def nastran_i1(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    dim_4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs an I1 section.

    Constructs an I1 section with the web's middle center at the origin ``(0, 0)``, with
    four parameters defining dimensions.

    Args:
        dim_1: Twice distance from web end to flange end
        dim_2: Thickness of web
        dim_3: Length of web (spacing between flanges)
        dim_4: Depth (y) of the I1-section
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        I1 section geometry

    Example:
        The following example creates an I1 cross-section with a depth of 5.0 and width
        of 1.75:

        .. plot::
            :include-source: True
            :caption: I1 section geometry

            from sectionproperties.pre.library import nastran_i1

            nastran_i1(dim_1=1.0, dim_2=0.75, dim_3=4.0, dim_4=5.0).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not dim_4 > dim_3:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    # construct the points and facets
    t = 0.5 * (dim_4 - dim_3)
    points = [
        (0.0, 0.0),
        (dim_1 + dim_2, 0.0),
        (dim_1 + dim_2, t),
        (0.5 * dim_1 + dim_2, t),
        (0.5 * dim_1 + dim_2, t + dim_3),
        (dim_1 + dim_2, t + dim_3),
        (dim_1 + dim_2, dim_4),
        (0.0, dim_4),
        (0.0, t + dim_3),
        (0.5 * dim_1, t + dim_3),
        (0.5 * dim_1, t),
        (0.0, t),
    ]

    geom = geometry.Geometry(geom=Polygon(points), material=material)

    c = (0.5 * (dim_1 + dim_2), 0.5 * dim_4)
    d = (0.5 * (dim_1 + dim_2), -0.5 * dim_4)
    e = (-0.5 * (dim_1 + dim_2), -0.5 * dim_4)
    f = (-0.5 * (dim_1 + dim_2), 0.5 * dim_4)
    geom.recovery_points = [c, d, e, f]

    return geom


def nastran_l(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    dim_4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs an L section.

    Constructs an L section with the intersection's center at the origin ``(0, 0)``,with
    four parameters defining dimensions.

    Args:
        dim_1: Width (x) of the L-section
        dim_2: Depth (y) of the L-section
        dim_3: Thickness of flange (horizontal portion)
        dim_4: Thickness of web (vertical portion)
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        L section geometry

    Example:
        The following example creates an L cross-section with a depth of 6.0 and width
        of 3.0:

        .. plot::
            :include-source: True
            :caption: L section geometry

            from sectionproperties.pre.library import nastran_l

            nastran_l(dim_1=3.0, dim_2=6.0, dim_3=0.375, dim_4=0.625).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not dim_4 < dim_1 or not dim_3 < dim_2:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    # construct the points and facets
    points = [
        (0, 0),
        (dim_1, 0),
        (dim_1, dim_3),
        (dim_4, dim_3),
        (dim_4, dim_2),
        (0, dim_2),
    ]

    geom = geometry.Geometry(geom=Polygon(points), material=material)

    c = (0.5 * dim_4, dim_2 - 0.5 * dim_3)
    d = (dim_1 - 0.5 * dim_4, -0.5 * dim_3)
    e = (-0.5 * dim_4, -0.5 * dim_3)
    f = (-0.5 * dim_4, dim_2 - 0.5 * dim_3)
    geom.recovery_points = [c, d, e, f]

    return geom


def nastran_rod(
    dim_1: float,
    n: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a circular rod section.

    Constructs a circular rod section with the center at the origin ``(0, 0)``, with one
    parameter defining dimension.

    Args:
        dim_1: Radius of the circular rod section
        n: Number of points discretising the circle
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Circular rod section geometry

    Example:
        The following example creates a circular rod with a radius of 3.0 and 50 points
        discretising the boundary:

        .. plot::
            :include-source: True
            :caption: Circular rod section geometry

            from sectionproperties.pre.library import nastran_rod

            nastran_rod(dim_1=3.0, n=50).plot_geometry()
    """
    # loop through each point on the circle
    d = 2.0 * dim_1
    points: list[tuple[float, float]] = []

    for i in range(n):
        # determine polar angle
        theta = i * 2 * np.pi * 1.0 / n

        # calculate location of the point
        x = 0.5 * d * np.cos(theta)
        y = 0.5 * d * np.sin(theta)

        # append the current point to the points list
        points.append((x, y))

    geom = geometry.Geometry(geom=Polygon(points), material=material)

    c = (0.0, dim_1)
    d_r = (dim_1, 0.0)
    e = (0.0, -dim_1)
    f = (-dim_1, 0.0)
    geom.recovery_points = [c, d_r, e, f]

    return geom


def nastran_tee(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    dim_4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a T section.

    Constructs a T section with the top flange's middle center at the origin ``(0, 0)``,
    with four parameters defining dimensions.

    Args:
        dim_1: Width (x) of top flange
        dim_2: Depth (y) of the T-section
        dim_3: Thickness of top flange
        dim_4: Thickness of web
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        T section geometry

    Example:
        The following example creates a T cross-section with a depth of 4.0 and width of
        3.0:

        .. plot::
            :include-source: True
            :caption: T section geometry

            from sectionproperties.pre.library import nastran_tee

            nastran_tee(dim_1=3.0, dim_2=4.0, dim_3=0.375, dim_4=0.25).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not dim_4 < dim_1 or not dim_3 < dim_2:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    d = dim_2
    b = dim_1
    t_f = dim_3
    t_w = dim_4
    r = 0.0
    n_r = 1

    points: list[tuple[float, float]] = []

    # add first two points
    points.append((b * 0.5 - t_w * 0.5, 0))
    points.append((b * 0.5 + t_w * 0.5, 0))

    # construct the top right radius
    pt = b * 0.5 + t_w * 0.5 + r, d - t_f - r
    points += sp_utils.draw_radius(pt=pt, r=r, theta=np.pi, n=n_r, ccw=False)

    # add next four points
    points.append((b, d - t_f))
    points.append((b, d))
    points.append((0, d))
    points.append((0, d - t_f))

    # construct the top left radius
    pt = b * 0.5 - t_w * 0.5 - r, d - t_f - r
    points += sp_utils.draw_radius(pt=pt, r=r, theta=0.5 * np.pi, n=n_r, ccw=False)

    geom = geometry.Geometry(geom=Polygon(points), material=material)

    c = (0.0, 0.5 * dim_3)
    d_r = (0.5 * dim_1, 0.5 * dim_3)
    e = (0.0, 0.5 * dim_3 - dim_2)
    f = (-0.5 * dim_1, 0.5 * dim_3)
    geom.recovery_points = [c, d_r, e, f]

    return geom


def nastran_tee1(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    dim_4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a T1 section.

    Constructs a T1 section with the right flange's middle center at the origin
    ``(0, 0)``, with four parameters defining dimensions.

    Args:
        dim_1: Depth (y) of T1-section
        dim_2: Length (x) of web
        dim_3: Thickness of right flange
        dim_4: Thickness of web
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        T1 section geometry

    Example:
        The following example creates a T1 cross-section with a depth of 3.0 and width
        of 3.875:

        .. plot::
            :include-source: True
            :caption: T1 section geometry

            from sectionproperties.pre.library import nastran_tee1

            nastran_tee1(dim_1=3.0, dim_2=3.5, dim_3=0.375, dim_4=0.25).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not dim_4 < dim_1:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    # construct the points and facets
    d1 = (dim_1 - dim_4) / 2.0
    points = [
        (0, 0),
        (dim_3, 0),
        (dim_3, dim_1),
        (0, dim_1),
        (0, d1 + dim_4),
        (-dim_2, d1 + dim_4),
        (-dim_2, d1),
        (0, d1),
    ]

    geom = geometry.Geometry(geom=Polygon(points), material=material)

    c = (0.5 * dim_3, 0)
    d = (0.5 * dim_3, -0.5 * dim_1)
    e = (-0.5 * dim_3 - dim_2, 0)
    f = (0.5 * dim_3, 0.5 * dim_1)
    geom.recovery_points = [c, d, e, f]

    return geom


def nastran_tee2(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    dim_4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a T2 section.

    Constructs a T2 section with the bottom flange's middle center at the origin
    ``(0, 0)``, with four parameters defining dimensions.

    Args:
        dim_1: Width (x) of T2-section
        dim_2: Depth (y) of T2-section
        dim_3: Thickness of bottom flange
        dim_4: Thickness of web
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        T2 section geometry

    Example:
        The following example creates a T2 cross-section with a depth of 4.0 and width
        of 3.0:

        .. plot::
            :include-source: True
            :caption: T2 section geometry

            from sectionproperties.pre.library import nastran_tee2

            nastran_tee2(dim_1=3.0, dim_2=4.0, dim_3=0.375, dim_4=0.5).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not dim_4 < dim_1 or not dim_3 < dim_2:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    # construct the points and facets
    d1 = 0.5 * (dim_1 - dim_4)
    points = [
        (0.0, 0.0),
        (dim_1, 0.0),
        (dim_1, dim_3),
        (dim_1 - d1, dim_3),
        (dim_1 - d1, dim_2),
        (d1, dim_2),
        (d1, dim_3),
        (0, dim_3),
    ]

    geom = geometry.Geometry(geom=Polygon(points), material=material)

    c = (0.5 * dim_4, dim_2 - 0.5 * dim_3)
    d = (0.5 * dim_1, -0.5 * dim_3)
    e = (-0.5 * dim_1, -0.5 * dim_3)
    f = (-0.5 * dim_4, dim_2 - 0.5 * dim_3)
    geom.recovery_points = [c, d, e, f]

    return geom


def nastran_tube(
    dim_1: float,
    dim_2: float,
    n: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a circular tube section.

    Constructs a circular tube section with the center at the origin ``(0, 0)``, with
    two parameters defining dimensions.

    Args:
        dim_1: Outer radius of the circular tube section
        dim_2: Inner radius of the circular tube section
        n: Number of points discretising the circle
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Circular tube section geometry

    Example:
        The following example creates a circular tube cross-section with an outer radius
        of 3.0 and an inner radius of 2.5:

        .. plot::
            :include-source: True
            :caption: Circular tube section geometry

            from sectionproperties.pre.library import nastran_tube

            nastran_tube(dim_1=3.0, dim_2=2.5, n=37).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not dim_2 < dim_1:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    d = 2.0 * dim_1
    t = dim_1 - dim_2
    points_inner: list[tuple[float, float]] = []
    points_outer: list[tuple[float, float]] = []

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
        points_outer.append((x_outer, y_outer))
        points_inner.append((x_inner, y_inner))

    exterior = geometry.Geometry(geom=Polygon(points_outer), material=material)
    interior = geometry.Geometry(geom=Polygon(points_inner), material=material)
    geom = exterior - interior

    c = (0.0, dim_1)
    d_r = (dim_1, 0.0)
    e = (0.0, -dim_1)
    f = (-dim_1, 0.0)
    geom.recovery_points = [c, d_r, e, f]

    return geom


def nastran_tube2(
    dim_1: float,
    dim_2: float,
    n: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a circular TUBE2 section.

    Constructs a circular TUBE2 section with the center at the origin ``(0, 0)``, with
    two parameters defining dimensions.

    Args:
        dim_1: Outer radius of the circular tube section
        dim_2: Thickness of wall
        n: Number of points discretising the circle
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        TUBE2 section geometry

    Example:
        The following example creates a circular TUBE2 cross-section with an outer
        radius of 3.0 and a wall thickness of 0.5:

        .. plot::
            :include-source: True
            :caption: TUBE2 section geometry

            from sectionproperties.pre.library import nastran_tube2

            nastran_tube2(dim_1=3.0, dim_2=0.5, n=37).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not dim_2 < dim_1:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    d = 2.0 * dim_1
    t = dim_2

    points_inner: list[tuple[float, float]] = []
    points_outer: list[tuple[float, float]] = []

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
        points_outer.append((x_outer, y_outer))
        points_inner.append((x_inner, y_inner))

    exterior = geometry.Geometry(geom=Polygon(points_outer), material=material)
    interior = geometry.Geometry(geom=Polygon(points_inner), material=material)
    geom = exterior - interior

    c = (0.0, dim_1)
    d_r = (dim_1, 0.0)
    e = (0.0, -dim_1)
    f = (-dim_1, 0.0)
    geom.recovery_points = [c, d_r, e, f]

    return geom


def nastran_zed(
    dim_1: float,
    dim_2: float,
    dim_3: float,
    dim_4: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a Z section.

    Constructs a Z section with the web's middle center at the origin ``(0, 0)``, with
    four parameters defining dimensions.

    Args:
        dim_1: Width (x) of horizontal members
        dim_2: Thickness of web
        dim_3: Spacing between horizontal members (length of web)
        dim_4: Depth (y) of Z-section
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Z section geometry

    Example:
        The following example creates a Z cross-section with a depth of 4.0 and width of
        1.125:

        .. plot::
            :include-source: True
            :caption: Z section geometry

            from sectionproperties.pre.library import nastran_zed

            nastran_zed(dim_1=1.125, dim_2=0.5, dim_3=3.5, dim_4=4.0).plot_geometry()
    """
    # Ensure dimensions are physically relevant
    if not dim_4 > dim_3:
        msg = "Invalid geometry specified."
        raise ValueError(msg)

    # construct the points and facets
    t = 0.5 * (dim_4 - dim_3)
    points = [
        (dim_1, 0.0),
        (2.0 * dim_1 + dim_2, 0.0),
        (2.0 * dim_1 + dim_2, t),
        (dim_1 + dim_2, t),
        (dim_1 + dim_2, dim_4),
        (0.0, dim_4),
        (0.0, dim_4 - t),
        (dim_1, dim_4 - t),
    ]

    geom = geometry.Geometry(geom=Polygon(points), material=material)

    c = (0.5 * dim_2, 0.5 * dim_4)
    d = (0.5 * dim_2 + dim_1, -0.5 * dim_4)
    e = (-0.5 * dim_2, -0.5 * dim_4)
    f = (-0.5 * dim_2 - dim_1, 0.5 * dim_4)

    geom.recovery_points = [c, d, e, f]

    return geom
