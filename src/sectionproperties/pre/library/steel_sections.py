"""Steel sections library."""

from __future__ import annotations

import numpy as np
from shapely import Polygon

import sectionproperties.pre.geometry as geometry
import sectionproperties.pre.library.utils as sp_utils
import sectionproperties.pre.pre as pre


def circular_hollow_section(
    d: float,
    t: float,
    n: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a circular hollow section (CHS).

    Constructs a circular hollow section (CHS) centered at the origin ``(0, 0)``, with
    diameter ``d`` and thickness ``t``, using ``n`` points to construct the inner and
    outer circles.

    Args:
        d: Outer diameter of the CHS
        t: Thickness of the CHS
        n: Number of points discretising the inner and outer circles
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        RuntimeError: If the geometry generation fails

    Returns:
        Circular hollow section geometry

    Example:
        The following example creates a CHS discretised with 64 points, with a diameter
        of 48 mm and thickness of 3.2mm:

        .. plot::
            :include-source: True
            :caption: Circular hollow section geometry

            from sectionproperties.pre.library import circular_hollow_section

            circular_hollow_section(d=48, t=3.2, n=64).plot_geometry()
    """
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

    inner_circle = Polygon(points_inner)
    outer_circle = Polygon(points_outer)
    poly_sub = outer_circle - inner_circle

    if isinstance(poly_sub, Polygon):
        return geometry.Geometry(geom=poly_sub, material=material)

    msg = "Geometry generation failed."
    raise RuntimeError(msg)


def elliptical_hollow_section(
    d_x: float,
    d_y: float,
    t: float,
    n: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs an elliptical hollow section (EHS).

    Constructs an elliptical hollow section (EHS) centered at the origin ``(0, 0)``,
    with outer horizontal diameter ``d_x``, outer vertical diameter ``d_y`` and
    thickness ``t``, using ``n`` points to construct the inner and outer ellipses.

    Args:
        d_x: Diameter of the ellipse in the x-dimension
        d_y: Diameter of the ellipse in the y-dimension
        t: Thickness of the EHS
        n: Number of points discretising the inner and outer ellipses
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        RuntimeError: If the geometry generation fails

    Returns:
        Elliptical hollow section geometry

    Example:
        The following example creates an EHS discretised with 30 points, with an outer
        horizontal diameter of 50 mm, outer vertical diameter of 25 mm and thickness of
        2 mm:

        .. plot::
            :include-source: True
            :caption: Elliptical hollow section geometry

            from sectionproperties.pre.library import elliptical_hollow_section

            elliptical_hollow_section(d_x=50, d_y=25, t=2, n=64).plot_geometry()
    """
    points_inner: list[tuple[float, float]] = []
    points_outer: list[tuple[float, float]] = []

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
        points_outer.append((x_outer, y_outer))
        points_inner.append((x_inner, y_inner))

    outer_polygon = Polygon(points_outer)
    inner_polygon = Polygon(points_inner)
    poly_sub = outer_polygon - inner_polygon

    if isinstance(poly_sub, Polygon):
        return geometry.Geometry(geom=poly_sub, material=material)

    msg = "Geometry generation failed."
    raise RuntimeError(msg)


def rectangular_hollow_section(
    d: float,
    b: float,
    t: float,
    r_out: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a rectangular hollow section (RHS).

    Constructs a rectangular hollow section (RHS) centered at ``(b/2, d/2)``, with depth
    ``d``, width ``b``, thickness ``t`` and outer radius ``r_out``, using ``n_r`` points
    to construct the inner and outer radii. If the outer radius is less than the
    thickness of the RHS, the inner radius is set to zero.

    Args:
        d: Depth of the RHS
        b: Width of the RHS
        t: Thickness of the RHS
        r_out: Outer radius of the RHS
        n_r: Number of points discretising the inner and outer radii
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        RuntimeError: If the geometry generation fails

    Returns:
        Rectangular hollow section geometry

    Example:
        The following example creates an RHS with a depth of 100 mm, a width of 50 mm,
        a thickness of 6 mm and an outer radius of 9 mm, using 8 points to discretise
        the inner and outer radii:

        .. plot::
            :include-source: True
            :caption: Rectangular hollow section geometry

            from sectionproperties.pre.library import rectangular_hollow_section

            rectangular_hollow_section(d=100, b=50, t=6, r_out=9, n_r=8).plot_geometry()
    """
    points_inner: list[tuple[float, float]] = []
    points_outer: list[tuple[float, float]] = []

    # calculate internal radius
    r_in = max(r_out - t, 0)

    # construct the outer radius points
    points_outer += sp_utils.draw_radius((r_out, r_out), r_out, np.pi, n_r)
    points_outer += sp_utils.draw_radius((b - r_out, r_out), r_out, 1.5 * np.pi, n_r)
    points_outer += sp_utils.draw_radius((b - r_out, d - r_out), r_out, 0, n_r)
    points_outer += sp_utils.draw_radius((r_out, d - r_out), r_out, 0.5 * np.pi, n_r)

    points_inner += sp_utils.draw_radius((t + r_in, t + r_in), r_in, np.pi, n_r)
    points_inner += sp_utils.draw_radius(
        (b - t - r_in, t + r_in),
        r_in,
        1.5 * np.pi,
        n_r,
    )
    points_inner += sp_utils.draw_radius((b - t - r_in, d - t - r_in), r_in, 0, n_r)
    points_inner += sp_utils.draw_radius(
        (t + r_in, d - t - r_in),
        r_in,
        0.5 * np.pi,
        n_r,
    )

    outer_polygon = Polygon(points_outer)
    inner_polygon = Polygon(points_inner)
    poly_sub = outer_polygon - inner_polygon

    if isinstance(poly_sub, Polygon):
        return geometry.Geometry(geom=poly_sub, material=material)

    msg = "Geometry generation failed."
    raise RuntimeError(msg)


def polygon_hollow_section(
    d: float,
    t: float,
    n_sides: int,
    r_in: float = 0.0,
    n_r: int = 1,
    rot: float = 0.0,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a regular hollow polygon section.

    Constructs a regular hollow polygon section centered at ``(0, 0)``, with a pitch
    circle diameter of bounding polygon ``d``, thickness ``t``, number of sides
    ``n_sides`` and an optional inner radius ``r_in``, using ``n_r`` points to construct
    the inner and outer radii (if radii is specified).

    Args:
        d: Pitch circle diameter of the outer bounding polygon (i.e. diameter of circle
            that passes through all vertices of the outer polygon)
        t: Thickness of the polygon section wall
        n_sides: Number of sides of the polygon
        r_in: Inner radius of the polygon corners. By default, if not specified, a
            polygon with no corner radii is generated. Defaults to ``0.0``.
        n_r: Number of points discretising the inner and outer radii, ignored if no
            inner radii is specified. Defaults to ``1``.
        rot: Initial counterclockwise rotation in degrees. By default bottom face is
            aligned with x axis. Defaults to ``0.0``.
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        ValueError: Number of sides in polygon must be greater than or equal to 3
        RuntimeError: If the geometry generation fails

    Returns:
        Regular hollow polygon section geometry

    Example:
        The following example creates an octagonal section (8 sides) with a diameter of
        200 mm, a thickness of 6 mm and an inner radius of 20 mm, using 12 points to
        discretise the inner and outer radii:

        .. plot::
            :include-source: True
            :caption: Rectangular hollow section geometry

            from sectionproperties.pre.library import polygon_hollow_section

            polygon_hollow_section(
                d=200, t=6, n_sides=8, r_in=20, n_r=12
            ).plot_geometry()
    """
    outer_points: list[tuple[float, float]] = []
    inner_points: list[tuple[float, float]] = []

    if n_sides < 3:
        msg = "n_sides required to be greater than 3 for polygon_section()"
        raise ValueError(msg)

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

    # calculate external radius, if r_in is zero, r_out also is zero
    if r_in == 0:
        r_out = 0.0
        n_r = 1
    else:
        r_out = r_in + t

    # equivalent side length of half the corner radii triangular segment
    c_out = r_out * (side_length_out / 2) / a_out
    c_in = r_in * (side_length_in / 2) / a_in

    # determine straight side length between corner radii (if present)
    side_length_straight_out = side_length_out - (2 * c_out)
    side_length_straight_in = side_length_in - (2 * c_in)

    # temp list for repeating geometry
    outer_base_points: list[tuple[float, float]] = []
    inner_base_points: list[tuple[float, float]] = []

    # start at bottom face, constructing one corner radii, then rotate by initial
    # rotation + alpha and repeat for n_side number of times to form full section
    # perimeter

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
        outer_base_points.append((x_outer, y_outer))
        inner_base_points.append((x_inner, y_inner))

    for i in range(n_sides):
        for point in outer_base_points:
            point_new = sp_utils.rotate(point, alpha * i + rot)
            outer_points.append(point_new)

        for point in inner_base_points:
            point_new = sp_utils.rotate(point, alpha * i + rot)
            inner_points.append(point_new)

    outer_polygon = Polygon(outer_points)
    inner_polygon = Polygon(inner_points)
    poly_sub = outer_polygon - inner_polygon

    if isinstance(poly_sub, Polygon):
        return geometry.Geometry(geom=poly_sub, material=material)

    msg = "Geometry generation failed."
    raise RuntimeError(msg)


def i_section(
    d: float,
    b: float,
    t_f: float,
    t_w: float,
    r: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs an I section.

    Constructs an I section centered at ``(b/2, d/2)``, with depth ``d``, width ``b``,
    flange thickness ``t_f``, web thickness ``t_w``, and root radius ``*, using ``n_r``
    points to construct the root radius.

    Args:
        d: Depth of the I section
        b: Width of the I section
        t_f: Flange thickness of the I section
        t_w: Web thickness of the I section
        r: Root radius of the I section
        n_r: Number of points discretising the root radius
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        I section geometry

    Example:
        The following example creates an I section with a depth of 203 mm, a width of
        133 mm, a flange thickness of 7.8 mm, a web thickness of 5.8 mm and a root
        radius of 8.9 mm, using 16 points to discretise the root radius:

        .. plot::
            :include-source: True
            :caption: I section geometry

            from sectionproperties.pre.library import i_section

            i_section(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=16).plot_geometry()
    """
    points: list[tuple[float, float]] = []

    # add first three points
    points.append((0, 0))
    points.append((b, 0))
    points.append((b, t_f))

    # construct the bottom right radius
    pt = b * 0.5 + t_w * 0.5 + r, t_f + r
    points += sp_utils.draw_radius(pt=pt, r=r, theta=1.5 * np.pi, n=n_r, ccw=False)

    # construct the top right radius
    pt = b * 0.5 + t_w * 0.5 + r, d - t_f - r
    points += sp_utils.draw_radius(pt=pt, r=r, theta=np.pi, n=n_r, ccw=False)

    # add the next four points
    points.append((b, d - t_f))
    points.append((b, d))
    points.append((0, d))
    points.append((0, d - t_f))

    # construct the top left radius
    pt = b * 0.5 - t_w * 0.5 - r, d - t_f - r
    points += sp_utils.draw_radius(pt=pt, r=r, theta=0.5 * np.pi, n=n_r, ccw=False)

    # construct the bottom left radius
    pt = b * 0.5 - t_w * 0.5 - r, t_f + r
    points += sp_utils.draw_radius(pt=pt, r=r, theta=0, n=n_r, ccw=False)

    # # add the last point
    points.append((0, t_f))

    i_section = Polygon(points)

    return geometry.Geometry(geom=i_section, material=material)


def mono_i_section(
    d: float,
    b_t: float,
    b_b: float,
    t_ft: float,
    t_fb: float,
    t_w: float,
    r: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a monosymmetric I section.

    Constructs a monosymmetric I section centered at ``(max(b_t, b_b)/2, d/2)``, with
    depth ``d``, top flange width ``b_t``, bottom flange width ``b_b``, top flange
    thickness ``t_ft``, top flange thickness ``t_fb``, web thickness ``t_w``, and root
    radius ``r``, using ``n_r`` points to construct the root radius.

    Args:
        d: Depth of the I section
        b_t: Top flange width
        b_b: Bottom flange width
        t_ft: Top flange thickness of the I section
        t_fb: Bottom flange thickness of the I section
        t_w: Web thickness of the I section
        r: Root radius of the I section
        n_r: Number of points discretising the root radius
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Monosymmetric I section geometry

    Example:
        The following example creates a monosymmetric I section with a depth of 200 mm,
        a top flange width of 50 mm, a top flange thickness of 12 mm, a bottom flange
        width of 130 mm, a bottom flange thickness of 8 mm, a web thickness of 6 mm and
        a root radius of 8 mm, using 16 points to discretise the root radius:

        .. plot::
            :include-source: True
            :caption: Monosymmetric I section geometry

            from sectionproperties.pre.library import mono_i_section

            mono_i_section(
                d=200, b_t=50, b_b=130, t_ft=12, t_fb=8, t_w=6, r=8, n_r=16
            ).plot_geometry()
    """
    points: list[tuple[float, float]] = []

    # calculate central axis
    x_central = max(b_t, b_b) * 0.5

    # add first three points
    points.append((x_central - b_b * 0.5, 0))
    points.append((x_central + b_b * 0.5, 0))
    points.append((x_central + b_b * 0.5, t_fb))

    # construct the bottom right radius
    pt = x_central + t_w * 0.5 + r, t_fb + r
    points += sp_utils.draw_radius(pt=pt, r=r, theta=1.5 * np.pi, n=n_r, ccw=False)

    # construct the top right radius
    pt = x_central + t_w * 0.5 + r, d - t_ft - r
    points += sp_utils.draw_radius(pt=pt, r=r, theta=np.pi, n=n_r, ccw=False)

    # add the next four points
    points.append((x_central + b_t * 0.5, d - t_ft))
    points.append((x_central + b_t * 0.5, d))
    points.append((x_central - b_t * 0.5, d))
    points.append((x_central - b_t * 0.5, d - t_ft))

    # construct the top left radius
    pt = x_central - t_w * 0.5 - r, d - t_ft - r
    points += sp_utils.draw_radius(pt=pt, r=r, theta=0.5 * np.pi, n=n_r, ccw=False)

    # construct the bottom left radius
    pt = x_central - t_w * 0.5 - r, t_fb + r
    points += sp_utils.draw_radius(pt=pt, r=r, theta=0, n=n_r, ccw=False)

    # add the last point
    points.append((x_central - b_b * 0.5, t_fb))

    polygon = Polygon(points)

    return geometry.Geometry(geom=polygon, material=material)


def tapered_flange_i_section(
    d: float,
    b: float,
    t_f: float,
    t_w: float,
    r_r: float,
    r_f: float,
    alpha: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a tapered flange I section.

    Constructs a tapered flange I section centered at ``(b/2, d/2)``, with depth ``d*``,
    mid-flange thickness ``t_f``, web thickness ``t_w``, root radius ``r_r``, flange
    radius ``r_f`` and flange angle ``alpha``, using ``n_r`` points to construct the
    radii.

    Args:
        d: Depth of the tapered flange I section
        b: Width of the tapered flange I section
        t_f: Mid-flange thickness of the tapered flange I section (measured at the point
            equidistant from the face of the web to the edge of the flange)
        t_w: Web thickness of the tapered flange I section
        r_r: Root radius of the tapered flange I section
        r_f: Flange radius of the tapered flange I section
        alpha: Flange angle of the tapered flange I section in degrees
        n_r: Number of points discretising the radii
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Tapered flange I section geometry

    Example:
        The following example creates a tapered flange I section with a depth of 588 mm,
        a width of 191 mm, a mid-flange thickness of 27.2 mm, a web thickness of 15.2
        mm, a root radius of 17.8 mm, a flange radius of 8.9 mm and a flange angle of
        8°, using 16 points to discretise the radii:

        .. plot::
            :include-source: True
            :caption: Tapered flange I section geometry

            from sectionproperties.pre.library import tapered_flange_i_section

            tapered_flange_i_section(
                d=588, b=191, t_f=27.2, t_w=15.2, r_r=17.8, r_f=8.9, alpha=8, n_r=16
            ).plot_geometry()
    """
    points: list[tuple[float, float]] = []

    # calculate alpha in radians
    alpha_rad = np.pi * alpha / 180

    # calculate the height of the flange toe and dimensions of the straight
    x1 = b * 0.25 - t_w * 0.25 - r_f * (1 - np.sin(alpha_rad))
    y1 = x1 * np.tan(alpha_rad)
    x2 = b * 0.25 - t_w * 0.25 - r_r * (1 - np.sin(alpha_rad))
    y2 = x2 * np.tan(alpha_rad)
    y_t = t_f - y1 - r_f * np.cos(alpha_rad)

    # add first two points
    points.append((0, 0))
    points.append((b, 0))

    # construct the bottom right flange toe radius
    if r_f == 0:
        points.append((b, y_t))
    else:
        for i in range(n_r):
            # determine polar angle
            theta = i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)

            # calculate the locations of the radius points
            x = b - r_f + r_f * np.cos(theta)
            y = y_t + r_f * np.sin(theta)

            # append the current points to the points list
            points.append((x, y))

    # construct the bottom right root radius
    if r_r == 0:
        points.append((b * 0.5 + t_w * 0.5, t_f + y2))
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
            points.append((x, y))

    # construct the top right root radius
    if r_r == 0:
        points.append((b * 0.5 + t_w * 0.5, d - t_f - y2))
    else:
        for i in range(n_r):
            # determine polar angle
            theta = np.pi - i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)

            # calculate the locations of the radius points
            x = b * 0.5 + t_w * 0.5 + r_r + r_r * np.cos(theta)
            y = d - t_f - y2 - r_r * np.cos(alpha_rad) + r_r * np.sin(theta)

            # append the current points to the points list
            points.append((x, y))

    # construct the top right flange toe radius
    if r_f == 0:
        points.append((b, d - y_t))
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
            points.append((x, y))

    # add the next two points
    points.append((b, d))
    points.append((0, d))

    # construct the top left flange toe radius
    if r_f == 0:
        points.append((0, d - y_t))
    else:
        for i in range(n_r):
            # determine polar angle
            theta = np.pi + (i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad))

            # calculate the locations of the radius points
            x = r_f + r_f * np.cos(theta)
            y = d - y_t + r_f * np.sin(theta)

            # append the current points to the points list
            points.append((x, y))

    # construct the top left root radius
    if r_r == 0:
        points.append((b * 0.5 - t_w * 0.5, d - t_f - y2))
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
            points.append((x, y))

    # construct the bottom left root radius
    if r_r == 0:
        points.append((b * 0.5 - t_w * 0.5, t_f + y2))
    else:
        for i in range(n_r):
            # determine polar angle
            theta = -i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)

            # calculate the locations of the radius points
            x = b * 0.5 - t_w * 0.5 - r_r + r_r * np.cos(theta)
            y = t_f + y2 + r_r * np.cos(alpha_rad) + r_r * np.sin(theta)

            # append the current points to the points list
            points.append((x, y))

    # construct the bottom left flange toe radius
    if r_f == 0:
        points.append((0, y_t))
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
            points.append((x, y))

    polygon = Polygon(points)

    return geometry.Geometry(geom=polygon, material=material)


def channel_section(
    d: float,
    b: float,
    t_f: float,
    t_w: float,
    r: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a parallel flange channel (PFC).

    Constructs a parallel flange channel (PFC) section with the bottom left corner at
    the origin ``(0, 0)``, with depth ``d``, width ``b``, flange thickness ``t_f``, web
    thickness ``t_w`` and root radius ``r``, using ``n_r`` points to construct the root
    radius.

    Args:
        d: Depth of the PFC section
        b: Width of the PFC section
        t_f: Flange thickness of the PFC section
        t_w: Web thickness of the PFC section
        r: Root radius of the PFC section
        n_r: Number of points discretising the root radius
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Parallel flange channel section geometry

    Example:
        The following example creates a PFC section with a depth of 250 mm, a width of
        90 mm, a flange thickness of 15 mm, a web thickness of 8 mm and a root radius of
        12 mm, using 8 points to discretise the root radius:

        .. plot::
            :include-source: True
            :caption: Parallel flange channel section geometry

            from sectionproperties.pre.library import channel_section

            channel_section(d=250, b=90, t_f=15, t_w=8, r=12, n_r=8).plot_geometry()
    """
    points: list[tuple[float, float]] = []

    # add first three points
    points.append((0, 0))
    points.append((b, 0))
    points.append((b, t_f))

    # construct the bottom right radius
    pt = t_w + r, t_f + r
    points += sp_utils.draw_radius(pt=pt, r=r, theta=1.5 * np.pi, n=n_r, ccw=False)

    # construct the top right radius
    pt = t_w + r, d - t_f - r
    points += sp_utils.draw_radius(pt=pt, r=r, theta=np.pi, n=n_r, ccw=False)

    # add last three points
    points.append((b, d - t_f))
    points.append((b, d))
    points.append((0, d))

    polygon = Polygon(points)

    return geometry.Geometry(geom=polygon, material=material)


def tapered_flange_channel(
    d: float,
    b: float,
    t_f: float,
    t_w: float,
    r_r: float,
    r_f: float,
    alpha: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a tapered flange channel section.

    Constructs a tapered flange channel section with the bottom left corner at the
    origin ``(0, 0)``, with depth ``d``, width *b*, mid-flange thickness ``t_f``, web
    thickness ``t_w``, root radius ``r_r``, flange radius ``r_f`` and flange angle
    ``alpha``, using ``n_r`` points to construct the radii.

    Args:
        d: Depth of the tapered flange channel section
        b: Width of the tapered flange channel section
        t_f: Mid-flange thickness of the tapered flange channel section (measured at the
            point equidistant from the face of the web to the edge of the flange)
        t_w: Web thickness of the tapered flange channel section
        r_r: Root radius of the tapered flange channel section
        r_f: Flange radius of the tapered flange channel section
        alpha: Flange angle of the tapered flange channel section in degrees
        n_r: Number of points discretising the radii
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Tapered flange channel section geometry

    Example:
        The following example creates a tapered flange channel section with a depth of
        10 in, a width of 3.5 in, a mid-flange thickness of 0.575 in, a web thickness of
        0.475 in, a root radius of 0.575 in, a flange radius of 0.4 in and a flange
        angle of 8°, using 16 points to discretise the radii:

        .. plot::
            :include-source: True
            :caption: Tapered flange channel section geometry

            from sectionproperties.pre.library import tapered_flange_channel

            tapered_flange_channel(
                d=10, b=3.5, t_f=0.575, t_w=0.475, r_r=0.575, r_f=0.4, alpha=8, n_r=16
            ).plot_geometry()
    """
    points: list[tuple[float, float]] = []

    # calculate alpha in radians
    alpha_rad = np.pi * alpha / 180

    # calculate the height of the flange toe and dimensions of the straight
    x1 = b * 0.5 - t_w * 0.5 - r_f * (1 - np.sin(alpha_rad))
    y1 = x1 * np.tan(alpha_rad)
    x2 = b * 0.5 - t_w * 0.5 - r_r * (1 - np.sin(alpha_rad))
    y2 = x2 * np.tan(alpha_rad)
    y_t = t_f - y1 - r_f * np.cos(alpha_rad)

    # add first two points
    points.append((0, 0))
    points.append((b, 0))

    # construct the bottom right flange toe radius
    if r_f == 0:
        points.append((b, y_t))
    else:
        for i in range(n_r):
            # determine polar angle
            theta = i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)

            # calculate the locations of the radius points
            x = b - r_f + r_f * np.cos(theta)
            y = y_t + r_f * np.sin(theta)

            # append the current points to the points list
            points.append((x, y))

    # construct the bottom right root radius
    if r_r == 0:
        points.append((t_w, t_f + y2))
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
            points.append((x, y))

    # construct the top right root radius
    if r_r == 0:
        points.append((t_w, d - t_f - y2))
    else:
        for i in range(n_r):
            # determine polar angle
            theta = np.pi - i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)

            # calculate the locations of the radius points
            x = t_w + r_r + r_r * np.cos(theta)
            y = d - t_f - y2 - r_r * np.cos(alpha_rad) + r_r * np.sin(theta)

            # append the current points to the points list
            points.append((x, y))

    # construct the top right flange toe radius
    if r_f == 0:
        points.append((b, d - y_t))
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
            points.append((x, y))

    # add the final two points
    points.append((b, d))
    points.append((0, d))

    polygon = Polygon(points)

    return geometry.Geometry(geom=polygon, material=material)


def tee_section(
    d: float,
    b: float,
    t_f: float,
    t_w: float,
    r: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a tee section.

    Constructs a tee section with the top left corner at ``(0, d)``, with depth ``d``,
    width ``b``, flange thickness ``t_f``, web thickness ``t_w`` and root radius ``r``,
    using ``n_r`` points to construct the root radius.

    Args:
        d: Depth of the tee section
        b: Width of the tee section
        t_f: Flange thickness of the tee section
        t_w: Web thickness of the tee section
        r: Root radius of the tee section
        n_r: Number of points discretising the root radius
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Tee section geometry

    Example:
        The following example creates a tee section with a depth of 200 mm, a width of
        100 mm, a flange thickness of 12 mm, a web thickness of 6 mm and a root radius
        of 8 mm, using 8 points to discretise the root radius:

        .. plot::
            :include-source: True
            :caption: Tee section geometry

            from sectionproperties.pre.library import tee_section

            tee_section(d=200, b=100, t_f=12, t_w=6, r=8, n_r=8).plot_geometry()
    """
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

    polygon = Polygon(points)

    return geometry.Geometry(geom=polygon, material=material)


def angle_section(
    d: float,
    b: float,
    t: float,
    r_r: float,
    r_t: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs an angle section.

    Constructs an angle section with the bottom left corner at the origin ``(0, 0)``,
    with depth ``d``, width ``b``, thickness ``t``, root radius ``r_r`` and toe radius
    ``r_t``, using ``n_r`` points to construct the radii.

    Args:
        d: Depth of the angle section
        b: Width of the angle section
        t: Thickness of the angle section
        r_r: Root radius of the angle section
        r_t: Toe radius of the angle section
        n_r: Number of points discretising the radii
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Angle section geometry

    Raises:
        ValueError: If the toe radius is larger than the thickness

    Example:
        The following example creates an angle section with a depth of 150 mm, a width
        of 100 mm, a thickness of 8 mm, a root radius of 12 mm and a toe radius of 5 mm,
        using 16 points to discretise the radii:

        .. plot::
            :include-source: True
            :caption: Angle section geometry

            from sectionproperties.pre.library import angle_section

            angle_section(d=150, b=100, t=8, r_r=12, r_t=5, n_r=16).plot_geometry()
    """
    if r_t > t:
        msg = "The radius of the toe (r_t) cannot be larger than the thickness (t)."
        raise ValueError(msg)

    points: list[tuple[float, float]] = []

    # add first two points, but don't add second if the toe radius equals the thickness
    points.append((0, 0))

    if not np.isclose(t, r_t):
        points.append((b, 0))

    # construct the bottom toe radius
    pt = b - r_t, t - r_t
    points += sp_utils.draw_radius(pt=pt, r=r_t, theta=0, n=n_r)

    # construct the root radius
    pt = t + r_r, t + r_r
    points += sp_utils.draw_radius(pt=pt, r=r_r, theta=1.5 * np.pi, n=n_r, ccw=False)

    # construct the top toe radius
    pt = t - r_t, d - r_t
    points += sp_utils.draw_radius(pt=pt, r=r_t, theta=0, n=n_r)

    # add the last point, only if the toe radius does not equal the thickness
    if not np.isclose(t, r_t):
        points.append((0, d))

    polygon = Polygon(points)

    return geometry.Geometry(geom=polygon, material=material)


def cee_section(
    d: float,
    b: float,
    l: float,
    t: float,
    r_out: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a cee section.

    Constructs a cee section (typical of cold-formed steel) with the bottom left corner
    at the origin ``(0, 0)``, with depth *d*, width ``b``, lip ``l``, thickness ``t``
    and outer radius ``r_out``, using ``n_r`` points to construct the radius. If the
    outer radius is less than the thickness of the Cee Section, the inner radius is set
    to zero.

    Args:
        d: Depth of the cee section
        b: Width of the cee section
        l: Lip of the cee section (which can be zero)
        t: Thickness of the cee section
        r_out: Outer radius of the cee section
        n_r: Number of points discretising the outer radius
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Cee section geometry

    Example:
        The following example creates a cee section with a depth of 125 mm, a width of
        30 mm, a lip of 30 mm, a thickness of 1.5 mm and an outer radius of 6 mm, using
        8 points to discretise the radius:

        .. plot::
            :include-source: True
            :caption: Cee section geometry

            from sectionproperties.pre.library import cee_section

            cee_section(d=125, b=50, l=30, t=1.5, r_out=6, n_r=8).plot_geometry()
    """
    points: list[tuple[float, float]] = []

    # calculate internal radius
    r_in = max(r_out - t, 0)

    # initialise lip variables
    r_out_l = 0.0
    r_in_l = 0.0

    # construct the outer bottom left radius
    points += sp_utils.draw_radius(pt=(r_out, r_out), r=r_out, theta=np.pi, n=n_r)

    # if the lip is longer than the outer radius (curve + straight section
    if l > r_out:
        # construct the outer bottom right radius
        points += sp_utils.draw_radius(
            pt=(b - r_out, r_out),
            r=r_out,
            theta=1.5 * np.pi,
            n=n_r,
        )

        # add next two points
        points.append((b, l))
        points.append((b - t, l))

        # construct the inner bottom right radius
        points += sp_utils.draw_radius(
            pt=(b - t - r_in, t + r_in),
            r=r_in,
            theta=0,
            n=n_r,
            ccw=False,
        )

    # if the lip is shorter than the outer radius (curve only)
    elif l > t and l <= r_out:
        # construct a smaller corner for bottom right if t < l < r_out
        r_out_l = l
        r_in_l = max(l - t, 0)
        points += sp_utils.draw_radius(
            pt=(b - r_out_l, r_out_l),
            r=r_out_l,
            theta=1.5 * np.pi,
            n=n_r,
        )
        points += sp_utils.draw_radius(
            pt=(b - t - r_in_l, t + r_in_l),
            r=r_in_l,
            theta=0,
            n=n_r,
            ccw=False,
        )

    # if the lip length is less than the section thickness (no lip)
    elif l <= t:
        # construct end as two points only
        points.append((b, 0))
        points.append((b, t))

    # construct the inner bottom left radius
    points += sp_utils.draw_radius(
        pt=(t + r_in, t + r_in),
        r=r_in,
        theta=1.5 * np.pi,
        n=n_r,
        ccw=False,
    )

    # construct the inner top left radius
    points += sp_utils.draw_radius(
        pt=(t + r_in, d - t - r_in),
        r=r_in,
        theta=np.pi,
        n=n_r,
        ccw=False,
    )

    # if the lip is longer than the outer radius (curve + straight section)
    if l > r_out:
        # construct the inner top right radius
        points += sp_utils.draw_radius(
            pt=(b - t - r_in, d - t - r_in),
            r=r_in,
            theta=0.5 * np.pi,
            n=n_r,
            ccw=False,
        )

        # add next two points
        points.append((b - t, d - l))
        points.append((b, d - l))

        # construct the outer top right radius
        points += sp_utils.draw_radius(
            pt=(b - r_out, d - r_out),
            r=r_out,
            theta=0,
            n=n_r,
        )

    # if the lip is shorter than the outer radius (curve only)
    elif l > t and l <= r_out:
        # construct a smaller corner for top right if t < l < r_out
        points += sp_utils.draw_radius(
            pt=(b - t - r_in_l, d - t - r_in_l),
            r=r_in_l,
            theta=0.5 * np.pi,
            n=n_r,
            ccw=False,
        )
        points += sp_utils.draw_radius(
            pt=(b - r_out_l, d - r_out_l),
            r=r_out_l,
            theta=0,
            n=n_r,
        )

    # if the lip length is less than the section thickness (no lip)
    elif l <= t:
        # construct end as two points only
        points.append((b, d - t))
        points.append((b, d))

    # construct the outer top left radius
    points += sp_utils.draw_radius(
        pt=(r_out, d - r_out),
        r=r_out,
        theta=0.5 * np.pi,
        n=n_r,
    )

    polygon = Polygon(points)

    return geometry.Geometry(geom=polygon, material=material)


def zed_section(
    d: float,
    b_l: float,
    b_r: float,
    l: float,
    t: float,
    r_out: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a zed section.

    Constructs a zed section with the bottom left corner at the origin ``(0, 0)``, with
    depth``*d``, left flange width ``b_l``, right flange width ``b_r``, lip ``l``,
    thickness ``t`` and outer radius ``r_out``, using ``n_r`` points to construct the
    radius. If the outer radius is less than the thickness of the Zed Section, the inner
    radius is set to zero.

    Args:
        d: Depth of the zed section
        b_l: Left flange width of the zed section
        b_r: Right flange width of the zed section
        l: Lip of the zed section (which can be zero)
        t: Thickness of the zed section
        r_out: Outer radius of the zed section
        n_r: Number of points discretising the outer radius
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Zed section geometry

    Example:
        The following example creates a zed section with a depth of 100 mm, a left
        flange width of 40 mm, a right flange width of 50 mm, a lip of 20 mm, a
        thickness of 1.2 mm and an outer radius of 5 mm, using 8 points to discretise
        the radius:

        .. plot::
            :include-source: True
            :caption: Zed section geometry

            from sectionproperties.pre.library import zed_section

            zed_section(
                d=100, b_l=40, b_r=50, l=20, t=1.2, r_out=5, n_r=8
            ).plot_geometry()
    """
    points: list[tuple[float, float]] = []

    # calculate internal radius
    r_in = max(r_out - t, 0)

    # initialise lip variables
    r_out_l = 0.0
    r_in_l = 0.0

    # construct the outer bottom left radius
    points += sp_utils.draw_radius(pt=(r_out, r_out), r=r_out, theta=np.pi, n=n_r)

    # if the lip is longer than the outer radius (curve + straight section
    if l > r_out:
        # construct the outer bottom right radius
        points += sp_utils.draw_radius(
            pt=(b_r - r_out, r_out),
            r=r_out,
            theta=1.5 * np.pi,
            n=n_r,
        )

        # add next two points
        points.append((b_r, l))
        points.append((b_r - t, l))

        # construct the inner bottom right radius
        points += sp_utils.draw_radius(
            pt=(b_r - t - r_in, t + r_in),
            r=r_in,
            theta=0,
            n=n_r,
            ccw=False,
        )

    # if the lip is shorter than the outer radius (curve only)
    elif l > t and l <= r_out:
        # construct a smaller corner for bottom right if t < l < r_out
        r_out_l = l
        r_in_l = max(l - t, 0)
        points += sp_utils.draw_radius(
            pt=(b_r - r_out_l, r_out_l),
            r=r_out_l,
            theta=1.5 * np.pi,
            n=n_r,
        )
        points += sp_utils.draw_radius(
            pt=(b_r - t - r_in_l, t + r_in_l),
            r=r_in_l,
            theta=0,
            n=n_r,
            ccw=False,
        )

    # if the lip length is less than the section thickness (no lip)
    elif l <= t:
        # construct end as two points only
        points.append((b_r, 0))
        points.append((b_r, t))

    # construct the inner bottom left radius
    points += sp_utils.draw_radius(
        pt=(t + r_in, t + r_in),
        r=r_in,
        theta=1.5 * np.pi,
        n=n_r,
        ccw=False,
    )

    # construct the outer top right radius
    points += sp_utils.draw_radius(pt=(t - r_out, d - r_out), r=r_out, theta=0, n=n_r)

    # if the lip is longer than the outer radius (curve + straight section
    if l > r_out:
        # construct the outer top left radius
        points += sp_utils.draw_radius(
            pt=(t - b_l + r_out, d - r_out),
            r=r_out,
            theta=0.5 * np.pi,
            n=n_r,
        )

        # add the next two points
        points.append((t - b_l, d - l))
        points.append((t - b_l + t, d - l))

        # construct the inner top left radius
        points += sp_utils.draw_radius(
            pt=(2 * t - b_l + r_in, d - t - r_in),
            r=r_in,
            theta=np.pi,
            n=n_r,
            ccw=False,
        )

    # if the lip is shorter than the outer radius (curve only)
    elif l > t and l <= r_out:
        # construct a smaller corner for top left if t < l < r_out
        points += sp_utils.draw_radius(
            pt=(t - b_l + r_out_l, d - r_out_l),
            r=r_out_l,
            theta=0.5 * np.pi,
            n=n_r,
        )
        points += sp_utils.draw_radius(
            pt=(2 * t - b_l + r_in_l, d - t - r_in_l),
            r=r_in_l,
            theta=np.pi,
            n=n_r,
            ccw=False,
        )

    # if the lip length is less than the section thickness (no lip)
    elif l <= t:
        # construct end as two points only
        points.append((t - b_l, d))
        points.append((t - b_l, d - t))

    # construct the inner top right radius
    points += sp_utils.draw_radius(
        pt=(-r_in, d - t - r_in),
        r=r_in,
        theta=0.5 * np.pi,
        n=n_r,
        ccw=False,
    )

    polygon = Polygon(points)

    return geometry.Geometry(geom=polygon, material=material)


def box_girder_section(
    d: float,
    b_t: float,
    b_b: float,
    t_ft: float,
    t_fb: float,
    t_w: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a box girder section.

    Constructs a box girder section centered at at ``(max(b_t, b_b)/2, d/2)``, with
    depth ``d``, top width ``b_t``, bottom width ``b_b``, top flange thickness ``t_ft``,
    bottom flange thickness ``t_fb`` and web thickness ``t_w``.

    Args:
        d: Depth of the box girder section
        b_t: Top width of the box girder section
        b_b: Bottom width of the box girder section
        t_ft: Top flange thickness of the box girder section
        t_fb: Bottom flange thickness of the box girder section
        t_w: Web thickness of the box girder section
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        RuntimeError: If the geometry generation fails

    Returns:
        Box girder section geometry

    Example:
        The following example creates a box girder section with a depth of 1200 mm, a
        top width of 1200 mm, a bottom width of 400 mm, a top flange thickness of 100
        mm, a bottom flange thickness of 80 mm and a web thickness of 50 mm:

        .. plot::
            :include-source: True
            :caption: Box girder section geometry

            from sectionproperties.pre.library import box_girder_section

            box_girder_section(
                d=1200, b_t=1200, b_b=400, t_ft=100, t_fb=80, t_w=50
            ).plot_geometry()
    """
    outer_points: list[tuple[float, float]] = []
    inner_points: list[tuple[float, float]] = []

    # calculate central axis
    x_c = max(b_t, b_b) * 0.5

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
    outer_points.append((x_c - 0.5 * b_b, 0))
    outer_points.append((x_c + 0.5 * b_b, 0))
    outer_points.append((x_c + 0.5 * b_t, d))
    outer_points.append((x_c - 0.5 * b_t, d))

    # add inner points
    inner_points.append((x_c - 0.5 * b_b - x_bot + web_x, t_fb))
    inner_points.append((x_c + 0.5 * b_b + x_bot - web_x, t_fb))
    inner_points.append((x_c + 0.5 * b_t + x_top - web_x, d - t_ft))
    inner_points.append((x_c - 0.5 * b_t - x_top + web_x, d - t_ft))

    outer_polygon = Polygon(outer_points)
    inner_polygon = Polygon(inner_points)
    poly_sub = outer_polygon - inner_polygon

    if isinstance(poly_sub, Polygon):
        return geometry.Geometry(geom=poly_sub, material=material)

    msg = "Geometry generation failed."
    raise RuntimeError(msg)


def bulb_section(
    d: float,
    b: float,
    t: float,
    r: float,
    n_r: int,
    d_b: float | None = None,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a bulb section.

    Constructs a bulb section with the bottom left corner at the point ``(-t / 2, 0)``,
    with depth ``d``, bulb depth ``d_b``, bulb width``*b``, web thickness ``t`` and
    radius ``r``, using ``n_r`` points to construct the radius.

    Args:
        d: Depth of the section
        b: Bulb width
        t: Web thickness
        r: Bulb radius
        d_b: Depth of the bulb (automatically calculated for standard sections, if
            provided the section may have sharp edges). Defaults to ``None``.
        n_r: Number of points discretising the radius
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Bulb section geometry

    Example:
        The following example creates a bulb section with a depth of 240 mm, a width of
        34 mm, a web thickness of 12 mm and a bulb radius of 16 mm, using 16 points to
        discretise the radius:

        .. plot::
            :include-source: True
            :caption: Bulb section geometry

            from sectionproperties.pre.library import bulb_section

            bulb_section(d=240, b=34, t=12, r=10, n_r=16).plot_geometry()
    """
    if d_b is None:
        d_b = r * np.cos(np.pi / 3) / np.cos(np.pi / 6) + r + b * np.tan(np.pi / 6)

    points: list[tuple[float, float]] = []

    # add first two points
    points.append((-t * 0.5, 0))
    points.append((t * 0.5, 0))

    # test of additional radius
    if d_b is not None:
        dc = r / np.sin(2 / 3 * np.pi / 2)
        ptb0 = (t * 0.5 + dc * np.cos(np.pi / 6), d - d_b - dc * np.cos(np.pi / 3))
        points += sp_utils.draw_radius(
            pt=ptb0,
            r=r,
            theta=np.pi,
            n=n_r,
            ccw=False,
            phi=np.pi / 3,
        )

    # end of test of additional radius
    ptb = b + t * 0.5 - r, d - r

    # build radius
    points += sp_utils.draw_radius(
        pt=ptb,
        r=r,
        theta=-np.pi * 1 / 3,
        n=n_r,
        ccw=True,
        phi=np.pi / 3,
    )
    points.pop()  # delete duplicate point
    points += sp_utils.draw_radius(pt=ptb, r=r, theta=0, n=n_r, ccw=True)

    # build the top part
    points.append((-t * 0.5, d))

    polygon = Polygon(points)

    return geometry.Geometry(geom=polygon, material=material)
