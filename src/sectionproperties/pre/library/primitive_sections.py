"""Primitive sections library."""

from __future__ import annotations

import numpy as np
from shapely import Polygon

import sectionproperties.pre.geometry as geometry
import sectionproperties.pre.library.utils as sp_utils
import sectionproperties.pre.pre as pre


def rectangular_section(
    d: float,
    b: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a rectangular section.

    Constructs a rectangular section with the bottom left corner at the origin
    ``(0, 0)``, with depth ``d`` and width ``b``.

    Args:
        d: Depth (``y``) of the rectangle
        b: Width (``x``) of the rectangle
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Rectangular section geometry

    Example:
        The following example creates a rectangular cross-section with a depth of 100 mm
        and width of 50 mm:

        .. plot::
            :include-source: True
            :caption: Rectangular section geometry

            from sectionproperties.pre.library import rectangular_section

            rectangular_section(d=100, b=50).plot_geometry()
    """
    points = [(0, 0), (b, 0), (b, d), (0, d)]
    rectangle = Polygon(points)

    return geometry.Geometry(geom=rectangle, material=material)


def circular_section(
    d: float,
    n: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a circular section.

    Constructs a solid circle centered at the origin ``(0, 0)`` with diameter ``d``,
    using ``n`` points to construct the circle.

    Args:
        d: Diameter of the circle
        n: Number of points discretising the circle
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Circular section geometry

    Example:
        The following example creates a circular geometry with a diameter of 50 mm with
        64 points:

        .. plot::
            :include-source: True
            :caption: Circular section geometry

            from sectionproperties.pre.library import circular_section

            circular_section(d=50, n=64).plot_geometry()
    """
    x_off, y_off = (0, 0)
    points: list[tuple[float, float]] = []

    # loop through each point on the circle
    for i in range(n):
        # determine polar angle
        theta = i * 2 * np.pi * 1.0 / n

        # calculate location of the point
        x = 0.5 * d * np.cos(theta) + x_off
        y = 0.5 * d * np.sin(theta) + y_off

        # append the current point to the points list
        points.append((x, y))

    circle = Polygon(points)

    return geometry.Geometry(geom=circle, material=material)


def circular_section_by_area(
    area: float,
    n: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    r"""Constructs a circular section defined by its area.

    Constructs a solid circle centered at the origin ``(0, 0`` defined by its ``area``,
    using ``n`` points to construct the circle.

    Args:
        area: Area of the circle
        n: Number of points discretising the circle
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Circular section geometry

    Example:
        The following example creates a circular geometry with an area of 310
        mm\ :sup:`2` with 32 points:

        .. plot::
            :include-source: True
            :caption: Circular section geometry

            from sectionproperties.pre.library import circular_section_by_area

            circular_section_by_area(area=310, n=32).plot_geometry()
    """
    s = 2 * np.sqrt(area / n) * np.sqrt(np.tan(np.pi / n))
    a = s / (2 * np.tan(np.pi / n))
    d = np.sqrt(a * a + (0.5 * s) * (0.5 * s)) * 2

    return circular_section(d=d, n=n, material=material)


def elliptical_section(
    d_x: float,
    d_y: float,
    n: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs an elliptical section.

    Constructs a solid ellipse centered at the origin ``(0, 0)`` with horizontal
    diameter ``d_x`` and vertical diameter ``d_y``, using ``n`` points to construct the
    ellipse.

    Args:
        d_x: Diameter of the ellipse in the x-dimension
        d_y: Diameter of the ellipse in the y-dimension
        n: Number of points discretising the ellipse
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Elliptical section geometry

    Example:
        The following example creates an elliptical cross-section with a horizontal
        diameter of 50 mm and a vertical diameter of 25 mm, with 40 points:

        .. plot::
            :include-source: True
            :caption: Elliptical section geometry

            from sectionproperties.pre.library import elliptical_section

            elliptical_section(d_x=50, d_y=25, n=40).plot_geometry()
    """
    points: list[tuple[float, float]] = []

    # loop through each point on the ellipse
    for i in range(n):
        # determine polar angle
        theta = i * 2 * np.pi * 1.0 / n

        # calculate location of the point
        x = 0.5 * d_x * np.cos(theta)
        y = 0.5 * d_y * np.sin(theta)

        # append the current point to the points list
        points.append((x, y))

    ellipse = Polygon(points)

    return geometry.Geometry(geom=ellipse, material=material)


def triangular_section(
    b: float,
    h: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a triangular section.

    Constructs a right angled triangle with points ``(0, 0)``, ``(b, 0)``, ``(0, h)``.

    Args:
        b: Base length of triangle
        h: Height of triangle
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Triangular section geometry

    Example:
        The following example creates a triangular cross-section with a base width of 10
        mm and height of 10 mm:

        .. plot::
            :include-source: True
            :caption: Triangular section geometry

            from sectionproperties.pre.library import triangular_section

            triangular_section(b=10, h=10).plot_geometry()
    """
    points = [(0, 0), (b, 0), (0, h)]
    triangle = Polygon(points)
    return geometry.Geometry(triangle, material)


def triangular_radius_section(
    b: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a triangular section with a radius.

    Constructs a right angled isosceles triangle with points``(0, 0)``, ``(b, 0)``,
    ``(0, h)`` and a concave radius on the hypotenuse.

    Args:
        b: Base length of triangle
        n_r: Number of points discretising the radius
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Triangular section geometry

    Example:
        The following example creates an isosceles triangular radius cross-section with
        a base width of 6 mm, using 16 points to construct the radius:

        .. plot::
            :include-source: True
            :caption: Triangular section geometry

            from sectionproperties.pre.library import triangular_radius_section

            triangular_radius_section(b=6, n_r=16).plot_geometry()
    """
    points = [(0.0, 0.0)]
    points += sp_utils.draw_radius(
        pt=(b, b), r=b, theta=3 * np.pi / 2, n=n_r, ccw=False
    )
    triangle = Polygon(points)

    return geometry.Geometry(geom=triangle, material=material)


def cruciform_section(
    d: float,
    b: float,
    t: float,
    r: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a cruciform section.

    Constructs a cruciform section centered at the origin ``(0, 0)``, with depth ``d``,
    width ``b``, thickness ``t`` and root radius ``r``, using ``n_r`` points to
    construct the root radius.

    Args:
        d: Depth of the cruciform section
        b: Width of the cruciform section
        t: Thickness of the cruciform section
        r: Root radius of the cruciform section
        n_r: Number of points discretising the root radius
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Cruciform section geometry

    Example:
        The following example creates a cruciform section with a depth of 250 mm, a
        width of 175 mm, a thickness of 12 mm and a root radius of 16 mm, using 16
        points to discretise the radius:

        .. plot::
            :include-source: True
            :caption: Cruciform section geometry

            from sectionproperties.pre.library import cruciform_section

            cruciform_section(d=250, b=175, t=12, r=16, n_r=16).plot_geometry()
    """
    points: list[tuple[float, float]] = []

    # add first two points
    points.append((-t * 0.5, -d * 0.5))
    points.append((t * 0.5, -d * 0.5))

    # construct the bottom right radius
    pt = 0.5 * t + r, -0.5 * t - r
    points += sp_utils.draw_radius(pt=pt, r=r, theta=np.pi, n=n_r, ccw=False)

    # add the next two points
    points.append((0.5 * b, -t * 0.5))
    points.append((0.5 * b, t * 0.5))

    # construct the top right radius
    pt = 0.5 * t + r, 0.5 * t + r
    points += sp_utils.draw_radius(pt=pt, r=r, theta=1.5 * np.pi, n=n_r, ccw=False)

    # add the next two points
    points.append((t * 0.5, 0.5 * d))
    points.append((-t * 0.5, 0.5 * d))

    # construct the top left radius
    pt = -0.5 * t - r, 0.5 * t + r
    points += sp_utils.draw_radius(pt=pt, r=r, theta=0, n=n_r, ccw=False)

    # add the next two points
    points.append((-0.5 * b, t * 0.5))
    points.append((-0.5 * b, -t * 0.5))

    # construct the bottom left radius
    pt = -0.5 * t - r, -0.5 * t - r
    points += sp_utils.draw_radius(pt=pt, r=r, theta=0.5 * np.pi, n=n_r, ccw=False)

    polygon = Polygon(points)

    return geometry.Geometry(geom=polygon, material=material)
