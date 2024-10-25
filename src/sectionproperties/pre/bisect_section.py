"""A number of functions required to bisect shapely geometries."""

from __future__ import annotations

from typing import TYPE_CHECKING

from shapely import GeometryCollection, LineString, Polygon

if TYPE_CHECKING:
    import numpy as np
    import numpy.typing as npt


def create_line_segment(
    point_on_line: tuple[float, float] | npt.NDArray[np.float64],
    vector: npt.NDArray[np.float64],
    bounds: tuple[float, float, float, float],
) -> LineString:
    """Create a line segment given a point, vector and bounds.

    Return a LineString of a line that contains ``point_on_line`` in the direction
    of ``vector`` bounded by ``bounds``.

    Args:
        point_on_line: Point on line
        vector: Line direction
        bounds: Min and max ordinates of geometry

    Returns:
        Line string defined as per docstring
    """
    p_x, p_y = point_on_line
    b_2 = max(bounds)
    b_1 = min(bounds)

    if abs(vector[0]) > 1e-12:  # Not a vertical line
        scale_factor_2 = (b_2 - p_x) / vector[0]
        y_2 = scale_factor_2 * vector[1] + p_y
        scale_factor_1 = (b_1 - p_x) / vector[0]
        y_1 = scale_factor_1 * vector[1] + p_y

        return LineString([(b_1, y_1), (b_2, y_2)])
    else:  # Vertical line
        scale_factor_2 = (b_2 - p_y) / vector[1]
        x_2 = scale_factor_2 * vector[0] + p_x
        scale_factor_1 = (b_1 - p_y) / vector[1]
        x_1 = scale_factor_1 * vector[0] + p_x

        return LineString([(x_1, b_1), (x_2, b_2)])


def group_top_and_bottom_polys(
    polys: GeometryCollection,
    line: LineString,
) -> tuple[list[Polygon], list[Polygon]]:
    """Sort polygons into groups above and below a line.

    Returns a tuple of two lists representing the list of polygons in ``polys`` on
    the "top" side of ``line`` and the list of polygons on the "bottom" side of the
    ``line`` after the original geometry has been split by ``line``.

    The 0th tuple element is the "top" polygons and the 1st element is the "bottom"
    polygons.

    In the event that ``line`` is a perfectly vertical line, the "top" polys are the
    polygons on the "right" of the ``line`` and the "bottom" polys are the polygons on
    the "left" of the ``line``.

    Args:
        polys: Collection of polygons to sort
        line: Line about which polygons have been split

    Raises:
        RuntimeError: If the split operation does not generate polygons

    Returns:
        List of polygons above and below the line
    """
    top_acc: list[Polygon] = []
    bot_acc: list[Polygon] = []

    for poly in polys.geoms:
        if not isinstance(poly, Polygon):
            msg = "Geometry split error."
            raise RuntimeError(msg)

        m, b = line_mx_plus_b(line)
        px, py = poly.representative_point().coords[0]

        if b is not None:  # Not a vertical line (special case)
            y_test = m * px + b
            if py < y_test:
                bot_acc.append(poly)
            elif py > y_test:
                top_acc.append(poly)

        else:  # The special case of vertical line
            lx, _ = line.coords[0]
            if px < lx:
                bot_acc.append(poly)
            elif px > lx:
                top_acc.append(poly)

    return top_acc, bot_acc


def line_mx_plus_b(
    line: LineString,
) -> tuple[float, float | None]:
    """Get the slope and y-intercept of a line.

    Returns a tuple representing the values of "m" and "b" from the definition of
    ``line`` as "y = mx + b".

    Args:
        line: Line to define

    Returns:
        Gradient and y-intercept of line, ``(1, None)`` if vertical
    """
    y2, y1 = line.coords[1][1], line.coords[0][1]
    x2, x1 = line.coords[1][0], line.coords[0][0]

    if x2 - x1 == 0:
        return (1, None)

    m_slope = (y2 - y1) / (x2 - x1)
    point_on_line = line.coords[0]
    p_x, p_y = point_on_line

    # solve line eqn for b given a known point on the line
    b_intercept = p_y - m_slope * p_x

    return (m_slope, b_intercept)
