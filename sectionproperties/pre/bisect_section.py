from typing import List, Optional, Tuple, Union

import numpy as np
import shapely.geometry as geom


def create_line_segment(
    point_on_line: Union[Tuple[float, float], np.ndarray],
    vector: np.ndarray,
    bounds: Tuple[float, float, float, float],
) -> geom.LineString:
    """Return a LineString of a line that contains ``point_on_line`` in the direction
    of ``vector`` bounded by ``bounds``.

    :param point_on_line: Point on line
    :param vector: Line direction
    :param bounds: Min and max ordinates of geometry

    :return: Line string defined as per docstring
    """

    p_x, p_y = point_on_line
    b_2 = max(bounds)
    b_1 = min(bounds)

    if abs(vector[0]) > 1e-12:  # Not a vertical line
        scale_factor_2 = (b_2 - p_x) / vector[0]
        y_2 = scale_factor_2 * vector[1] + p_y

        scale_factor_1 = (b_1 - p_x) / vector[0]
        y_1 = scale_factor_1 * vector[1] + p_y

        return geom.LineString([(b_1, y_1), (b_2, y_2)])
        
    else:  # Vertical line
        scale_factor_2 = (b_2 - p_y) / vector[1]
        x_2 = scale_factor_2 * vector[0] + p_x

        scale_factor_1 = (b_1 - p_y) / vector[1]
        x_1 = scale_factor_1 * vector[0] + p_x

        return geom.LineString([(x_1, b_1), (x_2, b_2)])


def group_top_and_bottom_polys(
    polys: geom.GeometryCollection,
    line: geom.LineString,
) -> Tuple[List[geom.Polygon], List[geom.Polygon]]:
    """Returns a tuple of two lists representing the list of polygons in ``polys`` on
    the "top" side of ``line`` and the list of polygons on the "bottom" side of the
    ``line`` after the original geometry has been split by ``line``.

    The 0th tuple element is the "top" polygons and the 1st element is the "bottom"
    polygons.

    In the event that ``line`` is a perfectly vertical line, the "top" polys are the
    polygons on the "right" of the ``line`` and the "bottom" polys are the polygons on
    the "left" of the ``line``.

    :param polys: Collection of polygons to sort
    :param line: Line about which polygons have been split

    :return: List of polygons above and below the line
    """

    top_acc = []
    bot_acc = []

    for poly in polys.geoms:
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
    line: geom.LineString,
) -> Tuple[float, Union[float, None]]:
    """Returns a tuple representing the values of "m" and "b" from the definition of
    ``line`` as "y = mx + b".

    :param line: Line to define

    :return: Gradient and y-intercept of line, ``(1, None)`` if vertical
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


def perp_mx_plus_b(
    m_slope: float,
    point_on_line: Tuple[float, float],
) -> Tuple[float, float]:
    """Returns a tuple representing the values of "m" and "b" from for a line that is
    perpendicular to ``m_slope`` and contains the ``point_on_line``, which represents
    an (x, y) coordinate.

    :param m_slope: Slope of the line
    :param point_on_line: A point on the line

    :return: Gradient and y-intercept of perpendicular line
    """

    m_perp = -1 / m_slope
    p_x, p_y = point_on_line
    b_intercept = p_y - m_perp * p_x

    return (m_perp, b_intercept)


def line_intersection(
    m_1: float,
    b_1: float,
    m_2: float,
    b_2: float,
) -> Optional[float]:
    """Returns a float representing the x-ordinate of the intersection point of the
    lines defined by y = m_1*x + b_1 and y = m_2*x + b_2.

    Returns None if the lines are parallel.

    :param m_1: Gradient of line 1
    :param b_1: y-intercept of line 1
    :param m_2: Gradient of line 2
    :param b_2: y-intercept of line 2

    :return: x-ordinate of intersection of lines
    """

    try:
        x = (b_2 - b_1) / (m_1 - m_2)
    except ZeroDivisionError:
        x = None
    
    return x


def sum_poly_areas(
    lop: List[geom.Polygon],
) -> float:
    """Returns a float representing the total area of all polygons in ``lop``, the list
    of polygons.

    :param lop: List of polygons to sum areas

    :return: Sum of areas
    """

    sum_acc = 0

    for poly in lop:
        sum_acc += poly.area
    
    return sum_acc
