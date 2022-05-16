import numpy as np
from typing import List


def draw_radius(
    pt: List[float],
    r: float,
    theta: float,
    n: int,
    ccw: bool = True,
    phi: float = np.pi * 0.5,
) -> List[List[float]]:
    """Adds a quarter radius of points to the points list - centered at point *pt*,
    with radius *r*, starting at angle *theta*, with *n* points. If r = 0, adds pt
    only. Phi describes the angle of rotation e.g. pi / 2 is a quarter circle.

    :param pt: Centre of radius *(x,y)*
    :type pt: List[float]
    :param float r: Radius
    :param float theta: Initial angle
    :param int n: Number of points
    :param bool ccw: Counter-clockwise rotation?
    :param float phi: Angle describing radius extent (default 90 degrees)

    :return: List of points
    :rtype: List[List[float]]
    """

    points = []

    if r == 0:
        points.append(pt)
        return points

    if ccw:
        mult = 1
    else:
        mult = -1

    # calculate radius of points
    for i in range(n):
        # determine angle
        t = theta + mult * i * 1.0 / max(1, n - 1) * phi

        x = pt[0] + r * np.cos(t)
        y = pt[1] + r * np.sin(t)

        points.append([x, y])

    return points


def rotate(point, angle: float):
    """
    Rotate a point counterclockwise by a given angle around origin [0, 0]

    :param list point: Point coordinates to be rotated
    :param float angle: Angle to rotate point coordinates
    :return: Coordinates of rotated point
    :rtype: list[float, float]
    """

    pt_x, pt_y = point

    c = np.cos(angle)
    s = np.sin(angle)

    new_x = c * pt_x - s * pt_y
    new_y = s * pt_x + c * pt_y

    return [new_x, new_y]
