"""Utilities used in the sections library."""

from __future__ import annotations

import numpy as np


def draw_radius(
    pt: tuple[float, float],
    r: float,
    theta: float,
    n: int,
    ccw: bool = True,
    phi: float = np.pi * 0.5,
) -> list[tuple[float, float]]:
    """Generates a list of points describing an arc.

    Generates an arc centered at point ``pt``, with radius ``r``, starting at angle
    ``theta``, with *n* points. If ``r=0``, adds ``pt`` only. ``phi`` describes the
    angle of rotation e.g. ``pi / 2`` is a quarter circle.

    Args:
        pt: Centre of radius (``x``, ``y``)
        r: Radius
        theta: Initial angle in radians
        n: Number of points
        ccw: If True, counter-clockwise rotation. Defaults to ``True``.
        phi: Angle describing radius extent in radians. Defaults to ``np.pi * 0.5``.

    Returns:
        List of points
    """
    points: list[tuple[float, float]] = []

    if r == 0:
        points.append(pt)
        return points

    mult = 1 if ccw else -1

    # calculate radius of points
    for i in range(n):
        # determine angle
        t = theta + mult * i * 1.0 / max(1, n - 1) * phi

        x = pt[0] + r * np.cos(t)
        y = pt[1] + r * np.sin(t)

        points.append((x, y))

    return points


def rotate(
    point: tuple[float, float],
    angle: float,
) -> tuple[float, float]:
    """Rotates a point counterclockwise by a given angle around origin ``(0, 0)``.

    Args:
        point: Point coordinates to be rotated
        angle: Angle to rotate point coordinates in radians

    Returns:
        Coordinates of rotated point
    """
    pt_x, pt_y = point

    c = np.cos(angle)
    s = np.sin(angle)

    new_x = c * pt_x - s * pt_y
    new_y = s * pt_x + c * pt_y

    return new_x, new_y
