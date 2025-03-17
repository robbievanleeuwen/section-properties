"""Timber sections library."""

from __future__ import annotations

import sectionproperties.pre.geometry as geometry
import sectionproperties.pre.library.primitive_sections as primitive_sections
import sectionproperties.pre.pre as pre


def clt_rectangular_section(
    d: list[float], lay_orient: list[pre.Material], b: float
) -> geometry.CompoundGeometry:
    """Constructs a timber rectangular section.

    Constructs a timber rectangular section of depth ``d`` and width ``b``.

    .. note::

    Args:
        d: Timber layer section thickness
        lay_orient: A list of materials for each layer from top to bottom
                    defined by the user.
        b: Timber section width

    Raises:
        ValueError: Geometry generation failed

    Returns:
        Timber rectangular section geometry

    Example:
        The following example creates a 120mm CLT cross-section.
    """
    layer_geom: list[geometry.Geometry] = []
    for idx in range(len(d)):
        di = float(d[idx])
        layer = lay_orient[idx]

        timb_mat = layer

        # create rectangular timber geometry
        layer = primitive_sections.rectangular_section(d=di, b=b, material=timb_mat)
        offset = -d[idx] * (idx + 1)
        layer = layer.shift_section(y_offset=offset)

        layer_geom.append(layer)

    # create compound geometry
    return geometry.CompoundGeometry(geoms=layer_geom)
