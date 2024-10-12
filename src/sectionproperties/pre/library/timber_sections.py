"""Timber sections library."""

from __future__ import annotations

import sectionproperties.pre.geometry as geometry
import sectionproperties.pre.library.primitive_sections as primitive_sections
import sectionproperties.pre.pre as pre


def clt_rectangular_section(
    d: list[float],
    lay_orient: list[int],
    b: float,
    timb_mat0: pre.Material = pre.DEFAULT_MATERIAL,
    timb_mat90: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.CompoundGeometry:
    """Constructs a timber rectangular section.

    Constructs a timber rectangular section of depth ``d`` and width ``b``.

    .. note::

    Args:
        d: Timber layer section thickness
        lay_orient: List of layer orientation
        b: Timber section width
        timb_mat0: Material object to assign to the timber area
                   parallel-to-grain
        timb_mat90: Material object to assign to the timber area,
                    perpendicular-to-grain

    Raises:
        ValueError: Geometry generation failed

    Returns:
        Timber rectangular section geometry

    Example:
        The following example creates a 120mm CLT cross-section.
    """

    layer_geom = list()
    for idx in range(len(d)):
        di = float(d[idx])
        layer = lay_orient[idx]

        if layer is int(0):
            timb_mat = timb_mat0
        else:
            timb_mat = timb_mat90

        # create rectangular timber geometry
        layer = primitive_sections.rectangular_section(d=di, b=b,
                                                       material=timb_mat)
        offset = -d[idx] * (idx + 1)
        layer = layer.shift_section(y_offset=offset)

        layer_geom.append(layer)

    # create compound geometry
    geom = geometry.CompoundGeometry(geoms=layer_geom)

    if isinstance(geom, geometry.CompoundGeometry):
        return geom
    else:
        raise ValueError("Timber section generation failed.")
