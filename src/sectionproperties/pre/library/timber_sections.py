"""Timber sections library."""

from __future__ import annotations

import sectionproperties.pre.geometry as geometry
import sectionproperties.pre.library.primitive_sections as primitive_sections
import sectionproperties.pre.pre as pre


def clt_rectangular_section(
    d: list[float],
    layer_mat: list[pre.Material],
    b: float,
) -> geometry.CompoundGeometry:
    """Constructs a CLT rectangular section.

    Constructs a CLT rectangular section, with layer depths ``d``, layer materials
    ``layer_mat``, and width ``b``.

    Args:
        d: CLT layer section thickness
        layer_mat: A list of timber materials for each layer (from top to bottom)
        b: CLT section width

    Returns:
        CLT rectangular section geometry

    Example:
        The following example creates a 120mm CLT cross-section:

        .. plot::
            :include-source: True
            :caption: 120mm CLT section geometry

            from sectionproperties.pre import Material
            from sectionproperties.pre.library import clt_rectangular_section
            from sectionproperties.analysis import Section

            timber0 = Material(
                name="Timber0",
                elastic_modulus=9.5e3,
                poissons_ratio=0.35,
                density=4.4e-7,
                yield_strength=5.5,
                color="burlywood",
            )
            timber90 = Material(
                name="Timber90",
                elastic_modulus=317,
                poissons_ratio=0.35,
                density=4.4e-7,
                yield_strength=5.5,
                color="orange",
            )

            geom = clt_rectangular_section(
                d=[40, 40, 40],
                layer_mat=[timber0, timber90, timber0],
                b=1000
            )

            geom.create_mesh(mesh_sizes=[0])  # a size of zero creates a coarse mesh
            Section(geometry=geom).plot_mesh()
    """
    layer_geom: list[geometry.Geometry] = []
    for idx in range(len(d)):
        di = float(d[idx])
        layer = layer_mat[idx]

        timb_mat = layer

        # create rectangular timber geometry
        layer = primitive_sections.rectangular_section(d=di, b=b, material=timb_mat)
        offset = -d[idx] * (idx + 1)
        layer = layer.shift_section(y_offset=offset)

        layer_geom.append(layer)

    # create compound geometry
    return geometry.CompoundGeometry(geoms=layer_geom)
