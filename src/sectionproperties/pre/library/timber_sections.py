"""Timber sections library."""

from __future__ import annotations

import numpy as np

import sectionproperties.pre.geometry as geometry
import sectionproperties.pre.library.primitive_sections as primitive_sections
import sectionproperties.pre.pre as pre


def timber_rectangular_section(
        d: float,
        b: float,
        timb_mat: pre.Material = pre.DEFAULT_MATERIAL
) -> geometry.Geometry:
        """Constructs a timber rectangular section.

        Constructs a timber rectangular section of depth ``d`` and width ``b``.

        .. note::

        Args:
            d: Timber section depth
            b: Timber section width
            timb_mat: Material object to assign to the timber area

        Raises:
            ValueError: Geometry generation failed

        Returns:
            Timber rectangular section geometry

        Example:
            The following example creates a 600mm deep x 300mm wide timber gluelaminated beam.

        .. plot::
            :include-source: True
            :caption: Timber rectangular section geometry

            from sectionproperties.pre import Material
            from sectionproperties.pre.library import timber_rectangular_section
            from sectionproperties.analysis import Section

            geom = timber_rectangular_section(500, 300)
            geom.create_mesh(mesh_sizes=[200])

            sec = Section(geometry=geom)
            sec.calculate_geometric_properties()

            ixx_c, iyy_c, ixy_c = sec.get_ic()
            print(f"Ixx = {ixx_c:.3e} mm4")
            print(f"Irec = {300 * 500**3 / 12:.3e} mm4")
            print(f"Iyy = {iyy_c:.3e} mm4")
            print(f"Irec = {500 * 300**3 / 12:.3e} mm4")
        """

        # create rectangular timber geometry
        geom = primitive_sections.rectangular_section(b=b, d=d, material=timb_mat)

        if isinstance(geom, geometry.Geometry):
            return geom
        else:
            raise ValueError("Timber section generation failed.")
