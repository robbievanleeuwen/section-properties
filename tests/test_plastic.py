"""Test plastic calculations."""

from __future__ import annotations

import pytest

import sectionproperties.pre.library.primitive_sections as sections
from sectionproperties.analysis.section import Section
from sectionproperties.pre.pre import Material


def test_rectangle():
    """Test plastic properties of a rectangle."""
    fy = 500
    e = 200e3
    b = 50
    d = 100

    steel = Material(
        name="Steel",
        elastic_modulus=e,
        poissons_ratio=0.3,
        yield_strength=fy,
        density=8.05e-6,
        color="grey",
    )

    sx = b * d * d / 4
    mp = sx * fy

    geom_mat = sections.rectangular_section(d=d, b=b, material=steel)
    geom_nomat = sections.rectangular_section(d=d, b=b)

    geom_mat.create_mesh(mesh_sizes=[2.5])
    geom_nomat.create_mesh(mesh_sizes=[2.5])

    sec_mat = Section(geometry=geom_mat)
    sec_nomat = Section(geometry=geom_nomat)

    sec_mat.calculate_geometric_properties()
    sec_mat.calculate_plastic_properties()

    sec_nomat.calculate_geometric_properties()
    sec_nomat.calculate_plastic_properties()

    assert sec_nomat.get_s()[0] == pytest.approx(sx)
    assert sec_mat.get_s()[0] == pytest.approx(mp)
    assert sec_mat.get_s()[0] / fy == sec_nomat.get_s()[0]
