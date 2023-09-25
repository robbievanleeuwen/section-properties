"""Test plastic calculations."""

from __future__ import annotations

import pytest

import sectionproperties.pre.geometry as sp_geom
import sectionproperties.pre.library.primitive_sections as sections
import sectionproperties.pre.library.steel_sections as steel_sections
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
    assert sec_mat.get_mp()[0] == pytest.approx(mp)
    assert sec_mat.get_mp()[0] / fy == sec_nomat.get_s()[0]


def test_plastic_centroid():
    """Test created in response to #114.

    Since the section being tested is a compound geometry with two different
    materials, this tests that the plastic centroid takes into account the
    correct "center" of the original section which is affected by EA of each
    of the constituent geometries.
    """
    steel = Material(
        name="Steel",
        elastic_modulus=200e3,
        poissons_ratio=0.3,
        density=7.85e-6,
        yield_strength=500,
        color="grey",
    )
    timber = Material(
        name="Timber",
        elastic_modulus=5e3,
        poissons_ratio=0.35,
        density=6.5e-7,
        yield_strength=20,
        color="burlywood",
    )

    # create 310UB40.4
    ub = steel_sections.i_section(
        d=304, b=165, t_f=10.2, t_w=6.1, r=11.4, n_r=8, material=steel
    )

    # create timber panel on top of the UB
    panel = sections.rectangular_section(d=50, b=600, material=timber)
    panel = panel.align_center(align_to=ub).align_to(other=ub, on="top")

    # merge the two sections into one geometry object
    geometry = sp_geom.CompoundGeometry(geoms=[ub, panel])

    # create a mesh - use a mesh size of 5 for the UB, 20 for the panel
    geometry.create_mesh(mesh_sizes=[100, 100])

    # create a Section object
    section = Section(geometry=geometry)

    # perform a geometric, warping and plastic analysis
    section.calculate_geometric_properties()
    section.calculate_plastic_properties()

    # Checking sections that were defined in test_geometry
    small_sq = sections.rectangular_section(d=100, b=75)
    small_hole = sections.rectangular_section(d=40, b=30).align_center(
        align_to=small_sq
    )
    nested_geom = (small_sq - small_hole) + small_hole
    nested_geom.create_mesh(mesh_sizes=[50])
    nested_sec = Section(geometry=nested_geom)
    overlay_geom = small_sq + small_hole
    overlay_geom.create_mesh(mesh_sizes=[50])
    overlay_sec = Section(geometry=overlay_geom)

    nested_sec.calculate_geometric_properties()
    nested_sec.calculate_plastic_properties()
    overlay_sec.calculate_geometric_properties()

    with pytest.warns(UserWarning):
        overlay_sec.calculate_plastic_properties()

    # section
    x_pc, y_pc = section.get_pc()
    assert x_pc == pytest.approx(82.5)
    assert y_pc == pytest.approx(250.360654576)

    # nested_sec
    x_pc, y_pc = nested_sec.get_pc()
    assert x_pc == pytest.approx(37.5)
    assert y_pc == pytest.approx(50)
