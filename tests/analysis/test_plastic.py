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
        d=304,
        b=165,
        t_f=10.2,
        t_w=6.1,
        r=11.4,
        n_r=8,
        material=steel,
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

    # perform a geometric and plastic analysis
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


def test_plastic_composite_simple():
    """Tests the plastic properties of a simple composite section."""
    mat1 = Material(
        name="mat1",
        elastic_modulus=10,
        poissons_ratio=0.3,
        density=1,
        yield_strength=50,
        color="grey",
    )

    mat2 = Material(
        name="mat2",
        elastic_modulus=20,
        poissons_ratio=0.3,
        density=1,
        yield_strength=80,
        color="black",
    )

    rect1 = sections.rectangular_section(d=100, b=50, material=mat1)
    rect2 = sections.rectangular_section(d=100, b=50, material=mat2).shift_section(
        x_offset=50
    )
    geom = rect1 + rect2
    geom.create_mesh(mesh_sizes=100)

    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()
    sec.calculate_plastic_properties()

    x_pc, y_pc = sec.get_pc()
    assert y_pc == pytest.approx(50)

    # calculate x_pc
    f_mat1 = 100 * 50 * 50  # force in mat1
    x_mat2 = f_mat1 / 80 / 100  # width of mat2 required to equal force in mat1
    x_left = (50 - x_mat2) / 2  # width of ma2 to left of pc
    assert x_pc == pytest.approx(x_left + 50)

    mp_xx, mp_yy = sec.get_mp()
    mp_xx_calc = 50 * 100**2 / 4 * (80 + 50)
    assert mp_xx == pytest.approx(mp_xx_calc)

    # calculate mp_yy
    mp_yy_calc = 0
    mp_yy_calc += f_mat1 * abs(25 - x_pc)  # mat 1 contribution
    # left mat 2 contribution
    mp_yy_calc += x_left * 100 * 80 * abs(50 + x_left / 2 - x_pc)
    # right mat 2 contribution
    mp_yy_calc += (50 - x_left) * 100 * 80 * abs(100 - (50 - x_left) / 2 - x_pc)
    assert mp_yy == pytest.approx(mp_yy_calc)

    # check principal calc with rotation 45 deg
    geom.rotate_section(angle=45)
    geom.create_mesh(mesh_sizes=100)
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()
    sec.calculate_plastic_properties()

    mp_11, mp_22 = sec.get_mp_p()
    assert mp_11 == pytest.approx(mp_xx_calc)
    assert mp_22 == pytest.approx(mp_yy_calc)


def test_plastic_composite_example():
    """Tests the composite example on the sectionproperties docs, refer issue #460.

    The steel section is simplified to have no root radii to allow hand calculation
    comparison.
    """
    # create the steel material
    steel = Material(
        name="Steel",
        elastic_modulus=200e3,
        poissons_ratio=0.3,
        density=7.85e-6,
        yield_strength=500,
        color="grey",
    )

    # create the timber material
    timber = Material(
        name="Timber",
        elastic_modulus=8e3,
        poissons_ratio=0.35,
        yield_strength=20,
        density=0.78e-6,
        color="burlywood",
    )

    # universal steel beam
    ub = steel_sections.i_section(
        d=304,
        b=165,
        t_f=10.2,
        t_w=6.1,
        r=0,
        n_r=2,
        material=steel,
    )

    # timber floor panel
    panel = sections.rectangular_section(d=100, b=600, material=timber)
    panel = panel.align_center(align_to=ub).align_to(other=ub, on="top")

    # combine geometry
    geom = ub + panel

    geom.create_mesh(mesh_sizes=[50, 500])
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()
    sec.calculate_plastic_properties()

    # get calculated values
    x_pc, y_pc = sec.get_pc()
    mp_xx, mp_yy = sec.get_mp()

    # hand calcs
    # 1) x-axis bending
    # calculate y_pc -> assume centroid in top flange
    f_timber = 600 * 100 * 20
    f_web = (304 - 2 * 10.2) * 6.1 * 500
    f_flange = 165 * 10.2 * 500

    # let y = distance from bottom of top flange to y_pc
    # A) bot = 165 * y * 500
    # B) top = 165 * (10.2 - y) * 500
    # C) (A) - (B) = f_timber - f_web - f_flange  = f_diff
    # 165 * 500 (2y - 10.2) = f_diff
    # y = (f_diff / (165 * 500) + 10.2) / 2
    f_diff = f_timber - f_web - f_flange
    y = (f_diff / (165 * 500) + 10.2) / 2
    assert y_pc == pytest.approx(304 - 10.2 + y)

    # calculate mp_xx
    mp_xx_calc = 0
    mp_xx_calc += f_timber * abs(304 + 50 - y_pc)  # timber contribution
    mp_xx_calc += (10.2 - y) * 165 * 500 * abs((10.2 - y) / 2)  # top flange p1
    mp_xx_calc += y * 165 * 500 * abs(y / 2)  # top flange p2
    mp_xx_calc += f_web * abs(10.2 + (304 - 2 * 10.2) / 2 - y_pc)  # web contribution
    mp_xx_calc += f_flange * abs(10.2 / 2 - y_pc)  # bot flange contribution
    assert mp_xx == pytest.approx(mp_xx_calc)

    # 2) y-axis bending
    assert x_pc == pytest.approx(165 / 2)

    # calculate mp_yy
    mp_yy_calc = 0
    mp_yy_calc += 2 * (10.2 * 165**2 / 4) * 500  # 2 flanges
    mp_yy_calc += (304 - 10.2 * 2) * 6.1**2 / 4 * 500  # web
    mp_yy_calc += 100 * 600**2 / 4 * 20  # timber
    assert mp_yy == pytest.approx(mp_yy_calc)
