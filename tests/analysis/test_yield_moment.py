"""Tests for the yield moment calculation."""

import numpy as np
import pytest

from sectionproperties.analysis import Section
from sectionproperties.post.stress_post import StressPost
from sectionproperties.pre import Material
from sectionproperties.pre.library import i_section, rectangular_section

STEEL = Material(
    name="Steel",
    elastic_modulus=200e3,
    poissons_ratio=0.3,
    yield_strength=500,
    density=7.85e-6,
    color="grey",
)


def test_get_without_analysis():
    """Test for raising an error if a geometric analysis has not been performed."""
    geom = rectangular_section(d=50, b=50, material=STEEL)
    geom.create_mesh(mesh_sizes=0, coarse=True)
    sec = Section(geometry=geom)

    with pytest.raises(RuntimeError, match="Conduct a geometric analysis."):
        sec.get_my()

    with pytest.raises(RuntimeError, match="Conduct a geometric analysis."):
        sec.get_my_p()


def test_non_composite():
    """Test for raising an error for non-composite analyses."""
    geom = rectangular_section(d=50, b=50)
    geom.create_mesh(mesh_sizes=0, coarse=True)
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()

    with pytest.raises(RuntimeError, match="Attempting to get a composite only"):
        sec.get_my()

    with pytest.raises(RuntimeError, match="Attempting to get a composite only"):
        sec.get_my_p()


def test_rectangle():
    """Test the yield moment of a simple rectangle."""
    geom = rectangular_section(d=100, b=50, material=STEEL)
    geom.create_mesh(mesh_sizes=0, coarse=True)
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()

    my_xx, my_yy = sec.get_my()
    my_11, my_22 = sec.get_my()
    my_xx_calc = 50 * 100 * 100 / 6 * 500
    my_yy_calc = 100 * 50 * 50 / 6 * 500

    # compare with theoretical values
    assert my_xx == pytest.approx(my_xx_calc)
    assert my_yy == pytest.approx(my_yy_calc)
    assert my_11 == pytest.approx(my_xx_calc)
    assert my_22 == pytest.approx(my_yy_calc)

    # compare with stress analysis
    stress = sec.calculate_stress(mxx=my_xx, myy=my_yy, m11=my_11, m22=my_22)
    assert check_yield_index(stress, "sig_zz_mxx") == pytest.approx(1.0)
    assert check_yield_index(stress, "sig_zz_myy") == pytest.approx(1.0)
    assert check_yield_index(stress, "sig_zz_m11") == pytest.approx(1.0)
    assert check_yield_index(stress, "sig_zz_m22") == pytest.approx(1.0)

    # check shape factors
    sec.calculate_plastic_properties()
    sf_xx, _, sf_yy, _ = sec.get_sf()
    sf_11, _, sf_22, _ = sec.get_sf_p()
    assert sf_xx == pytest.approx(1.5)
    assert sf_yy == pytest.approx(1.5)
    assert sf_11 == pytest.approx(1.5)
    assert sf_22 == pytest.approx(1.5)


def test_rectangle_rotated():
    """Test the yield moment of a simple rotated rectangle."""
    geom = rectangular_section(d=100, b=50, material=STEEL)
    geom.rotate_section(angle=45.0)
    geom.create_mesh(mesh_sizes=0, coarse=True)
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()

    my_11, my_22 = sec.get_my()
    my_xx_calc = 50 * 100 * 100 / 6 * 500
    my_yy_calc = 100 * 50 * 50 / 6 * 500

    # compare with theoretical values
    assert my_11 == pytest.approx(my_xx_calc)
    assert my_22 == pytest.approx(my_yy_calc)

    # compare with stress analysis
    stress = sec.calculate_stress(m11=my_11, m22=my_22)
    assert check_yield_index(stress, "sig_zz_m11") == pytest.approx(1.0)
    assert check_yield_index(stress, "sig_zz_m22") == pytest.approx(1.0)

    # check shape factors
    sec.calculate_plastic_properties()
    sf_11, _, sf_22, _ = sec.get_sf_p()
    assert sf_11 == pytest.approx(1.5)
    assert sf_22 == pytest.approx(1.5)


def test_isection():
    """Test the yield moment of an isection."""
    geom = i_section(d=200, b=100, t_f=10, t_w=5, r=12, n_r=8, material=STEEL)
    geom.create_mesh(mesh_sizes=0, coarse=True)
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()

    my_xx, my_yy = sec.get_my()
    my_11, my_22 = sec.get_my()
    ezxx, _, ezyy, _ = sec.get_ez()
    my_xx_calc = ezxx / 200e3 * 500
    my_yy_calc = ezyy / 200e3 * 500

    # compare with theoretical values
    assert my_xx == pytest.approx(my_xx_calc)
    assert my_yy == pytest.approx(my_yy_calc)
    assert my_11 == pytest.approx(my_xx_calc)
    assert my_22 == pytest.approx(my_yy_calc)

    # compare with stress analysis
    stress = sec.calculate_stress(mxx=my_xx, myy=my_yy, m11=my_11, m22=my_22)
    assert check_yield_index(stress, "sig_zz_mxx") == pytest.approx(1.0)
    assert check_yield_index(stress, "sig_zz_myy") == pytest.approx(1.0)
    assert check_yield_index(stress, "sig_zz_m11") == pytest.approx(1.0)
    assert check_yield_index(stress, "sig_zz_m22") == pytest.approx(1.0)

    # check that shape factors with a geometric-only analysis match the composite
    sec.calculate_plastic_properties()
    sf_xx_c, _, sf_yy_c, _ = sec.get_sf()

    geom = i_section(d=200, b=100, t_f=10, t_w=5, r=12, n_r=8)
    geom.create_mesh(mesh_sizes=0, coarse=True)
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()
    sec.calculate_plastic_properties()

    sf_xx_plus, sf_xx_minus, sf_yy_plus, sf_yy_minus = sec.get_sf()
    assert sf_xx_c == pytest.approx(sf_xx_plus)
    assert sf_xx_c == pytest.approx(sf_xx_minus)
    assert sf_yy_c == pytest.approx(sf_yy_plus)
    assert sf_yy_c == pytest.approx(sf_yy_minus)


def test_rectangle_composite():
    """Test the yield moment of a composite rectangular section."""
    mat1 = Material("a", 2, 0, 2, 1, "b")
    mat2 = Material("b", 1, 0, 1, 1, "r")

    rect1 = rectangular_section(d=20, b=20, material=mat1)
    rect2 = rectangular_section(d=20, b=20, material=mat2).align_to(rect1, "top")
    rect3 = rectangular_section(d=20, b=20, material=mat1).align_to(rect2, "top")
    geom = rect1 + rect2 + rect3
    geom.create_mesh(mesh_sizes=0, coarse=True)
    sec = Section(geom)
    sec.calculate_geometric_properties()

    my_xx, my_yy = sec.get_my()
    my_11, my_22 = sec.get_my()
    eixx_calc = (20 * 20**3 / 12) + 2 * 2 * ((20 * 20**3 / 12) + (20 * 20 * 20**2))
    eiyy_calc = 5 * (20 * 20**3 / 12)
    # myield = f * Ieff / y
    my_xx_calc = 2 * (eixx_calc / 2) / 30
    my_yy_calc = 1 * eiyy_calc / 10

    # compare with theoretical values
    assert my_xx == pytest.approx(my_xx_calc)
    assert my_yy == pytest.approx(my_yy_calc)
    assert my_11 == pytest.approx(my_xx_calc)
    assert my_22 == pytest.approx(my_yy_calc)

    # compare with stress analysis
    stress = sec.calculate_stress(mxx=my_xx, myy=my_yy, m11=my_11, m22=my_22)
    assert check_yield_index(stress, "sig_zz_mxx") == pytest.approx(1.0)
    assert check_yield_index(stress, "sig_zz_myy") == pytest.approx(1.0)
    assert check_yield_index(stress, "sig_zz_m11") == pytest.approx(1.0)
    assert check_yield_index(stress, "sig_zz_m22") == pytest.approx(1.0)


def test_composite_example():
    """Tests the composite example from the docs, compares to stress analysis."""
    timber = Material(
        name="Timber",
        elastic_modulus=8e3,
        poissons_ratio=0.35,
        yield_strength=20,
        density=0.78e-6,
        color="burlywood",
    )

    ub = i_section(d=304, b=165, t_f=10.2, t_w=6.1, r=11.4, n_r=8, material=STEEL)
    panel = rectangular_section(d=100, b=600, material=timber)
    panel = panel.align_center(align_to=ub).align_to(other=ub, on="top")
    geom = ub + panel
    geom.create_mesh(mesh_sizes=[10, 500])
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()

    # compare with stress analysis
    my_xx, my_yy = sec.get_my()
    my_11, my_22 = sec.get_my()
    stress = sec.calculate_stress(mxx=my_xx, myy=my_yy, m11=my_11, m22=my_22)
    assert check_yield_index(stress, "sig_zz_mxx") == pytest.approx(1.0)
    assert check_yield_index(stress, "sig_zz_myy") == pytest.approx(1.0)
    assert check_yield_index(stress, "sig_zz_m11") == pytest.approx(1.0)
    assert check_yield_index(stress, "sig_zz_m22") == pytest.approx(1.0)


def test_yield_internal():
    """Test where the internal material yields first."""
    mat1 = Material(
        name="Material_1",
        elastic_modulus=0.75,
        poissons_ratio=0,
        yield_strength=1,
        density=1,
        color="gold",
    )
    mat2 = Material(
        name="Material_2",
        elastic_modulus=1,
        poissons_ratio=0,
        yield_strength=1,
        density=1,
        color="blue",
    )
    mat3 = Material(
        name="Material 3",
        elastic_modulus=3,
        poissons_ratio=0,
        yield_strength=1,
        density=1,
        color="red",
    )

    sq1 = rectangular_section(d=100, b=100, material=mat1).align_center()
    sq2 = rectangular_section(d=75, b=75, material=mat2).align_center()
    sq3 = rectangular_section(d=50, b=50, material=mat3).align_center()
    hole = rectangular_section(d=25, b=25).align_center()
    compound = (sq1 - sq2) + (sq2 - sq3) + (sq3 - hole)
    compound.create_mesh(10)
    sec = Section(compound)
    sec.calculate_geometric_properties()

    # compare with stress analysis
    my_xx, my_yy = sec.get_my()
    my_11, my_22 = sec.get_my()
    stress = sec.calculate_stress(mxx=my_xx, myy=my_yy, m11=my_11, m22=my_22)
    assert check_yield_index(stress, "sig_zz_mxx") == pytest.approx(1.0)
    assert check_yield_index(stress, "sig_zz_myy") == pytest.approx(1.0)
    assert check_yield_index(stress, "sig_zz_m11") == pytest.approx(1.0)
    assert check_yield_index(stress, "sig_zz_m22") == pytest.approx(1.0)


def check_yield_index(stress: StressPost, key: str) -> float:
    """Returns the largest yield index.

    Given the output from StressPost object and a dict key representing the type of
    stress, returns the largest yield index.
    """
    yield_index = 0.0

    for idx, mat_dict in enumerate(stress.get_stress()):
        fy = stress.material_groups[idx].material.yield_strength
        sigs = mat_dict[key]
        sig_max = max(np.absolute(sigs))

        yield_index = max(yield_index, sig_max / fy)

    return yield_index
