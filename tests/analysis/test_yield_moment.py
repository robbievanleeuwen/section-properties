"""Tests for the yield moment calculation."""

import pytest

from sectionproperties.analysis import Section
from sectionproperties.pre import Material
from sectionproperties.pre.library import rectangular_section

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
    geom = rectangular_section(d=50, b=50)
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

    my_xx_calc = 50 * 100 * 100 / 4 * 500
    my_yy_calc = 100 * 50 * 50 / 4 * 500

    assert my_xx == pytest.approx(my_xx_calc)
    assert my_yy == pytest.approx(my_yy_calc)
    assert my_11 == pytest.approx(my_xx_calc)
    assert my_22 == pytest.approx(my_yy_calc)


# TODO: add more tests
