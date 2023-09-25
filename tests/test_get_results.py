"""Validation of cross-section property get methods."""

from __future__ import annotations

import pytest

from sectionproperties.analysis import Section
from sectionproperties.pre import Material
from sectionproperties.pre.library import rectangular_section


# sections setup
geom_no_mat = rectangular_section(d=1, b=1)
dummy_mat = Material(
    name="test",
    elastic_modulus=5,
    poissons_ratio=0,
    yield_strength=3,
    density=2,
    color="w",
)
geom_mat = rectangular_section(d=1, b=1, material=dummy_mat)
geom_no_mat.create_mesh(mesh_sizes=0)
geom_mat.create_mesh(mesh_sizes=0)
rect_no_mat = Section(geom_no_mat)
rect_mat = Section(geom_mat)


def test_is_composite():
    """Check whether section is composite or not."""
    assert not rect_no_mat.is_composite()
    assert rect_mat.is_composite()


def test_get_e_ref():
    """Check get_e_ref results."""
    test_mat = Material(
        name="test",
        elastic_modulus=10,
        poissons_ratio=0,
        yield_strength=1,
        density=1,
        color="w",
    )

    assert rect_no_mat.get_e_ref(1) == 1
    assert rect_no_mat.get_e_ref(100) == 100
    assert rect_no_mat.get_e_ref(dummy_mat) == 5
    assert rect_no_mat.get_e_ref(test_mat) == 10


def test_no_analysis():
    """Check errors when no analysis has been conducted.

    RuntimeError = incorrect analysis type (geometric only vs. composite), takes
    precedence over...
    AssertionError = relevant analysis has been conducted
    """
    # check area
    with pytest.raises(AssertionError):
        rect_no_mat.get_area()

    # check perimeter
    with pytest.raises(AssertionError):
        rect_no_mat.get_perimeter()

    # check mass
    with pytest.raises(RuntimeError):
        rect_no_mat.get_mass()

    with pytest.raises(AssertionError):
        rect_mat.get_mass()

    # check ea
    with pytest.raises(RuntimeError):
        rect_no_mat.get_mass()

    with pytest.raises(AssertionError):
        rect_mat.get_mass()

    # check q
    with pytest.raises(AssertionError):
        rect_no_mat.get_q()

    with pytest.raises(RuntimeError):
        rect_mat.get_q()

    # check eq
    with pytest.raises(RuntimeError):
        rect_no_mat.get_eq()

    with pytest.raises(AssertionError):
        rect_mat.get_eq()

    # check ig
    with pytest.raises(AssertionError):
        rect_no_mat.get_ig()

    with pytest.raises(RuntimeError):
        rect_mat.get_ig()

    # check eig
    with pytest.raises(RuntimeError):
        rect_no_mat.get_eig()

    with pytest.raises(AssertionError):
        rect_mat.get_eig()


def test_get_geometric_only():
    """Check errors and results when a geometric analysis has been conducted.

    RuntimeError = incorrect analysis type (geometric only vs. composite)
    """
    rect_no_mat.calculate_geometric_properties()
    rect_mat.calculate_geometric_properties()

    # check area
    assert rect_no_mat.get_area() == pytest.approx(1.0)
    assert rect_mat.get_area() == pytest.approx(1.0)

    # check perimeter
    assert rect_no_mat.get_perimeter() == pytest.approx(4.0)
    assert rect_mat.get_perimeter() == pytest.approx(4.0)

    # check mass
    with pytest.raises(RuntimeError):
        rect_no_mat.get_mass()

    assert rect_mat.get_mass() == pytest.approx(2.0)

    # check ea
    with pytest.raises(RuntimeError):
        rect_no_mat.get_ea()

    assert rect_mat.get_ea() == pytest.approx(5.0)
    assert rect_mat.get_ea(e_ref=2) == pytest.approx(2.5)
    assert rect_mat.get_ea(e_ref=dummy_mat) == pytest.approx(1.0)

    # check q
    qx, qy = rect_no_mat.get_q()
    assert qx == pytest.approx(0.5)
    assert qy == pytest.approx(0.5)

    with pytest.raises(RuntimeError):
        rect_mat.get_q()

    # check eq
    with pytest.raises(RuntimeError):
        rect_no_mat.get_eq()

    eqx, eqy = rect_mat.get_eq()
    assert eqx == pytest.approx(2.5)
    assert eqy == pytest.approx(2.5)
    eqx, eqy = rect_mat.get_eq(e_ref=2)
    assert eqx == pytest.approx(1.25)
    assert eqy == pytest.approx(1.25)
    eqx, eqy = rect_mat.get_eq(e_ref=dummy_mat)
    assert eqx == pytest.approx(0.5)
    assert eqy == pytest.approx(0.5)

    # check ig
    igxx, igyy, igxy = rect_no_mat.get_ig()
    assert igxx == pytest.approx(1 / 3.0)
    assert igyy == pytest.approx(1 / 3.0)
    assert igxy == pytest.approx(0.25)

    with pytest.raises(RuntimeError):
        rect_mat.get_ig()

    # check eig
    with pytest.raises(RuntimeError):
        rect_no_mat.get_eig()

    eigxx, eigyy, eigxy = rect_mat.get_eig()
    assert eigxx == pytest.approx(5.0 / 3)
    assert eigyy == pytest.approx(5.0 / 3)
    assert eigxy == pytest.approx(1.25)
    eigxx, eigyy, eigxy = rect_mat.get_eig(e_ref=2)
    assert eigxx == pytest.approx(5.0 / 6)
    assert eigyy == pytest.approx(5.0 / 6)
    assert eigxy == pytest.approx(0.625)
    eigxx, eigyy, eigxy = rect_mat.get_eig(e_ref=dummy_mat)
    assert eigxx == pytest.approx(1 / 3.0)
    assert eigyy == pytest.approx(1 / 3.0)
    assert eigxy == pytest.approx(0.25)
