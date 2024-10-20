"""Tests for the concrete sections library."""

from __future__ import annotations

import numpy as np
import pytest
import pytest_check as check

import sectionproperties.pre.library.concrete_sections as cs
import sectionproperties.pre.library.primitive_sections as ps
import sectionproperties.pre.pre as pre

r_tol = 1e-6


# material setup
@pytest.fixture
def get_materials() -> tuple[pre.Material, pre.Material]:
    """Creates a concrete and steel material.

    Returns:
        Material objects
    """
    conc_mat = pre.Material(
        name="Concrete",
        elastic_modulus=32.8e3,
        poissons_ratio=0.2,
        density=2.4e-6,
        yield_strength=40,
        color="lightgrey",
    )

    steel_mat = pre.Material(
        name="Steel",
        elastic_modulus=200e3,
        poissons_ratio=0.3,
        density=7.85e-6,
        yield_strength=500,
        color="grey",
    )

    return conc_mat, steel_mat


def test_concrete_rectangular_section(get_materials):
    """Tests the concrete_rectangular_section() method."""
    conc_mat, steel_mat = get_materials

    rect = cs.concrete_rectangular_section(
        d=600,
        b=300,
        dia_top=16,
        area_top=200,
        n_top=3,
        c_top=30,
        dia_bot=20,
        area_bot=310,
        n_bot=3,
        c_bot=30,
        n_circle=16,
        conc_mat=conc_mat,
        steel_mat=steel_mat,
    )

    # check geometry is created correctly
    conc_area = 0
    steel_area = 0

    for geom in rect.geoms:
        if geom.material == conc_mat:
            conc_area += geom.calculate_area()
        elif geom.material == steel_mat:
            steel_area += geom.calculate_area()

    net_area = 600 * 300
    actual_steel_area = 3 * (200 + 310)

    # check areas
    check.almost_equal(conc_area, net_area - actual_steel_area, rel=r_tol)
    check.almost_equal(steel_area, actual_steel_area, rel=r_tol)


def test_concrete_tee_section(get_materials):
    """Tests the concrete_tee_section() method."""
    conc_mat, steel_mat = get_materials

    rect = cs.concrete_tee_section(
        d=900,
        b=300,
        d_f=200,
        b_f=1200,
        dia_top=20,
        area_top=310,
        n_top=6,
        c_top=30,
        dia_bot=24,
        area_bot=450,
        n_bot=3,
        c_bot=30,
        n_circle=16,
        conc_mat=conc_mat,
        steel_mat=steel_mat,
    )

    # check geometry is created correctly
    conc_area = 0
    steel_area = 0

    for geom in rect.geoms:
        if geom.material == conc_mat:
            conc_area += geom.calculate_area()
        elif geom.material == steel_mat:
            steel_area += geom.calculate_area()

    net_area = 700 * 300 + 1200 * 200
    actual_steel_area = 6 * 310 + 3 * 450

    # check areas
    check.almost_equal(conc_area, net_area - actual_steel_area, rel=r_tol)
    check.almost_equal(steel_area, actual_steel_area, rel=r_tol)


def test_concrete_circular_section(get_materials):
    """Tests the concrete_circular_section() method."""
    conc_mat, steel_mat = get_materials

    rect = cs.concrete_circular_section(
        d=600,
        area_conc=np.pi * 600 * 600 / 4,
        n_conc=64,
        dia_bar=20,
        area_bar=310,
        n_bar=8,
        cover=45,
        n_circle=16,
        conc_mat=conc_mat,
        steel_mat=steel_mat,
    )

    # check geometry is created correctly
    conc_area = 0
    steel_area = 0

    for geom in rect.geoms:
        if geom.material == conc_mat:
            conc_area += geom.calculate_area()
        elif geom.material == steel_mat:
            steel_area += geom.calculate_area()

    net_area = np.pi * 600 * 600 / 4
    actual_steel_area = 8 * 310

    # check areas
    check.almost_equal(conc_area, net_area - actual_steel_area, rel=r_tol)
    check.almost_equal(steel_area, actual_steel_area, rel=r_tol)


def test_concrete_column_section(get_materials):
    """Tests the concrete_column_section() method."""
    conc_mat, steel_mat = get_materials

    concrete = pre.Material(
        name="Concrete",
        elastic_modulus=30.1e3,
        poissons_ratio=0.2,
        yield_strength=32,
        density=2.4e-6,
        color="lightgrey",
    )
    steel = pre.Material(
        name="Steel",
        elastic_modulus=200e3,
        poissons_ratio=0.3,
        yield_strength=500,
        density=7.85e-6,
        color="grey",
    )

    geometry = cs.concrete_column_section(
        d=600,
        b=300,
        dia_bar=40,
        area_bar=500,
        n_x=3,
        n_y=6,
        cover=40,
        n_circle=4,
        steel_mat=steel,
        conc_mat=concrete,
        filled=False,
    )  # NOTE: Bar diam and Bar area do not match. This is intentional.
    geometry.create_mesh(mesh_sizes=[500])

    # check geometry is created correctly
    conc_area = 0
    steel_area = 0

    for geom in geometry.geoms:
        if geom.material == concrete:
            conc_area += geom.calculate_area()
        elif geom.material == steel:
            steel_area += geom.calculate_area()

    net_area = 300 * 600
    actual_steel_area = 14 * 500.0

    check.almost_equal(conc_area, net_area - actual_steel_area, rel=r_tol)
    check.almost_equal(steel_area, actual_steel_area, rel=r_tol)

    bar_centroids = [tuple(geom.geom.centroid.coords[0]) for geom in geometry.geoms[1:]]

    from collections import Counter

    x_coords = Counter(round(coord[0], 0) for coord in bar_centroids)
    y_coords = Counter(round(coord[1], 0) for coord in bar_centroids)

    # Validate that we have 14 bars with the correct x-coordinates
    check.equal(x_coords.get(60), 6)
    check.equal(x_coords.get(150), 2)
    check.equal(x_coords.get(240), 6)

    # Validate that we have 14 bars with the correct y-coordinates
    check.equal(y_coords.get(60), 3)
    check.equal(y_coords.get(156), 2)
    check.equal(y_coords.get(252), 2)
    check.equal(y_coords.get(348), 2)
    check.equal(y_coords.get(444), 2)
    check.equal(y_coords.get(540), 3)


def test_add_bar():
    """Tests the add_bar() method."""
    rect = ps.rectangular_section(b=400, d=600)
    steel = pre.Material(
        name="Steel",
        elastic_modulus=200e3,
        poissons_ratio=0.3,
        yield_strength=500,
        density=7.85e-6,
        color="grey",
    )
    rect = cs.add_bar(geometry=rect, area=500, x=100, y=100, n=4, material=steel)
    rect = cs.add_bar(geometry=rect, area=500, x=200, y=200, n=4, material=steel)
    rect_area = rect.geoms[0].geom.area
    steel_area = 2 * 500.0
    check.almost_equal(rect_area, 400 * 600 - steel_area)

    bar_1 = rect.geoms[1]
    bar_2 = rect.geoms[2]
    check.almost_equal(bar_1.calculate_centroid(), (100, 100))
    check.almost_equal(bar_2.calculate_centroid(), (200, 200))
