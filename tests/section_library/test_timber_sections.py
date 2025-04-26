"""Tests for the timber sections library."""

from __future__ import annotations

import pytest
import pytest_check as check

import sectionproperties.pre.library.timber_sections as ts
import sectionproperties.pre.pre as pre

r_tol = 1e-6


# material setup
@pytest.fixture
def get_materials() -> tuple[pre.Material, pre.Material]:
    """Creates a timber material parallel and perpendicular-to-grain.

    Returns:
        Material objects
    """
    timb_mat0 = pre.Material(
        name="Timber E0",
        elastic_modulus=9.5e3,
        poissons_ratio=0.35,
        density=4.4e-7,
        yield_strength=5.5,
        color="burlywood",
    )

    timb_mat90 = pre.Material(
        name="Timber90",
        elastic_modulus=317,
        poissons_ratio=0.35,
        density=4.4e-7,
        yield_strength=5.5,
        color="orange",
    )

    return timb_mat0, timb_mat90


def test_timber_clt_rectangular_section(get_materials):
    """Tests the timber clt_rectangular_section() method."""
    timber0, timber90 = get_materials

    rect = ts.clt_rectangular_section(
        d=[40, 40, 40], layer_mat=[timber0, timber90, timber0], b=1000
    )

    # check geometry is created correctly
    timb0_area = 0
    timb90_area = 0

    for geom in rect.geoms:
        if geom.material == timber0:
            timb0_area += geom.calculate_area()
        elif geom.material == timber90:
            timb90_area += geom.calculate_area()

    actual_timb0_area = 2 * 40 * 1000
    actual_timb90_area = 40 * 1000

    # check areas
    check.almost_equal(timb0_area, actual_timb0_area)
    check.almost_equal(timb90_area, actual_timb90_area)
