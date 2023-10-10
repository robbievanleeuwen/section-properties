"""Benchmark tests for sectionproperties."""

from __future__ import annotations

import pytest
import pytest_benchmark

from sectionproperties.pre import Material
from sectionproperties.pre.library import (
    circular_hollow_section,
    concrete_column_section,
    rectangular_section,
)


@pytest.fixture
def conc_steel() -> tuple[Material, Material]:
    """Creates a concrete and steel material object.

    Returns:
        Concrete and steel
    """
    concrete = Material(
        name="Concrete",
        elastic_modulus=30.1e3,
        poissons_ratio=0.2,
        yield_strength=32,
        density=2.4e-6,
        color="lightgrey",
    )
    steel = Material(
        name="Steel",
        elastic_modulus=200e3,
        poissons_ratio=0.3,
        yield_strength=500,
        density=7.85e-6,
        color="grey",
    )

    return concrete, steel


def test_create_rectangle_geometry(benchmark):
    """Benchmark test for creating rectangular geometry."""
    benchmark(rectangular_section, d=100, b=50)


def test_create_chs_geometry(benchmark):
    """Benchmark test for creating CHS geometry."""
    benchmark(circular_hollow_section, d=100, t=3, n=128)


def test_create_concrete_geometry(benchmark, conc_steel):
    """Benchmark test for creating concrete geometry."""
    concrete, steel = conc_steel

    benchmark(
        concrete_column_section,
        d=600,
        b=300,
        dia_bar=25,
        area_bar=500,
        n_x=3,
        n_y=6,
        cover=35,
        n_circle=4,
        filled=False,
        conc_mat=concrete,
        steel_mat=steel,
    )
