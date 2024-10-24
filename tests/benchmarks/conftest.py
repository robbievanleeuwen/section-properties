"""pytest benchmark configurations."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest

from sectionproperties.pre import Geometry, Material
from sectionproperties.pre.library import (
    circular_hollow_section,
    circular_section,
    concrete_column_section,
    rectangular_section,
)

if TYPE_CHECKING:
    from collections.abc import Callable


@pytest.fixture
def concrete() -> Material:
    """Creates a concrete material object.

    Returns:
        Concrete
    """
    return Material(
        name="Concrete",
        elastic_modulus=30.1e3,
        poissons_ratio=0.2,
        yield_strength=32,
        density=2.4e-6,
        color="lightgrey",
    )


@pytest.fixture
def steel() -> Material:
    """Creates a steel material object.

    Returns:
        Steel
    """
    return Material(
        name="Steel",
        elastic_modulus=200e3,
        poissons_ratio=0.3,
        yield_strength=500,
        density=7.85e-6,
        color="grey",
    )


@pytest.fixture
def rect_geom() -> Geometry:
    """Creates a rectangular geometry.

    Returns:
        Geometry
    """
    return rectangular_section(d=100, b=50)


@pytest.fixture
def chs_geom() -> Geometry:
    """Creates a rectangular geometry.

    Returns:
        Geometry
    """
    return circular_hollow_section(d=100, t=3, n=128)


@pytest.fixture
def concrete_column_with_hole(concrete, steel) -> Callable:
    """Creates a concrete column with a hole at its centre.

    Args:
        concrete: Concrete material
        steel: Steel material

    Returns:
        Generator function
    """

    def _generate_geom() -> Geometry:
        geom = concrete_column_section(
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
        hole = circular_section(d=100, n=32).shift_section(x_offset=150, y_offset=300)

        return geom - hole

    return _generate_geom


@pytest.fixture
def analysis_geometry() -> Callable:
    """Create a geometry to be used for analysis.

    Returns:
        Generator function
    """

    def _generate_geom(num_elements: int) -> Geometry:
        mat_a = Material("a", 1, 0, 1, 1, color="b")
        mat_b = Material("b", 10, 0, 1, 1, color="g")
        mat_c = Material("c", 5, 0, 1, 1, color="r")
        mat_d = Material("d", 2, 0, 1, 1, color="y")

        a = rectangular_section(20, 20, mat_a)
        b = rectangular_section(20, 20, mat_b).align_to(a, "right")
        c = rectangular_section(20, 20, mat_c).align_to(a, "top")
        d = rectangular_section(20, 20, mat_d).align_to(a, "top").align_to(a, "right")
        geom = a + b + c + d
        mesh_area = geom.calculate_area() / num_elements * 1.6
        return geom.create_mesh([mesh_area])

    return _generate_geom
