"""Tests for secionpropertiers.post.stress_post."""

from __future__ import annotations

import platform

import pytest

from sectionproperties.analysis import Section
from sectionproperties.pre import Material
from sectionproperties.pre.library import rectangular_section

linux_only = pytest.mark.skipif(
    platform.system() != "Linux",
    reason="Only test plotting on Linux",
)


@pytest.fixture
def example_section() -> tuple[Section, Material]:
    """Creates an example section with geometric properties.

    Returns:
        Section and material
    """
    mat_a = Material("a", 1, 0.2, 2, 1, "k")
    mat_b = Material("b", 2, 0.2, 2, 1, "r")
    rect_a = rectangular_section(d=1, b=1, material=mat_a)
    rect_b = rectangular_section(d=1, b=1, material=mat_b).align_to(rect_a, "right")
    geom = rect_a + rect_b
    geom.create_mesh(mesh_sizes=[0])

    return Section(geometry=geom), mat_b


@linux_only
def test_stress_plot(example_section):
    """Tests the plot_stress() method."""
    sec, mat_b = example_section
    sec.calculate_geometric_properties()
    stress = sec.calculate_stress(n=1)
    stress.plot_stress(stress="n_zz", render=False)
    stress.plot_stress(stress="n_zz", stress_limits=[0, 1], render=False)
    stress.plot_stress(stress="n_zz", material_list=[mat_b], render=False)


@linux_only
def test_stress_plot_constant():
    """Tests the plot_stress() method with a constant stress."""
    geom = rectangular_section(d=1, b=1)
    geom.create_mesh(mesh_sizes=[0])
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()
    sec.calculate_stress(n=1).plot_stress(stress="n_zz", render=False)


@linux_only
def test_stress_vector_plot(example_section):
    """Tests the plot_stress_vector() method."""
    sec, _ = example_section
    sec.calculate_geometric_properties()
    sec.calculate_warping_properties()
    sec.calculate_stress(mzz=1).plot_stress_vector(stress="zxy", render=False)


@linux_only
def test_plot_mohrs_circles(example_section):
    """Tests the plot_mohrs_circles() method."""
    sec, _ = example_section
    sec.calculate_geometric_properties()
    stress = sec.calculate_stress(mxx=1)
    stress.plot_mohrs_circles(x=1, y=0.5, render=False)

    with pytest.raises(ValueError, match="is not within mesh"):
        stress.plot_mohrs_circles(x=2, y=2, render=False)
