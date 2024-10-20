"""Tests for secionpropertiers.post.post."""

from __future__ import annotations

import platform

import matplotlib.pyplot as plt
import pytest
import pytest_check as check

from sectionproperties.analysis import Section
from sectionproperties.pre import Material
from sectionproperties.pre.library import rectangular_section

linux_only = pytest.mark.skipif(
    platform.system() != "Linux",
    reason="Only test plotting on Linux",
)


@pytest.fixture
def example_section() -> Section:
    """Creates an example section with geometric properties.

    Returns:
        Section
    """
    geom = rectangular_section(d=1, b=1)
    geom.create_mesh(mesh_sizes=0)

    return Section(geometry=geom)


def test_as_dict(example_section):
    """Test as_dict method."""
    sec = example_section
    sec.calculate_geometric_properties()
    res = sec.section_props.asdict()

    check.almost_equal(res["area"], sec.get_area())


def test_elastic_centroid_error(example_section):
    """Test the ValueError when calculating elastic centroids."""
    sec = example_section

    with pytest.raises(RuntimeError, match="Calculate geometric properties first"):
        sec.section_props.calculate_elastic_centroid()


def test_centroidal_properties_error(example_section):
    """Test the ValueError when calculating centroidal_properties."""
    sec = example_section

    with pytest.raises(RuntimeError, match="Calculate geometric properties first"):
        sec.section_props.calculate_centroidal_properties(
            node_list=sec.mesh["vertices"],
        )


@linux_only
def test_save_plot(example_section, tmp_path):
    """Tests saving a plot."""
    sec = example_section
    d = tmp_path / "sub"
    d.mkdir()

    sec.plot_mesh(filename=d / "fig.png")
    plt.close("all")


@linux_only
def test_supplied_axis(example_section):
    """Tests supplying an axis to a plot."""
    sec = example_section
    _, ax = plt.subplots()

    sec.plot_mesh(ax=ax, render=False)
    plt.close("all")
    sec.plot_mesh(nrows=2, axis_index=1, render=False)
    plt.close("all")

    with pytest.raises(ValueError, match="is not compatible"):
        sec.plot_mesh(nrows=2, ncols=2, axis_index=5, render=False)

    plt.close("all")


@linux_only
def test_plot_centroids(example_section):
    """Tests plotting centroids."""
    sec = example_section
    sec.calculate_geometric_properties()
    sec.plot_centroids(render=False)
    plt.close("all")


def test_print_results(example_section):
    """Tests printing results."""
    sec = example_section
    sec.display_results()
    sec.calculate_geometric_properties()
    sec.calculate_warping_properties()
    sec.calculate_plastic_properties()
    sec.display_results()

    mat = Material("a", 1, 0.2, 2, 1, "k")
    geom_mat = rectangular_section(d=1, b=1, material=mat)
    geom_mat.create_mesh(mesh_sizes=0)
    sec_mat = Section(geometry=geom_mat)
    sec_mat.display_results()
    sec_mat.calculate_geometric_properties()
    sec_mat.calculate_warping_properties()
    sec_mat.calculate_plastic_properties()
    sec_mat.display_results()
