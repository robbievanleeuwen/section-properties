"""Benchmark tests for sectionproperties analysis."""

import pytest

from sectionproperties.analysis import Section


@pytest.mark.benchmark_suite
@pytest.mark.benchmark_analysis
@pytest.mark.parametrize("elements", [50, 500, 5000])
def test_create_section(benchmark, analysis_geometry, elements):
    """Benchmark test for creating a Section object."""
    geom = analysis_geometry(elements)

    def create_section():
        Section(geometry=geom)

    benchmark(create_section)


@pytest.mark.benchmark_suite
@pytest.mark.benchmark_analysis
@pytest.mark.parametrize("elements", [50, 500, 5000])
def test_geometric_analysis(benchmark, analysis_geometry, elements):
    """Benchmark test for conducting a geometric analysis."""
    geom = analysis_geometry(elements)
    sec = Section(geometry=geom)

    def geometric_analysis():
        sec.calculate_geometric_properties()

    benchmark(geometric_analysis)


@pytest.mark.benchmark_suite
@pytest.mark.benchmark_analysis
def test_plastic_analysis(benchmark, analysis_geometry):
    """Benchmark test for conducting a plastic analysis.

    Note that a plastic analysis is mesh-independent.
    """
    geom = analysis_geometry(1)
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()

    def plastic_analysis():
        sec.calculate_plastic_properties()

    benchmark(plastic_analysis)


@pytest.mark.benchmark_suite
@pytest.mark.benchmark_analysis
@pytest.mark.parametrize("elements", [50, 500, 5000])
def test_warping_analysis(benchmark, analysis_geometry, elements):
    """Benchmark test for conducting a warping analysis."""
    geom = analysis_geometry(elements)
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()

    def warping_analysis():
        sec.calculate_warping_properties()

    benchmark(warping_analysis)


@pytest.mark.benchmark_suite
@pytest.mark.benchmark_analysis
@pytest.mark.parametrize("elements", [50, 500, 5000])
def test_frame_analysis(benchmark, analysis_geometry, elements):
    """Benchmark test for conducting a frame analysis."""
    geom = analysis_geometry(elements)
    sec = Section(geometry=geom)

    def frame_analysis():
        sec.calculate_frame_properties()

    benchmark(frame_analysis)


@pytest.mark.benchmark_suite
@pytest.mark.benchmark_analysis
@pytest.mark.parametrize("elements", [50, 500, 5000])
def test_stress_analysis(benchmark, analysis_geometry, elements):
    """Benchmark test for conducting a stress analysis."""
    geom = analysis_geometry(elements)
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()
    sec.calculate_warping_properties()

    def stress_analysis():
        sec.calculate_stress(n=1, vx=1, vy=1, mxx=1, m22=1, mzz=1)

    benchmark(stress_analysis)
