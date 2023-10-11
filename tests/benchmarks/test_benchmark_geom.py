"""Benchmark tests for sectionproperties geometry creation."""

import pytest

from sectionproperties.pre.library import circular_hollow_section, rectangular_section


@pytest.mark.benchmark_suite
@pytest.mark.benchmark_geom
def test_create_simple_geometry(benchmark):
    """Benchmark test for creating rectangular geometry."""
    benchmark(rectangular_section, d=100, b=50)


@pytest.mark.benchmark_suite
@pytest.mark.benchmark_geom
def test_create_intermediate_geometry(benchmark):
    """Benchmark test for creating CHS geometry."""
    benchmark(circular_hollow_section, d=100, t=3, n=128)


@pytest.mark.benchmark_suite
@pytest.mark.benchmark_geom
def test_create_complex_geometry(benchmark, concrete_column_with_hole):
    """Benchmark test for creating concrete geometry."""
    benchmark(concrete_column_with_hole)
