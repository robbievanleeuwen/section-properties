"""Benchmark tests for sectionproperties mesh creation."""

import pytest


@pytest.mark.benchmark_suite
@pytest.mark.benchmark_mesh
@pytest.mark.parametrize("ms", [0.0, 50.0, 5.0])
def test_create_simple_mesh(benchmark, rect_geom, ms):
    """Benchmark test for creating a mesh for a rectangular geometry."""
    geom = rect_geom
    benchmark(geom.create_mesh, ms)


@pytest.mark.benchmark_suite
@pytest.mark.benchmark_mesh
@pytest.mark.parametrize("ms", [0.0, 1.0, 0.3])
def test_create_intermediate_mesh(benchmark, chs_geom, ms):
    """Benchmark test for creating a mesh for a CHS geometry."""
    geom = chs_geom
    benchmark(geom.create_mesh, ms)


@pytest.mark.benchmark_suite
@pytest.mark.benchmark_mesh
@pytest.mark.parametrize("ms", [0.0, 100.0, 20.0])
def test_create_complex_mesh(benchmark, concrete_column_with_hole, ms):
    """Benchmark test for creating a mesh for a concrete geometry."""
    geom = concrete_column_with_hole()
    benchmark(geom.create_mesh, ms)
