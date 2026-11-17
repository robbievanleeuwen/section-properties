"""Tests for the direct warping solver that fixes a degree of freedom.

``solver.WarpingSolver`` must give the same solution as the Lagrangian multiplier
method it replaces, with and without PARDISO.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING

import numpy as np
import pytest

import sectionproperties.analysis.solver as solver
from sectionproperties.analysis.section import Section
from sectionproperties.pre import Material
from sectionproperties.pre.library import (
    circular_hollow_section,
    i_section,
    rectangular_section,
)

if TYPE_CHECKING:
    from collections.abc import Callable

STEEL = Material(
    name="Steel",
    elastic_modulus=200e3,
    poissons_ratio=0.3,
    yield_strength=500,
    density=7.85e-6,
    color="grey",
)
TIMBER = Material(
    name="Timber",
    elastic_modulus=8e3,
    poissons_ratio=0.35,
    yield_strength=20,
    density=6.5e-7,
    color="burlywood",
)


def i_beam() -> Section:
    """Doubly symmetric I-section."""
    geom = i_section(d=200, b=100, t_f=10, t_w=6, r=8, n_r=8)
    geom.create_mesh(mesh_sizes=10)

    return Section(geometry=geom)


def composite() -> Section:
    """Timber beam with a steel plate on top, not symmetric about x."""
    timber = rectangular_section(d=200, b=80, material=TIMBER)
    plate = rectangular_section(d=10, b=120, material=STEEL)
    plate = plate.align_center(align_to=timber).align_to(other=timber, on="top")
    geom = timber + plate
    geom.create_mesh(mesh_sizes=[20, 5])

    return Section(geometry=geom)


def chs() -> Section:
    """Circular hollow section, a closed and multiply connected mesh."""
    geom = circular_hollow_section(d=100, t=5, n=32)
    geom.create_mesh(mesh_sizes=5)

    return Section(geometry=geom)


SECTIONS = [i_beam, composite, chs]


@pytest.fixture(params=[False, True], ids=["superlu", "pardiso"])
def pardiso(request: pytest.FixtureRequest, monkeypatch: pytest.MonkeyPatch) -> bool:
    """Runs a test with and without PARDISO, skipping PARDISO when not installed."""
    use_pardiso = bool(request.param)

    if use_pardiso and not solver.HAS_PARDISO:
        pytest.skip("pypardiso is not installed")

    monkeypatch.setattr(solver, "HAS_PARDISO", use_pardiso)

    return use_pardiso


@pytest.mark.usefixtures("pardiso")
@pytest.mark.parametrize("build", SECTIONS)
def test_matches_lagrangian(build: Callable[[], Section]) -> None:
    """The warping function equals the Lagrangian multiplier solution."""
    section = build()
    row, col, data, f_torsion, c_torsion = section.assemble_torsion_entries()
    k_lg, f_lg = section.assemble_torsion()

    expected = solver.solve_direct_lagrange(k_lg=k_lg, f=f_lg)
    actual = solver.WarpingSolver(row=row, col=col, data=data, c=c_torsion).solve(
        f=f_torsion
    )

    np.testing.assert_allclose(
        actual, expected, rtol=0, atol=1e-9 * np.abs(expected).max()
    )


@pytest.mark.usefixtures("pardiso")
@pytest.mark.parametrize("build", SECTIONS)
def test_satisfies_constraint(build: Callable[[], Section]) -> None:
    """Every solution satisfies the constraint that removes the constant."""
    section = build()
    row, col, data, f_torsion, c_torsion = section.assemble_torsion_entries()
    warping_solver = solver.WarpingSolver(row=row, col=col, data=data, c=c_torsion)

    for f in (f_torsion, -2.5 * f_torsion):
        u = warping_solver.solve(f=f)
        assert abs(c_torsion @ u) <= 1e-12 * np.abs(c_torsion).sum() * np.abs(u).max()


@pytest.mark.usefixtures("pardiso")
def test_several_right_hand_sides() -> None:
    """Solving several right hand sides together equals solving them one by one."""
    section = composite()
    row, col, data, f_torsion, c_torsion = section.assemble_torsion_entries()
    warping_solver = solver.WarpingSolver(row=row, col=col, data=data, c=c_torsion)

    rhs = np.column_stack((f_torsion, -2.5 * f_torsion, 0.5 * f_torsion))
    together = warping_solver.solve(f=rhs)

    for idx in range(rhs.shape[1]):
        np.testing.assert_allclose(
            together[:, idx],
            warping_solver.solve(f=rhs[:, idx]),
            rtol=0,
            atol=1e-12 * np.abs(together).max(),
        )


@pytest.mark.parametrize("build", SECTIONS)
def test_warping_properties_match_cgs(build: Callable[[], Section]) -> None:
    """A full warping analysis agrees with the iterative Lagrangian solver."""
    direct = build()
    direct.calculate_geometric_properties()
    direct.calculate_warping_properties(solver_type="direct")

    iterative = build()
    iterative.calculate_geometric_properties()
    iterative.calculate_warping_properties(solver_type="cgs")

    for name in ("j", "gamma", "x_se", "y_se", "a_sx", "a_sy"):
        actual = getattr(direct.section_props, name)
        expected = getattr(iterative.section_props, name)
        assert actual == pytest.approx(expected, rel=1e-5, abs=1e-6), name


def test_frame_analysis_matches_warping_analysis() -> None:
    """The frame analysis uses the same solver and gives the same torsion constant."""
    section = i_beam()
    section.calculate_geometric_properties()
    section.calculate_warping_properties()

    frame = i_beam()
    frame.calculate_frame_properties()

    assert frame.get_j() == pytest.approx(section.get_j(), rel=1e-12)


def test_disjoint_mesh_fixes_one_dof_per_region() -> None:
    """A mesh of two separate regions still solves, with one dof fixed in each."""
    left = rectangular_section(d=100, b=20)
    right = rectangular_section(d=100, b=20).shift_section(x_offset=50)
    geom = left + right
    geom.create_mesh(mesh_sizes=20)
    section = Section(geometry=geom)

    row, col, data, f_torsion, c_torsion = section.assemble_torsion_entries()
    warping_solver = solver.WarpingSolver(row=row, col=col, data=data, c=c_torsion)
    omega = warping_solver.solve(f=f_torsion)

    assert len(warping_solver.pinned) == 2
    assert np.all(np.isfinite(omega))

    section.calculate_geometric_properties()

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        section.calculate_warping_properties()

    assert any("disjoint" in str(warning.message) for warning in caught)
    assert np.isfinite(section.get_j())
