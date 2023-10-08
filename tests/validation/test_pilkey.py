"""Validation tests from Pilkey.

Bibtex reference:
@book{PilkeyWalterD2002AaDo,
    author = {Pilkey, Walter D},
    address = {Newark},
    booktitle = {Analysis and Design of Elastic Beams},
    edition = {First},
    isbn = {0471381527},
    language = {eng},
    publisher = {Wiley},
    title = {Analysis and Design of Elastic Beams: Computational Methods},
    year = {2002},
}
"""

from __future__ import annotations

from typing import Callable

import numpy as np
import pytest
import pytest_check as check
from shapely import Polygon

from sectionproperties.analysis import Section
from sectionproperties.pre import Geometry, Material
from sectionproperties.pre.library import (
    circular_hollow_section,
    circular_section_by_area,
    elliptical_section,
    rectangular_section,
)


@pytest.fixture
def figure_1_8() -> Callable:
    """Defines the geometry in Figure 1.8 (page 18).

    Returns:
        Section object
    """

    def _generate_sec(a: float, t: float, ms: float) -> Section:
        points = [
            [-0.5 * t, -2 * a],
            [0.5 * t, -2 * a],
            [0.5 * t, -0.5 * t],
            [a, -0.5 * t],
            [a, 0.5 * t],
            [-0.5 * t, 0.5 * t],
        ]
        geom = Geometry(Polygon(points))
        geom.create_mesh(mesh_sizes=ms)
        sec = Section(geometry=geom)
        sec.calculate_geometric_properties()
        return sec

    return _generate_sec


@pytest.fixture
def al_cu() -> tuple[Material, Material]:
    """Defines the Aluminium and Copper materials.

    Aluminium: E = 10.4e6; nu = 0.3
    Copper: E = 18.5e6; nu = 0.3

    Returns:
        Aluminium and Copper material
    """
    al = Material("Al", 10.4e6, 0.3, 1, 1, "r")
    cu = Material("Cu", 18.5e6, 0.3, 1, 1, "b")

    return al, cu


def test_example_1_1(figure_1_8):
    """Validation test for Example 1.1 (page 18).

    Note results calculated based on on thin-wall assumptions.
    """
    tol = 1e-3
    stress_tol = 5e-2  # note stress inaccuracy due to thin wall assumption in Pilkey

    a = 1
    t = 0.1
    ms = 1e-3
    sec = figure_1_8(a=a, t=t, ms=ms)

    # check section properties
    check.almost_equal(sec.section_props.area, 0.3, rel=tol)
    check.almost_equal(sec.section_props.qx, -0.1999, rel=tol)
    check.almost_equal(sec.section_props.qy, 0.0499, rel=tol)
    check.almost_equal(sec.section_props.cx, 0.1663, rel=tol)
    check.almost_equal(sec.section_props.cy, -0.6663, rel=tol)
    check.almost_equal(sec.section_props.ixx_c, 0.1336, rel=tol)
    check.almost_equal(sec.section_props.iyy_c, 0.0252, rel=tol)
    check.almost_equal(sec.section_props.ixy_c, 0.0332, rel=tol)
    check.almost_equal(sec.section_props.phi, -15.7589, rel=tol)
    check.almost_equal(sec.section_props.i11_c, 0.1430, rel=tol)
    check.almost_equal(sec.section_props.i22_c, 0.0158, rel=2 * tol)

    # check stresses
    pts = [(a, 0), (0, 0), (0, -2 * a)]
    sigs = sec.get_stress_at_points(pts=pts, mxx=-1)
    check.almost_equal(sigs[0][0], 5, rel=stress_tol)
    check.almost_equal(sigs[1][0], -10, rel=stress_tol)
    check.almost_equal(sigs[2][0], 12.5, rel=stress_tol)


@pytest.mark.parametrize("d", [1.0, 50.0, 500.0])
def test_example_5_1(d: float):
    """Validation test for Example 5.1 (page 174).

    Note higher tolerance due to circle discretisation.

    J = Ixx + Iyy
    J = pi * d**4 / 32
    tau_zx = -Mzz * y / J
    tau_zy = Mzz * x / J
    CHS -> J = pi / 32 * (do**4 - di**4)
    """
    tol = 1e-3

    # create solid geometry
    area = np.pi * d**2 / 4
    ms = area / 1000.0
    geom_solid = circular_section_by_area(area=area, n=64)
    geom_solid.create_mesh(mesh_sizes=[ms])
    sec_solid = Section(geometry=geom_solid)

    # conduct solid analysis
    sec_solid.calculate_geometric_properties()
    sec_solid.calculate_warping_properties()

    # check J results
    j_numerical = sec_solid.get_j()
    check.almost_equal(j_numerical, sum(sec_solid.get_ic()[:2]), rel=tol)
    check.almost_equal(j_numerical, np.pi * d**4 / 32, rel=tol)

    # conduct stress analysis
    pts = [(d / 2, 0), (0, d / 2)]
    mzz = d * d
    sigs = sec_solid.get_stress_at_points(pts=pts, mzz=mzz)

    # check stress results
    check.almost_equal(sigs[0][1], 0, abs=tol)
    check.almost_equal(sigs[0][2], mzz * d / 2 / j_numerical, rel=tol)
    check.almost_equal(sigs[1][1], -mzz * d / 2 / j_numerical, rel=tol)
    check.almost_equal(sigs[1][2], 0, abs=tol)

    # create CHS geometry
    geom_chs = circular_hollow_section(d=d, t=0.05 * d, n=64)
    geom_chs.create_mesh(mesh_sizes=[ms / 2])
    sec_chs = Section(geometry=geom_chs)

    # conduct solid analysis
    sec_chs.calculate_geometric_properties()
    sec_chs.calculate_frame_properties()

    # check J results
    d_in = d - 0.1 * d
    check.almost_equal(sec_chs.get_j(), np.pi / 32 * (d**4 - d_in**4), rel=1e-2)


@pytest.mark.parametrize(("a", "b"), [(50, 75), (75, 50), (10, 50)])
def test_example_5_3(a: float, b: float):
    """Validation test for Example 5.3 (page 176).

    Note higher tolerance due to circle discretisation.

    J = pi * a**3 * b**3 / (a**2 + b**2)
    tau_zx = -2 * Mzz * y / (a * b**3 * pi)
    tau_zy = 2 * Mzz * x / (a**3 * b * pi)
    """
    tol = 1e-2

    # create geometry
    geom = elliptical_section(d_x=2 * a, d_y=2 * b, n=64)
    ms = geom.calculate_area() / 1000.0
    geom.create_mesh(mesh_sizes=[ms])
    sec = Section(geometry=geom)

    # conduct analysis
    sec.calculate_geometric_properties()
    sec.calculate_warping_properties()

    # check J results
    j_theory = np.pi * a**3 * b**3 / (a**2 + b**2)
    check.almost_equal(sec.get_j(), j_theory, rel=tol)

    # conduct stress analysis
    pts = [(a, 0), (0, b * 0.75)]
    mzz = a * b
    sigs = sec.get_stress_at_points(pts=pts, mzz=mzz)

    # check stress results
    factor_x = -2 / (a * b**3 * np.pi)
    factor_y = 2 / (a**3 * b * np.pi)
    check.almost_equal(sigs[0][1], 0, abs=tol)
    check.almost_equal(sigs[0][2], factor_y * mzz * a, rel=2 * tol)
    check.almost_equal(sigs[1][1], factor_x * mzz * b * 0.75, rel=2 * tol)
    check.almost_equal(sigs[1][2], 0, abs=tol)


def test_example_5_6():
    """Validation test for Example 5.6 (page 184).

    J = h**4 / (15 * sqrt(3))
    """
    tol = 1e-6

    a = 1  # base length
    h = a * np.sqrt(3) / 2

    # create geometry
    geom = Geometry(Polygon([(-a / 2, 0), (a / 2, 0), (0, h)]))
    ms = geom.calculate_area() / 1000.0
    geom.create_mesh(mesh_sizes=[ms])
    sec = Section(geometry=geom)

    # conduct analysis and check j
    sec.calculate_frame_properties()
    check.almost_equal(sec.get_j(), h**4 / (15 * np.sqrt(3)), rel=tol)


def test_example_5_7():
    """Validation test for Example 5.7 (page 189).

    8 x 0.25 CHS + same with slit. Note higher tolerance due to circle discretisation.

    J_chs = 91.49
    J_slit = 0.13
    tau_chs = 0.0437
    tau_slit = 1.92
    """
    tol = 5e-2

    # create chs
    chs = circular_hollow_section(d=8, t=0.25, n=64)
    ms = chs.calculate_area() / 1000.0
    chs.create_mesh(mesh_sizes=[ms])
    sec_chs = Section(geometry=chs)

    # conduct analysis and check results
    sec_chs.calculate_geometric_properties()
    sec_chs.calculate_warping_properties()
    check.almost_equal(sec_chs.get_j(), 91.49, rel=tol)
    sigs = sec_chs.get_stress_at_points(pts=[(0, 4)], mzz=1)
    check.almost_equal(sigs[0][1], -0.0437, rel=tol)

    # create chs with slit
    slit_chs = circular_hollow_section(d=8, t=0.25, n=64)
    slit_chs.points[0] = (4, -1e-3)
    slit_chs.points[64] = (3.75, 1e-3)
    slit_chs = Geometry(Polygon(slit_chs.points))
    slit_chs.create_mesh(mesh_sizes=[ms])
    sec_slit_chs = Section(geometry=slit_chs)

    # conduct analysis and check results
    sec_slit_chs.calculate_geometric_properties()
    sec_slit_chs.calculate_warping_properties()
    check.almost_equal(sec_slit_chs.get_j(), 0.13, rel=tol)
    sigs = sec_slit_chs.get_stress_at_points(pts=[(0, 4)], mzz=1)
    check.almost_equal(sigs[0][1], -1.92, rel=tol)


def test_example_5_13(al_cu):
    """Validation test for Example 5.13 (page 220).

    2 off 15 x 2 bars (1 Aluminium, 1 Copper)

    EA_al = 83.36
    J_al = 106.22
    """
    al, cu = al_cu
    tol = 1e-3

    # create geometry
    geom_al = rectangular_section(d=2, b=15, material=al)
    geom_cu = rectangular_section(d=2, b=15, material=cu).align_to(geom_al, "right")
    geom = geom_al + geom_cu
    geom.create_mesh(mesh_sizes=[0.1])
    sec = Section(geometry=geom)

    # analysis and results
    sec.calculate_frame_properties()
    check.almost_equal(sec.get_ea(e_ref=al), 83.36, rel=tol)
    check.almost_equal(sec.get_ej(e_ref=al), 106.22, rel=tol)


def test_example_5_14(figure_1_8):
    """Validation test for Example 5.14 (page 221).

    a = 1, t = 0.1

    J = 0.00099
    tau = 98e3 from mzz=1e3
    """
    tol = 5e-3
    stress_tol = 5e-2

    # define geometry
    sec = figure_1_8(a=1, t=0.1, ms=1e-3)

    # analysis and results
    sec.calculate_warping_properties()
    # check stress at top mid-point along edge
    sigs = sec.get_stress_at_points(pts=[(0.5, -0.05)], mzz=1e3)
    check.almost_equal(sec.get_j(), 0.00099, rel=tol)
    check.almost_equal(sigs[0][1], 98e3, rel=stress_tol)


# TODO:
# EXAMPLE 6.1
# EXAMPLE 6.2
# EXAMPLE 6.3
# EXAMPLE 6.4
# EXAMPLE 6.5
# EXAMPLE 6.6
# EXAMPLE 6.7
# EXAMPLE B.1
# EXAMPLE B.2
# EXAMPLE B.3
# EXAMPLE B.4
# EXAMPLE B.5
# EXAMPLE B.7
# EXAMPLE B.8
