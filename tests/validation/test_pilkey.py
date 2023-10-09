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
def figure_6_7() -> Callable:
    """Defines the geometry in Figure 6.7 (page 253).

    Returns:
        Section object
    """

    def _generate_sec(c: float, ms: float, mat: Material | None = None) -> Section:
        points = [
            [-0.05 * c, -2 * c - 0.05 * c],
            [c, -2 * c - 0.05 * c],
            [c, -2 * c + 0.05 * c],
            [0.05 * c, -2 * c + 0.05 * c],
            [0.05 * c, -0.05 * c],
            [c, -0.05 * c],
            [c, 0.05 * c],
            [-0.05 * c, 0.05 * c],
        ]

        if mat:
            geom = Geometry(Polygon(points), material=mat)
        else:
            geom = Geometry(Polygon(points))

        geom.create_mesh(mesh_sizes=ms)
        sec = Section(geometry=geom)
        sec.calculate_geometric_properties()
        return sec

    return _generate_sec


@pytest.fixture
def figure_6_8() -> Callable:
    """Defines the geometry in Figure 6.8 (page 263).

    Returns:
        Section object
    """

    def _generate_sec(d: float, t: float) -> Section:
        geom = circular_hollow_section(d=d, t=t, n=128)
        ms = geom.calculate_area() / 500.0
        geom.points[0] = (d / 2, -1e-3)
        geom.points[128] = (d / 2 - t, 1e-3)
        geom = Geometry(Polygon(geom.points))
        geom = geom.rotate_section(angle=-90)
        geom.create_mesh(mesh_sizes=[ms])
        sec = Section(geometry=geom)
        sec.calculate_geometric_properties()

        return sec

    return _generate_sec


@pytest.fixture
def figure_6_9() -> Callable:
    """Defines the geometry in Figure 6.9 (page 264).

    Returns:
        Section object
    """

    def _generate_sec(beta: float) -> Section:
        geom = rectangular_section(d=1, b=2).align_center(align_to=(0, 0))
        geom.rotate_section(angle=beta)
        geom.create_mesh(mesh_sizes=geom.calculate_area() / 1000.0)
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


@pytest.mark.parametrize("d", [1.0, 50.0, 150.0])
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
    ms = area / 500.0
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
    ms = geom.calculate_area() / 500.0
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


def test_example_5_7(figure_6_8):
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
    ms = chs.calculate_area() / 500.0
    chs.create_mesh(mesh_sizes=[ms])
    sec_chs = Section(geometry=chs)

    # conduct analysis and check results
    sec_chs.calculate_geometric_properties()
    sec_chs.calculate_warping_properties()
    check.almost_equal(sec_chs.get_j(), 91.49, rel=tol)
    sigs = sec_chs.get_stress_at_points(pts=[(0, 4)], mzz=1)
    check.almost_equal(sigs[0][1], -0.0437, rel=tol)

    # create chs with slit
    sec_slit_chs = figure_6_8(d=8, t=0.25)

    # conduct analysis and check results
    sec_slit_chs.calculate_warping_properties()
    check.almost_equal(sec_slit_chs.get_j(), 0.13, rel=tol)
    sigs = sec_slit_chs.get_stress_at_points(pts=[(0, 4)], mzz=1)
    check.almost_equal(sigs[0][1], -1.92, rel=tol)


def test_example_5_13(al_cu: tuple[Material, Material]):
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


def test_example_6_1(figure_1_8):
    """Validation test for Example 6.1 (page 233).

    a = 10, t = 1

    tau_min = P / (12 * a)
    tau_max = 25 * P / (36 * a)
    """
    tol = 1e-2

    # define geometry
    sec = figure_1_8(a=10, t=1, ms=0.1)

    # analysis and results
    sec.calculate_warping_properties()
    sigs = sec.calculate_stress(vy=-1).get_stress()[0]
    vx = sigs["sig_zx_vy"]
    vy = sigs["sig_zy_vy"]

    # check stresses
    v_min = np.sqrt(max(vx) ** 2 + max(vy) ** 2)  # assume occurs at same spot
    check.almost_equal(v_min, 1 / (12 * 10), rel=5 * tol)
    check.almost_equal(min(vy), -25 * 1 / (36 * 10), rel=tol)  # vx insignificant


def test_example_6_2(figure_6_7):
    """Validation test for Example 6.2 (page 253)."""
    tol = 1e-3

    # define geometry
    mat = Material("a", 1, 0.3, 1, 1, "k")
    sec = figure_6_7(c=1.0, ms=1e-3, mat=mat)

    # analysis and results
    sec.calculate_warping_properties()
    x_c, y_c = sec.get_c()
    x_sc, y_sc = sec.get_sc()
    area = sec.get_area()
    a_sx, a_sy = sec.get_eas()
    check.almost_equal(x_c, 0.24937, rel=tol)
    check.almost_equal(y_c, -1.0, rel=tol)
    check.almost_equal(x_sc - x_c, -0.62056, rel=tol)
    check.almost_equal(y_sc - y_c, 0.0, abs=tol)
    check.almost_equal(area / a_sx, 3.09621, rel=tol)
    check.almost_equal(area / a_sy, 2.34102, rel=tol)
    check.almost_equal(sec.get_ej(), 0.00133, rel=tol)


@pytest.mark.parametrize(
    ("nu", "result"), [(0.0, -0.62055), (0.333, -0.62056), (0.5, -0.62057)]
)
def test_example_6_3(figure_6_7, nu: float, result: float):
    """Validation test for Example 6.3 (page 256)."""
    tol = 1e-3

    # define geometry and perform analysis
    mat = Material("a", 1, nu, 1, 1, "k")
    sec = figure_6_7(c=1.0, ms=1e-3, mat=mat)
    sec.calculate_warping_properties()

    # check results
    x_c, _ = sec.get_c()
    x_sc, _ = sec.get_sc()
    x_sct, _ = sec.get_sc()
    check.almost_equal(x_sc - x_c, result, rel=tol)
    check.almost_equal(x_sct - x_c, -0.62055, rel=tol)


def test_example_6_4(figure_6_8):
    """Validation test for Example 6.4 (page 263)."""
    tol = 1e-3

    # define geometry and perform analysis
    sec = figure_6_8(d=17.25, t=1.25)
    sec.calculate_warping_properties()

    # check results
    _, y_sc = sec.get_sc()
    area = sec.get_area()
    a_sx, a_sy = sec.get_as()
    check.almost_equal(y_sc, 15.90306, rel=tol)
    check.almost_equal(area / a_sx, 5.93977, rel=tol)
    check.almost_equal(area / a_sy, 1.98015, rel=tol)
    check.almost_equal(sec.get_j(), 32.23967, rel=2 * tol)


@pytest.mark.parametrize(("nu", "nu_idx"), [(0.0, 0), (0.3, 1), (0.5, 2)])
@pytest.mark.parametrize(
    ("beta", "beta_idx"), [(0.0, 0), (30.0, 1), (45.0, 2), (90.0, 3)]
)
def test_example_6_5(nu: float, nu_idx: int, beta: float, beta_idx: int):
    """Validation test for Example 6.5 (page 264).

    Note that the table entry for beta = 60.0 should read beta = 45.0 (note that a_x =
    a_y in the results), error in text.

    Note that the a_xy values for beta = 30.0 out by a factor of 10, error in text.
    """
    results = [
        # nu = 0.0
        [
            {"a_x": 1.2, "a_y": 1.2, "a_xy": 0.0},  # beta = 0.0
            {"a_x": 1.2, "a_y": 1.2, "a_xy": 0.0},  # beta = 30.0
            {"a_x": 1.2, "a_y": 1.2, "a_xy": 0.0},  # beta = 45.0
            {"a_x": 1.2, "a_y": 1.2, "a_xy": 0.0},  # beta = 90.0
        ],
        # nu = 0.3
        [
            {"a_x": 1.20056, "a_y": 1.27479, "a_xy": 0.0},  # beta = 0.0
            {"a_x": 1.21912, "a_y": 1.25624, "a_xy": -0.032142},  # beta = 30.0
            {"a_x": 1.23768, "a_y": 1.23768, "a_xy": -0.03711},  # beta = 45.0
            {"a_x": 1.27479, "a_y": 1.20056, "a_xy": 0.0},  # beta = 90.0
        ],
        # nu = 0.5
        [
            {"a_x": 1.20118, "a_y": 1.35605, "a_xy": 0.0},  # beta = 0.0
            {"a_x": 1.23990, "a_y": 1.31733, "a_xy": -0.067062},  # beta = 30.0
            {"a_x": 1.27861, "a_y": 1.27861, "a_xy": -0.07744},  # beta = 45.0
            {"a_x": 1.35605, "a_y": 1.20118, "a_xy": 0.0},  # beta = 90.0
        ],
    ]

    # set tol and get result
    tol = 1e-3
    result = results[nu_idx][beta_idx]

    # define materials
    mat = Material("a", 1, nu, 1, 1, "w")

    # define geometry and perform analysis
    geom = rectangular_section(d=1.0, b=2.0, material=mat).rotate_section(angle=beta)
    geom.create_mesh(mesh_sizes=[0.01])
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()
    sec.calculate_warping_properties()

    # check results
    area = sec.get_area()
    a_sx, a_sy = sec.get_eas()
    a_sxy = sec.section_props.a_sxy
    assert a_sxy

    check.almost_equal(area / a_sx, result["a_x"], rel=tol)
    check.almost_equal(area / a_sy, result["a_y"], rel=tol)

    # check if a_xy is zero, adjust tolerance
    if np.isclose(result["a_xy"], 0):
        check.almost_equal(area / a_sxy, result["a_xy"], abs=tol)
    else:
        check.almost_equal(area / a_sxy, result["a_xy"], rel=tol)


@pytest.mark.parametrize(("nu", "idx"), [(0.0, 0), (0.3, 1), (0.5, 2)])
def test_example_6_6(nu: float, idx: int):
    """Validation test for Example 6.6 (page 267)."""
    results = [
        # nu = 0.0
        {
            "a_x": 1.313417,
            "a_y": 1.169516,
            "a_xy": -0.063785,
            "x_sc": -0.069460,
            "y_sc": -0.072013,
            "phi": -4.860108,
            "a_11": 1.323148,
            "a_22": 1.159780,
        },
        # nu = 0.3
        {
            "a_x": 1.465389,
            "a_y": 1.170704,
            "a_xy": -0.076064,
            "x_sc": -0.071690,
            "y_sc": -0.090135,
            "phi": -4.860108,
            "a_11": 1.476062,
            "a_22": 1.159977,
        },
        # nu = 0.5
        {
            "a_x": 1.630496,
            "a_y": 1.171995,
            "a_xy": -0.089406,
            "x_sc": -0.072681,
            "y_sc": -0.098190,
            "phi": -4.860108,
            "a_11": 1.642190,
            "a_22": 1.160191,
        },
    ]

    # set tol and get result
    tol = 1e-3
    result = results[idx]

    # define materials
    mat = Material("a", 1, nu, 1, 1, "w")

    # define geometry and perform analysis
    geom = Geometry(Polygon([(-0.5, -1), (0.5, -1), (0.5, 2), (-0.5, 1)]), material=mat)
    geom.create_mesh(mesh_sizes=[0.01])
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()
    sec.calculate_warping_properties()

    # check results
    area = sec.get_area()
    c_x, c_y = sec.get_c()
    x_sc, y_sc = sec.get_sc()
    a_sx, a_sy = sec.get_eas()
    a_s11, a_s22 = sec.get_eas_p()
    a_sxy = sec.section_props.a_sxy
    assert a_sxy

    check.almost_equal(area / a_sx, result["a_x"], rel=tol)
    check.almost_equal(area / a_sy, result["a_y"], rel=tol)
    check.almost_equal(area / a_sxy, result["a_xy"], rel=tol)
    check.almost_equal(x_sc - c_x, result["x_sc"], rel=tol)
    check.almost_equal(y_sc - c_y, result["y_sc"], rel=tol)
    check.almost_equal(sec.get_phi(), result["phi"], rel=tol)
    check.almost_equal(area / a_s11, result["a_11"], rel=tol)
    check.almost_equal(area / a_s22, result["a_22"], rel=tol)


@pytest.mark.parametrize(("nu", "idx"), [(0.0, 0), (0.3, 1), (0.5, 2)])
@pytest.mark.parametrize("sym", [True, False])
def test_example_6_7(nu: float, idx: int, sym: bool):
    """Validation test for Example 6.7 (page 268)."""
    results_sym = [
        # nu = 0.0
        {
            "j": 6.539189,
            "a_x": 2.301059,
            "a_y": 2.301059,
            "a_xy": 0.014969,
            "x_sc": -2.464750,
            "y_sc": -2.464750,
            "phi": -135.0,  # converted to degrees
        },
        # nu = 0.3
        {
            "j": 6.539189,
            "a_x": 2.302935,
            "a_y": 2.302935,
            "a_xy": 0.016782,
            "x_sc": -2.464971,
            "y_sc": -2.464971,
            "phi": -135.0,  # converted to degrees
        },
        # nu = 0.5
        {
            "j": 6.539189,
            "a_x": 2.304973,
            "a_y": 2.304973,
            "a_xy": 0.018753,
            "x_sc": -2.465070,
            "y_sc": -2.465070,
            "phi": -135.0,  # converted to degrees
        },
    ]

    results_unsym = [
        # nu = 0.0
        {
            "j": 8.211705,
            "a_x": 3.058207,
            "a_y": 1.898375,
            "a_xy": 0.039510,
            "x_sc": -1.997641,
            "y_sc": -4.423925,
            "phi": 0.430248 * 180 / np.pi - 180.0,  # converted to degrees
        },
        # nu = 0.3
        {
            "j": 8.211705,
            "a_x": 3.061764,
            "a_y": 1.899138,
            "a_xy": 0.041123,
            "x_sc": -1.997605,
            "y_sc": -4.424391,
            "phi": 0.430248 * 180 / np.pi - 180.0,  # converted to degrees
        },
        # nu = 0.5
        {
            "j": 8.211705,
            "a_x": 3.065628,
            "a_y": 1.899967,
            "a_xy": 0.042876,
            "x_sc": -1.997589,
            "y_sc": -4.424592,
            "phi": 0.430248 * 180 / np.pi - 180.0,  # converted to degrees
        },
    ]

    # set tol and get result
    tol = 1e-3

    if sym:
        result = results_sym[idx]
    else:
        result = results_unsym[idx]

    # define materials
    mat = Material("a", 1, nu, 1, 1, "w")

    # define geometry and perform analysis
    if sym:
        up1 = (1, 10.5)
        up2 = (0, 10.5)
    else:
        up1 = (1, 15.5)
        up2 = (0, 15.5)

    geom = Geometry(
        Polygon([(0, 0), (10.5, 0), (10.5, 1), (1, 1), up1, up2]), material=mat
    )
    geom.create_mesh(mesh_sizes=[0.1])
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()
    sec.calculate_warping_properties()

    # check results
    area = sec.get_area()
    c_x, c_y = sec.get_c()
    x_sc, y_sc = sec.get_sc()
    a_sx, a_sy = sec.get_eas()
    a_sxy = sec.section_props.a_sxy
    assert a_sxy

    check.almost_equal(sec.get_ej(), result["j"], rel=tol)
    check.almost_equal(area / a_sx, result["a_x"], rel=tol)
    check.almost_equal(area / a_sy, result["a_y"], rel=tol)
    check.almost_equal(area / a_sxy, result["a_xy"], rel=0.02)  # within 2% error
    check.almost_equal(x_sc - c_x, result["x_sc"], rel=tol)
    check.almost_equal(y_sc - c_y, result["y_sc"], rel=tol)
    check.almost_equal(sec.get_phi(), result["phi"], rel=tol)


# TODO:
# EXAMPLE B.1
# EXAMPLE B.2
# EXAMPLE B.3
# EXAMPLE B.4
# EXAMPLE B.5
# EXAMPLE B.7
# EXAMPLE B.8
