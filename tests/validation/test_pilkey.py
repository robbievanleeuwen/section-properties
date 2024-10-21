"""Validation tests from Pilkey.

BibTeX reference:
@book{Pilkey,
    author = {Pilkey, Walter D},
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

from typing import TYPE_CHECKING

import numpy as np
import pytest
import pytest_check as check
from shapely import LineString, Polygon, buffer

from sectionproperties.analysis import Section
from sectionproperties.pre import Geometry, Material
from sectionproperties.pre.library import (
    circular_hollow_section,
    circular_section_by_area,
    elliptical_hollow_section,
    elliptical_section,
    rectangular_section,
)

if TYPE_CHECKING:
    from collections.abc import Callable


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

    def _generate_sec(
        d: float,
        t: float,
        mat: Material | None = None,
        shift: tuple[float, float] | None = None,
    ) -> Section:
        geom = circular_hollow_section(d=d, t=t, n=128)
        ms = geom.calculate_area() / 500.0
        geom.points[0] = (d / 2, -1e-3)
        geom.points[128] = (d / 2 - t, 1e-3)

        if mat:
            geom = Geometry(Polygon(geom.points), material=mat)
        else:
            geom = Geometry(Polygon(geom.points))

        geom = geom.rotate_section(angle=-90)

        if shift:
            geom = geom.shift_section(x_offset=shift[0], y_offset=shift[1])

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

    def _generate_sec(beta: float, mat: Material) -> Section:
        geom = rectangular_section(d=1.0, b=2.0, material=mat).rotate_section(
            angle=beta,
        )
        geom.create_mesh(mesh_sizes=0.01)
        sec = Section(geometry=geom)
        sec.calculate_geometric_properties()
        sec.calculate_warping_properties()

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
    """Validation test for Example 5.13 - Composite Rectangular Strip (page 220).

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
    ("nu", "result"),
    [(0.0, -0.62055), (0.333, -0.62056), (0.5, -0.62057)],
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


@pytest.mark.parametrize(("nu", "nu_idx"), [(0.0, 0), (0.3, 1), (0.5, 2)])
@pytest.mark.parametrize(
    ("beta", "beta_idx"),
    [(0.0, 0), (30.0, 1), (45.0, 2), (90.0, 3)],
)
def test_example_6_5(figure_6_9, nu: float, nu_idx: int, beta: float, beta_idx: int):
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
    sec = figure_6_9(beta=beta, mat=mat)

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
    result = results_sym[idx] if sym else results_unsym[idx]

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
        Polygon([(0, 0), (10.5, 0), (10.5, 1), (1, 1), up1, up2]),
        material=mat,
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


def test_example_b_1():
    """Validation test for Example B.1 - Closed Elliptical Tube (page 434)."""
    # define input values
    a_mid = 8  # x-distance from centre to mid-line of EHS
    b_mid = 5  # y-distance from centre to mid-line of EHS
    t = 1  # thickness of EHS
    mat = Material("a", 2.1e8, 0.33333, 1, 1, "w")
    tol = 2e-2  # note differences in spline approach used by Pilkey (not an ellipse)

    # create geometry, mesh, section, and run analysis
    geom = elliptical_hollow_section(
        d_x=2 * a_mid + t,
        d_y=2 * b_mid + t,
        t=t,
        n=128,
        material=mat,
    )
    geom.create_mesh(mesh_sizes=[0.1])
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()
    sec.calculate_warping_properties()

    # conduct checks
    ixx, iyy, ixy = sec.get_eic(e_ref=mat)
    zxx, _, zyy, _ = sec.get_ez(e_ref=mat)
    rx, ry = sec.get_rc()
    a_sx, a_sy = sec.get_eas(e_ref=mat)
    a_sxy = sec.section_props.a_sxy
    assert a_sxy
    check.almost_equal(sec.get_area(), 41.38626, rel=tol)
    check.almost_equal(ixx, 580.42697, rel=tol)
    check.almost_equal(iyy, 1180.33120, rel=tol)
    check.almost_equal(ixy, 0, abs=1e-3)
    check.almost_equal(zxx, 105.54076, rel=tol)
    check.almost_equal(zyy, 138.86247, rel=tol)
    check.almost_equal(rx, 3.74495, rel=tol)
    check.almost_equal(ry, 5.34040, rel=tol)
    check.almost_equal(sec.get_phi(), -90.0, rel=tol)
    check.almost_equal(sec.get_area() / a_sx, 1.51457, rel=tol)
    check.almost_equal(sec.get_area() / a_sy, 3.05985, rel=tol)
    check.almost_equal(sec.get_ea() / a_sxy, 0.0, abs=1e-3)
    check.almost_equal(sec.get_ej(e_ref=mat), 1537.38165, rel=tol)
    # gamma is within 7%, possible spline error amplifies this?
    check.almost_equal(sec.get_egamma(e_ref=mat), 451.90976, rel=7e-2)


def test_example_b_2():
    """Validation test for Example B.2 - Symmetric Channel Section (page 437)."""
    # define input values
    b = 8
    h = 18
    t = 1
    mat = Material("a", 2.1e8, 0.33333, 1, 1, "w")
    tol = 1e-5
    warp_tol = 1e-3
    zero_tol = 1e-3

    # create geometry, mesh, section, and run analysis
    geom = Geometry(
        Polygon(
            [
                (-0.5 * t, -0.5 * h - 0.5 * t),
                (b, -0.5 * h - 0.5 * t),
                (b, -0.5 * h + 0.5 * t),
                (0.5 * t, -0.5 * h + 0.5 * t),
                (0.5 * t, 0.5 * h - 0.5 * t),
                (b, 0.5 * h - 0.5 * t),
                (b, 0.5 * h + 0.5 * t),
                (-0.5 * t, 0.5 * h + 0.5 * t),
            ],
        ),
        material=mat,
    )
    geom.create_mesh(mesh_sizes=[0.1])
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()
    sec.calculate_warping_properties()

    # conduct checks
    qx, qy = sec.get_eq(e_ref=mat)
    cx, cy = sec.get_c()
    ixx_g, iyy_g, ixy_g = sec.get_eig(e_ref=mat)
    ixx_c, iyy_c, ixy_c = sec.get_eic(e_ref=mat)
    zxx_p, zxx_m, zyy_p, zyy_m = sec.get_ez(e_ref=mat)
    rx, ry = sec.get_rc()
    x_sc, y_sc = sec.get_sc()
    x_sct, y_sct = sec.get_sc_t()
    a_sx, a_sy = sec.get_eas(e_ref=mat)
    a_sxy = sec.section_props.a_sxy
    assert a_sxy
    check.almost_equal(sec.get_area(), 34.0, rel=tol)
    check.almost_equal(qx, 0.0, abs=zero_tol)
    check.almost_equal(qy, 63.75, rel=tol)
    check.almost_equal(cx, 1.875, rel=tol)
    check.almost_equal(cy, 0, abs=zero_tol)
    check.almost_equal(x_sc, -2.86769, rel=warp_tol)
    check.almost_equal(y_sc, 0, abs=zero_tol)
    check.almost_equal(x_sct, -2.86759, rel=warp_tol)
    check.almost_equal(y_sct, 0, abs=zero_tol)
    # textbook erroneously printed 1 instead of 7, ixx_g should equal ixx_c
    check.almost_equal(ixx_g, 1787.83333, rel=tol)
    check.almost_equal(iyy_g, 342.83333, rel=tol)
    check.almost_equal(ixy_g, 0.0, abs=zero_tol)
    check.almost_equal(ixx_c, 1787.83333, rel=tol)
    check.almost_equal(iyy_c, 223.30208, rel=tol)
    check.almost_equal(ixy_c, 0.0, abs=zero_tol)
    check.almost_equal(min(zxx_p, zxx_m), 188.19298, rel=tol)
    check.almost_equal(min(zyy_p, zyy_m), 36.45748, rel=tol)
    check.almost_equal(rx, 7.25144, rel=tol)
    check.almost_equal(ry, 2.56275, rel=tol)
    check.almost_equal(sec.get_phi(), 0.0, abs=zero_tol)
    check.almost_equal(sec.get_area() / a_sx, 3.40789, rel=warp_tol)
    check.almost_equal(sec.get_area() / a_sy, 2.15337, rel=warp_tol)
    check.almost_equal(sec.get_ea() / a_sxy, 0.0, abs=zero_tol)
    check.almost_equal(sec.get_ej(e_ref=mat), 11.28862, rel=warp_tol)
    check.almost_equal(sec.get_egamma(e_ref=mat), 12763.15184, rel=warp_tol)


def test_example_b_3():
    """Validation test for Example B.3 - L Section Without Symmetry (page 441)."""
    # define input values
    b = 5.625
    h = 7.625
    t = 0.75
    mat = Material("a", 2.1e8, 0.33333, 1, 1, "w")
    tol = 1e-5
    warp_tol = 1e-3

    # create geometry, mesh, section, and run analysis
    geom = Geometry(
        Polygon(
            [
                (-0.5 * t, -0.5 * t),
                (b, -0.5 * t),
                (b, 0.5 * t),
                (0.5 * t, 0.5 * t),
                (0.5 * t, h),
                (-0.5 * t, h),
            ],
        ),
        material=mat,
    )
    geom.create_mesh(mesh_sizes=[0.05])
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()
    sec.calculate_warping_properties()

    # conduct checks
    qx, qy = sec.get_eq(e_ref=mat)
    cx, cy = sec.get_c()
    ixx_c, iyy_c, ixy_c = sec.get_eic(e_ref=mat)
    x_sc, y_sc = sec.get_sc()
    i11_c, i22_c = sec.get_eip(e_ref=mat)
    check.almost_equal(sec.get_area(), 9.93750, rel=tol)
    check.almost_equal(qx, 21.75000, rel=tol)
    check.almost_equal(qy, 11.81250, rel=tol)
    check.almost_equal(cx, 1.18868, rel=tol)
    check.almost_equal(cy, 2.18868, rel=tol)
    check.almost_equal(x_sc, 0.00313, rel=10 * warp_tol)  # lack of precision in digits
    check.almost_equal(y_sc, 0.05735, rel=warp_tol)
    check.almost_equal(ixx_c, 63.42455, rel=tol)
    check.almost_equal(iyy_c, 30.72142, rel=tol)
    check.almost_equal(ixy_c, -25.85377, rel=tol)
    check.almost_equal(sec.get_phi(), 28.84407 - 180.0, rel=tol)
    check.almost_equal(i11_c, 77.66369, rel=tol)
    check.almost_equal(i22_c, 16.48228, rel=tol)


def test_example_b_4(figure_6_8):
    """Validation test for Example B.4 - Open Circular Cross Section (page 444)."""
    # define input values
    r_mid = 8
    t = 1.25
    mat = Material("a", 2.1e8, 0.33333, 1, 1, "w")
    tol = 1e-3
    warp_tol = 1e-3
    zero_tol = 1e-2

    # define geometry and perform analysis
    sec = figure_6_8(d=2 * r_mid + t, t=t, mat=mat, shift=(0.0, 8.0))
    sec.calculate_warping_properties()

    # conduct checks
    qx, qy = sec.get_eq(e_ref=mat)
    cx, cy = sec.get_c()
    ixx_g, iyy_g, ixy_g = sec.get_eig(e_ref=mat)
    ixx_c, iyy_c, ixy_c = sec.get_eic(e_ref=mat)
    zxx_p, zxx_m, zyy_p, zyy_m = sec.get_ez(e_ref=mat)
    rx, ry = sec.get_rc()
    x_sc, y_sc = sec.get_sc()
    x_sct, y_sct = sec.get_sc_t()
    a_sx, a_sy = sec.get_eas(e_ref=mat)
    a_sxy = sec.section_props.a_sxy
    assert a_sxy
    check.almost_equal(sec.get_area(), 62.83182, rel=tol)
    check.almost_equal(qx, 502.65457, rel=tol)
    check.almost_equal(qy, 0.0, abs=zero_tol)
    check.almost_equal(cx, 0.0, abs=zero_tol)
    check.almost_equal(cy, 8.0, rel=tol)
    check.almost_equal(x_sc, 0.0, abs=zero_tol)
    check.almost_equal(y_sc, 23.90306, rel=warp_tol)
    check.almost_equal(x_sct, 0.0, abs=zero_tol)
    check.almost_equal(y_sct, 23.90282, rel=warp_tol)
    check.almost_equal(ixx_g, 6044.12555, rel=tol)
    check.almost_equal(iyy_g, 2022.88907, rel=tol)
    check.almost_equal(ixy_g, 0.0, abs=zero_tol)
    check.almost_equal(ixx_c, 2022.88890, rel=tol)
    check.almost_equal(iyy_c, 2022.88907, rel=tol)
    check.almost_equal(ixy_c, 0.0, abs=0.1)
    check.almost_equal(min(zxx_p, zxx_m), 234.53784, rel=tol)
    check.almost_equal(min(zyy_p, zyy_m), 234.53789, rel=tol)
    check.almost_equal(rx, 5.67409, rel=tol)
    check.almost_equal(ry, 5.67409, rel=tol)
    check.almost_equal(sec.get_area() / a_sx, 5.93977, rel=warp_tol)
    check.almost_equal(sec.get_area() / a_sy, 1.98015, rel=warp_tol)
    check.almost_equal(sec.get_ea() / a_sxy, 0.0, abs=zero_tol)
    check.almost_equal(sec.get_ej(e_ref=mat), 32.23967, rel=2 * warp_tol)
    check.almost_equal(sec.get_egamma(e_ref=mat), 331651.29223, rel=2 * warp_tol)


def test_example_b_5():
    """Validation test for Example B.5 - Welded Hat Section (page 445)."""
    # define input values
    mat = Material("a", 2.1e8, 0.33333, 1, 1, "w")
    tol = 5e-2  # note differences in modelling and mesh (pilkey uses splines)
    zero_tol = 1e-3

    # define geometry and perform analysis
    bot_plate = rectangular_section(d=0.75, b=31, material=mat).shift_section(
        y_offset=-0.75,
    )
    top_plate_l = rectangular_section(d=0.75, b=6, material=mat)
    top_plate_r = rectangular_section(d=0.75, b=6, material=mat).shift_section(
        x_offset=25,
    )

    sep = 1e-3  # separation between hat and bottom plate to simulate lack of weld
    # hat centreline
    line_pts = [
        (6, 0.375),
        (8, 0.375 + sep),
        (12, 6),
        (19, 6),
        (23, 0.375 + sep),
        (25, 0.375),
    ]
    line = LineString(coordinates=line_pts)
    # buffer hat centreline
    poly = buffer(geometry=line, distance=0.375, cap_style="flat", join_style="mitre")
    hat = Geometry(geom=poly, material=mat)

    # we must correct the ends of the hat to make sure they are orthogonal to x-y axes
    hat.points[4] = (25, 0.75)
    hat.points[5] = (25, 0)
    hat.points[10] = (6, 0)
    hat.points[11] = (6, 0.75)

    geom = bot_plate + hat + top_plate_l + top_plate_r
    geom.holes = [geom.holes[1]]  # remove erroneous holes
    geom.create_mesh(mesh_sizes=[0.1])
    sec = Section(geometry=geom)

    # perform analysis
    sec.calculate_geometric_properties()
    sec.calculate_warping_properties()

    # conduct checks
    qx, qy = sec.get_eq(e_ref=mat)
    cx, cy = sec.get_c()
    ixx_g, iyy_g, ixy_g = sec.get_eig(e_ref=mat)
    ixx_c, iyy_c, ixy_c = sec.get_eic(e_ref=mat)
    zxx_p, zxx_m, zyy_p, zyy_m = sec.get_ez(e_ref=mat)
    rx, ry = sec.get_rc()
    x_sc, y_sc = sec.get_sc()
    a_sx, a_sy = sec.get_eas(e_ref=mat)
    a_sxy = sec.section_props.a_sxy
    assert a_sxy
    check.almost_equal(sec.get_area(), 50.10830, rel=tol)
    check.almost_equal(qx, 57.90772, rel=tol)
    check.almost_equal(qy, 776.70074, rel=tol)
    check.almost_equal(cx, 15.50044, rel=tol)
    check.almost_equal(cy, 1.15565, rel=tol)
    check.almost_equal(x_sc, 15.50097, rel=tol)
    check.almost_equal(y_sc, 2.29918, rel=tol)
    check.almost_equal(ixx_g, 316.10331, rel=tol)
    check.almost_equal(iyy_g, 15874.79940, rel=tol)
    check.almost_equal(ixy_g, 897.57735, rel=tol)
    check.almost_equal(ixx_c, 249.18219, rel=tol)
    check.almost_equal(iyy_c, 3835.59608, rel=tol)
    check.almost_equal(ixy_c, 0.0, abs=1e-1)
    check.almost_equal(min(zxx_p, zxx_m), 47.74201, rel=tol)
    check.almost_equal(min(zyy_p, zyy_m), 247.45078, rel=tol)
    check.almost_equal(rx, 2.22999, rel=tol)
    check.almost_equal(ry, 8.74906, rel=tol)
    check.almost_equal(sec.get_phi(), -90, rel=tol)
    check.almost_equal(sec.get_area() / a_sx, 1.47731, rel=tol)
    check.almost_equal(sec.get_area() / a_sy, 9.98545, rel=tol)
    check.almost_equal(sec.get_ea() / a_sxy, 0.0, abs=zero_tol)
    check.almost_equal(sec.get_ej(e_ref=mat), 401.99583, rel=tol)
    check.almost_equal(sec.get_egamma(e_ref=mat), 1856.54586, rel=tol)


def test_example_b_7():
    """Validation test for Example B.7 - Circular Arc (page 451)."""
    # define input values
    mat = Material("a", 2.1e8, 0.33333, 1, 1, "w")
    r = 16  # radius
    t = 0.5  # thickness
    alpha = 2 * np.pi / 3  # arc angle
    n = 128  # number of points defining arc
    tol = 1e-4
    zero_tol = 1e-3

    # define points on the arc
    pts = []
    for idx in range(n):
        theta = -alpha / 2 + idx / (n - 1) * alpha
        x = r * np.sin(theta)
        y = r * np.cos(theta)
        pts.append((x, y))

    # create a centreline defining the arc
    line = LineString(coordinates=pts)
    # buffer the arc by the thickness
    poly = poly = buffer(
        geometry=line,
        distance=0.5 * t,
        cap_style="flat",
        join_style="mitre",
    )
    # create geometry, mesh, section and perform analysis
    geom = Geometry(geom=poly, material=mat)
    geom.create_mesh(mesh_sizes=[0.1])
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()
    sec.calculate_warping_properties()

    # conduct checks
    qx, qy = sec.get_eq(e_ref=mat)
    cx, cy = sec.get_c()
    ixx_g, iyy_g, ixy_g = sec.get_eig(e_ref=mat)
    ixx_c, iyy_c, ixy_c = sec.get_eic(e_ref=mat)
    zxx_p, zxx_m, zyy_p, zyy_m = sec.get_ez(e_ref=mat)
    rx, ry = sec.get_rc()
    x_sc, y_sc = sec.get_sc()
    a_sx, a_sy = sec.get_eas(e_ref=mat)
    a_sxy = sec.section_props.a_sxy
    assert a_sxy
    check.almost_equal(sec.get_area(), 16.75516, rel=tol)
    check.almost_equal(qx, 221.72054, rel=tol)
    check.almost_equal(qy, 0.0, abs=zero_tol)
    check.almost_equal(cx, 0.0, abs=zero_tol)
    check.almost_equal(cy, 13.23297, rel=tol)
    check.almost_equal(x_sc, 0.0, abs=zero_tol)
    check.almost_equal(y_sc, 17.83662, rel=tol)
    check.almost_equal(ixx_g, 3032.21070, rel=tol)
    check.almost_equal(iyy_g, 1258.15764, rel=tol)
    check.almost_equal(ixy_g, 0.0, abs=zero_tol)
    check.almost_equal(ixx_c, 98.18931, rel=tol)
    check.almost_equal(iyy_c, 1258.15764, rel=tol)
    check.almost_equal(ixy_c, 0.0, abs=zero_tol)
    check.almost_equal(min(zxx_p, zxx_m), 18.32584, rel=5 * tol)
    check.almost_equal(min(zyy_p, zyy_m), 89.40279, rel=tol)
    check.almost_equal(rx, 2.42079, rel=tol)
    check.almost_equal(ry, 8.66549, rel=tol)
    check.almost_equal(sec.get_phi(), -90, rel=tol)
    check.almost_equal(sec.get_area() / a_sx, 1.50823, rel=tol)
    check.almost_equal(sec.get_area() / a_sy, 4.60034, rel=tol)
    check.almost_equal(sec.get_ea() / a_sxy, 0.0, abs=zero_tol)
    check.almost_equal(sec.get_ej(e_ref=mat), 1.38355, rel=tol)
    check.almost_equal(sec.get_egamma(e_ref=mat), 1046.49221, rel=tol)
