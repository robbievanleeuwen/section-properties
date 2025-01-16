"""Validation of cross-section properties for a rectangular section."""

from __future__ import annotations

import pytest
import pytest_check as check

from sectionproperties.analysis import Section
from sectionproperties.pre.library import rectangular_section

# constants
tol = 1e-6
warp_tol = 1e-4
zero_tol = 5e-3


@pytest.fixture
def rect_section() -> Section:
    """Creates a 100x50 rectangular section with a 10 units square mesh.

    Returns:
        Section object
    """
    geom = rectangular_section(d=100, b=50)
    geom.create_mesh(mesh_sizes=10)
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()
    return sec


def test_rectangular_section_geometric(rect_section):
    """Test rectangular section geometric properties.

    Geometric properties calculated from first principles.
    """
    sec = rect_section

    check.almost_equal(sec.section_props.area, 100 * 50, rel=tol)
    check.almost_equal(sec.section_props.perimeter, 2 * 100 + 2 * 50, rel=tol)
    check.almost_equal(sec.section_props.mass, 1 * 100 * 50, rel=tol)
    check.almost_equal(sec.section_props.ea, 1 * 100 * 50, rel=tol)
    check.almost_equal(sec.section_props.qx, 100 * 50 * 50, rel=tol)
    check.almost_equal(sec.section_props.qy, 100 * 50 * 25, rel=tol)
    check.almost_equal(sec.section_props.ixx_g, 50 * 100**3 / 3, rel=tol)
    check.almost_equal(sec.section_props.iyy_g, 100 * 50**3 / 3, rel=tol)
    check.almost_equal(sec.section_props.ixy_g, 100 * 50 * 50 * 25, rel=tol)
    check.almost_equal(sec.section_props.cx, 50 / 2, rel=tol)
    check.almost_equal(sec.section_props.cy, 100 / 2, rel=tol)
    check.almost_equal(sec.section_props.ixx_c, 50 * 100**3 / 12, rel=tol)
    check.almost_equal(sec.section_props.iyy_c, 100 * 50**3 / 12, rel=tol)
    check.almost_equal(sec.section_props.ixy_c, 0, abs=zero_tol)
    check.almost_equal(sec.section_props.zxx_plus, 50 * 100**2 / 6, rel=tol)
    check.almost_equal(sec.section_props.zxx_minus, 50 * 100**2 / 6, rel=tol)
    check.almost_equal(sec.section_props.zyy_plus, 100 * 50**2 / 6, rel=tol)
    check.almost_equal(sec.section_props.zyy_minus, 100 * 50**2 / 6, rel=tol)
    check.almost_equal(
        sec.section_props.rx_c,
        (50 * 100**3 / 12 / 100 / 50) ** 0.5,
        rel=tol,
    )
    check.almost_equal(
        sec.section_props.ry_c,
        (100 * 50**3 / 12 / 100 / 50) ** 0.5,
        rel=tol,
    )
    check.almost_equal(sec.section_props.i11_c, (50 * 100**3 / 12), rel=tol)
    check.almost_equal(sec.section_props.i22_c, (100 * 50**3 / 12), rel=tol)
    check.almost_equal(sec.section_props.phi, 0, abs=zero_tol)
    check.almost_equal(sec.section_props.z11_plus, 50 * 100**2 / 6, rel=tol)
    check.almost_equal(sec.section_props.z11_minus, 50 * 100**2 / 6, rel=tol)
    check.almost_equal(sec.section_props.z22_plus, 100 * 50**2 / 6, rel=tol)
    check.almost_equal(sec.section_props.z22_minus, 100 * 50**2 / 6, rel=tol)
    check.almost_equal(
        sec.section_props.r11_c,
        (50 * 100**3 / 12 / 100 / 50) ** 0.5,
        rel=tol,
    )
    check.almost_equal(
        sec.section_props.r22_c,
        (100 * 50**3 / 12 / 100 / 50) ** 0.5,
        rel=tol,
    )


def test_rectangular_section_warping(rect_section):
    """Test rectangular section warping properties.

    Several non-trivial warping results were obtained from sectionproperties v3.0.2
    with a refined mesh [mesh_sizes=1.0], (previously inaccurate numerical results were
    from Strand7).
    """
    sec = rect_section
    sec.calculate_warping_properties()

    x_se, y_se = sec.get_sc()
    x11_se, y22_se = sec.get_sc_p()
    x_st, y_st = sec.get_sc_t()

    check.almost_equal(sec.section_props.j, 2.858521e06, rel=warp_tol)
    check.almost_equal(sec.section_props.gamma, 3.175417e08, rel=warp_tol)
    check.almost_equal(x_se, 50 / 2, rel=tol)
    check.almost_equal(y_se, 100 / 2, rel=tol)
    check.almost_equal(x11_se, 0, abs=zero_tol)
    check.almost_equal(y22_se, 0, abs=zero_tol)
    check.almost_equal(x_st, 50 / 2, rel=tol)
    check.almost_equal(y_st, 100 / 2, rel=tol)
    check.almost_equal(sec.section_props.a_sx, 4.166667e03, rel=warp_tol)
    check.almost_equal(sec.section_props.a_sy, 4.166667e03, rel=warp_tol)
    check.almost_equal(sec.section_props.a_s11, 4.166667e03, rel=warp_tol)
    check.almost_equal(sec.section_props.a_s22, 4.166667e03, rel=warp_tol)
    check.almost_equal(sec.section_props.beta_x_plus, 0, abs=zero_tol)
    check.almost_equal(sec.section_props.beta_x_minus, 0, abs=zero_tol)
    check.almost_equal(sec.section_props.beta_y_plus, 0, abs=zero_tol)
    check.almost_equal(sec.section_props.beta_y_minus, 0, abs=zero_tol)
    check.almost_equal(sec.section_props.beta_11_plus, 0, abs=zero_tol)
    check.almost_equal(sec.section_props.beta_11_minus, 0, abs=zero_tol)
    check.almost_equal(sec.section_props.beta_22_plus, 0, abs=zero_tol)
    check.almost_equal(sec.section_props.beta_22_minus, 0, abs=zero_tol)

    # check cgs method
    sec.calculate_warping_properties(solver_type="cgs")

    x_se, y_se = sec.get_sc()
    x11_se, y22_se = sec.get_sc_p()
    x_st, y_st = sec.get_sc_t()

    check.almost_equal(sec.section_props.j, 2.858521e06, rel=warp_tol)
    check.almost_equal(sec.section_props.gamma, 3.175417e08, rel=warp_tol)
    check.almost_equal(x_se, 50 / 2, rel=tol)
    check.almost_equal(y_se, 100 / 2, rel=tol)
    check.almost_equal(x11_se, 0, abs=zero_tol)
    check.almost_equal(y22_se, 0, abs=zero_tol)
    check.almost_equal(x_st, 50 / 2, rel=tol)
    check.almost_equal(y_st, 100 / 2, rel=tol)
    check.almost_equal(sec.section_props.a_sx, 4.166667e03, rel=warp_tol)
    check.almost_equal(sec.section_props.a_sy, 4.166667e03, rel=warp_tol)
    check.almost_equal(sec.section_props.a_s11, 4.166667e03, rel=warp_tol)
    check.almost_equal(sec.section_props.a_s22, 4.166667e03, rel=warp_tol)
    check.almost_equal(sec.section_props.beta_x_plus, 0, abs=zero_tol)
    check.almost_equal(sec.section_props.beta_x_minus, 0, abs=zero_tol)
    check.almost_equal(sec.section_props.beta_y_plus, 0, abs=zero_tol)
    check.almost_equal(sec.section_props.beta_y_minus, 0, abs=zero_tol)
    check.almost_equal(sec.section_props.beta_11_plus, 0, abs=zero_tol)
    check.almost_equal(sec.section_props.beta_11_minus, 0, abs=zero_tol)
    check.almost_equal(sec.section_props.beta_22_plus, 0, abs=zero_tol)
    check.almost_equal(sec.section_props.beta_22_minus, 0, abs=zero_tol)


def test_rectangular_section_plastic(rect_section):
    """Test rectangular section plastic properties.

    Plastic properties calculated from first principles.
    """
    sec = rect_section
    sec.calculate_plastic_properties()

    x_pc, y_pc = sec.get_pc()
    x11_pc, y22_pc = sec.get_pc_p()

    check.almost_equal(x_pc, 50 / 2, rel=tol)
    check.almost_equal(y_pc, 100 / 2, rel=tol)
    check.almost_equal(x11_pc, 50 / 2, rel=tol)
    check.almost_equal(y22_pc, 100 / 2, rel=tol)
    check.almost_equal(sec.section_props.sxx, 50 * 100**2 / 4, rel=tol)
    check.almost_equal(sec.section_props.syy, 100 * 50**2 / 4, rel=tol)
    check.almost_equal(sec.section_props.s11, 50 * 100**2 / 4, rel=tol)
    check.almost_equal(sec.section_props.s22, 100 * 50**2 / 4, rel=tol)
    check.almost_equal(sec.section_props.sf_xx_plus, 1.5, rel=tol)
    check.almost_equal(sec.section_props.sf_xx_minus, 1.5, rel=tol)
    check.almost_equal(sec.section_props.sf_yy_plus, 1.5, rel=tol)
    check.almost_equal(sec.section_props.sf_yy_minus, 1.5, rel=tol)
    check.almost_equal(sec.section_props.sf_11_plus, 1.5, rel=tol)
    check.almost_equal(sec.section_props.sf_11_minus, 1.5, rel=tol)
    check.almost_equal(sec.section_props.sf_22_plus, 1.5, rel=tol)
    check.almost_equal(sec.section_props.sf_22_minus, 1.5, rel=tol)
