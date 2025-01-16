"""Validation of cross-section properties for an angle section."""

from __future__ import annotations

import pytest
import pytest_check as check

from sectionproperties.analysis import Section
from sectionproperties.pre.library import angle_section

# constants
tol = 1e-6
warp_tol = 1e-4


@pytest.fixture
def ang_section() -> Section:
    """Creates a 150x90x12 unequal angle section with a 5 units square mesh.

    Returns:
        Section object
    """
    geom = angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
    geom.create_mesh(mesh_sizes=5)
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()
    return sec


def test_angle_section_geometric(ang_section):
    """Test angle section geometric properties.

    Geometric properties calculated from sectionproperties v3.0.2 with a refined mesh
    [mesh_sizes=0.5].
    """
    sec = ang_section

    check.almost_equal(sec.section_props.area, 2.747059e03, rel=tol)
    check.almost_equal(sec.section_props.perimeter, 4.713501e02, rel=tol)
    check.almost_equal(sec.section_props.mass, 2.747059e03, rel=tol)
    check.almost_equal(sec.section_props.ea, 2.747059e03, rel=tol)
    check.almost_equal(sec.section_props.qx, 1.400486e05, rel=tol)
    check.almost_equal(sec.section_props.qy, 5.830033e04, rel=tol)
    check.almost_equal(sec.section_props.ixx_g, 1.342632e07, rel=tol)
    check.almost_equal(sec.section_props.iyy_g, 2.955753e06, rel=tol)
    check.almost_equal(sec.section_props.ixy_g, 1.086603e06, rel=tol)
    check.almost_equal(sec.section_props.cx, 2.122282e01, rel=tol)
    check.almost_equal(sec.section_props.cy, 5.098127e01, rel=tol)
    check.almost_equal(sec.section_props.ixx_c, 6.286470e06, rel=tol)
    check.almost_equal(sec.section_props.iyy_c, 1.718455e06, rel=tol)
    check.almost_equal(sec.section_props.ixy_c, -1.885622e06, rel=tol)
    check.almost_equal(sec.section_props.zxx_plus, 6.348769e04, rel=tol)
    check.almost_equal(sec.section_props.zxx_minus, 1.233094e05, rel=tol)
    check.almost_equal(sec.section_props.zyy_plus, 2.498584e04, rel=tol)
    check.almost_equal(sec.section_props.zyy_minus, 8.097207e04, rel=tol)
    check.almost_equal(sec.section_props.rx_c, 4.783761e01, rel=tol)
    check.almost_equal(sec.section_props.ry_c, 2.501124e01, rel=tol)
    check.almost_equal(sec.section_props.i11_c, 6.964263e06, rel=tol)
    check.almost_equal(sec.section_props.i22_c, 1.040662e06, rel=tol)
    check.almost_equal(sec.section_props.phi, -1.602289e02, rel=tol)
    check.almost_equal(sec.section_props.z11_plus, 9.775662e04, rel=tol)
    check.almost_equal(sec.section_props.z11_minus, 6.939239e04, rel=tol)
    check.almost_equal(sec.section_props.z22_plus, 2.796211e04, rel=tol)
    check.almost_equal(sec.section_props.z22_minus, 2.076613e04, rel=tol)
    check.almost_equal(sec.section_props.r11_c, 5.035048e01, rel=tol)
    check.almost_equal(sec.section_props.r22_c, 1.946350e01, rel=tol)


def test_angle_section_warping(ang_section):
    """Test angle section warping properties.

    Warping properties calculated from sectionproperties v3.0.2 with a refined mesh
    [mesh_sizes=0.5].
    """
    sec = ang_section
    sec.calculate_warping_properties()

    x_se, y_se = sec.get_sc()
    x11_se, y22_se = sec.get_sc_p()
    x_st, y_st = sec.get_sc_t()

    check.almost_equal(sec.section_props.j, 1.354614e05, rel=1.5 * warp_tol)
    check.almost_equal(sec.section_props.gamma, 1.622210e08, rel=warp_tol)
    check.almost_equal(x_se, 6.124892e00, rel=warp_tol)
    check.almost_equal(y_se, 8.126346e00, rel=warp_tol)
    check.almost_equal(x11_se, 2.870420e01, rel=warp_tol)
    check.almost_equal(y22_se, 3.522160e01, rel=warp_tol)
    check.almost_equal(x_st, 6.124892e00, rel=warp_tol)
    check.almost_equal(y_st, 8.126346e00, rel=warp_tol)
    check.almost_equal(sec.section_props.a_sx, 8.586193e02, rel=warp_tol)
    check.almost_equal(sec.section_props.a_sy, 1.539971e03, rel=warp_tol)
    check.almost_equal(sec.section_props.a_s11, 8.855853e02, rel=warp_tol)
    check.almost_equal(sec.section_props.a_s22, 1.460224e03, rel=warp_tol)
    check.almost_equal(sec.section_props.beta_x_plus, -1.061144e02, rel=warp_tol)
    check.almost_equal(sec.section_props.beta_x_minus, 1.061144e02, rel=warp_tol)
    check.almost_equal(sec.section_props.beta_y_plus, -5.774378e01, rel=warp_tol)
    check.almost_equal(sec.section_props.beta_y_minus, 5.774378e01, rel=warp_tol)
    check.almost_equal(sec.section_props.beta_11_plus, 8.547674e01, rel=warp_tol)
    check.almost_equal(sec.section_props.beta_11_minus, -8.547674e01, rel=warp_tol)
    check.almost_equal(sec.section_props.beta_22_plus, 1.419115e02, rel=warp_tol)
    check.almost_equal(sec.section_props.beta_22_minus, -1.419115e02, rel=warp_tol)


def test_angle_section_plastic(ang_section):
    """Test angle section plastic properties.

    Warping properties calculated from sectionproperties v3.0.2 with a refined mesh
    [mesh_sizes=0.5].
    """
    sec = ang_section
    sec.calculate_plastic_properties()

    x_pc, y_pc = sec.get_pc()
    x11_pc, y22_pc = sec.get_pc_p()

    check.almost_equal(x_pc, 9.159481e00, rel=tol)
    check.almost_equal(y_pc, 3.507843e01, rel=tol)
    check.almost_equal(x11_pc, 2.311371e01, rel=tol)
    check.almost_equal(y22_pc, 4.123001e01, rel=tol)
    check.almost_equal(sec.section_props.sxx, 1.135392e05, rel=tol)
    check.almost_equal(sec.section_props.syy, 4.572265e04, rel=tol)
    check.almost_equal(sec.section_props.s11, 1.210275e05, rel=tol)
    check.almost_equal(sec.section_props.s22, 4.376054e04, rel=tol)
    check.almost_equal(sec.section_props.sf_xx_plus, 1.788366e00, rel=tol)
    check.almost_equal(sec.section_props.sf_xx_minus, 9.207672e-01, rel=tol)
    check.almost_equal(sec.section_props.sf_yy_plus, 1.829943e00, rel=tol)
    check.almost_equal(sec.section_props.sf_yy_minus, 5.646718e-01, rel=tol)
    check.almost_equal(sec.section_props.sf_11_plus, 1.238049e00, rel=tol)
    check.almost_equal(sec.section_props.sf_11_minus, 1.744103e00, rel=tol)
    check.almost_equal(sec.section_props.sf_22_plus, 1.564994e00, rel=tol)
    check.almost_equal(sec.section_props.sf_22_minus, 2.107303e00, rel=tol)
