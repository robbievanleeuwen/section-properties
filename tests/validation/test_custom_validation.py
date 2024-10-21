"""Validation of cross-section properties for a custom section."""

from __future__ import annotations

import pytest
import pytest_check as check
from shapely import Polygon

from sectionproperties.analysis import Section
from sectionproperties.pre import Geometry

# constants
tol = 1e-6
plastic_tol = 1e-5
warp_tol = 1e-3


@pytest.fixture
def custom_section() -> Section:
    """Creates a custom section with a 5 units square mesh.

    Returns:
        Section object
    """
    points = [
        [-10, 0],
        [110, 0],
        [100, 10],
        [55, 10],
        [55, 90],
        [100, 90],
        [110, 100],
        [110, 110],
        [-10, 110],
        [-10, 100],
        [0, 90],
        [45, 90],
        [45, 10],
        [-10, 10],
    ]
    geom = Geometry(Polygon(points))
    geom.create_mesh(mesh_sizes=5)
    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()
    return sec


def test_custom_section_geometric(custom_section):
    """Test custom section geometric properties.

    Geometric properties calculated from sectionproperties v3.0.2 with a refined mesh
    [mesh_sizes=0.5].
    """
    sec = custom_section

    check.almost_equal(sec.section_props.area, 4250, rel=tol)
    check.almost_equal(sec.section_props.perimeter, 6.624264e02, rel=tol)
    check.almost_equal(sec.section_props.mass, 4250, rel=tol)
    check.almost_equal(sec.section_props.ea, 4250, rel=tol)
    check.almost_equal(sec.section_props.qx, 2.763333e05, rel=tol)
    check.almost_equal(sec.section_props.qy, 2.096667e05, rel=tol)
    check.almost_equal(sec.section_props.ixx_g, 2.567250e07, rel=tol)
    check.almost_equal(sec.section_props.iyy_g, 1.418583e07, rel=tol)
    check.almost_equal(sec.section_props.ixy_g, 1.379792e07, rel=tol)
    check.almost_equal(sec.section_props.cx, 4.933333e01, rel=tol)
    check.almost_equal(sec.section_props.cy, 6.501961e01, rel=tol)
    check.almost_equal(sec.section_props.ixx_c, 7.705415e06, rel=tol)
    check.almost_equal(sec.section_props.iyy_c, 3.842278e06, rel=tol)
    check.almost_equal(sec.section_props.ixy_c, 1.654722e05, rel=tol)
    check.almost_equal(sec.section_props.zxx_plus, 1.713061e05, rel=tol)
    check.almost_equal(sec.section_props.zxx_minus, 1.185091e05, rel=tol)
    check.almost_equal(sec.section_props.zyy_plus, 6.333425e04, rel=tol)
    check.almost_equal(sec.section_props.zyy_minus, 6.475749e04, rel=tol)
    check.almost_equal(sec.section_props.rx_c, 4.257979e01, rel=tol)
    check.almost_equal(sec.section_props.ry_c, 3.006768e01, rel=tol)
    check.almost_equal(sec.section_props.i11_c, 7.712490e06, rel=tol)
    check.almost_equal(sec.section_props.i22_c, 3.835203e06, rel=tol)
    check.almost_equal(sec.section_props.phi, -2.448209e00, rel=tol)
    check.almost_equal(sec.section_props.z11_plus, 1.622630e05, rel=tol)
    check.almost_equal(sec.section_props.z11_minus, 1.142680e05, rel=tol)
    check.almost_equal(sec.section_props.z22_plus, 6.050295e04, rel=tol)
    check.almost_equal(sec.section_props.z22_minus, 6.266613e04, rel=tol)
    check.almost_equal(sec.section_props.r11_c, 4.259934e01, rel=tol)
    check.almost_equal(sec.section_props.r22_c, 3.003998e01, rel=tol)


def test_custom_section_warping(custom_section):
    """Test custom section warping properties.

    Warping properties calculated from sectionproperties v3.0.2 with a refined mesh
    [mesh_sizes=0.5].
    """
    sec = custom_section
    sec.calculate_warping_properties()

    x_se, y_se = sec.get_sc()
    x11_se, y22_se = sec.get_sc_p()
    x_st, y_st = sec.get_sc_t()

    check.almost_equal(sec.section_props.j, 3.473083e05, rel=2 * warp_tol)
    check.almost_equal(sec.section_props.gamma, 7.535654e09, rel=warp_tol)
    check.almost_equal(x_se, 5.137486e01, rel=warp_tol)
    check.almost_equal(y_se, 6.795580e01, rel=warp_tol)
    check.almost_equal(x11_se, 1.914242e00, rel=2 * warp_tol)
    check.almost_equal(y22_se, 3.020718e00, rel=2 * warp_tol)
    check.almost_equal(x_st, 5.137486e01, rel=warp_tol)
    check.almost_equal(y_st, 6.795580e01, rel=warp_tol)
    check.almost_equal(sec.section_props.a_sx, 2.952086e03, rel=warp_tol)
    check.almost_equal(sec.section_props.a_sy, 9.545730e02, rel=2 * warp_tol)
    check.almost_equal(sec.section_props.a_s11, 2.943977e03, rel=warp_tol)
    check.almost_equal(sec.section_props.a_s22, 9.554239e02, rel=2 * warp_tol)
    check.almost_equal(sec.section_props.beta_x_plus, 2.530083e01, rel=warp_tol)
    check.almost_equal(sec.section_props.beta_x_minus, -2.530083e01, rel=warp_tol)
    check.almost_equal(sec.section_props.beta_y_plus, 5.645168e00, rel=warp_tol)
    check.almost_equal(sec.section_props.beta_y_minus, -5.645168e00, rel=warp_tol)
    check.almost_equal(sec.section_props.beta_11_plus, 2.546759e01, rel=warp_tol)
    check.almost_equal(sec.section_props.beta_11_minus, -2.546759e01, rel=warp_tol)
    check.almost_equal(sec.section_props.beta_22_plus, 3.724649e00, rel=2 * warp_tol)
    check.almost_equal(sec.section_props.beta_22_minus, -3.724649e00, rel=2 * warp_tol)


def test_custom_section_plastic(custom_section):
    """Test custom section plastic properties.

    Plastic properties calculated from sectionproperties v3.0.2 with a refined mesh
    [mesh_sizes=0.5].
    """
    sec = custom_section
    sec.calculate_plastic_properties()

    x_pc, y_pc = sec.get_pc()
    x11_pc, y22_pc = sec.get_pc_p()

    check.almost_equal(x_pc, 4.977273e01, rel=plastic_tol)
    check.almost_equal(y_pc, 9.172040e01, rel=plastic_tol)
    check.almost_equal(x11_pc, 5.133714e01, rel=plastic_tol)
    check.almost_equal(y22_pc, 9.158984e01, rel=plastic_tol)
    check.almost_equal(sec.section_props.sxx, 1.531971e05, rel=plastic_tol)
    check.almost_equal(sec.section_props.syy, 1.014943e05, rel=plastic_tol)
    check.almost_equal(sec.section_props.s11, 1.533463e05, rel=plastic_tol)
    check.almost_equal(sec.section_props.s22, 1.015010e05, rel=plastic_tol)
    check.almost_equal(sec.section_props.sf_xx_plus, 8.942884e-01, rel=plastic_tol)
    check.almost_equal(sec.section_props.sf_xx_minus, 1.292703e00, rel=plastic_tol)
    check.almost_equal(sec.section_props.sf_yy_plus, 1.602519e00, rel=plastic_tol)
    check.almost_equal(sec.section_props.sf_yy_minus, 1.567298e00, rel=plastic_tol)
    check.almost_equal(sec.section_props.sf_11_plus, 9.450478e-01, rel=plastic_tol)
    check.almost_equal(sec.section_props.sf_11_minus, 1.341988e00, rel=plastic_tol)
    check.almost_equal(sec.section_props.sf_22_plus, 1.677621e00, rel=plastic_tol)
    check.almost_equal(sec.section_props.sf_22_minus, 1.619711e00, rel=plastic_tol)
