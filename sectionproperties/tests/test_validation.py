import pytest_check as check
from shapely.geometry import Polygon
from sectionproperties.pre.geometry import Geometry
import sectionproperties.pre.library.primitive_sections as sections
import sectionproperties.pre.library.steel_sections as steel_sections
from sectionproperties.analysis.section import Section
from sectionproperties.tests.helper_functions import validate_properties


# Setup for angle section
angle = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
angle.create_mesh(mesh_sizes=2.5)
angle_section = Section(angle)
angle_section.calculate_geometric_properties()
angle_section.calculate_plastic_properties()
angle_section.calculate_warping_properties()


def test_angle_all_properties():
    check.almost_equal(angle_section.section_props.area, 2747.059)
    check.almost_equal(angle_section.section_props.perimeter, 471.3501)
    check.almost_equal(angle_section.section_props.cx, 2.122282e1)
    check.almost_equal(angle_section.section_props.cy, 5.098127e1)
    check.almost_equal(angle_section.section_props.ixx_g, 1.342632e7)
    check.almost_equal(angle_section.section_props.iyy_g, 2.955753e6)
    check.almost_equal(angle_section.section_props.ixy_g, 1.086603e6)
    check.almost_equal(angle_section.section_props.ixx_c, 6.286470e6)
    check.almost_equal(angle_section.section_props.iyy_c, 1.718455e6)
    check.almost_equal(angle_section.section_props.ixy_c, -1.885622e6)
    check.almost_equal(angle_section.section_props.zxx_plus, 6.348769e4)
    check.almost_equal(angle_section.section_props.zxx_minus, 1.233094e5)
    check.almost_equal(angle_section.section_props.zyy_plus, 2.498584e04)
    check.almost_equal(angle_section.section_props.zyy_minus, 8.097207e4)
    check.almost_equal(angle_section.section_props.rx_c, 4.783761e1)
    check.almost_equal(angle_section.section_props.ry_c, 2.501124e1)
    check.almost_equal(angle_section.section_props.i11_c, 6.964263e6)
    check.almost_equal(angle_section.section_props.i22_c, 1.040662e6)
    check.almost_equal(angle_section.section_props.phi, -1.602289e2)
    check.almost_equal(angle_section.section_props.z11_plus, 9.775662e4)
    check.almost_equal(angle_section.section_props.z11_minus, 6.939239e4)
    check.almost_equal(angle_section.section_props.z22_plus, 2.796211e4)
    check.almost_equal(angle_section.section_props.z22_minus, 2.076613e4)
    check.almost_equal(angle_section.section_props.r11_c, 5.035048e1)
    check.almost_equal(angle_section.section_props.r22_c, 1.946350e1)
    check.almost_equal(angle_section.section_props.sxx, 1.135392e5)
    check.almost_equal(
        angle_section.section_props.syy, 4.572267e4
    )  # Altered from 4.572269e4
    check.almost_equal(angle_section.section_props.sf_xx_plus, 1.788366)
    check.almost_equal(angle_section.section_props.sf_xx_minus, 9.207672e-1)
    check.almost_equal(
        angle_section.section_props.sf_yy_plus, 1.829943
    )  # Altered from 1.829944
    check.almost_equal(
        angle_section.section_props.sf_yy_minus, 5.646721e-1
    )  # Altered from 5.646723e-1
    check.almost_equal(angle_section.section_props.s11, 1.210275e5)
    check.almost_equal(angle_section.section_props.s22, 4.376054e4)
    check.almost_equal(angle_section.section_props.sf_11_plus, 1.238049)
    check.almost_equal(angle_section.section_props.sf_11_minus, 1.744103)
    check.almost_equal(angle_section.section_props.sf_22_plus, 1.564994)
    check.almost_equal(angle_section.section_props.sf_22_minus, 2.107303)
    check.almost_equal(
        angle_section.section_props.j, 1.354663e5
    )  # Altered from 1.354663e5
    check.almost_equal(angle_section.section_props.gamma, 162220735.49)
    check.almost_equal(angle_section.section_props.A_s11, 8.855951e2)
    check.almost_equal(angle_section.section_props.A_s22, 1.460240e3)
    check.almost_equal(angle_section.section_props.x11_se, 2.870404e1)
    check.almost_equal(angle_section.section_props.y22_se, 3.522141e1)


# Setup custom section
custom_geom_points = [
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
custom_geom = Geometry(Polygon(custom_geom_points))
custom_geom.create_mesh(mesh_sizes=5)
custom_section = Section(custom_geom)
custom_section.calculate_geometric_properties()
custom_section.calculate_plastic_properties()
custom_section.calculate_warping_properties()


def test_custom_section_all_properties():
    check.almost_equal(custom_section.section_props.area, 4250)
    check.almost_equal(custom_section.section_props.cx, 4.933333e1)
    check.almost_equal(custom_section.section_props.cy, 6.501961e1)
    check.almost_equal(custom_section.section_props.ixx_g, 2.567250e7)
    check.almost_equal(custom_section.section_props.iyy_g, 1.418583e7)
    check.almost_equal(custom_section.section_props.ixy_g, 1.379792e7)
    check.almost_equal(custom_section.section_props.ixx_c, 7.705415e6)
    check.almost_equal(custom_section.section_props.iyy_c, 3.842278e6)
    check.almost_equal(custom_section.section_props.ixy_c, 1.654722e5)
    check.almost_equal(custom_section.section_props.zxx_plus, 1.713061e5)
    check.almost_equal(custom_section.section_props.zxx_minus, 1.185091e5)
    check.almost_equal(custom_section.section_props.zyy_plus, 6.333425e4)
    check.almost_equal(custom_section.section_props.zyy_minus, 6.475749e4)
    check.almost_equal(custom_section.section_props.rx_c, 4.257979e01)
    check.almost_equal(custom_section.section_props.ry_c, 3.006768e1)
    check.almost_equal(custom_section.section_props.phi, -2.448209)
    check.almost_equal(custom_section.section_props.i11_c, 7.712490e6)
    check.almost_equal(custom_section.section_props.i22_c, 3.835203e6)
    check.almost_equal(custom_section.section_props.z11_plus, 1.622630e5)
    check.almost_equal(custom_section.section_props.z11_minus, 1.142680e5)
    check.almost_equal(custom_section.section_props.z22_plus, 6.050295e4)
    check.almost_equal(custom_section.section_props.z22_minus, 6.266613e4)
    check.almost_equal(custom_section.section_props.r11_c, 4.259934e1)
    check.almost_equal(custom_section.section_props.r22_c, 3.003998e1)
    check.almost_equal(custom_section.section_props.sxx, 1.531971e5)
    check.almost_equal(custom_section.section_props.syy, 1.014943e5)
    check.almost_equal(
        custom_section.get_sf(), (8.942884e-01, 1.292703, 1.602519, 1.567298)
    )
    check.almost_equal(custom_section.section_props.s11, 1.533463e5)
    check.almost_equal(custom_section.section_props.s22, 1.015010e5)
    check.almost_equal(custom_section.section_props.sf_11_plus, 9.450478e-1)
    check.almost_equal(custom_section.section_props.sf_11_minus, 1.341988)
    check.almost_equal(custom_section.section_props.sf_22_plus, 1.677621)
    check.almost_equal(custom_section.section_props.sf_22_minus, 1.619711)
    check.almost_equal(custom_section.section_props.j, 3.477399e5)
    check.almost_equal(custom_section.section_props.gamma, 7.532929e9)
    check.almost_equal(custom_section.section_props.A_s11, 2.945692e3)
    check.almost_equal(custom_section.section_props.A_s22, 9.564143e2)
    check.almost_equal(custom_section.section_props.x11_se, 1.916270)
    check.almost_equal(custom_section.section_props.y22_se, 3.017570)
