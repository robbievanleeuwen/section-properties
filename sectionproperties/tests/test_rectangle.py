import pytest_check as check

# import unittest
import sectionproperties.pre.library.primitive_sections as primitive_sections
import sectionproperties.pre.pre as pre
from sectionproperties.analysis.section import Section
from sectionproperties.tests.helper_functions import validate_properties

import sectionproperties.analysis.section as file

# Rectangle section setup
rectangle_geometry = primitive_sections.rectangular_section(b=50, d=100)
rectangle_geometry.create_mesh(mesh_sizes=100)
rectangle_section = Section(rectangle_geometry)
rectangle_section.calculate_geometric_properties()
rectangle_section.calculate_warping_properties()
rectangle_section.calculate_plastic_properties()

tol = 1e-6
warp_tol = 1e-4


def test_rectangular_section_geometric():
    check.almost_equal(rectangle_section.section_props.area, 100 * 50, rel=tol)
    check.almost_equal(
        rectangle_section.section_props.perimeter, 2 * 100 + 2 * 50, rel=tol
    )
    check.almost_equal(rectangle_section.section_props.mass, 1 * 100 * 50, rel=tol)
    check.almost_equal(rectangle_section.section_props.ea, 1 * 100 * 50, rel=tol)
    check.almost_equal(rectangle_section.section_props.qx, 100 * 50 * 50, rel=tol)
    check.almost_equal(rectangle_section.section_props.qy, 100 * 50 * 25, rel=tol)
    check.almost_equal(
        rectangle_section.section_props.ixx_g, 50 * 100 ** 3 / 3, rel=tol
    )
    check.almost_equal(
        rectangle_section.section_props.iyy_g, 100 * 50 ** 3 / 3, rel=tol
    )
    check.almost_equal(
        rectangle_section.section_props.ixy_g, 100 * 50 * 50 * 25, rel=tol
    )
    check.almost_equal(rectangle_section.section_props.cx, 50 / 2, rel=tol)
    check.almost_equal(rectangle_section.section_props.cy, 100 / 2, rel=tol)
    check.almost_equal(
        rectangle_section.section_props.ixx_c, 50 * 100 ** 3 / 12, rel=tol
    )
    check.almost_equal(
        rectangle_section.section_props.iyy_c, 100 * 50 ** 3 / 12, rel=tol
    )
    check.almost_equal(rectangle_section.section_props.ixy_c, 0, abs=tol)
    check.almost_equal(
        rectangle_section.section_props.zxx_plus, 50 * 100 ** 2 / 6, rel=tol
    )
    check.almost_equal(
        rectangle_section.section_props.zxx_minus, 50 * 100 ** 2 / 6, rel=tol
    )
    check.almost_equal(
        rectangle_section.section_props.zyy_plus, 100 * 50 ** 2 / 6, rel=tol
    )
    check.almost_equal(
        rectangle_section.section_props.zyy_minus, 100 * 50 ** 2 / 6, rel=tol
    )
    check.almost_equal(
        rectangle_section.section_props.rx_c,
        (50 * 100 ** 3 / 12 / 100 / 50) ** 0.5,
        rel=tol,
    )
    check.almost_equal(
        rectangle_section.section_props.ry_c,
        (100 * 50 ** 3 / 12 / 100 / 50) ** 0.5,
        rel=tol,
    )
    check.almost_equal(
        rectangle_section.section_props.i11_c, (50 * 100 ** 3 / 12), rel=tol
    )
    check.almost_equal(
        rectangle_section.section_props.i22_c, (100 * 50 ** 3 / 12), rel=tol
    )
    check.almost_equal(rectangle_section.section_props.phi, 0, rel=tol)
    check.almost_equal(
        rectangle_section.section_props.z11_plus, 50 * 100 ** 2 / 6, rel=tol
    )
    check.almost_equal(
        rectangle_section.section_props.z11_minus, 50 * 100 ** 2 / 6, rel=tol
    )
    check.almost_equal(
        rectangle_section.section_props.z22_plus, 100 * 50 ** 2 / 6, rel=tol
    )
    check.almost_equal(
        rectangle_section.section_props.z22_minus, 100 * 50 ** 2 / 6, rel=tol
    )
    check.almost_equal(
        rectangle_section.section_props.r11_c,
        (50 * 100 ** 3 / 12 / 100 / 50) ** 0.5,
        rel=tol,
    )
    check.almost_equal(
        rectangle_section.section_props.r22_c,
        (100 * 50 ** 3 / 12 / 100 / 50) ** 0.5,
        rel=tol,
    )


def test_rectangular_section_plastic():
    check.almost_equal(rectangle_section.get_pc(), (50 / 2, 100 / 2))
    check.almost_equal(rectangle_section.get_pc_p(), (50 / 2, 100 / 2))
    check.almost_equal(
        rectangle_section.get_s(), (50 * 100 ** 2 / 4, 100 * 50 ** 2 / 4)
    )
    check.almost_equal(
        rectangle_section.get_sp(), (50 * 100 ** 2 / 4, 100 * 50 ** 2 / 4)
    )
    check.almost_equal(rectangle_section.get_sf(), (1.5, 1.5, 1.5, 1.5))
    check.almost_equal(rectangle_section.get_sf_p(), (1.5, 1.5, 1.5, 1.5))


def test_rectangular_section_warping():
    check.almost_equal(
        rectangle_section.section_props.j, 2861002, rel=2 * warp_tol
    )  # roark's
    # check.almost_equal(rectangle_section.section_props.j, 2.85852e6, rel=2e-5) #st7
    check.almost_equal(
        rectangle_section.section_props.j, 2.861326e06, rel=tol
    )  # main branch
    check.almost_equal(
        rectangle_section.section_props.gamma, 3.177234e08, rel=tol
    )  # main branch
    check.almost_equal(rectangle_section.get_sc(), (50 / 2, 100 / 2), rel=warp_tol)
    check.almost_equal(
        rectangle_section.get_sc_p(), (-4.103589e-04, 1.164891e-03), rel=tol
    )
    check.almost_equal(rectangle_section.get_sc_t(), (50 / 2, 100 / 2), rel=warp_tol)
    check.almost_equal(rectangle_section.get_As(), (4.168418e03, 4.166821e03), rel=tol)
    check.almost_equal(
        rectangle_section.get_As_p(), (4.168418e03, 4.166821e03), rel=tol
    )
