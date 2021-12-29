import pytest_check as check
import numpy as np
import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import Section


r_tol = 1e-3


def test_rectangular_perimeter():
    rect = sections.rectangular_section(d=500, b=300)
    rect.create_mesh([200])
    section = Section(rect)
    section.calculate_geometric_properties()
    assert section.get_perimeter() == 2 * (500 + 300)


def test_i_section():
    i_section = sections.i_section(d=308, b=305, t_f=15.4, t_w=9.9, r=16.5, n_r=16)
    i_section.create_mesh([100])
    section = Section(i_section)
    section.calculate_geometric_properties()
    perim = (
        (2 * 305)
        + (4 * 15.4)
        + 2 * (305 - 9.9 - 2 * 16.5)
        + (2 * np.pi * 16.5)
        + 2 * (308 - 2 * 15.4 - 2 * 16.5)
    )
    check.almost_equal(section.get_perimeter(), perim, rel=r_tol)


def test_box_girder_perimeter():
    box_girder = sections.box_girder_section(
        d=400, b_t=700, b_b=100, t_ft=20, t_fb=20, t_w=12
    )
    box_girder.create_mesh([100])
    section = Section(box_girder)
    section.calculate_geometric_properties()
    assert section.get_perimeter() == 700 + 100 + 2 * 500


def test_custom_geometry_perimeter():
    points = [[0, 0], [5, 0], [11, 8], [3, 2], [0, 2]]
    facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 0]]
    control_points = [[5, 5]]
    custom = sections.Geometry.from_points(points, facets, control_points, holes=None)
    custom.create_mesh([100])
    section = Section(custom)
    section.calculate_geometric_properties()
    assert section.get_perimeter() == 5 + 2 + 10 + 10 + 3


def test_compound_rectangular_perimeter():
    rect1 = sections.rectangular_section(d=100, b=100)
    rect2 = sections.rectangular_section(d=50, b=50).align_to(rect1, on="top")
    rect3 = sections.rectangular_section(d=50, b=100).align_to(rect1, on="right")
    rect4 = (
        sections.rectangular_section(d=50, b=50)
        .align_to(rect3, on="bottom")
        .align_to(rect3, on="right", inner=True)
    )
    geom = rect1 + rect2 + rect3 + rect4
    geom.create_mesh([100])
    section = Section(geom)
    section.calculate_geometric_properties()
    assert section.get_perimeter() == (
        150 + 50 + 50 + 100 + 100 + 50 + 50 + 50 + 50 + 150
    )


def test_compound_rectangular_isection_perimeter1():
    d = 300
    b = 150
    tf = 10
    tw = 6
    r = 12
    b_p = 250
    t_p = 16
    ub = sections.i_section(d=d, b=b, t_f=tf, t_w=tw, r=r, n_r=16)
    plate = (
        sections.rectangular_section(b=b_p, d=t_p)
        .align_center(ub)
        .align_to(ub, on="top")
    )
    geom = ub + plate
    geom.create_mesh([100])
    section = Section(geom)
    section.calculate_geometric_properties()
    perim = (
        b
        + (4 * tf)
        + 2 * (b - tw - 2 * r)
        + (2 * np.pi * r)
        + 2 * (d - 2 * tf - 2 * r)
        + (b_p - b)
        + (2 * t_p)
        + b_p
    )
    check.almost_equal(section.get_perimeter(), perim, rel=r_tol)


def test_compound_rectangular_isection_perimeter2():
    i_section = sections.i_section(d=308, b=305, t_f=15.4, t_w=9.9, r=16.5, n_r=16)
    rect1 = (
        sections.rectangular_section(d=330, b=16)
        .align_center(i_section)
        .align_to(i_section, "left")
    )
    rect2 = (
        sections.rectangular_section(d=330, b=16)
        .align_center(i_section)
        .align_to(i_section, "right")
    )
    geom = i_section + rect1 + rect2
    geom.create_mesh([100])
    section = Section(geom)
    section.calculate_geometric_properties()
    assert section.get_perimeter() == 2 * 330 + 4 * 16 + 2 * 305 + 2 * (330 - 308)


def test_compound_rhs_isection_perimeter():
    d = 200
    b = 150
    t = 9
    r = 15
    b_p = 250
    t_p = 16
    rhs = sections.rectangular_hollow_section(d=d, b=b, t=t, r_out=r, n_r=16)
    plate1 = (
        sections.rectangular_section(b=b_p, d=t_p)
        .align_center(rhs)
        .align_to(rhs, on="top")
    )
    plate2 = (
        sections.rectangular_section(b=b_p, d=t_p)
        .align_center(rhs)
        .align_to(rhs, on="bottom")
    )
    geom = rhs + plate1 + plate2
    geom.create_mesh([100])
    section = Section(geom)
    section.calculate_geometric_properties()
    perim = (
        (2 * b_p)
        + (4 * t_p)
        + 2 * (b_p - b + 2 * r)
        + (2 * np.pi * r)
        + 2 * (d - 2 * r)
    )
    check.almost_equal(section.get_perimeter(), perim, rel=r_tol)
