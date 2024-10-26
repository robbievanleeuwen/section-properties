"""Test cases for the offset perimeter method."""

from __future__ import annotations

import numpy as np
import pytest_check as check
from shapely import Polygon

import sectionproperties.pre.library.primitive_sections as sections
import sectionproperties.pre.library.steel_sections as steel_sections
from sectionproperties.analysis.section import Section
from sectionproperties.pre.geometry import Geometry, check_geometry_overlaps

r_tol = 1e-3


def test_rectangular_offset():
    """Tests simple rectangular offsets."""
    # exterior negative offset
    rect = sections.rectangular_section(d=500, b=300)
    rect = rect.offset_perimeter(amount=-10, where="exterior")
    rect.create_mesh(mesh_sizes=[0])
    section = Section(geometry=rect)
    section.calculate_geometric_properties()
    area = 480 * 280
    check.almost_equal(section.get_area(), area, rel=r_tol)

    # exterior positive offset
    rect = sections.rectangular_section(d=500, b=300)
    rect = rect.offset_perimeter(amount=10, where="exterior")
    rect.create_mesh(mesh_sizes=[0])
    section = Section(geometry=rect)
    section.calculate_geometric_properties()
    area = 520 * 320 - (20 * 20 - np.pi * 10 * 10)
    check.almost_equal(section.get_area(), area, rel=r_tol)


def test_box_offset():
    """Tests offsets for RHS sections."""
    # exterior negative offset
    box = steel_sections.rectangular_hollow_section(d=200, b=100, t=10, r_out=0, n_r=1)
    box = box.offset_perimeter(amount=-5, where="exterior")
    box.create_mesh(mesh_sizes=[0])
    section = Section(geometry=box)
    section.calculate_geometric_properties()
    area = 190 * 90 - 180 * 80
    check.almost_equal(section.get_area(), area, rel=r_tol)

    # exterior positve offset
    box = steel_sections.rectangular_hollow_section(d=200, b=100, t=10, r_out=0, n_r=1)
    box = box.offset_perimeter(amount=5, where="exterior")
    box.create_mesh(mesh_sizes=[0])
    section = Section(geometry=box)
    section.calculate_geometric_properties()
    area = 210 * 110 - (10 * 10 - np.pi * 5 * 5) - 180 * 80
    check.almost_equal(section.get_area(), area, rel=r_tol)

    # interior negative offset
    box = steel_sections.rectangular_hollow_section(d=200, b=100, t=10, r_out=0, n_r=1)
    box = box.offset_perimeter(amount=-5, where="interior")
    box.create_mesh(mesh_sizes=[0])
    section = Section(geometry=box)
    section.calculate_geometric_properties()
    area = 200 * 100 - 170 * 70
    check.almost_equal(section.get_area(), area, rel=r_tol)

    # interior positive offset
    box = steel_sections.rectangular_hollow_section(d=200, b=100, t=10, r_out=0, n_r=1)
    box = box.offset_perimeter(amount=5, where="interior")
    box.create_mesh(mesh_sizes=[0])
    section = Section(geometry=box)
    section.calculate_geometric_properties()
    area = 200 * 100 - 190 * 90 + (10 * 10 - np.pi * 5 * 5)
    check.almost_equal(section.get_area(), area, rel=r_tol)

    # all negative offset
    box = steel_sections.rectangular_hollow_section(d=200, b=100, t=10, r_out=0, n_r=1)
    box = box.offset_perimeter(amount=-2.5, where="all")
    box.create_mesh(mesh_sizes=[50])
    section = Section(geometry=box)
    section.calculate_geometric_properties()
    area = 195 * 95 - 185 * 85 + (5 * 5 - np.pi * 2.5 * 2.5)
    check.almost_equal(section.get_area(), area, rel=r_tol)

    # all positive offset
    box = steel_sections.rectangular_hollow_section(d=200, b=100, t=10, r_out=0, n_r=1)
    box = box.offset_perimeter(amount=5, where="all")
    box.create_mesh(mesh_sizes=[0])
    section = Section(geometry=box)
    section.calculate_geometric_properties()
    area = 210 * 110 - (10 * 10 - np.pi * 5 * 5) - 170 * 70
    check.almost_equal(section.get_area(), area, rel=r_tol)


def test_compound_rectangular_offset():
    """Tests offsets for compound rectangular sections."""
    rect1 = sections.rectangular_section(d=50, b=50)
    rect2 = sections.rectangular_section(d=50, b=50).align_to(rect1, "right")
    geom = rect1 + rect2

    # shrink
    geom_neg = geom.offset_perimeter(amount=-5, where="exterior")
    geom_neg.create_mesh(mesh_sizes=[0])
    section = Section(geometry=geom_neg)
    section.calculate_geometric_properties()
    area = 90 * 40
    check.almost_equal(section.get_area(), area, rel=r_tol)

    # balloon case (two rectangles)
    geom = rect1 + rect2
    geom = geom.offset_perimeter(amount=5, where="exterior")

    # ensure there are no overlaps
    assert not check_geometry_overlaps([g.geom for g in geom.geoms])

    # calculate area
    geom.create_mesh(mesh_sizes=[0])
    section = Section(geometry=geom)
    section.calculate_geometric_properties()
    area = 100 * 50 + 2 * (5 * 100 + 5 * 50) + np.pi * 5**2
    check.almost_equal(section.get_area(), area, rel=r_tol)

    # balloon case (three rectangles)
    geom = rect1 + rect2 + rect1.align_to(rect2, "right")
    geom = geom.offset_perimeter(amount=5, where="exterior")

    # ensure there are no overlaps
    assert not check_geometry_overlaps([g.geom for g in geom.geoms])

    # calculate area
    geom.create_mesh(mesh_sizes=[0])
    section = Section(geometry=geom)
    section.calculate_geometric_properties()
    area = 150 * 50 + 2 * (5 * 150 + 5 * 50) + np.pi * 5**2
    check.almost_equal(section.get_area(), area, rel=r_tol)


def test_compound_rectangular_isection_offset_corrode():
    """Tests offsets for a complex compound section."""
    d = 300
    b = 150
    tf = 10
    tw = 8
    r = 12
    b_p = 250
    t_p = 16
    ub = steel_sections.i_section(d=d, b=b, t_f=tf, t_w=tw, r=r, n_r=16)
    plate = (
        sections.rectangular_section(b=b_p, d=t_p)
        .align_center(align_to=ub)
        .align_to(other=ub, on="top")
    )
    geom_test = ub + plate
    geom_test = geom_test.offset_perimeter(amount=-2, where="exterior")
    geom_test.create_mesh(mesh_sizes=[0])
    section_test = Section(geometry=geom_test)
    section_test.calculate_geometric_properties()

    ub_corroded = steel_sections.mono_i_section(
        d=298,
        b_t=146,
        b_b=146,
        t_ft=8,
        t_fb=6,
        t_w=4,
        r=14,
        n_r=16,
    )
    plate_corroded1 = (
        sections.rectangular_section(b=146, d=2)
        .align_center(align_to=ub_corroded)
        .align_to(other=ub_corroded, on="top")
    )
    plate_corroded2 = (
        sections.rectangular_section(b=246, d=12)
        .align_center(align_to=ub_corroded)
        .align_to(other=plate_corroded1, on="top")
    )
    rad_l = (
        draw_radius(2, 8)
        .align_to(other=plate_corroded1, on="left")
        .align_to(other=plate_corroded2, on="bottom")
    )
    rad_r = (
        draw_radius(2, 8)
        .mirror_section("y", (2, 0))
        .align_to(other=plate_corroded1, on="right")
        .align_to(other=plate_corroded2, on="bottom")
    )
    geom_corroded = ub_corroded + plate_corroded1 + plate_corroded2 + rad_l + rad_r
    geom_corroded.create_mesh(mesh_sizes=[0])
    section_corroded = Section(geometry=geom_corroded)
    section_corroded.calculate_geometric_properties()

    check.almost_equal(section_test.get_area(), section_corroded.get_area(), rel=r_tol)


def test_compound_stiffened_isection():
    """Tests full erosion of a geometry.

    Tests that plates 1 and 2 can be eroded to nothing and a valid Section can
    still be generated without errors.
    """
    uc = steel_sections.i_section(d=400, b=400, t_f=25, t_w=25, r=30, n_r=8)
    plate1 = (
        sections.rectangular_section(b=500, d=10)
        .align_center(align_to=uc)
        .align_to(other=uc, on="top")
    )
    plate2 = (
        sections.rectangular_section(b=500, d=10)
        .align_center(align_to=uc)
        .align_to(other=uc, on="bottom")
    )
    geom = uc + plate1 + plate2

    new_geom = geom.offset_perimeter(-9)
    new_geom.create_mesh(mesh_sizes=[0])
    Section(geometry=new_geom)

    new_geom = geom.offset_perimeter(-10)
    new_geom.create_mesh(mesh_sizes=[0])
    Section(geometry=new_geom)

    new_geom = geom.offset_perimeter(-11)
    new_geom.create_mesh(mesh_sizes=[0])
    Section(geometry=new_geom)


def draw_radius(r: float, n: int) -> Geometry:
    """Draws a radius.

    Args:
        r: Radius
        n: Number of points

    Returns:
        Geometry object describing the radius
    """
    points = []

    # calculate radius of points
    for i in range(n):
        # determine angle
        t = i * 1.0 / max(1, n - 1) * np.pi * 0.5

        x = r * np.cos(t)
        y = r * np.sin(t)
        points.append([x, y])

    points.append([r, r])

    return Geometry(Polygon(points))
