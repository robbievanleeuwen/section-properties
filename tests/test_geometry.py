"""Tests for various geometry methods and classes."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest
from shapely import (
    GeometryCollection,
    LineString,
    MultiPolygon,
    Point,
    Polygon,
    box,
    wkt,
)

import sectionproperties.pre.geometry as sp_geom
import sectionproperties.pre.library.nastran_sections as nastran_sections
import sectionproperties.pre.library.primitive_sections as primitive_sections
import sectionproperties.pre.library.steel_sections as steel_sections
from sectionproperties.analysis.section import Section
from sectionproperties.pre.pre import Material


big_sq = primitive_sections.rectangular_section(d=300, b=250)
small_sq = primitive_sections.rectangular_section(d=100, b=75)
small_hole = primitive_sections.rectangular_section(d=40, b=30).align_center(
    align_to=small_sq
)
i_sec = steel_sections.i_section(d=200, b=100, t_f=20, t_w=10, r=12, n_r=12)
small_sq_w_hole = small_sq - small_hole
composite = (
    big_sq
    + small_sq_w_hole.align_to(other=big_sq, on="top", inner=True).align_to(
        other=big_sq, on="top"
    )
    + i_sec.align_to(other=big_sq, on="bottom", inner=True).align_to(
        other=big_sq, on="right"
    )
)
composite.create_mesh(mesh_sizes=[200])
comp_sec = Section(geometry=composite)
comp_sec.calculate_geometric_properties()
comp_sec.calculate_plastic_properties()

# Subtractive modelling
nested_geom = (small_sq - small_hole) + small_hole
nested_geom.create_mesh(mesh_sizes=[50])
nested_sec = Section(geometry=nested_geom)

# Overlapped modelling
overlay_geom = small_sq + small_hole
overlay_geom.create_mesh(mesh_sizes=[50])
overlay_sec = Section(geometry=overlay_geom)

steel = Material(
    name="steel",
    elastic_modulus=200e3,
    poissons_ratio=0.3,
    density=7.85e-6,
    yield_strength=400,
    color="grey",
)


def test_material_persistence():
    """Tests material persistence.

    Test ensures that the material attribute gets transformed through all of the
    Geometry transformation methods, each which returns a new Geometry object. The
    material assignment should persist through all of the transformations.
    """
    big_sq.material = steel
    new_geom = (
        big_sq.align_to(other=small_sq, on="left", inner=False)
        .align_center()
        .rotate_section(angle=23)
        .mirror_section(axis="y")
        .offset_perimeter(amount=1)
    )
    assert new_geom.material == steel


def test_for_incidental_holes():
    """Tests for incidental hole creation.

    One hole in the geometry was explicitly created through subtraction. Another hole in
    the geometry was created accidentally by sticking a I Section up against a
    rectangle. There should be two holes created after .compile_geometry()
    """
    assert len(composite.holes) == 2
    assert len(nested_geom.holes) == 0


def test__sub__():
    """Tests the geometry subtraction method."""
    small_hole.material = steel
    top_left = (
        small_hole.align_to(other=big_sq, on="left")
        .align_to(other=big_sq, on="top")
        .shift_section(x_offset=20, y_offset=-20)
    )
    top_right = top_left.shift_section(x_offset=200)

    compound = big_sq - top_left
    compound = compound + top_left
    compound = compound - top_right
    compound = compound + top_right
    compound = compound - small_hole

    assert len(compound.control_points) == 3
    assert len(compound.holes) == 1


def test_geometry_from_points():
    """Tests the Geometry.from_points() method.

    Geometry.from_points() tests a shape with exactly one exterior and an arbitrary
    number of interiors being built from the legacy points, facets, holes,
    control_points interface of sectionproperties.
    """
    exterior = [(-6, 10), (6, 10), (6, -10), (-6, -10)]
    interior1 = [(-4, 8), (4, 8), (4, 4), (-4, 4)]
    interior2 = [(-4, -8), (4, -8), (4, -4), (-4, -4)]
    points = exterior + interior1 + interior2
    facets = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (8, 9),
        (9, 10),
        (10, 11),
        (11, 7),
    ]
    control_points = [(0.0, 0.0)]
    holes = [(0, 6), (0, -6)]
    material = Material(
        name="mat1",
        elastic_modulus=2,
        poissons_ratio=0.3,
        density=1e-6,
        yield_strength=5,
        color="grey",
    )
    new_geom = sp_geom.Geometry.from_points(
        points=points,
        facets=facets,
        control_points=control_points,
        holes=holes,
        material=material,
    )
    poly = "POLYGON ((6 10, 6 -10, -6 -10, -6 10, 6 10), (-4 4, 4 4, 4 8, -4 8, -4 4), "
    poly += "(4 -8, 4 -4, -4 -4, -4 -8, 4 -8))"
    wkt_test_geom = wkt.loads(data=poly)
    assert (new_geom.geom - wkt_test_geom) == Polygon()
    assert (new_geom.material) == material


def test_compound_geometry_from_points():
    """Tests the CompoundGeometry.from_points() method.

    CompoundGeometry.from_points() tests a shape with an arbitrary number of exteriors
    and an arbitrary number of interiors being built from the legacy points, facets
     holes, control_points interface of sectionproperties.
    """
    a = 1
    b = 2
    t = 0.1

    # build the lists of points, facets, holes and control points
    points = [
        (-t / 2, -2 * a),
        (t / 2, -2 * a),
        (t / 2, -t / 2),
        (a, -t / 2),
        (a, t / 2),
        (-t / 2, t / 2),
        (-b / 2, -2 * a),
        (b / 2, -2 * a),
        (b / 2, -2 * a - t),
        (-b / 2, -2 * a - t),
    ]
    facets = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 4),
        (4, 5),
        (5, 0),
        (6, 7),
        (7, 8),
        (8, 9),
        (9, 6),
    ]
    control_points = [(0, 0), (0, -2 * a - t / 2)]
    mat1 = Material(
        name="mat1",
        elastic_modulus=2,
        poissons_ratio=0.3,
        density=1e-6,
        yield_strength=5,
        color="grey",
    )
    mat2 = Material(
        name="mat2",
        elastic_modulus=5,
        poissons_ratio=0.2,
        density=2e-6,
        yield_strength=10,
        color="blue",
    )
    materials = [mat1, mat2]
    new_geom = sp_geom.CompoundGeometry.from_points(
        points=points, facets=facets, control_points=control_points, materials=materials
    )
    poly = "MULTIPOLYGON (((-0.05 -2, 0.05 -2, 0.05 -0.05, 1 -0.05, 1 0.05, -0.05 0.05,"
    poly += " -0.05 -2)), ((-1 -2, 1 -2, 1 -2.1, -1 -2.1, -1 -2)))"
    wkt_test_geom = wkt.loads(data=poly)
    assert (new_geom.geom - wkt_test_geom) == Polygon()

    # test materials
    for idx, geom in enumerate(new_geom.geoms):
        assert geom.material == materials[idx]


def test_multi_nested_compound_geometry_from_points():
    """Testing a multi-nested section.

    This section contains three nested materials in concentric square rings with a hole
    going through the center of the whole section. This test confirms that the section
    can be successfully built using .from_points, that the control_points and hole nodes
    persist in the right locations, and that the plastic section calculation raises a
    warning because the nested regions overlap.
    """
    points = [
        (-50.0, 50.0),
        (50.0, 50.0),
        (50.0, -50.0),
        (-50.0, -50.0),
        (37.5, -37.5),
        (37.5, 37.5),
        (-37.5, 37.5),
        (-37.5, -37.5),
        (25.0, -25.0),
        (25.0, 25.0),
        (-25.0, 25.0),
        (-25.0, -25.0),
        (12.5, -12.5),
        (12.5, 12.5),
        (-12.5, 12.5),
        (-12.5, -12.5),
    ]
    facets = [
        (0, 1),
        (1, 2),
        (2, 3),
        (3, 0),
        (4, 5),
        (5, 6),
        (6, 7),
        (7, 4),
        (8, 9),
        (9, 10),
        (10, 11),
        (11, 8),
        (12, 13),
        (13, 14),
        (14, 15),
        (15, 12),
    ]
    control_points = [(-43.75, 0.0), (-31.25, 0.0), (-18.75, 0.0)]
    holes = [(0.0, 0.0)]
    mat1 = Material(
        name="mat1",
        elastic_modulus=2,
        poissons_ratio=0.3,
        density=1e-6,
        yield_strength=5,
        color="grey",
    )
    mat2 = Material(
        name="mat2",
        elastic_modulus=5,
        poissons_ratio=0.2,
        density=2e-6,
        yield_strength=10,
        color="blue",
    )
    mat3 = Material(
        name="mat3",
        elastic_modulus=1,
        poissons_ratio=0.25,
        density=1.5e-6,
        yield_strength=3,
        color="green",
    )
    materials = [mat1, mat2, mat3]
    nested_compound = sp_geom.CompoundGeometry.from_points(
        points=points,
        facets=facets,
        control_points=control_points,
        holes=holes,
        materials=materials,
    )
    poly = "MULTIPOLYGON (((50 50, 50 -50, -50 -50, -50 50, 50 50), (12.5 12.5, -12.5 "
    poly += "12.5, -12.5 -12.5, 12.5 -12.5, 12.5 12.5)), ((-37.5 -37.5, -37.5 37.5, "
    poly += "37.5 37.5, 37.5 -37.5, -37.5 -37.5), (12.5 12.5, -12.5 12.5, -12.5 -12.5, "
    poly += "12.5 -12.5, 12.5 12.5)), ((-25 -25, -25 25, 25 25, 25 -25, -25 -25), "
    poly += "(12.5 12.5, -12.5 12.5, -12.5 -12.5, 12.5 -12.5, 12.5 12.5)))"
    wkt_test_geom = wkt.loads(data=poly)
    assert (nested_compound.geom - wkt_test_geom) == Polygon()

    assert nested_compound.control_points == [
        (-43.75, 0.0),
        (-31.25, 0.0),
        (-18.75, 0.0),
    ]
    assert nested_compound.holes == [(0, 0)]

    # test materials
    for idx, geom in enumerate(nested_compound.geoms):
        assert geom.material == materials[idx]

    # Section contains overlapping geometries which will result in potentially incorrect
    # plastic properties calculation (depends on user intent and geometry).
    # Test to ensure a warning is raised about this to notify the user.
    nested_compound.create_mesh([25, 30, 35])
    nested_compound_sec = Section(nested_compound)
    nested_compound_sec.calculate_geometric_properties()
    with pytest.warns(UserWarning):
        nested_compound_sec.calculate_plastic_properties()


def test_geometry_from_dxf():
    """Tests loading geometry from a .dxf file."""
    section_holes_dxf = Path(__file__).parent.absolute() / "section_holes.dxf"
    # print(section_holes_dxf)
    poly = "POLYGON ((-0.338658834889 -0.395177702895, -0.338658834889 29.092318216393,"
    poly += " 31.962257588776 29.092318216393, 31.962257588776 -0.395177702895, "
    poly += "-0.338658834889 -0.395177702895), (16.684315862478 2.382629883704, "
    poly += "29.683030851053 2.382629883704, 29.683030851053 24.355800152063, "
    poly += "16.684315862478 24.355800152063, 16.684315862478 2.382629883704), "
    poly += "(1.548825807288 3.344178663681, 14.547540795863 3.344178663681, "
    poly += "14.547540795863 27.382898163101, 1.548825807288 27.382898163101, "
    poly += "1.548825807288 3.344178663681))"
    assert sp_geom.Geometry.from_dxf(section_holes_dxf).geom.wkt == poly


@pytest.mark.skipif(sys.version_info >= (3, 10), reason="requires python 3.9")
def test_geometry_from_3dm_file_simple():
    """Tests loading geometry from a simple .3dm file."""
    section = Path(__file__).parent.absolute() / "3in x 2in.3dm"
    exp = Polygon([(0, 0), (0, 3), (2, 3), (2, 0), (0, 0)])
    test = sp_geom.Geometry.from_3dm(filepath=section)
    assert (test.geom - exp).is_empty


@pytest.mark.skipif(sys.version_info >= (3, 10), reason="requires python 3.9")
def test_geometry_from_3dm_file_complex():
    """Tests loading geometry from a complex .3dm file."""
    section_3dm = Path(__file__).parent.absolute() / "complex_shape.3dm"
    section_wkt = Path(__file__).parent.absolute() / "complex_shape.txt"
    with open(section_wkt) as file:
        wkt_str = file.readlines()
    exp = wkt.loads(wkt_str[0])
    test = sp_geom.Geometry.from_3dm(filepath=section_3dm)
    assert (test.geom - exp).is_empty


@pytest.mark.skipif(sys.version_info >= (3, 10), reason="requires python 3.9")
def test_geometry_from_3dm_file_compound():
    """Tests loading compound geometry from a .3dm file."""
    section_3dm = Path(__file__).parent.absolute() / "compound_shape.3dm"
    section_wkt = Path(__file__).parent.absolute() / "compound_shape.txt"
    with open(section_wkt) as file:
        wkt_str = file.readlines()
    exp = [wkt.loads(wkt_str[0]), wkt.loads(wkt_str[1])]
    test = sp_geom.CompoundGeometry.from_3dm(filepath=section_3dm)
    assert (MultiPolygon([ii.geom for ii in test.geoms]) - MultiPolygon(exp)).is_empty


@pytest.mark.skipif(sys.version_info >= (3, 10), reason="requires python 3.9")
def test_geometry_from_3dm_encode():
    """Tests loading compound geometry from a .json file."""
    section_3dm = Path(__file__).parent.absolute() / "rhino_data.json"
    with open(section_3dm) as file:
        brep_encoded = json.load(file)
    exp = Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])
    test = sp_geom.Geometry.from_rhino_encoding(r3dm_brep=brep_encoded)
    assert (test.geom - exp).is_empty


def test_shift_points():
    """Tests the shift_points() method."""
    assymetrical_chan = nastran_sections.nastran_chan(
        dim_1=75, dim_2=200, dim_3=8, dim_4=16
    ).shift_points(point_idxs=1, dy=-10)
    assert (
        assymetrical_chan.geom.wkt
        == "POLYGON ((0 0, 75 -10, 75 16, 8 16, 8 184, 75 184, 75 200, 0 200, 0 0))"
    )


def test_mirror_section():
    """Tests the mirror_section() method."""
    assymetrical_chan = nastran_sections.nastran_chan(
        dim_1=75, dim_2=200, dim_3=8, dim_4=16
    ).shift_points(point_idxs=1, dy=-10)
    assert (
        assymetrical_chan.mirror_section(axis="x").geom.wkt
        == "POLYGON ((0 190, 75 200, 75 174, 8 174, 8 6, 75 6, 75 -10, 0 -10, 0 190))"
    )
    assert (
        assymetrical_chan.mirror_section(axis="y").geom.wkt
        == "POLYGON ((75 0, 0 -10, 0 16, 67 16, 67 184, 0 184, 0 200, 75 200, 75 0))"
    )
    p = "POLYGON ((100 0, 25 -10, 25 16, 92 16, 92 184, 25 184, 25 200, 100 200, "
    p += "100 0))"
    assert (
        assymetrical_chan.mirror_section(axis="y", mirror_point=(50, 50)).geom.wkt == p
    )
    p = "POLYGON ((0 100, 75 110, 75 84, 8 84, 8 -84, 75 -84, 75 -100, 0 -100, 0 100))"
    assert (
        assymetrical_chan.mirror_section(axis="x", mirror_point=(50, 50)).geom.wkt == p
    )


def test_filter_non_polygons():
    """Tests the filter_non_polygons() method."""
    point1 = (0, 0)
    point2 = (1, 1)
    point3 = (1, 0)
    line = LineString([point1, point2])
    poly = Polygon([point1, point2, point3])
    multi_poly = MultiPolygon([poly, poly])
    collection = GeometryCollection([poly, Point(point1), line])
    assert sp_geom.filter_non_polygons(input_geom=poly) == poly
    assert sp_geom.filter_non_polygons(input_geom=multi_poly) == multi_poly
    assert sp_geom.filter_non_polygons(input_geom=Point(point1)) == Polygon()
    assert sp_geom.filter_non_polygons(input_geom=line) == Polygon()
    assert sp_geom.filter_non_polygons(input_geom=collection) == poly


def test_round_polygon_vertices():
    """Tests the round_polygon_vertices() method."""
    big_box = box(0, 0, 200, 200)
    bottom_box = box(10.00001, 10.000001, 50.100, 50.2)
    upper_box = box(120.000011, 120.000032, 169.999987, 170.0001)
    test_shape = big_box - bottom_box - upper_box
    poly = "POLYGON ((0 200, 200 200, 200 0, 0 0, 0 200), (10 50, 10 10, 50 10, 50 50, "
    poly += "10 50), (170 170, 120 170, 120 120, 170 120, 170 170))"
    assert test_shape.wkt != poly
    test_shape_rounded = sp_geom.round_polygon_vertices(poly=test_shape, tol=0)
    assert test_shape_rounded.wkt == poly


def test_check_geometry_overlaps():
    """Tests the geometry overlap checker."""
    big_sq = primitive_sections.rectangular_section(d=300, b=250)
    small_sq = primitive_sections.rectangular_section(d=100, b=75)
    small_hole = primitive_sections.rectangular_section(d=40, b=30).align_center(
        align_to=small_sq
    )

    rect = primitive_sections.rectangular_section(d=50, b=50)
    circ = primitive_sections.circular_section(d=50, n=32).shift_section(
        x_offset=125, y_offset=25
    )

    assert sp_geom.check_geometry_overlaps(lop=[small_sq.geom, small_hole.geom])
    assert sp_geom.check_geometry_overlaps(lop=[small_sq.geom, small_sq.geom])
    assert not sp_geom.check_geometry_overlaps(
        lop=[big_sq.geom, small_sq.shift_section(x_offset=270).geom]
    )
    assert sp_geom.check_geometry_overlaps(
        lop=[big_sq.geom, small_sq.shift_section(x_offset=200, y_offset=150).geom]
    )

    assert not sp_geom.check_geometry_overlaps(lop=[rect.geom, circ.geom])


def test_check_geometry_disjoint():
    """Tests the geometry disjoint checker."""
    rect = primitive_sections.rectangular_section(d=50, b=50)
    circ = primitive_sections.circular_section(d=50, n=32).shift_section(
        x_offset=125, y_offset=25
    )

    rect2 = primitive_sections.rectangular_section(d=50, b=50).shift_section(
        x_offset=50
    )
    assert not sp_geom.check_geometry_disjoint(lop=[rect.geom, rect2.geom])

    rect3 = primitive_sections.rectangular_section(d=25, b=25).shift_section(
        y_offset=50
    )
    assert not sp_geom.check_geometry_disjoint(lop=[rect.geom, rect2.geom, rect3.geom])

    assert sp_geom.check_geometry_disjoint(lop=[rect.geom, circ.geom])

    rect2 = primitive_sections.rectangular_section(d=50, b=50).shift_section(
        x_offset=50
    )
    assert not sp_geom.check_geometry_disjoint(lop=[rect.geom, rect2.geom])


def test_warping_disjoint_warning():
    """Tests that the warning occurs when analysing (warping) a disjoint geometry."""
    rect = primitive_sections.rectangular_section(d=50, b=50)
    circ = primitive_sections.circular_section(d=50, n=32).shift_section(
        x_offset=125, y_offset=25
    )
    geom = (rect + circ).create_mesh(mesh_sizes=[0])

    sec = Section(geometry=geom)
    sec.calculate_geometric_properties()
    with pytest.warns(UserWarning):
        sec.calculate_warping_properties()


def test_align_center():
    """Tests the align_center() method."""
    rect = primitive_sections.rectangular_section(d=200, b=70)
    circ = primitive_sections.circular_section(d=200, n=30)
    rect = rect.rotate_section(angle=-45, rot_point=(0, 0))
    rect_point = rect.points[1]
    circ = circ.align_center(align_to=rect_point)
    circ_x, circ_y = circ.calculate_centroid()
    assert pytest.approx(circ_x) == 49.497474683057995
    assert pytest.approx(circ_y) == -49.49747468305799

    circ = circ.align_center()
    circ_x, circ_y = circ.calculate_centroid()
    assert pytest.approx(circ_x) == 0
    assert pytest.approx(circ_y) == 0

    circ = circ.align_center(align_to=rect)
    circ_x, circ_y = circ.calculate_centroid()
    assert pytest.approx(circ_x) == 95.45941546018399
    assert pytest.approx(circ_y) == 45.961940777125974
