"""Tests for various geometry methods and classes."""

from __future__ import annotations

import json
import platform
from pathlib import Path

import pytest
import pytest_check as check
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
from sectionproperties.pre import CompoundGeometry, Geometry
from sectionproperties.pre.pre import Material


@pytest.fixture
def big_square() -> Geometry:
    """Generates a 300 x 250 rectangular geometry object.

    Returns:
        Geometry
    """
    return primitive_sections.rectangular_section(d=300, b=250)


@pytest.fixture
def small_square() -> Geometry:
    """Generates a 100 x 75 rectangular geometry object.

    Returns:
        Geometry
    """
    return primitive_sections.rectangular_section(d=100, b=75)


@pytest.fixture
def unit_square() -> Geometry:
    """Generates a unit square geometry object.

    Returns:
        Geometry
    """
    return Geometry(Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]))


@pytest.fixture
def unit_square_compound() -> CompoundGeometry:
    """Generates a unit square geometry object.

    Returns:
        Geometry
    """
    left = Geometry(Polygon([(0, 0), (0.5, 0), (0.5, 1), (0, 1)]))
    right = Geometry(Polygon([(0.5, 0), (1, 0), (1, 1), (0.5, 1)]))

    return left + right


@pytest.fixture
def small_hole(small_square: Geometry) -> Geometry:
    """Generates a 100 x 75 rectangular geometry object.

    Args:
        small_square: Small square fixture

    Returns:
        Geometry
    """
    small_sq = small_square

    return primitive_sections.rectangular_section(d=40, b=30).align_center(
        align_to=small_sq,
    )


@pytest.fixture
def small_square_with_hole(small_square: Geometry, small_hole: Geometry) -> Geometry:
    """Generates a 100 x 75 rectangular geometry object with 40 x 30 hole.

    Args:
        small_square: Small square fixture
        small_hole: Small hole fixture

    Returns:
        Geometry
    """
    small_sq = small_square
    small_hl = small_hole

    return small_sq - small_hl


@pytest.fixture
def i_section_geom() -> Geometry:
    """Generates a 200 deep i-section geometry object.

    Returns:
        Geometry
    """
    return steel_sections.i_section(d=200, b=100, t_f=20, t_w=10, r=12, n_r=12)


@pytest.fixture
def composite_geom(
    big_square: Geometry,
    small_square_with_hole: Geometry,
    i_section_geom: Geometry,
) -> Geometry:
    """Generates a composite geometry object.

    Args:
        big_square: Big square fixture
        small_square_with_hole: Small square with hole fixture
        i_section_geom: I-section fixture

    Returns:
        Geometry
    """
    big_sq = big_square
    small_sq_w_hole = small_square_with_hole
    i_sec = i_section_geom

    composite = (
        big_sq
        + small_sq_w_hole.align_to(other=big_sq, on="top", inner=True).align_to(
            other=big_sq,
            on="top",
        )
        + i_sec.align_to(other=big_sq, on="bottom", inner=True).align_to(
            other=big_sq,
            on="right",
        )
    )

    return composite


@pytest.fixture
def nested_geometry(small_square: Geometry, small_hole: Geometry) -> Geometry:
    """Generates a nested geometry object.

    Args:
        small_square: Small square fixture
        small_hole: Small hole fixture

    Returns:
        Geometry
    """
    small_sq = small_square
    small_hl = small_hole

    return (small_sq - small_hl) + small_hl


@pytest.fixture
def steel_material() -> Material:
    """Generates a steel material.

    Returns:
        Material
    """
    return Material(
        name="steel",
        elastic_modulus=200e3,
        poissons_ratio=0.3,
        density=7.85e-6,
        yield_strength=400,
        color="grey",
    )


def test_material_persistence(
    big_square: Geometry,
    small_square: Geometry,
    steel_material: Material,
):
    """Tests material persistence.

    Test ensures that the material attribute gets transformed through all of the
    Geometry transformation methods, each which returns a new Geometry object. The
    material assignment should persist through all of the transformations.
    """
    big_sq = big_square
    small_sq = small_square
    steel = steel_material

    big_sq.material = steel
    new_geom = (
        big_sq.align_to(other=small_sq, on="left", inner=False)
        .align_center()
        .rotate_section(angle=23)
        .mirror_section(axis="y")
        .offset_perimeter(amount=1)
    )
    assert new_geom.material == steel


def test_for_incidental_holes(composite_geom: Geometry, nested_geometry: Geometry):
    """Tests for incidental hole creation.

    One hole in the geometry was explicitly created through subtraction. Another hole in
    the geometry was created accidentally by sticking a I Section up against a
    rectangle. There should be two holes created after .compile_geometry()
    """
    composite = composite_geom
    nested_geom = nested_geometry

    assert len(composite.holes) == 2
    assert len(nested_geom.holes) == 0


def test__sub__(steel_material: Material, big_square: Geometry, small_hole: Geometry):
    """Tests the geometry subtraction method."""
    steel = steel_material
    big_sq = big_square
    small_hl = small_hole

    small_hl.material = steel
    top_left = (
        small_hl.align_to(other=big_sq, on="left")
        .align_to(other=big_sq, on="top")
        .shift_section(x_offset=20, y_offset=-20)
    )
    top_right = top_left.shift_section(x_offset=200)

    compound = big_sq - top_left
    compound = compound + top_left
    compound = compound - top_right
    compound = compound + top_right
    compound = compound - small_hl

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
        points=points,
        facets=facets,
        control_points=control_points,
        materials=materials,
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


def test_geometry_from_3dm_file_simple():
    """Tests loading geometry from a simple .3dm file."""
    section = Path(__file__).parent.absolute() / "3in x 2in.3dm"
    exp = Polygon([(0, 0), (0, 3), (2, 3), (2, 0), (0, 0)])
    test = sp_geom.Geometry.from_3dm(filepath=section)
    assert (test.geom - exp).is_empty


def test_geometry_from_3dm_file_complex():
    """Tests loading geometry from a complex .3dm file."""
    section_3dm = Path(__file__).parent.absolute() / "complex_shape.3dm"
    section_wkt = Path(__file__).parent.absolute() / "complex_shape.txt"
    with open(section_wkt) as file:
        wkt_str = file.readlines()
    exp = wkt.loads(wkt_str[0])
    test = sp_geom.Geometry.from_3dm(filepath=section_3dm)
    assert (test.geom - exp).is_empty


def test_geometry_from_3dm_file_compound():
    """Tests loading compound geometry from a .3dm file."""
    section_3dm = Path(__file__).parent.absolute() / "compound_shape.3dm"
    section_wkt = Path(__file__).parent.absolute() / "compound_shape.txt"
    with open(section_wkt) as file:
        wkt_str = file.readlines()
    exp = [wkt.loads(wkt_str[0]), wkt.loads(wkt_str[1])]
    test = sp_geom.CompoundGeometry.from_3dm(filepath=section_3dm)
    assert (MultiPolygon([ii.geom for ii in test.geoms]) - MultiPolygon(exp)).is_empty


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
        dim_1=75,
        dim_2=200,
        dim_3=8,
        dim_4=16,
    ).shift_points(point_idxs=1, dy=-10)
    assert (
        assymetrical_chan.geom.wkt
        == "POLYGON ((0 0, 75 -10, 75 16, 8 16, 8 184, 75 184, 75 200, 0 200, 0 0))"
    )


def test_mirror_section():
    """Tests the mirror_section() method."""
    assymetrical_chan = nastran_sections.nastran_chan(
        dim_1=75,
        dim_2=200,
        dim_3=8,
        dim_4=16,
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
        align_to=small_sq,
    )

    rect = primitive_sections.rectangular_section(d=50, b=50)
    circ = primitive_sections.circular_section(d=50, n=32).shift_section(
        x_offset=125,
        y_offset=25,
    )

    assert sp_geom.check_geometry_overlaps(lop=[small_sq.geom, small_hole.geom])
    assert sp_geom.check_geometry_overlaps(lop=[small_sq.geom, small_sq.geom])
    assert not sp_geom.check_geometry_overlaps(
        lop=[big_sq.geom, small_sq.shift_section(x_offset=270).geom],
    )
    assert sp_geom.check_geometry_overlaps(
        lop=[big_sq.geom, small_sq.shift_section(x_offset=200, y_offset=150).geom],
    )

    assert not sp_geom.check_geometry_overlaps(lop=[rect.geom, circ.geom])


def test_check_geometry_disjoint():
    """Tests the geometry disjoint checker."""
    rect = primitive_sections.rectangular_section(d=50, b=50)
    circ = primitive_sections.circular_section(d=50, n=32).shift_section(
        x_offset=125,
        y_offset=25,
    )

    rect2 = primitive_sections.rectangular_section(d=50, b=50).shift_section(
        x_offset=50,
    )
    assert not sp_geom.check_geometry_disjoint(lop=[rect.geom, rect2.geom])

    rect3 = primitive_sections.rectangular_section(d=25, b=25).shift_section(
        y_offset=50,
    )
    assert not sp_geom.check_geometry_disjoint(lop=[rect.geom, rect2.geom, rect3.geom])

    assert sp_geom.check_geometry_disjoint(lop=[rect.geom, circ.geom])

    rect2 = primitive_sections.rectangular_section(d=50, b=50).shift_section(
        x_offset=50,
    )
    assert not sp_geom.check_geometry_disjoint(lop=[rect.geom, rect2.geom])


def test_warping_disjoint_warning():
    """Tests that the warning occurs when analysing (warping) a disjoint geometry."""
    rect = primitive_sections.rectangular_section(d=50, b=50)
    circ = primitive_sections.circular_section(d=50, n=32).shift_section(
        x_offset=125,
        y_offset=25,
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


def test_geom_obj_value_error():
    """Tests the ValueErrors when creating a Geometry object."""
    # create some points lists
    pts = [(0, 0), (1, 0), (1, 1), (0, 1)]
    pts_2 = [(0, 1), (1, 1), (1, 2), (0, 2)]
    multi_poly = MultiPolygon([Polygon(pts), Polygon(pts_2)])

    with pytest.raises(ValueError, match="Use CompoundGeometry"):
        Geometry(geom=multi_poly)

    with pytest.raises(ValueError, match="Argument is not a valid shapely.Polygon"):
        Geometry(geom=pts)


def test_repr_svg(unit_square, unit_square_compound, capfd: pytest.CaptureFixture[str]):
    """Tests the rep_svg() method."""
    geom = unit_square
    geom._repr_svg_()

    out, _ = capfd.readouterr()
    assert "sectionproperties.pre.geometry.Geometry" in out
    assert "Material: default" in out

    geom = unit_square_compound
    geom._repr_svg_()

    out, _ = capfd.readouterr()
    assert "sectionproperties.pre.geometry.CompoundGeometry" in out
    assert "Materials incl.: ['default']" in out


def test_assign_control_point(unit_square):
    """Tests the assign_control_point() method."""
    geom = unit_square
    geom = geom.assign_control_point((0.25, 0.25))

    assert geom.control_points == [(0.25, 0.25)]


def test_from_points_value_error():
    """Tests the from_points() ValueError."""
    pts = [(0, 0), (1, 0), (1, 1), (0, 1)]
    fcts = [(0, 1), (1, 2), (2, 3), (3, 0)]
    ctrl_pts = [(0.24, 0.24), (0.56, 0.76)]

    with pytest.raises(ValueError, match="Control points for Geometry instances "):
        Geometry.from_points(pts, fcts, ctrl_pts)


def test_create_mesh_value_error(unit_square):
    """Tests the create_mesh() ValueError."""
    geom = unit_square

    with pytest.raises(ValueError, match="for a Geometry must be either"):
        geom.create_mesh([1, 2])


def test_align_to_point(unit_square):
    """Tests the align_to() method with a point."""
    geom = unit_square
    geom = geom.align_to(other=(5, 0), on="right")

    assert geom.points == [(5, 0), (6, 0), (6, 1), (5, 1)]


def test_align_centre_value_error(unit_square):
    """Tests the align_centre() ValueError."""
    geom = unit_square

    with pytest.raises(ValueError, match="align_to must be either a Geometry object"):
        geom.align_center(align_to=[1])


def test_rotate_control_point():
    """Tests the align_centre() ValueError."""
    geom = Geometry(
        Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
        control_points=(0.25, 0.25),
    )
    geom = geom.rotate_section(angle=180, rot_point=(0, 0))

    check.almost_equal(geom.control_points[0][0], -0.25)
    check.almost_equal(geom.control_points[0][1], -0.25)

    geom = Geometry(
        Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
        control_points=(0.25, 0.25),
    )
    geom = geom.rotate_section(angle=180, rot_point="center")

    check.almost_equal(geom.control_points[0][0], 0.75)
    check.almost_equal(geom.control_points[0][1], 0.75)


def test_split_section(unit_square):
    """Tests the split_section() method."""
    point_i = (0.5, 0.5)
    point_j = (1, 0.5)
    vector = (0.5, 1)

    top_geoms, bot_geoms = unit_square.split_section(point_i, point_j)

    assert top_geoms[0].geom.equals(Polygon([(0, 0.5), (1, 0.5), (1, 1), (0, 1)]))
    assert bot_geoms[0].geom.equals(Polygon([(0, 0), (1, 0), (1, 0.5), (0, 0.5)]))

    top_geoms, bot_geoms = unit_square.split_section(point_i, vector)

    assert top_geoms[0].geom.equals(Polygon([(0.5, 0), (1, 0), (1, 1), (0.5, 1)]))
    assert bot_geoms[0].geom.equals(Polygon([(0, 0), (0.5, 0), (0.5, 1), (0, 1)]))

    with pytest.raises(ValueError, match="Either a second point"):
        unit_square.split_section(point_i)


@pytest.mark.skipif(platform.system() != "Linux", reason="Only test plotting on Linux")
def test_plot_geometry(small_square_with_hole):
    """Tests plot geometry executes without errors."""
    geom = small_square_with_hole
    geom.plot_geometry(render=False)
    assert True


def test_recovery_points(unit_square):
    """Tests the recovery points."""
    geom = unit_square
    nastran = nastran_sections.nastran_bar(1, 1)

    assert len(geom.recovery_points) == 0
    assert nastran.recovery_points == [
        (0.5, 0.5),
        (0.5, -0.5),
        (-0.5, -0.5),
        (-0.5, 0.5),
    ]


def test_union(unit_square):
    """Tests the union operator."""
    geom_1 = unit_square
    geom_2 = unit_square.shift_section(0.5, 0.5)
    geom = geom_1 | geom_2

    assert geom.geom.equals(
        Polygon(
            [
                (0, 0),
                (1, 0),
                (1, 0.5),
                (1.5, 0.5),
                (1.5, 1.5),
                (0.5, 1.5),
                (0.5, 1),
                (0, 1),
            ],
        ),
    )


def test_xor(unit_square):
    """Tests the xor operator."""
    geom_1 = unit_square
    geom_2 = unit_square.shift_section(0.5, 0.5)
    geom = geom_1 ^ geom_2

    # check area is as expected
    geom.create_mesh(0)
    sec = Section(geom)
    sec.calculate_geometric_properties()
    check.almost_equal(sec.get_area(), 1.5)


def test_compound_from_points_errors():
    """Tests the error raised by CompoundGeometry.from_points()."""
    pts = [(0, 0), (0.5, 0), (0.5, 1), (0, 1), (0.5, 0), (1, 0), (1, 1), (0.5, 1)]
    fcts = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4)]
    mat1 = Material("a", 1, 0, 2, 1, "k")
    mat2 = Material("b", 2, 0.3, 5, 2, "r")

    with pytest.raises(ValueError, match="Materials cannot be assigned"):
        CompoundGeometry.from_points(
            points=pts,
            facets=fcts,
            control_points=None,
            materials=[mat1, mat2],
        )

    with pytest.raises(ValueError, match="the number of materials in the list"):
        CompoundGeometry.from_points(
            points=pts,
            facets=fcts,
            control_points=[(0.25, 0.25), (0.75, 0.75)],
            materials=[mat1],
        )

    with pytest.raises(ValueError, match="The number of exterior regions"):
        CompoundGeometry.from_points(
            points=pts,
            facets=fcts,
            control_points=[(0.25, 0.25)],
        )


def test_compound_mesh_as_float(unit_square_compound):
    """Tests mesh size given as a float."""
    geom1 = unit_square_compound
    geom2 = unit_square_compound
    geom1.create_mesh(0.1)
    geom2.create_mesh([0.1, 0.1])

    assert len(geom1.mesh["vertices"]) == len(geom2.mesh["vertices"])


def test_compound_rotate(unit_square_compound):
    """Tests the rotate_section method for CompoundGeometry."""
    geom = unit_square_compound
    geom_rot = geom.rotate_section(180)

    assert geom.geoms[0].geom.equals(geom_rot.geoms[1].geom)

    geom = unit_square_compound
    geom_rot = geom.rotate_section(180, (0, 0))

    assert geom_rot.geoms[0].geom.equals(
        Polygon([(0, 0), (-0.5, 0), (-0.5, -1), (0, -1)]),
    )


def test_compound_mirror(unit_square_compound):
    """Tests the mirror method for CompoundGeometry."""
    geom = unit_square_compound
    geom_mir = geom.mirror_section("y", (0, 0))

    assert geom_mir.geoms[0].geom.equals(
        Polygon([(0, 0), (-0.5, 0), (-0.5, 1), (0, 1)]),
    )


def test_compound_align_center(unit_square_compound):
    """Tests the align_center method for CompoundGeometry."""
    pts = [(1, 1), (2, 1), (2, 2), (1, 2)]
    geom_orig = Geometry(Polygon(pts))
    geom = unit_square_compound
    geom_shift = geom.align_center(geom_orig)

    assert geom_shift.geoms[0].geom.equals(
        Polygon([(1, 1), (1.5, 1), (1.5, 2), (1, 2)]),
    )

    geom = unit_square_compound
    geom_shift = geom.align_center((1.5, 1.5))
    assert geom_shift.geoms[0].geom.equals(
        Polygon([(1, 1), (1.5, 1), (1.5, 2), (1, 2)]),
    )

    with pytest.raises(ValueError, match="align_to must be either a Geometry object"):
        geom_shift = geom.align_center(1.5)
