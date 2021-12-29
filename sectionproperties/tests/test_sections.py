import pathlib
import pytest

from sectionproperties.pre.sections import *
from sectionproperties.pre.nastran_sections import *
from sectionproperties.analysis.cross_section import Section
from sectionproperties.pre.pre import Material
from sectionproperties.pre.rhino import load_3dm, load_brep_encoding
from shapely.geometry import (
    Polygon,
    MultiPolygon,
    LineString,
    Point,
    GeometryCollection,
)
from shapely import wkt
import json

big_sq = rectangular_section(d=300, b=250)
small_sq = rectangular_section(d=100, b=75)
small_hole = rectangular_section(d=40, b=30)
i_sec = i_section(d=200, b=100, t_f=20, t_w=10, r=12, n_r=12)

small_sq = small_sq - small_hole.align_center(small_sq)
composite = (
    big_sq
    + small_sq.align_to(big_sq, on="top", inner=True).align_to(big_sq, on="top")
    + i_sec.align_to(big_sq, on="bottom", inner=True).align_to(big_sq, on="right")
)
composite.compile_geometry()
composite.create_mesh([200])
comp_sec = Section(composite)
comp_sec.calculate_geometric_properties()
comp_sec.calculate_plastic_properties()


def test_material_persistence():
    # Test ensures that the material attribute gets transformed
    # through all of the Geometry transformation methods, each which
    # returns a new Geometry object.
    # The material assignment should persist through all of the
    # transformations
    steel = Material("steel", 200e3, 0.3, 7.85e-6, 400, "grey")
    big_sq.material = steel
    new_geom = (
        big_sq.align_to(small_sq, on="left", inner=False)
        .align_center()
        .rotate_section(23)
        .mirror_section(axis="y")
        .offset_perimeter(amount=1)
    )
    new_geom.material == steel


def test_for_incidental_holes():
    # One hole in the geometry was explicitly created through subtraction
    # Another hole in the geometry was created accidentally by sticking
    # a I-section up against a rectangle.
    # There should be two holes created after .compile_geometry()
    assert len(composite.holes) == 2


def test_geometry_from_points():
    # Geometry.from_points() tests a shape with exactly one exterior
    # and an arbitrary number of interiors being built from the legacy
    # points, facets, holes, control_points interface of sectionproperties
    exterior = [[-6, 10], [6, 10], [6, -10], [-6, -10]]
    interior1 = [[-4, 8], [4, 8], [4, 4], [-4, 4]]
    interior2 = [[-4, -8], [4, -8], [4, -4], [-4, -4]]
    points = exterior + interior1 + interior2
    facets = [
        [0, 1],
        [1, 2],
        [2, 3],
        [3, 0],
        [4, 5],
        [5, 6],
        [6, 7],
        [7, 4],
        [8, 9],
        [9, 10],
        [10, 11],
        [11, 7],
    ]
    control_points = [[0, 0]]
    holes = [[0, 6], [0, -6]]
    new_geom = Geometry.from_points(
        points=points, facets=facets, control_points=control_points, holes=holes
    )
    assert (
        new_geom.geom.wkt
        == "POLYGON ((6 10, 6 -10, -6 -10, -6 10, 6 10), (-4 4, 4 4, 4 8, -4 8, -4 4), (4 -8, 4 -4, -4 -4, -4 -8, 4 -8))"
    )  # Note, the order of point coordinates can change based on geometry operations, e.g. +/- and how many; seems to be a shapely issue


def test_compound_geometry_from_points():
    # CompoundGeometry.from_points() tests a shape with an arbitrary
    # number of exteriors and an arbitrary number of interiors being
    # built from the legacy
    # points, facets, holes, control_points interface of sectionproperties
    a = 1
    b = 2
    t = 0.1

    # build the lists of points, facets, holes and control points
    points = [
        [-t / 2, -2 * a],
        [t / 2, -2 * a],
        [t / 2, -t / 2],
        [a, -t / 2],
        [a, t / 2],
        [-t / 2, t / 2],
        [-b / 2, -2 * a],
        [b / 2, -2 * a],
        [b / 2, -2 * a - t],
        [-b / 2, -2 * a - t],
    ]
    facets = [
        [0, 1],
        [1, 2],
        [2, 3],
        [3, 4],
        [4, 5],
        [5, 0],
        [6, 7],
        [7, 8],
        [8, 9],
        [9, 6],
    ]
    control_points = [[0, 0], [0, -2 * a - t / 2]]
    new_geom = CompoundGeometry.from_points(points, facets, control_points)
    assert (
        new_geom.geom.wkt
        == "MULTIPOLYGON (((-0.05 -2, 0.05 -2, 0.05 -0.05, 1 -0.05, 1 0.05, -0.05 0.05, -0.05 -2)), ((-1 -2, 1 -2, 1 -2.1, -1 -2.1, -1 -2)))"
    )


def test_nested_compound_geometry_from_points():
    """
    Tests a nested compound geometry can be built .from_points, that the control_points
    and hole nodes persist in the right locations, and that ...
    """
    points = [
        [-50.0, 50.0],
        [50.0, 50.0],
        [50.0, -50.0],
        [-50.0, -50.0],
        [37.5, -37.5],
        [37.5, 37.5],
        [-37.5, 37.5],
        [-37.5, -37.5],
        [25.0, -25.0],
        [25.0, 25.0],
        [-25.0, 25.0],
        [-25.0, -25.0],
        [12.5, -12.5],
        [12.5, 12.5],
        [-12.5, 12.5],
        [-12.5, -12.5],
    ]
    facets = [
        [0, 1],
        [1, 2],
        [2, 3],
        [3, 0],
        [4, 5],
        [5, 6],
        [6, 7],
        [7, 4],
        [8, 9],
        [9, 10],
        [10, 11],
        [11, 8],
        [12, 13],
        [13, 14],
        [14, 15],
        [15, 12],
    ]
    control_points = [[-43.75, 0.0], [-31.25, 0.0], [-18.75, 0.0]]
    holes = [[0, 0]]
    nested_compound = CompoundGeometry.from_points(
        points=points, facets=facets, control_points=control_points, holes=holes
    )
    assert (
        nested_compound.geom.wkt
        == "MULTIPOLYGON (((50 50, 50 -50, -50 -50, -50 50, 50 50), (12.5 12.5, -12.5 12.5, -12.5 -12.5, 12.5 -12.5, 12.5 12.5)), ((-37.5 -37.5, -37.5 37.5, 37.5 37.5, 37.5 -37.5, -37.5 -37.5), (12.5 12.5, -12.5 12.5, -12.5 -12.5, 12.5 -12.5, 12.5 12.5)), ((-25 -25, -25 25, 25 25, 25 -25, -25 -25), (12.5 12.5, -12.5 12.5, -12.5 -12.5, 12.5 -12.5, 12.5 12.5)))"
    )
    assert nested_compound.control_points == [
        (-43.75, 0.0),
        (-31.25, 0.0),
        (-18.75, 0.0),
    ]
    assert nested_compound.holes == [(0, 0), (0, 0), (0, 0)]


def test_geometry_from_dxf():
    section_holes_dxf = (
        pathlib.Path.cwd() / "sectionproperties" / "tests" / "section_holes.dxf"
    )
    assert Geometry.from_dxf(section_holes_dxf).geom.wkt == (
        "POLYGON ((-0.3386588348890669 -0.3951777028951984, "
        "-0.3386588348890669 29.09231821639347, 31.96225758877617 29.09231821639347, "
        "31.96225758877617 -0.3951777028951984, -0.3386588348890669 -0.3951777028951984), "
        "(16.68431586247806 2.382629883704458, 29.68303085105337 2.382629883704458, "
        "29.68303085105337 24.35580015206329, 16.68431586247806 24.35580015206329, "
        "16.68431586247806 2.382629883704458), (1.548825807287628 3.34417866368126, "
        "14.54754079586294 3.34417866368126, 14.54754079586294 27.38289816310137, "
        "1.548825807287628 27.38289816310137, 1.548825807287628 3.34417866368126))"
    )


def test_plastic_centroid():
    ## Test created in response to #114
    # Since the section being tested is a compound geometry with two different
    # materials, this tests that the plastic centroid takes into account the
    # correct "center" of the original section which is affected by EA of each
    # of the constituent geometries.

    steel = Material(
        name="Steel",
        elastic_modulus=200e3,
        poissons_ratio=0.3,
        density=7.85e-6,
        yield_strength=500,
        color="grey",
    )
    timber = Material(
        name="Timber",
        elastic_modulus=5e3,
        poissons_ratio=0.35,
        density=6.5e-7,
        yield_strength=20,
        color="burlywood",
    )

    # create 310UB40.4
    ub = i_section(d=304, b=165, t_f=10.2, t_w=6.1, r=11.4, n_r=8, material=steel)

    # create timber panel on top of the UB
    panel = rectangular_section(d=50, b=600, material=timber)
    panel = panel.align_center(ub).align_to(ub, on="top")

    # merge the two sections into one geometry object
    geometry = CompoundGeometry([ub, panel])

    # create a mesh - use a mesh size of 5 for the UB, 20 for the panel
    geometry.create_mesh(mesh_sizes=[100, 100])

    # create a Section object
    section = Section(geometry)

    # perform a geometric, warping and plastic analysis
    section.calculate_geometric_properties()
    section.calculate_plastic_properties()

    x_pc, y_pc = section.get_pc()
    assert x_pc == pytest.approx(82.5)
    assert y_pc == pytest.approx(250.360654576)


def test_geometry_from_3dm_file_simple():
    section = pathlib.Path.cwd() / "sectionproperties" / "tests" / "3in x 2in.3dm"
    exp = Polygon([(0, 0), (0, 3), (2, 3), (2, 0), (0, 0)])
    test = Geometry.from_3dm(section)
    assert (test.geom - exp).is_empty


def test_geometry_from_3dm_file_complex():
    section_3dm = (
        pathlib.Path.cwd() / "sectionproperties" / "tests" / "complex_shape.3dm"
    )
    section_wkt = (
        pathlib.Path.cwd() / "sectionproperties" / "tests" / "complex_shape.txt"
    )
    with open(section_wkt) as file:
        wkt_str = file.readlines()
    exp = wkt.loads(wkt_str[0])
    test = Geometry.from_3dm(section_3dm)
    assert (test.geom - exp).is_empty


def test_geometry_from_3dm_file_compound():
    section_3dm = (
        pathlib.Path.cwd() / "sectionproperties" / "tests" / "compound_shape.3dm"
    )
    section_wkt = (
        pathlib.Path.cwd() / "sectionproperties" / "tests" / "compound_shape.txt"
    )
    with open(section_wkt) as file:
        wkt_str = file.readlines()
    exp = [wkt.loads(wkt_str[0]), wkt.loads(wkt_str[1])]
    test = CompoundGeometry.from_3dm(section_3dm)
    assert (MultiPolygon([ii.geom for ii in test.geoms]) - MultiPolygon(exp)).is_empty


def test_geometry_from_3dm_encode():
    section_3dm = pathlib.Path.cwd() / "sectionproperties" / "tests" / "rhino_data.json"
    with open(section_3dm) as file:
        brep_encoded = json.load(file)
    exp = Polygon([(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)])
    test = Geometry.from_rhino_encoding(brep_encoded)
    assert (test.geom - exp).is_empty


def test_shift_points():
    assymetrical_chan = nastran_chan(75, 200, 8, 16).shift_points(1, dy=-10)
    assert (
        assymetrical_chan.geom.wkt
        == "POLYGON ((0 0, 75 -10, 75 16, 8 16, 8 184, 75 184, 75 200, 0 200, 0 0))"
    )


def test_mirror_section():
    assymetrical_chan = nastran_chan(75, 200, 8, 16).shift_points(1, dy=-10)
    assert (
        assymetrical_chan.mirror_section(axis="x").geom.wkt
        == "POLYGON ((0 190, 75 200, 75 174, 8 174, 8 6, 75 6, 75 -10, 0 -10, 0 190))"
    )
    assert (
        assymetrical_chan.mirror_section(axis="y").geom.wkt
        == "POLYGON ((75 0, 0 -10, 0 16, 67 16, 67 184, 0 184, 0 200, 75 200, 75 0))"
    )
    assert (
        assymetrical_chan.mirror_section(axis="y", mirror_point=[50, 50]).geom.wkt
        == "POLYGON ((100 0, 25 -10, 25 16, 92 16, 92 184, 25 184, 25 200, 100 200, 100 0))"
    )
    assert (
        assymetrical_chan.mirror_section(axis="x", mirror_point=[50, 50]).geom.wkt
        == "POLYGON ((0 100, 75 110, 75 84, 8 84, 8 -84, 75 -84, 75 -100, 0 -100, 0 100))"
    )


def test_filter_non_polygons():
    point1 = Point([0, 0])
    point2 = Point([1, 1])
    point3 = Point([1, 0])
    line = LineString([point1, point2])
    poly = Polygon([point1, point2, point3])
    multi_poly = MultiPolygon([poly, poly])
    collection = GeometryCollection([poly, point1, line])
    out = filter_non_polygons(collection)
    assert filter_non_polygons(poly) == poly
    assert filter_non_polygons(multi_poly) == multi_poly
    assert filter_non_polygons(point1) == Polygon()
    assert filter_non_polygons(line) == Polygon()
    assert filter_non_polygons(collection) == poly
