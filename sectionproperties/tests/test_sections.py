import pathlib

from sectionproperties.pre.sections import *
from sectionproperties.analysis.cross_section import Section
from sectionproperties.pre.pre import Material

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
    steel = Material("steel", 200e3, 0.3, 400, "grey")
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
    holes = [[0, 6], [0, -6]]
    new_geom = Geometry.from_points(points, facets, holes, control_points=[])
    assert (
        new_geom.geom.wkt
        == 'POLYGON ((6 10, 6 -10, -6 -10, -6 10, 6 10), (-4 4, 4 4, 4 8, -4 8, -4 4), (4 -8, 4 -4, -4 -4, -4 -8, 4 -8))'
    )


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
    holes = []
    control_points = [[0, 0], [0, -2 * a - t / 2]]
    new_geom = CompoundGeometry.from_points(points, facets, holes, control_points)
    assert (
        new_geom.geom.wkt
        == "MULTIPOLYGON (((-0.05 -2, 0.05 -2, 0.05 -0.05, 1 -0.05, 1 0.05, -0.05 0.05, -0.05 -2)), ((-1 -2, 1 -2, 1 -2.1, -1 -2.1, -1 -2)))"
    )


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
