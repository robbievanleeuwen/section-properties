from sectionproperties.pre.sections import *
from sectionproperties.analysis.cross_section import Section
from sectionproperties.pre.pre import Material

big_sq = rectangular_section(d=300, b= 250)
small_sq = rectangular_section(d=100, b=75)
small_hole = rectangular_section(d=40, b=30)
i_sec = i_section(d=200, b= 100, t_f=20, t_w=10, r=12, n_r=12)

small_sq = small_sq - small_hole.align_center(small_sq)
composite = big_sq + small_sq.align_top(big_sq, inner=True).align_right(big_sq) + i_sec.align_bottom(big_sq, inner=True).align_right(big_sq)
composite.compile_geometry()
composite.create_mesh(200)
comp_sec = Section(composite)
comp_sec.calculate_geometric_properties()
comp_sec.calculate_plastic_properties()

def test_material_persistence():
    # Test ensures that the material attribute gets transformed
    # through all of the Geometry transformation methods, each which
    # returns a new Geometry object.
    # The material assignment should persist through all of the
    # transformations
    steel = Material("steel", 200e3, 0.3, 400)
    big_sq.material = steel
    new_geom = (big_sq
        .align_left(small_sq)
        .align_right(small_sq)
        .align_top(small_hole)
        .align_bottom(small_sq)
        .align_center()
        .rotate_section(23)
        .mirror_section(axis='y')
        .offset_section_perimeter(amount=1)
    )
    new_geom.material == steel

def test_for_incidental_holes():
    # One hole in the geometry was explicitly created through subtraction
    # Another hole in the geometry was created accidentally by sticking
    # a I-section up against a rectangle.
    # There should be two holes created after .compile_geometry()
    assert len(composite.holes) == 2

def test_geometry_from_points():
    exterior = [[-6, 10], [6, 10], [6, -10], [-6, -10]]
    interior1 = [[-4, 8], [4, 8], [4, 4], [-4, 4]]
    interior2 = [[-4, -8], [4, -8], [4, -4], [-4, -4]]
    points = exterior + interior1 + interior2
    facets = [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7,4], [8,9], [9,10], [10, 11], [11,7]]
    holes = [[0, 6], [0, -6]]
    new_geom = Geometry.from_points(points, facets, holes, control_points=[])
    assert new_geom.geom.wkt == 'POLYGON ((-6 10, 6 10, 6 -10, -6 -10, -6 10), (-4 8, -4 4, 4 4, 4 8, -4 8), (-4 -8, 4 -8, 4 -4, -4 -4, -4 -8))'
    
 
