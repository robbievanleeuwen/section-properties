import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection

afloat = st.floats(min_value=0)

@given(
    d=st.floats(min_value=0.001),
    b=st.floats(min_value=0.001),
    shift=st.lists(afloat, min_size=2, max_size=2)
)
def test_rectangle_geom(d, b, shift):
    '''
    This function tests the rectangular geometry.
    The pre module is currently undergoing an overhaul, so this is really just
    serving as the simplest point of entry for me to learn pytest/hypothesis.
    '''
    xsect = sections.RectangularSection(d, b, shift=shift)

    assert len(xsect.points) == 4
    assert len(xsect.facets) == 4
    assert len(xsect.control_points) == 1
    assert xsect.control_points[0][0] == b/2 + shift[0]
    assert xsect.control_points[0][1] == d/2 + shift[1]



@given(
    d=st.floats(min_value=0.1, max_value=999),
    b=st.floats(min_value=0.1, max_value=999)
)
def test_rectangle_mesh(d, b):
    '''
    This function tests the rectangular mesh.
    The pre module is currently undergoing an overhaul, so this is really just
    serving as the simplest point of entry for me to learn pytest/hypothesis.
    '''
    xsect = sections.RectangularSection(d, b)
    mesh = xsect.create_mesh([min(d,b)/50])


@given(
    d=st.floats(min_value=.01, max_value=10),
    b=st.floats(min_value=.01, max_value=10)
)
def test_rectangle_mesh_size(d, b):
    '''
    This function verifies that the ValueError is raised when mesh size is
    set too small.
    '''
    with pytest.raises(ValueError):
        xsect = sections.RectangularSection(d, b)
        mesh = xsect.create_mesh([0.1/251])



@given(
    geom = st.builds(
        sections.RectangularSection,
        d = st.floats(min_value=0.1, max_value=10),
        b = st.floats(min_value=0.1, max_value=10)
    )
)
@settings(deadline=1000)
def test_cross_section(geom):
    '''
    This function tests a basic CrossSection object, instead of 
    just Geometry and Mesh from the prior funtions. It builds a 
    RectangularSection Geometry, then passes that to the
    CrossSection Class for analysis.
    '''
    b = geom.points[2][0]
    d = geom.points[2][1]
    mesh = geom.create_mesh([min(b,d)/5])

    xsect = CrossSection(geom, mesh)
    assert xsect.geometry == geom
    assert xsect.num_nodes == len(mesh.points)
    assert len(xsect.elements) == len(mesh.elements)

    xsect.calculate_geometric_properties()
    xsect.calculate_frame_properties()
    xsect.calculate_warping_properties()
    xsect.calculate_geometric_properties()