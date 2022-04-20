"""
This file tests a couple of distinct examples  from
'Aircraft Structures,' by Peery. These cases have 
known results, and the output from SectionProperties 
is compared for accuracy. These examples represent a 
more rigourous 'proof' against a 'real' problem. 
Only results that have values in the reference material 
are tested here.

BibTeX Entry for reference:
@Book{Peery,
    title = {Aircraft Structures},
    author = {David J. Peery},
    organization = {Pensylvania State University},
    publisher = {McGraw-Hill Book Company},
    year = {1950},
    edition = {First},
    ISBN = {978-0486485805}
}
"""
import pytest
import pytest_check as check

from typing import Tuple

from sectionproperties.pre.library import nastran_sections
from sectionproperties.analysis.section import Section


## Classes
class Z_Section:
    """
    This is basically just a fixture for testing purposes.
    It's called by the actual pytest fixtures to generate
    the Z-sections for analysis.

    We have this class for fixtures, just to have
    a method for the load application, and simpler fixtures,
    along with the same base for multiple Z_Sections.
    """

    def __init__(self, DIM1, DIM2, DIM3, DIM4, shift, m, name):
        # Setup the analysis, and calculate properties
        base_geom = nastran_sections.nastran_zed(DIM1, DIM2, DIM3, DIM4)
        self.geom = base_geom.shift_section(*shift)
        self.geom = self.geom.create_mesh(mesh_sizes=[m])
        self.xsect = Section(self.geom)
        self.xsect.calculate_geometric_properties()
        # This plotting code was just for verifying the section offsets.
        # ax = self.xsect.plot_centroids(pause=False, render=False)
        # ax.grid(1, which='both', linestyle=':')
        # fig = ax.get_figure()
        # fig.savefig(f'{name}_geom.png')

    def apply_load(self, v):
        """
        This method applies the suplied load to the section.
        v is a list-like with the first entry being Mxx, and
        second entry Myy.
        """
        self.xsect.calculate_warping_properties()
        self.stress = self.xsect.calculate_stress(Mxx=v[0], Myy=v[1])


## Utility
def get_node(nodes, coord) -> Tuple[int, tuple]:
    """
    This function will loop over the node list provided,
    finding the index of the coordinates you want.
    Returns the index in the nodes list, and the coords.
    """
    for index, var in enumerate(nodes):
        if all(var == coord):
            return index, var
        else:
            continue

    raise ValueError(f"No node found with coordinates: {coord}")


## Fixtures
@pytest.fixture
def PeeryEx6_2_1():
    """
    Example 1 in Sec. 6.2 (Symmetric Bending)
    This is a symmetric I-section with no lateral supports,
    undergoing pure unidirectional cantilever bending.
    Note that units here are **inches**, to match the text.
    """
    name = "Peery_6.2.1"
    geom = nastran_sections.nastran_i(6, 3, 3, 1, 1, 1)
    geom = geom.shift_section(0, -3)
    geom = geom.create_mesh([0.25])
    xsect = Section(geom)
    xsect.calculate_geometric_properties()
    # This plotting code was just for verifying the section offsets.
    # ax = xsect.plot_centroids(pause=False, render=False)
    # ax.grid(1, which='both', linestyle=':')
    # fig = ax.get_figure()
    # fig.savefig(f'{name}_geom.png')

    return geom, xsect


@pytest.fixture
def PeeryEx7_2_1():
    """
    Example 1 in Sec. 7.2. (Unsymmetric Bending)
    This is an unsymmetric Z-section with no lateral supports.
    Note that units here are **inches**, to match the text.
    """
    return Z_Section(
        DIM1=4, DIM2=2, DIM3=8, DIM4=12, shift=[-5, -6], m=0.25, name="Peery_7.2.1"
    )


## Tests
def test_symmetric_ixx(PeeryEx6_2_1):
    # Directly from the example, we know that
    # the 2nd moment of inertia resisting bending is.
    _geom, xsect = PeeryEx6_2_1

    check.almost_equal(xsect.section_props.ixx_g, 43.3, rel=1e-3)


def test_symmetric_fb(PeeryEx6_2_1):
    "Max bending stress on the section."
    _geom, xsect = PeeryEx6_2_1
    # Defined in the text
    moment = 8e5
    y = 3
    I = xsect.section_props.ixx_g
    xsect.calculate_warping_properties()
    stress = xsect.calculate_stress(Mxx=moment)

    # The number quoted in the book. (Peery rounds this to the hundreds)
    # 55400 = 55427.3
    perfect_result = 55427.3

    # The number from the textbook equation
    computed_result = moment * y / I
    check.almost_equal(perfect_result, computed_result, rel=1e-3)

    # The max stress, computed through FEA on our mesh.
    numerical_result = max(stress.get_stress()[0]["sig_zz"])
    check.almost_equal(numerical_result, perfect_result, rel=1e-3)


def test_unsymmetric_ixx(PeeryEx7_2_1):
    # Directly from the example, we know what
    # the section properties should be.
    xsect = PeeryEx7_2_1.xsect
    check.almost_equal(xsect.section_props.ixx_g, 693.3, rel=1e-3)


def test_unsymmetric_iyy(PeeryEx7_2_1):
    # Directly from the example, we know what
    # the section properties should be.
    xsect = PeeryEx7_2_1.xsect
    check.almost_equal(xsect.section_props.iyy_g, 173.3, rel=1e-3)


def test_unsymmetric_ixy(PeeryEx7_2_1):
    # Directly from the example, we know what
    # the section properties should be.
    xsect = PeeryEx7_2_1.xsect
    check.almost_equal(xsect.section_props.ixy_g, -240, rel=1e-3)


def test_unsymmetric_i11(PeeryEx7_2_1):
    # Directly from the example, we know what
    # the section properties should be.
    xsect = PeeryEx7_2_1.xsect
    check.almost_equal(xsect.section_props.i11_c, 787, rel=1e-3)


def test_unsymmetric_i22(PeeryEx7_2_1):
    # Directly from the example, we know what
    # the section properties should be.
    xsect = PeeryEx7_2_1.xsect
    check.almost_equal(xsect.section_props.i22_c, 79.5, rel=1e-3)


def test_fb_C(PeeryEx7_2_1):
    """Check the stress at point C."""
    # Load from the text
    v = [-1e5, 1e4]
    # Coordinates of point C
    C = (1, 6)
    # The answer in the example
    # For this point, Peery rounds to the tens place,
    # thus -2380 is the exact number written in the book
    # but -2384 is the answer computed from his values.
    perfect_result = -2384
    # The simplified textbook equation
    text_result = round(-494 * 1 + -315 * 6)
    nodes = PeeryEx7_2_1.xsect.mesh_nodes
    assert len(nodes > 0)
    index, _ = get_node(nodes, C)
    _ = PeeryEx7_2_1.apply_load(v)
    computed_result = PeeryEx7_2_1.stress.get_stress()[0]["sig_zz"][index]

    check.almost_equal(text_result, perfect_result)
    check.almost_equal(computed_result, perfect_result, rel=1e-3)


def test_fb_B(PeeryEx7_2_1):
    """Check the stress at point B."""
    # Load from the text
    v = [-1e5, 1e4]
    # Coordinates of point B
    B = (-5, 6)
    # The answer in the example
    perfect_result = 580
    # The sipmlified textbook equation
    text_result = round(-494 * -5 + -315 * 6)
    nodes = PeeryEx7_2_1.xsect.mesh_nodes
    index, _ = get_node(nodes, B)
    _ = PeeryEx7_2_1.apply_load(v)
    computed_result = PeeryEx7_2_1.stress.get_stress()[0]["sig_zz"][index]

    check.almost_equal(text_result, perfect_result)
    check.almost_equal(computed_result, perfect_result, rel=1e-3)


def test_fb_A(PeeryEx7_2_1):
    """Check the stress at point A."""
    # Load from the text
    v = [-1e5, 1e4]
    # Coordinates of point A
    A = (-5, 4)
    # The answer in the example
    perfect_result = 1210
    # The simplified textbook equation
    text_result = round(-494 * -5 + -315 * 4)
    nodes = PeeryEx7_2_1.xsect.mesh_nodes
    index, _ = get_node(nodes, A)
    _ = PeeryEx7_2_1.apply_load(v)
    computed_result = PeeryEx7_2_1.stress.get_stress()[0]["sig_zz"][index]

    check.almost_equal(text_result, perfect_result)
    check.almost_equal(computed_result, perfect_result, rel=1e-3)
