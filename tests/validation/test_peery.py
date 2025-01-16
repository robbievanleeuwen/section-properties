"""Validation tests from Peery.

This file tests a couple of distinct examples  from 'Aircraft Structures,' by Peery.
These cases have known results, and the output from SectionProperties is compared for
accuracy. These examples represent a more rigourous 'proof' against a 'real' problem.
Only results that have values in the reference material are tested here.

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

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest
import pytest_check as check

from sectionproperties.analysis.section import Section
from sectionproperties.pre.library import nastran_sections

if TYPE_CHECKING:
    from collections.abc import Callable


class ZSection:
    """Class for a Z shaped section.

    This is basically just a fixture for testing purposes. It's called by the actual
    pytest fixtures to generate the Z-sections for analysis.

    We have this class for fixtures, just to havea method for the load application,
    and simpler fixtures, along with the same base for multiple ZSections.
    """

    def __init__(
        self,
        dim1: float,
        dim2: float,
        dim3: float,
        dim4: float,
        shift: tuple[float, float],
        m: float,
    ) -> None:
        """Inits the ZSection class.

        Args:
            dim1: Zed section dimension 1
            dim2: Zed section dimension 2
            dim3: Zed section dimension 3
            dim4: Zed section dimension 4
            shift: Coordinates to shift section
            m: Mesh size
        """
        # Setup the analysis, and calculate properties
        base_geom = nastran_sections.nastran_zed(
            dim_1=dim1,
            dim_2=dim2,
            dim_3=dim3,
            dim_4=dim4,
        )
        self.geom = base_geom.shift_section(*shift)
        self.geom = self.geom.create_mesh(mesh_sizes=[m])
        self.xsect = Section(self.geom)
        self.xsect.calculate_geometric_properties()
        # This plotting code was just for verifying the section offsets.
        # ax = self.xsect.plot_centroids(pause=False, render=False)
        # ax.grid(1, which='both', linestyle=':')
        # fig = ax.get_figure()
        # fig.savefig(f'geom.png')

    def apply_load(
        self,
        v: tuple[float, float],
    ) -> None:
        """Applies a load to the section.

        This method applies the suplied load to the section. ``v`` is a tuple with the
        first entry being ``mxx``, and second entry ``myy``.

        Args:
            v: Load to apply
        """
        self.xsect.calculate_warping_properties()
        self.stress = self.xsect.calculate_stress(mxx=v[0], myy=v[1])


def get_node(
    nodes: list[tuple[float, float]],
    coord: tuple[float, float],
) -> tuple[int, tuple[float, float]]:
    """Finds a node given coordinates.

    This function will loop over the node list provided, finding the index of the
    coordinates you want. Returns the index in the nodes list, and the coords.

    Args:
        nodes: List of nodes
        coord: Coordinate to find

    Raises:
        RuntimeError: Node cannot be found

    Returns:
        Node index and coorindate
    """
    for index, var in enumerate(nodes):
        if var[0] == coord[0] and var[1] == coord[1]:
            return index, var
        else:
            continue

    msg = f"No node found with coordinates: {coord}"
    raise RuntimeError(msg)


@pytest.fixture
def peery_ex_6_2_1() -> Section:
    """Ex 6.2.1.

    Example 1 in Sec. 6.2 (Symmetric Bending) This is a symmetric I-section with no
    lateral supports, undergoing pure unidirectional cantilever bending. Note that units
    here are **inches**, to match the text.

    Returns:
        Section object
    """
    geom = nastran_sections.nastran_i(
        dim_1=6,
        dim_2=3,
        dim_3=3,
        dim_4=1,
        dim_5=1,
        dim_6=1,
    )
    geom = geom.shift_section(0, -3)
    geom = geom.create_mesh([0.25])
    xsect = Section(geom)
    xsect.calculate_geometric_properties()

    return xsect


@pytest.fixture
def peery_ex_7_2_1() -> ZSection:
    """Example 7.2.1.

    Example 1 in Sec. 7.2. (Unsymmetric Bending) This is an unsymmetric Z-section with
    no lateral supports. Note that units here are **inches**, to match the text.

    Returns:
        Geometry object
    """
    return ZSection(dim1=4, dim2=2, dim3=8, dim4=12, shift=(-5, -6), m=0.25)


def test_symmetric_ixx(peery_ex_6_2_1: Callable) -> None:
    """Tests Example 6.2.1 ixx.

    Args:
        peery_ex_6_2_1: peery_ex_6_2_1 test fixture
    """
    # Directly from the example, we know that
    # the 2nd moment of inertia resisting bending is.
    xsect = peery_ex_6_2_1

    check.almost_equal(xsect.section_props.ixx_g, 43.3, rel=1e-3)


def test_symmetric_fb(peery_ex_6_2_1: Callable) -> None:
    """Max bending stress on the section.

    Args:
        peery_ex_6_2_1: peery_ex_6_2_1 test fixture
    """
    xsect = peery_ex_6_2_1

    # Defined in the text
    moment = 8e5
    y = 3
    i = xsect.section_props.ixx_g
    xsect.calculate_warping_properties()
    stress = xsect.calculate_stress(mxx=moment)

    # The number quoted in the book. (Peery rounds this to the hundreds)
    # 55400 = 55427.3
    perfect_result = 55427.3

    # The number from the textbook equation
    computed_result = moment * y / i
    check.almost_equal(perfect_result, computed_result, rel=1e-3)

    # The max stress, computed through FEA on our mesh.
    numerical_result = max(stress.get_stress()[0]["sig_zz"])
    check.almost_equal(numerical_result, perfect_result, rel=1e-3)


def test_unsymmetric_ixx(peery_ex_7_2_1: Callable) -> None:
    """Tests Example 7.2.1 ixx.

    Args:
        peery_ex_7_2_1: peery_ex_7_2_1 test fixture
    """
    # Directly from the example, we know what
    # the section properties should be.
    xsect = peery_ex_7_2_1.xsect
    check.almost_equal(xsect.section_props.ixx_g, 693.3, rel=1e-3)


def test_unsymmetric_iyy(peery_ex_7_2_1: Callable) -> None:
    """Tests Example 7.2.1 iyy.

    Args:
        peery_ex_7_2_1: peery_ex_7_2_1 test fixture
    """
    # Directly from the example, we know what
    # the section properties should be.
    xsect = peery_ex_7_2_1.xsect
    check.almost_equal(xsect.section_props.iyy_g, 173.3, rel=1e-3)


def test_unsymmetric_ixy(peery_ex_7_2_1: Callable) -> None:
    """Tests Example 7.2.1 ixy.

    Args:
        peery_ex_7_2_1: peery_ex_7_2_1 test fixture
    """
    # Directly from the example, we know what
    # the section properties should be.
    xsect = peery_ex_7_2_1.xsect
    check.almost_equal(xsect.section_props.ixy_g, -240, rel=1e-3)


def test_unsymmetric_i11(peery_ex_7_2_1: Callable) -> None:
    """Tests Example 7.2.1 i11.

    Args:
        peery_ex_7_2_1: peery_ex_7_2_1 test fixture
    """
    # Directly from the example, we know what
    # the section properties should be.
    xsect = peery_ex_7_2_1.xsect
    check.almost_equal(xsect.section_props.i11_c, 787, rel=1e-3)


def test_unsymmetric_i22(peery_ex_7_2_1: Callable) -> None:
    """Tests Example 7.2.1 i22.

    Args:
        peery_ex_7_2_1: peery_ex_7_2_1 test fixture
    """
    # Directly from the example, we know what
    # the section properties should be.
    xsect = peery_ex_7_2_1.xsect
    check.almost_equal(xsect.section_props.i22_c, 79.5, rel=1e-3)


def test_fb_c(peery_ex_7_2_1: Callable) -> None:
    """Check the stress at point C.

    Args:
        peery_ex_7_2_1: peery_ex_7_2_1 test fixture
    """
    # Load from the text
    v = [-1e5, 1e4]
    # Coordinates of point C
    c = (1, 6)
    # The answer in the example
    # For this point, Peery rounds to the tens place,
    # thus -2380 is the exact number written in the book
    # but -2384 is the answer computed from his values.
    perfect_result = -2384
    # The simplified textbook equation
    text_result = round(-494 * 1 + -315 * 6)
    nodes = peery_ex_7_2_1.xsect.mesh_nodes
    assert len(nodes > 0)
    index, _ = get_node(nodes, c)
    _ = peery_ex_7_2_1.apply_load(v)
    computed_result = peery_ex_7_2_1.stress.get_stress()[0]["sig_zz"][index]

    check.almost_equal(text_result, perfect_result)
    check.almost_equal(computed_result, perfect_result, rel=1e-3)


def test_fb_b(peery_ex_7_2_1: Callable) -> None:
    """Check the stress at point B.

    Args:
        peery_ex_7_2_1: peery_ex_7_2_1 test fixture
    """
    # Load from the text
    v = [-1e5, 1e4]
    # Coordinates of point B
    b = (-5, 6)
    # The answer in the example
    perfect_result = 580
    # The sipmlified textbook equation
    text_result = round(-494 * -5 + -315 * 6)
    nodes = peery_ex_7_2_1.xsect.mesh_nodes
    index, _ = get_node(nodes, b)
    _ = peery_ex_7_2_1.apply_load(v)
    computed_result = peery_ex_7_2_1.stress.get_stress()[0]["sig_zz"][index]

    check.almost_equal(text_result, perfect_result)
    check.almost_equal(computed_result, perfect_result, rel=1e-3)


def test_fb_a(peery_ex_7_2_1: Callable) -> None:
    """Check the stress at point A.

    Args:
        peery_ex_7_2_1: peery_ex_7_2_1 test fixture
    """
    # Load from the text
    v = [-1e5, 1e4]
    # Coordinates of point A
    a = (-5, 4)
    # The answer in the example
    perfect_result = 1210
    # The simplified textbook equation
    text_result = round(-494 * -5 + -315 * 4)
    nodes = peery_ex_7_2_1.xsect.mesh_nodes
    index, _ = get_node(nodes, a)
    _ = peery_ex_7_2_1.apply_load(v)
    computed_result = peery_ex_7_2_1.stress.get_stress()[0]["sig_zz"][index]

    check.almost_equal(text_result, perfect_result)
    check.almost_equal(computed_result, perfect_result, rel=1e-3)
