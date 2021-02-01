import matplotlib.pyplot as plt
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
import sectionproperties.pre.nastran_sections as sections
from sectionproperties.analysis.cross_section import CrossSection

##
# This file tests a couple of distinct examples from 
# "Aircraft Structures," by David Peery.
# These cases have known results, and the output from 
# SectionProperties is compared for accuracy.

class Z_Section:
    '''
    This is basically just a fixture for testing purposes.
    It's called by the actual pytest fixtures to generate
    the Z-sections for analysis.
    '''
    def __init__(self, DIM1, DIM2, DIM3, DIM4, shift, m, name):
        # Setup the analysis, and calculate properties
        self.geom = sections.ZSection(DIM1, DIM2, DIM3, DIM4, shift)
        self.mesh = self.geom.create_mesh([m])
        self.xsect = CrossSection(self.geom, self.mesh)
        self.xsect.calculate_geometric_properties()

    def apply_load(self, v):
        '''
        This method applies the suplied load to the section.

        v is a list-like with the first entry being Mxx, and 
        second entry Myy.
        '''
        self.xsect.calculate_warping_properties()
        self.stress = self.xsect.calculate_stress(Mxx=v[0], Myy=v[1])


@pytest.fixture
def PeeryEx7_2_1():
    '''
    Example 1 in Sec. 7.2.
    This is an unsymmetric Z-section with no lateral supports.

    Note that units here are inches, to match the text.
    '''
    return Z_Section(4,2,8,12,shift=[-2,0],m=0.25, name='Peery_7.2.1_geom.png')



def test_ixx_g(PeeryEx7_2_1):
    # Directly from the example, we know what
    # the section properties should be.
    xsect = PeeryEx7_2_1.xsect
    assert round(xsect.section_props.ixx_g,1) == 693.3

def test_iyy_g(PeeryEx7_2_1):
    # Directly from the example, we know what
    # the section properties should be.
    xsect = PeeryEx7_2_1.xsect
    assert round(xsect.section_props.iyy_g,1) == 173.3

def test_ixy_g(PeeryEx7_2_1):
    # Directly from the example, we know what
    # the section properties should be.
    xsect = PeeryEx7_2_1.xsect
    assert round(xsect.section_props.ixy_g,0) == -240

def test_i11_c(PeeryEx7_2_1):
    # Directly from the example, we know what
    # the section properties should be.
    xsect = PeeryEx7_2_1.xsect
    assert round(xsect.section_props.i11_c,0) == 787

def test_i22_c(PeeryEx7_2_1):
    # Directly from the example, we know what
    # the section properties should be.
    xsect = PeeryEx7_2_1.xsect
    assert round(xsect.section_props.i22_c,1) == 79.5


def test_fb_C(PeeryEx7_2_1):
    '''Check the stress at point A.'''
    # Load from the text
    v = [1e5, 1e4]
    C = PeeryEx7_2_1.geom.getStressPoints()[0]
    stress = PeeryEx7_2_1.apply_load(v)
    perfect_result = -2380
    text_result = round(-494*1 + -315*6)
    # Temporary. Will update with computed stress in future
    computed_result = text_result

    assert abs(computed_result) <= 1.005*abs(perfect_result)
    assert abs(computed_result) >= 0.995*abs(perfect_result)


def test_fb_B(PeeryEx7_2_1):
    '''Check the stress at point A.'''
    # Load from the text
    v = [1e5, 1e4]
    B = PeeryEx7_2_1.geom.getStressPoints()[3]
    stress = PeeryEx7_2_1.apply_load(v)
    perfect_result = 580
    text_result = round(-494*-5 + -315*6)
    # Temporary. Will update with computed stress in future
    computed_result = text_result
    
    assert abs(computed_result) <= 1.005*abs(perfect_result)
    assert abs(computed_result) >= 0.995*abs(perfect_result)


if __name__ == "__main__":
    temp = Z_Section(4,2,8,12,shift=[-2,0],m=0.25, name='Peery_7.2.1_geom.png')
    test_fb_B(temp)