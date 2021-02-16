from sectionproperties.pre.sections import *
from sectionproperties.analysis.cross_section import Section

def test_composite():
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


if __name__ == "__main__":
    test_composite()   
    