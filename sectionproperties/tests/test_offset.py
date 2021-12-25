import pytest_check as check
import numpy as np
import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import Section


r_tol = 1e-3

def test_rectangular_offset():
    # exterior negative offset
    rect = sections.rectangular_section(d=500, b=300)
    rect = rect.offset_perimeter(amount=-10, where='exterior')
    rect.create_mesh([200])
    section = Section(rect)
    section.calculate_geometric_properties()
    area = 480 * 280
    check.almost_equal(section.get_area(), area, rel=r_tol)

    # exterior positive offset
    rect = sections.rectangular_section(d=500, b=300)
    rect = rect.offset_perimeter(amount=10, where='exterior')
    rect.create_mesh([200])
    section = Section(rect)
    section.calculate_geometric_properties()
    area = 520 * 320 - (20 * 20 - np.pi * 10 * 10)
    check.almost_equal(section.get_area(), area, rel=r_tol)

def test_box_offset():
    # exterior negative offset
    box = sections.rectangular_hollow_section(d=200, b=100, t=10, r_out=0, n_r=1)
    box = box.offset_perimeter(amount=-5, where='exterior')
    box.create_mesh([50])
    section = Section(box)
    section.calculate_geometric_properties()
    area = 190 * 90 - 180 * 80
    check.almost_equal(section.get_area(), area, rel=r_tol)

    # exterior positve offset
    box = sections.rectangular_hollow_section(d=200, b=100, t=10, r_out=0, n_r=1)
    box = box.offset_perimeter(amount=5, where='exterior')
    box.create_mesh([50])
    section = Section(box)
    section.calculate_geometric_properties()
    area = 210 * 110 - (10 * 10 - np.pi * 5 * 5) - 180 * 80
    check.almost_equal(section.get_area(), area, rel=r_tol)

    # interior negative offset
    box = sections.rectangular_hollow_section(d=200, b=100, t=10, r_out=0, n_r=1)
    box = box.offset_perimeter(amount=-5, where='interior')
    box.create_mesh([50])
    section = Section(box)
    section.calculate_geometric_properties()
    area = 200 * 100 - 170 * 70
    check.almost_equal(section.get_area(), area, rel=r_tol)

    # interior positive offset
    box = sections.rectangular_hollow_section(d=200, b=100, t=10, r_out=0, n_r=1)
    box = box.offset_perimeter(amount=5, where='interior')
    box.create_mesh([50])
    section = Section(box)
    section.calculate_geometric_properties()
    area = 200 * 100 - 190 * 90 + (10 * 10 - np.pi * 5 * 5)
    check.almost_equal(section.get_area(), area, rel=r_tol)

    # all negative offset
    box = sections.rectangular_hollow_section(d=200, b=100, t=10, r_out=0, n_r=1)
    box = box.offset_perimeter(amount=-2.5, where='all')
    box.create_mesh([50])
    section = Section(box)
    section.calculate_geometric_properties()
    area = 195 * 95 - 185 * 85 + (5 * 5 - np.pi * 2.5 * 2.5)
    check.almost_equal(section.get_area(), area, rel=r_tol)

    # all positive offset
    box = sections.rectangular_hollow_section(d=200, b=100, t=10, r_out=0, n_r=1)
    box = box.offset_perimeter(amount=5, where='all')
    box.create_mesh([50])
    section = Section(box)
    section.calculate_geometric_properties()
    area = 210 * 110 - (10 * 10 - np.pi * 5 * 5) - 170 * 70
    check.almost_equal(section.get_area(), area, rel=r_tol)

def test_compound_rectangular_offset():
    # This fails to generate the correct offset geometry
    rect1 = sections.rectangular_section(d=50, b=50)
    rect2 = sections.rectangular_section(d=50, b=50).align_to(rect1, "right")
    geom = rect1 + rect2
    geom = geom.offset_perimeter(amount=-5, where='exterior')
    geom.create_mesh([50])
    section = Section(geom)
    section.plot_mesh()

def test_compound_rectangular_isection_offset():
    # This fails to generate the correct offset geometry
    d = 300
    b = 150
    tf = 10
    tw = 8
    r = 12
    b_p = 250
    t_p = 16
    ub = sections.i_section(d=d, b=b, t_f=tf, t_w=tw, r=r, n_r=16)
    plate = sections.rectangular_section(b=b_p, d=t_p).align_center(ub).align_to(ub, on="top")
    geom = ub + plate
    geom = geom.offset_perimeter(amount=-2, where='exterior')
    geom.create_mesh([100])
    section = Section(geom)
    section.plot_mesh()
