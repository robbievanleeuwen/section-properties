import pytest
import math
import sectionproperties.pre.library.primitive_sections as primitive_sections
from sectionproperties.analysis.section import Section

# define geometry
rect = primitive_sections.rectangular_section(b=50, d=100)
rect.create_mesh(mesh_sizes=0)  # coarse mesh
sec = Section(rect)


# test stress runtime errors
def test_stress_runtime_errors():
    # check runtime error with no analysis performed
    with pytest.raises(RuntimeError):
        sec.calculate_stress()

    sec.calculate_geometric_properties()

    # check runtime errors with shear/torsion applied, no warping analysis
    with pytest.raises(RuntimeError):
        sec.calculate_stress(Vx=1)
        sec.get_stress_at_points(pts=[[10, 10]], Vx=1)

    with pytest.raises(RuntimeError):
        sec.calculate_stress(Vy=1)
        sec.get_stress_at_points(pts=[[10, 10]], Vy=1)

    with pytest.raises(RuntimeError):
        sec.calculate_stress(Mzz=1)
        sec.get_stress_at_points(pts=[[10, 10]], Mzz=1)

    # check no runtime errors with no shear/torsion applied
    sec.calculate_stress(N=1)
    sec.calculate_stress(Mxx=1)
    sec.calculate_stress(Myy=1)
    sec.calculate_stress(M11=1)
    sec.calculate_stress(M22=1)

    sec.calculate_warping_properties()

    # check no runtime errors after warping analysis
    sec.calculate_stress(N=1)
    sec.calculate_stress(Mxx=1)
    sec.calculate_stress(Myy=1)
    sec.calculate_stress(M11=1)
    sec.calculate_stress(M22=1)
    sec.calculate_stress(Vx=1)
    sec.calculate_stress(Vy=1)
    sec.calculate_stress(Mzz=1)
    sec.get_stress_at_points(pts=[[10, 10]], Mzz=1)


def test_rectangle():
    Mxx = 7
    Sy = 50.0 * 100.0**2 / 6.0
    sig_max = Mxx / Sy
    (sig_0, sig_1, sig_2) = sec.get_stress_at_points(
        pts=[[25, 50], [25, 75], [25, 100]], Mxx=Mxx
    )
    assert sig_0 == pytest.approx((0, 0, 0))
    assert sig_1 == pytest.approx((sig_max / 2.0, 0, 0))
    assert sig_2 == pytest.approx((sig_max, 0, 0))


def test_rotated_rectangle():
    b = 50
    d = 100
    angle = math.atan(100 / 50)
    cx = b / 2 * math.cos(angle) - d / 2 * math.sin(angle)
    cy = b / 2 * math.sin(angle) + d / 2 * math.cos(angle)
    Sy = b * d / 6.0 * cy
    Mxx = 7
    sig_max = Mxx / Sy
    rot_rect = (
        primitive_sections.rectangular_section(b=b, d=d)
        .shift_section(-b / 2, -d / 2)
        .rotate_section(angle, use_radians=True)
    )
    rot_rect.create_mesh(mesh_sizes=0)  # coarse mesh
    rot_sec = Section(rot_rect)
    rot_sec.calculate_geometric_properties()
    rot_sec.calculate_warping_properties()
    (sig_0, sig_1) = rot_sec.get_stress_at_points(pts=[[cx, 0], [cx, cy]], Mxx=Mxx)
    assert sig_1 == pytest.approx((sig_max, 0, 0))
