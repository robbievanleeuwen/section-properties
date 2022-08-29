import pytest
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
