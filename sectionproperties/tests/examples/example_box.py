import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection

geometry = sections.Rhs(d=100, b=150, t=6, r_out=15, n_r=8)
geometry.plot_geometry(pause=False)
mesh = geometry.create_mesh(mesh_sizes=[5])

section = CrossSection(geometry, mesh)
section.display_mesh_info()
section.plot_mesh()
section.calculate_geometric_properties(time_info=True)
section.calculate_warping_properties(time_info=True)
# section.calculate_plastic_properties(time_info=True, verbose=True)
stress_result = section.calculate_stress(
    N=1, Vx=2, Vy=3, Mxx=4, Myy=5, M11=6, M22=7, Mzz=8,  time_info=True)

section.plot_centroids()
section.display_results()
stress_result.plot_stress_m_zz()
stress_result.plot_stress_n_zz()
# stress_result.plot_stress_vm()
