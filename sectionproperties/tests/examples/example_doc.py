import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection

geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
mesh = geometry.create_mesh(mesh_sizes=[2.5])
section = CrossSection(geometry, mesh)

section.calculate_geometric_properties()
section.calculate_warping_properties()
stress_post = section.calculate_stress(N=50e3, Mxx=-5e6, M22=2.5e6, Mzz=0.5e6,
                                       Vx=10e3, Vy=5e3)

stress_post.plot_stress_vm()
