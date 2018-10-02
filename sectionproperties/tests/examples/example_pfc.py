import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection

pfc1 = sections.PfcSection(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8)
pfc2 = sections.PfcSection(d=150, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8)
pfc1.mirror_section(axis='y', mirror_point=[0, 0])
geometry = sections.MergedSection([pfc1, pfc2])
geometry.rotate_section(angle=30)
geometry.clean_geometry(verbose=True)
geometry.plot_geometry()
mesh = geometry.create_mesh(mesh_sizes=[5, 4])

section = CrossSection(geometry, mesh)
section.display_mesh_info()
section.plot_mesh()
section.calculate_geometric_properties(time_info=True)
section.calculate_warping_properties(time_info=True)
# section.calculate_plastic_properties()
stress_result = section.calculate_stress(N=1e3, Vy=3e3, Mxx=1e6, Mzz=5e5,
                                         time_info=True)

section.plot_centroids()
section.display_results()
stress_result.plot_stress_n_zz()
stress_result.plot_vector_mzz_zxy()
