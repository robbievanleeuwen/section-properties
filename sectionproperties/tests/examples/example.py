import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection

isection = sections.ISection(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8)
box = sections.RhsSection(d=100, b=150, t=6, r_out=15, n_r=8,
                          shift=[-8.5, 203])

geometry = sections.MergedSection([isection, box])
geometry.plot_geometry()

mesh = geometry.create_mesh(mesh_sizes=[5, 2.5])

section = CrossSection(geometry, mesh)
section.display_mesh_info()
section.plot_mesh()

section.calculate_geometric_properties(time_info=True)
section.calculate_warping_properties(time_info=True)
section.plot_centroids()

section.display_results()
