import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection

a = 1
t = 0.1

points = [[-t / 2, -2 * a], [t / 2, -2 * a], [t / 2, -t / 2], [a, -t / 2],
          [a, t / 2], [-t / 2, t / 2]]
facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]
holes = []
control_points = [[0, 0]]

geometry = sections.CustomSection(points, facets, holes, control_points)
# geometry.plot_geometry()
mesh = geometry.create_mesh(mesh_sizes=[0.0005])

section = CrossSection(geometry, mesh)
section.display_mesh_info()
# section.plot_mesh()
section.calculate_geometric_properties(time_info=True)
# section.calculate_warping_properties(time_info=True)
section.calculate_plastic_properties(time_info=True)
# stress_result = section.calculate_stress(N=1, Mxx=1e6, time_info=True)

# section.plot_centroids()
section.display_results()
# stress_result.plot_stress_n_zz()
# stress_result.plot_stress_mxx_zz()
