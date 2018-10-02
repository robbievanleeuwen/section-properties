import sectionproperties.pre.sections as sections
from sectionproperties.pre.pre import Material
from sectionproperties.analysis.cross_section import CrossSection

# initialise lists holding geometry, mesh sizes and materials
geometries = []
mesh_sizes = []
materials = []

# create materials
concrete = Material(name='Concrete', elastic_modulus=30000,
                    poissons_ratio=0.2, yield_strength=30,
                    color='lightgrey')
steel = Material(name='Steel', elastic_modulus=200000,
                 poissons_ratio=0.3, yield_strength=500,
                 color='grey')

# create concrete cross-section
d = 600
w = 1200
t = 100

points = [[0, 0], [w, 0], [w, d], [0, d],
          [t, t], [w - t, t], [w - t, d - t], [t, d - t]]
facets = [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4]]
holes = [[w / 2, d / 2]]
control_points = [[t / 2, d / 2]]
box = sections.CustomSection(points, facets, holes, control_points)
geometries.append(box)
mesh_sizes.append(200)
materials.append(concrete)

# create steel bars
diameter = 24
cover = 30
n_bars = 8
bar_mesh = 10
spacing = (w - 2 * cover - diameter) / (n_bars - 1)

for i in range(n_bars):
    shift = [cover + diameter / 2 + i * spacing, cover + diameter / 2]
    bar = sections.CircularSection(d=diameter, n=16, shift=shift)

    geometries.append(bar)
    mesh_sizes.append(bar_mesh)
    materials.append(steel)

geometry = sections.MergedSection(geometries)
geometry.plot_geometry()
mesh = geometry.create_mesh(mesh_sizes=mesh_sizes)

section = CrossSection(geometry, mesh, materials)
section.display_mesh_info()
section.plot_mesh(materials=True, alpha=0.5)

# section.calculate_geometric_properties(time_info=True)
# section.calculate_warping_properties(time_info=True)
# section.display_results()
# section.plot_centroids()
#
# stress_result = section.calculate_stress(Mzz=100e6, time_info=True)
# stress_result.plot_vector_mzz_zxy()
