import sectionproperties.pre.sections as sections
from sectionproperties.pre.pre import Material
from sectionproperties.analysis.cross_section import CrossSection

# concrete variables
depth = 600
width = 300
concrete_mesh = 500
concrete_e = 30000
concrete_nu = 0.2
concrete_yield = 30

# steel variables
diameter = 20
cover = 30
n_bars = 3
bar_mesh = 50
steel_e = 200000
steel_nu = 0.3
steel_yield = 500

# initialise lists holding geometry, mesh sizes and materials
geometries = []
mesh_sizes = []
materials = []

# create materials
concrete = Material(name='Concrete', elastic_modulus=concrete_e,
                    poissons_ratio=concrete_nu, yield_strength=concrete_yield,
                    color='lightgrey')
steel = Material(name='Steel', elastic_modulus=steel_e,
                 poissons_ratio=steel_nu, yield_strength=steel_yield,
                 color='grey')

# create concrete cross-section
geometries.append(sections.RectangularSection(d=depth, b=width))
mesh_sizes.append(concrete_mesh)
materials.append(concrete)

# create steel bars
spacing = (width - 2 * cover - diameter) / (n_bars - 1)

for i in range(n_bars):
    shift = [cover + diameter / 2 + i * spacing, cover + diameter / 2]
    bar = sections.CircularSection(d=diameter, n=16, shift=shift)

    geometries.append(bar)
    mesh_sizes.append(bar_mesh)
    materials.append(steel)

geometry = sections.MergedSection(geometries)
geometry.plot_geometry()
mesh = geometry.create_mesh(mesh_sizes=mesh_sizes)

section = CrossSection(geometry, mesh, materials, time_info=True)
section.display_mesh_info()
section.plot_mesh(materials=True, alpha=0.5)

section.calculate_geometric_properties(time_info=True)
# section.display_results()
# section.calculate_warping_properties(time_info=True)
section.calculate_plastic_properties(time_info=True, verbose=True)
# stress_result = section.calculate_stress(N=1e3, Vy=3e3, Mxx=1e6, Mzz=5e5,
#                                          time_info=True)

section.plot_centroids()
# section.display_results()
# stress_result.plot_stress_n_zz()
# stress_result.plot_vector_mzz_zxy(pause=False)
# stress_result.plot_vector_vy_zxy(pause=False)
# stress_result.plot_vector_zxy()
