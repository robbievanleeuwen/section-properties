import sectionproperties.pre.sections as sections
from sectionproperties.pre.pre import Material
from sectionproperties.analysis.cross_section import Section

# create material properties
steel = Material(name='Steel', elastic_modulus=200e3, poissons_ratio=0.3,
                 yield_strength=500, color='grey')
timber = Material(name='Timber', elastic_modulus=8e3, poissons_ratio=0.35,
                  yield_strength=20, color='burlywood')

# create 310UB40.4
ub = sections.i_section(d=304, b=165, t_f=10.2, t_w=6.1, r=11.4, n_r=8, material=steel)

# create timber panel on top of the UB
panel = sections.rectangular_section(d=50, b=600, material=timber)
panel = panel.align_center(ub).align_to(ub, on="top")

# merge the two sections into one geometry object
geometry = sections.CompoundGeometry([ub, panel])

# create a mesh - use a mesh size of 5 for the UB, 20 for the panel
geometry.create_mesh(mesh_sizes=[5, 20])

# create a Section object
section = Section(geometry, time_info=True)
section.display_mesh_info()  # display the mesh information

# plot the mesh with coloured materials and a line transparency of 0.6
section.plot_mesh(materials=True, alpha=0.6)

# perform a geometric, warping and plastic analysis
section.calculate_geometric_properties()
section.calculate_warping_properties()
section.calculate_plastic_properties(verbose=True)

# perform a stress analysis with N = 100 kN, Mxx = 120 kN.m and Vy = 75 kN
stress_post = section.calculate_stress(N=-100e3, Mxx=-120e6, Vy=-75e3)

# print the results to the terminal
section.display_results()

# plot the centroids
section.plot_centroids()

stress_post.plot_stress_n_zz(pause=False)  # plot the axial stress
stress_post.plot_stress_m_zz(pause=False)  # plot the bending stress
stress_post.plot_stress_v_zxy()  # plot the shear stress
