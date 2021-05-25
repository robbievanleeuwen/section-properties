"""
.. _ref_ex_composite:

Creating a Composite Cross-Section
----------------------------------

Create a section of mixed materials.

The following example demonstrates how to create a composite cross-section by assigning
different material properties to various regions of the mesh. A steel 310UB40.4 is modelled
with a 50Dx600W timber panel placed on its top flange.

The geometry and mesh are plotted, and the mesh information printed to the terminal
before the analysis is carried out. All types of cross-section analyses are carried
out, with an axial force, bending moment and shear force applied during the stress
analysis. Once the analysis is complete, the cross-section properties are printed
to the terminal and a plot of the centroids and cross-section stresses generated
"""

# sphinx_gallery_thumbnail_number = 2

import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection
from sectionproperties.pre.pre import Material


# %%
# create material properties
steel = Material(
    name='Steel', elastic_modulus=200e3, poissons_ratio=0.3, yield_strength=500, color='grey'
)
timber = Material(
    name='Timber', elastic_modulus=8e3, poissons_ratio=0.35, yield_strength=20, color='burlywood'
)

# %%
# create 310UB40.4
ub = sections.ISection(d=304, b=165, t_f=10.2, t_w=6.1, r=11.4, n_r=8)

# create timber panel on top of the UB
panel = sections.RectangularSection(d=50, b=600, shift=[-217.5, 304])

# merge the two sections into one geometry object
geometry = sections.MergedSection([ub, panel])
geometry.clean_geometry()
geometry.plot_geometry()

# create a mesh - use a mesh size of 5 for the UB, 20 for the panel
mesh = geometry.create_mesh(mesh_sizes=[5, 20])

# create a CrossSection object - take care to list the materials in the same order as entered into
# the MergedSection
section = CrossSection(geometry, mesh, materials=[steel, timber])
section.display_mesh_info()
section.plot_mesh(materials=True, alpha=0.5)

# %%
# perform a geometric, warping and plastic analysis
section.calculate_geometric_properties(time_info=True)
section.calculate_warping_properties(time_info=True)
section.calculate_plastic_properties(time_info=True, verbose=True)

# perform a stress analysis with N = 100 kN, Mxx = 120 kN.m and Vy = 75 kN
stress_post = section.calculate_stress(N=-100e3, Mxx=-120e6, Vy=-75e3, time_info=True)

# %%
# print the results to the terminal
section.display_results()

# %%
# plot the centroids
section.plot_centroids()

# %%
# plot the axial stress
stress_post.plot_stress_n_zz(pause=False)

# %%
# plot the bending stress
stress_post.plot_stress_m_zz(pause=False)
# %%
# plot the shear stress
stress_post.plot_stress_v_zxy()
