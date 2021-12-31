r"""
.. _ref_ex_composite:

Creating a Composite Section
----------------------------

Create a section of mixed materials.

The following example demonstrates how to create a composite cross-section by assigning
different material properties to various regions of the mesh. A steel 310UB40.4 is modelled
with a 50Dx600W timber panel placed on its top flange.

The geometry and mesh are plotted, and the mesh information printed to the terminal
before the analysis is carried out. All types of cross-section analyses are carried
out, with an axial force, bending moment and shear force applied during the stress
analysis. Once the analysis is complete, the cross-section properties are printed
to the terminal and a plot of the centroids and cross-section stresses generated.
"""

# sphinx_gallery_thumbnail_number = 2

import sectionproperties.pre.library.standard_sections as sections
import sectionproperties.pre.library.steel_sections as steel_sections
from sectionproperties.pre.geometry import CompoundGeometry
from sectionproperties.pre.pre import Material
from sectionproperties.analysis.section import Section

# %%
# Create material properties
steel = Material(name='Steel', elastic_modulus=200e3, poissons_ratio=0.3,
                 yield_strength=500, density=8.05e-6 ,color='grey')
timber = Material(name='Timber', elastic_modulus=8e3, poissons_ratio=0.35,
                  yield_strength=20, density=0.78e-6, color='burlywood')

# %%
# Create 310UB40.4
ub = steel_sections.i_section(d=304, b=165, t_f=10.2, t_w=6.1, r=11.4, n_r=8, material=steel)

# %%
# Create timber panel on top of the UB
panel = sections.rectangular_section(d=50, b=600, material=timber)
panel = panel.align_center(ub).align_to(ub, on="top")

# %%
# Merge the two sections into one geometry object
geometry = sections.CompoundGeometry([ub, panel])

# %%
# Create a mesh and a Section object. For the mesh use a mesh size of 5 for
# the UB, 20 for the panel
geometry.create_mesh(mesh_sizes=[5, 20])
section = Section(geometry, time_info=True)
section.display_mesh_info()  # display the mesh information

# %%
# Plot the mesh with coloured materials and a line transparency of 0.6
section.plot_mesh(materials=True, alpha=0.6)

# %%
# Perform a geometric, warping and plastic analysis
section.calculate_geometric_properties()
section.calculate_warping_properties()
section.calculate_plastic_properties(verbose=True)

# %%
# Perform a stress analysis with N = 100 kN, Mxx = 120 kN.m and Vy = 75 kN
stress_post = section.calculate_stress(N=-100e3, Mxx=-120e6, Vy=-75e3)

# %%
# Print the results to the terminal
section.display_results()

# %%
# Plot the centroids
section.plot_centroids()

# %%
# Plot the axial stress
stress_post.plot_stress_n_zz(pause=False)

# %%
# Plot the bending stress
stress_post.plot_stress_m_zz(pause=False)

# %%
# Plot the shear stress
stress_post.plot_stress_v_zxy()