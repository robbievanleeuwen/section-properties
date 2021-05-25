"""
.. _ref_ex_mirror_rot:

Mirroring and Rotating Geometry
-------------------------------

Mirror and rotate a cross section.

The following example demonstrates how geometry objects can be mirrored and
rotated. A 200PFC and 150PFC are placed back-to-back by using the
:func:`~sectionproperties.pre.sections.Geometry.mirror_section` method and are
rotated counter-clockwise by 30 degrees by using the
:func:`~sectionproperties.pre.sections.Geometry.rotate_section` method. The
geometry is cleaned to ensure there are no overlapping facets along the junction
between the two PFCs. A geometric, warping and plastic analysis is then carried out.

The geometry and mesh are plotted, and the mesh information printed to the terminal
before the analysis is carried out. Detailed time information is printed to the
terminal during the cross-section analysis stage and iteration information printed
for the plastic analysis. Once the analysis is complete, a plot of the various
calculated centroids is generated
"""

# sphinx_gallery_thumbnail_number = 1

import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection


# %%
# create a 200PFC and a 150PFC
pfc1 = sections.PfcSection(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8)
pfc2 = sections.PfcSection(d=150, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8, shift=[0, 26.5])

# %%
# mirror the 200 PFC about the y-axis
pfc1.mirror_section(axis='y', mirror_point=[0, 0])

# %%
# merge the pfc sections
geometry = sections.MergedSection([pfc1, pfc2])

# %%
# rotate the geometry counter-clockwise by 30 degrees
geometry.rotate_section(angle=30)

# %%
# clean the geometry - print cleaning information to the terminal
geometry.clean_geometry(verbose=True)
geometry.plot_geometry()

# %%
# create a mesh - use a mesh size of 5 for the 200PFC and 4 for the 150PFC
mesh = geometry.create_mesh(mesh_sizes=[5, 4])

# %%
# create a CrossSection object
section = CrossSection(geometry, mesh)
section.display_mesh_info()
section.plot_mesh()

# %%
# perform a geometric, warping and plastic analysis, displaying the time info and the iteration
# info for the plastic analysis
section.calculate_geometric_properties(time_info=True)
section.calculate_warping_properties(time_info=True)
section.calculate_plastic_properties(time_info=True, verbose=True)

# %%
# plot the centroids
section.plot_centroids()
