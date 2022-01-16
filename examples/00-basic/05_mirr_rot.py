r"""
.. _ref_ex_mirror_rot:

Mirroring and Rotating Geometry
-------------------------------

Mirror and rotate a cross section.

The following example demonstrates how geometry objects can be mirrored and
rotated. A 200PFC and 150PFC are placed back-to-back by using the
:func:`~sectionproperties.pre.geometry.Geometry.mirror_section` method and are
rotated counter-clockwise by 30 degrees by using the
:func:`~sectionproperties.pre.geometry.Geometry.rotate_section` method. The
geometry is cleaned to ensure there are no overlapping facets along the junction
between the two PFCs. A geometric, warping and plastic analysis is then carried out.

The geometry and mesh are plotted, and the mesh information printed to the terminal
before the analysis is carried out. Detailed time information is printed to the
terminal during the cross-section analysis stage and iteration information printed
for the plastic analysis. Once the analysis is complete, a plot of the various
calculated centroids is generated.
"""

# sphinx_gallery_thumbnail_number = 1

import sectionproperties.pre.library.steel_sections as steel_sections
from sectionproperties.analysis.section import Section

# %%
# Create a 200PFC and a 150PFC
pfc1 = steel_sections.channel_section(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8)
pfc2 = steel_sections.channel_section(
    d=150, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8
).shift_section(0, 26.5)

# %%
# Mirror the 200 PFC about the y-axis
pfc1 = pfc1.mirror_section(axis="y", mirror_point=[0, 0])

# %%
# Merge the pfc sections
geometry = ((pfc1 - pfc2) | pfc1) + pfc2

# %%
# Rotate the geometry counter-clockwise by 30 degrees
geometry = geometry.rotate_section(angle=30)
geometry.plot_geometry()

# %%
# Create a mesh and section. For the mesh, use a mesh size of 5 for the 200PFC
# and 4 for the 150PFC
geometry.create_mesh(mesh_sizes=[5, 4])

section = Section(geometry, time_info=True)
section.display_mesh_info()  # display the mesh information
section.plot_mesh()  # plot the generated mesh

# %%
# Perform a geometric, warping and plastic analysis, displaying the time info
# and the iteration info for the plastic analysis
section.calculate_geometric_properties()
section.calculate_warping_properties()
section.calculate_plastic_properties(verbose=True)

section.plot_centroids()
