r"""
.. _ref_ex_merged:

Creating a Merged Section
-------------------------

Merge two sections together into a single larger section.

The following example demonstrates how to merge multiple geometry objects into
a single geometry object. A 150x100x6 RHS is modelled with a solid 50x50 triangular
section on its top and a 100x100x6 EA section on its right side. The three geometry
objects are merged together using the :class:`~sectionproperties.pre.sections.MergedSection`
class. The order of the geometry objects in the list that is passed into the constructor of the
:class:`~sectionproperties.pre.sections.MergedSection` class is important, as this same
order relates to specifying mesh sizes and material properties.

Once the geometry has been merged, it is vital to clean the geometry to remove
any artefacts that may impede the meshing algorithm. A mesh is created with a mesh
size of 2.5 mm\ :sup:`2` for the RHS (first in ``section_list``), 5 mm\ :sup:`2` for the triangle
(second in ``section_list``) and 3 mm\ :sup:`2` for the angle (last in ``section_list``).

The geometry and mesh are plotted, and the mesh information printed to the terminal
before the analysis is carried out. Detailed time information is printed to the
terminal during the cross-section analysis stage. Once the analysis is complete,
the centroids are plotted
"""

# sphinx_gallery_thumbnail_number = 1

import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection


# %%
# create a 150x100x6 RHS
rhs = sections.Rhs(d=150, b=100, t=6, r_out=15, n_r=8)

# %%
# create a triangular section on top of the RHS
points = [[0, 0], [50, 0], [25, 50]]
facets = [[0, 1], [1, 2], [2, 0]]
holes = []
control_points = [[25, 25]]
triangle = sections.CustomSection(points, facets, holes, control_points, shift=[25, 150])

# %%
# create a 100x100x6 EA on the right of the RHS
angle = sections.AngleSection(d=100, b=100, t=6, r_r=8, r_t=5, n_r=8, shift=[100, 25])

# %%
# create a list of the sections to be merged
section_list = [rhs, triangle, angle]

# %%
# merge the three sections into one geometry object
geometry = sections.MergedSection(section_list)

# %%
# clean the geometry - print cleaning information to the terminal
geometry.clean_geometry(verbose=True)
geometry.plot_geometry()  # plot the geometry

# %%
# create a mesh - use a mesh size of 2.5 for the RHS, 5 for the triangle and 3 for the angle
mesh = geometry.create_mesh(mesh_sizes=[2.5, 5, 3])

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
