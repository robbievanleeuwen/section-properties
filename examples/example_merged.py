r"""
.. _ref_ex_merged:

Creating a Built-Up Section
---------------------------

Merge two sections together into a single larger section.

The following example demonstrates how to combine multiple geometry objects into
a single geometry object. A 150x100x6 RHS is modelled with a solid 50x50 triangular
section on its top and a 100x100x6 angle section on its right side.
The three geometry objects are combined together as a :class:`~sectionproperties.pre.sections.CompoundGeometry`
object using the `+` operator.

To manipulate individual geometries into the final shape, there are a variety of
methods available to move and align. This example uses `.align_center()`, `.align_to()`,
and `.shift_section()`.

The geometry and mesh are plotted, and the mesh information printed to the terminal
before the analysis is carried out. Detailed time information is printed to the
terminal during the cross-section analysis stage. Once the analysis is complete,
the centroids are plotted.
"""

# sphinx_gallery_thumbnail_number = 1

import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import Section

# %%
# Create a 150x100x6 RHS
rhs = sections.rectangular_hollow_section(d=150, b=100, t=6, r_out=15, n_r=8)

# %%
# Create a triangular section from points
# only the points are needed because they are in sequential order
# and represent only a single, contiguous region
points = [[0, 0], [50, 0], [25, 50]]
triangle = sections.Geometry.from_points(points)
triangle = triangle.align_center(rhs).align_to(rhs, on="top")

# %%
# Create a 100x100x6 angle and position it on the right of the RHS
angle = sections.angle_section(d=100, b=100, t=6, r_r=8, r_t=5, n_r=8)
angle = angle.shift_section(x_offset=100, y_offset=25)

# %%
# Combine the sections into a CompoundGeometry with `+` operator
geometry = rhs + triangle + angle
geometry.plot_geometry()  # plot the geometry

# %%
# Create a mesh and section. For the mesh, use a mesh size of 2.5 for 
# the RHS, 5 for the triangle and 3 for the angle.
geometry.create_mesh(mesh_sizes=[2.5, 5, 3])

section = Section(geometry, time_info=True)
section.display_mesh_info()  # display the mesh information
section.plot_mesh()  # plot the generated mesh

# %%
# Perform a geometric, warping and plastic analysis, displaying the time info
# and the iteration info for the plastic analysis
section.calculate_geometric_properties()
section.calculate_warping_properties()
section.calculate_plastic_properties(verbose=True)

# plot the centroids
section.plot_centroids()
