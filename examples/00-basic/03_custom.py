"""
.. _ref_ex_custom:

Creating Custom Geometry
------------------------


Calculate section properties of a user-defined section from points and facets.

The following example demonstrates how geometry objects can be created from a
list of points, facets, holes and control points. An straight angle section with
a plate at its base is created from a list of points and facets. The bottom plate
is assigned a separate control point meaning two discrete regions are created.
Creating separate regions allows the user to control the mesh size in each region
and assign material properties to different regions. The geometry is cleaned to
remove the overlapping facet at the junction of the angle and the plate. A
geometric, warping and plastic analysis is then carried out.

The geometry and mesh are plotted before the analysis is carried out. Once the
analysis is complete, a plot of the various calculated centroids is generated.
"""

# sphinx_gallery_thumbnail_number = 2

import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import Section

# %%
# Define parameters for the angle section
a = 1
b = 2
t = 0.1

# %%
# Build the lists of points, facets, holes and control points
points = [[-t/2, -2*a], [t/2, -2*a], [t/2, -t/2], [a, -t/2], [a, t/2],
          [-t/2, t/2], [-b/2, -2*a], [b/2, -2*a], [b/2, -2*a-t],
          [-b/2, -2*a-t]]
facets = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], [6, 7], [7, 8],
          [8, 9], [9, 6]]
holes = []
control_points = [[0, 0], [0, -2 * a - t / 2]]

# %%
# Because we have two separate geometry regions (as indicated by our control_points)
# we create a CompoundGeometry from points
geometry = sections.CompoundGeometry.from_points(points, facets, control_points, holes)

# %%
# Create the mesh and section. For the mesh, use a smaller refinement for the angle region.
geometry.create_mesh(mesh_sizes=[0.0005, 0.001])

section = Section(geometry)
section.plot_mesh()  # plot the generated mesh

# %%
# Perform a geometric, warping and plastic analysis
section.calculate_geometric_properties()
section.calculate_warping_properties()
section.calculate_plastic_properties()

section.plot_centroids()
