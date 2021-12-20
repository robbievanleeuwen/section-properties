import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import Section

# create a 150x100x6 RHS
rhs = sections.rectangular_hollow_section(d=150, b=100, t=6, r_out=15, n_r=8)

# create a triangular section from points
# only the points are needed because they are in sequential order
# and represent only a single, contiguous region
points = [[0, 0], [50, 0], [25, 50]]
triangle = sections.Geometry.from_points(points)
triangle = triangle.align_center(rhs).align_to(rhs, on="top")

# create a 100x100x6 angle and position it on the right of the RHS
angle = sections.angle_section(d=100, b=100, t=6, r_r=8, r_t=5, n_r=8)
angle = angle.shift_section(x_offset=100, y_offset=25)

# combine the sections into a CompoundGeometry with `+` operator
geometry = rhs + triangle + angle
geometry.plot_geometry()  # plot the geometry

# create a mesh - use a mesh size of 2.5 for the RHS, 5 for the triangle and
# 3 for the angle.
geometry.create_mesh(mesh_sizes=[2.5, 5, 3])

# create a Section object
section = Section(geometry, time_info=True)
section.display_mesh_info()  # display the mesh information
section.plot_mesh()  # plot the generated mesh

# perform a geometric, warping and plastic analysis, displaying the time info
# and the iteration info for the plastic analysis
section.calculate_geometric_properties()
section.calculate_warping_properties()
section.calculate_plastic_properties(verbose=True)

# plot the centroids
section.plot_centroids()
