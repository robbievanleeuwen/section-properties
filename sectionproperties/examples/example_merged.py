import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection

# create a 150x100x6 RHS
rhs = sections.Rhs(d=150, b=100, t=6, r_out=15, n_r=8)

# create a triangular section on top of the RHS
points = [[0, 0], [50, 0], [25, 50]]
facets = [[0, 1], [1, 2], [2, 0]]
holes = []
control_points = [[25, 25]]
triangle = sections.CustomSection(points, facets, holes, control_points, shift=[25, 150])

# create a 100x100x6 EA on the right of the RHS
angle = sections.AngleSection(d=100, b=100, t=6, r_r=8, r_t=5, n_r=8, shift=[100, 25])

# create a list of the sections to be merged
section_list = [rhs, triangle, angle]

# merge the three sections into one geometry object
geometry = sections.MergedSection(section_list)

# clean the geometry - print cleaning information to the terminal
geometry.clean_geometry(verbose=True)
geometry.plot_geometry()  # plot the geometry

# create a mesh - use a mesh size of 2.5 for the RHS, 5 for the triangle and 3 for the angle
mesh = geometry.create_mesh(mesh_sizes=[2.5, 5, 3])

# create a CrossSection object
section = CrossSection(geometry, mesh)
section.display_mesh_info()  # display the mesh information
section.plot_mesh()  # plot the generated mesh

# perform a geometric, warping and plastic anaylsis, displaying the time info and the iteration
# info for the plastic analysis
section.calculate_geometric_properties(time_info=True)
section.calculate_warping_properties(time_info=True)
section.calculate_plastic_properties(time_info=True, verbose=True)

# plot the centroids
section.plot_centroids()
