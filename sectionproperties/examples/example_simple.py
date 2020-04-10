import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection

# create a 50 diameter circle discretised by 64 points
geometry = sections.CircularSection(d=50, n=64)
geometry.plot_geometry()  # plot the geometry

# create a mesh with a mesh size of 2.5
mesh = geometry.create_mesh(mesh_sizes=[2.5])

section = CrossSection(geometry, mesh)  # create a CrossSection object
section.display_mesh_info()  # display the mesh information
section.plot_mesh()  # plot the generated mesh

# perform a geometric, warping and plastic anaylsis, displaying the time info
section.calculate_geometric_properties(time_info=True)
section.calculate_warping_properties(time_info=True)
section.calculate_plastic_properties(time_info=True)

# print the results to the terminal
section.display_results()

# get the second moments of area and the torsion constant
(ixx_c, iyy_c, ixy_c) = section.get_ic()
j = section.get_j()

# print the sum of the second moments of area and the torsion constant
print("Ixx + Iyy = {0:.3f}".format(ixx_c + iyy_c))
print("J = {0:.3f}".format(j))
