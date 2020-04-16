import sectionproperties.pre.nastran_sections as nsections
from sectionproperties.analysis.cross_section import CrossSection

# create a HAT1 section
geometry = nsections.HAT1Section(DIM1=4.0, DIM2=2.0, DIM3=1.5, DIM4=0.1875, DIM5=0.375)
geometry.plot_geometry()  # plot the geometry

# create a mesh with a maximum elemental area of 0.005
mesh = geometry.create_mesh(mesh_sizes=[0.005])

section = CrossSection(geometry, mesh)  # create a CrossSection object
section.display_mesh_info()  # display the mesh information
section.plot_mesh()  # plot the generated mesh`

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
