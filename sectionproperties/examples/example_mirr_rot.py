import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection

# create a 200PFC and a 150PFC
pfc1 = sections.PfcSection(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8)
pfc2 = sections.PfcSection(d=150, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8, shift=[0, 26.5])

# mirror the 200 PFC about the y-axis
pfc1.mirror_section(axis='y', mirror_point=[0, 0])

# merge the pfc sections
geometry = sections.MergedSection([pfc1, pfc2])

# rotate the geometry counter-clockwise by 30 degrees
geometry.rotate_section(angle=30)

# clean the geometry - print cleaning information to the terminal
geometry.clean_geometry(verbose=True)
geometry.plot_geometry()  # plot the geometry

# create a mesh - use a mesh size of 5 for the 200PFC and 4 for the 150PFC
mesh = geometry.create_mesh(mesh_sizes=[5, 4])

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
