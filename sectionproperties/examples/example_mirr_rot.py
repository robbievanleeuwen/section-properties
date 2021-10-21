import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import Section

# create a 200PFC and a 150PFC
pfc1 = sections.channel_section(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8)
pfc2 = sections.channel_section(
    d=150, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8).shift_section(0, 26.5)

# mirror the 200 PFC about the y-axis
pfc1 = pfc1.mirror_section(axis='y', mirror_point=[0, 0])

# merge the pfc sections
geometry = ((pfc1 - pfc2) | pfc1) + pfc2

# rotate the geometry counter-clockwise by 30 degrees
geometry = geometry.rotate_section(angle=30)

# plot the geometry
geometry.plot_geometry()  

# create a mesh - use a mesh size of 5 for the 200PFC and 4 for the 150PFC
geometry.create_mesh(mesh_sizes=[5, 4])

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
