import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import Section

rect = sections.rectangular_section(d=500, b=300)
box_girder = sections.box_girder_section(d=400, b_t=500, b_b=200, t_ft=12, t_fb=12, t_w=6)

points = [[0,0], [5,0], [11,8], [0,2]]
facets = [[0,1], [1,2], [2,3], [3,0]]
control_points = [[5,5]]
custom = sections.Geometry.from_points(points, facets, control_points, holes=None)

box_girder.plot_geometry()
