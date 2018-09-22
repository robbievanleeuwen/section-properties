import sectionproperties.pre.sections as sections
import sectionproperties.pre.pre as pre

points = [[0, 0], [0.5, 0.5], [2, 2], [3, 3]]
facets = [[1, 0], [0, 2], [0, 3]]
holes = []
control_points = []

geometry = sections.CustomSection(points, facets, holes, control_points)
geometry = pre.GeometryCleaner(geometry).clean_geometry()
print(geometry.facets)
geometry.plot_geometry()
