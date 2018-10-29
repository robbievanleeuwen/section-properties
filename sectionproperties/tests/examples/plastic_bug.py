import sectionproperties.pre.sections as sections
from sectionproperties.pre.pre import Material
from sectionproperties.analysis.cross_section import CrossSection

steel = Material(name='Steel', elastic_modulus=200e3, poissons_ratio=0.3,
                 yield_strength=500, color='grey')
timber = Material(name='Timber', elastic_modulus=8e3, poissons_ratio=0.35,
                  yield_strength=20, color='burlywood')

geometry = sections.RectangularSection(d=100, b=50)
geometry.control_points.insert(0, [-5, -5])
geometry.plot_geometry()

mesh = geometry.create_mesh(mesh_sizes=[1, 20])

section = CrossSection(geometry, mesh, [steel, timber])
section.plot_mesh(materials=True)
