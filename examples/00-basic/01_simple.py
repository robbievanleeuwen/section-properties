r"""
.. _ref_ex_simple:

Simple Example
--------------

Calculate section properties of a circle.

The following example calculates the geometric, warping and plastic properties
of a 50 mm diameter circle. The circle is discretised with 64 points and a mesh
size of 2.5 mm\ :sup:`2`.

The geometry and mesh are plotted, and the mesh information printed to the terminal
before the analysis is carried out. Detailed time information is printed to the
terminal during the cross-section analysis stage. Once the analysis is complete,
the cross-section properties are printed to the terminal. The centroidal
axis second moments of area and torsion constant are saved to variables and it
is shown that, for a circle, the torsion constant is equal to the sum of the
second moments of area.
"""

# sphinx_gallery_thumbnail_number = 1
import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection


# %%
# Create a 50 diameter circle discretised by 64 points
geometry = sections.CircularSection(d=50, n=64)
geometry.plot_geometry()

# %%
# Create a mesh with a mesh size of 2.5 and display information about it
mesh = geometry.create_mesh(mesh_sizes=[2.5])

section = CrossSection(geometry, mesh)
section.display_mesh_info()
section.plot_mesh()

# %%
# Perform a geometric, warping and plastic analysis
section.calculate_geometric_properties(time_info=True)
section.calculate_warping_properties(time_info=True)
section.calculate_plastic_properties(time_info=True)

# %%
# Print the results to the terminal
section.display_results()

# %%
# Get and print the second moments of area and the torsion constant
(ixx_c, iyy_c, ixy_c) = section.get_ic()
j = section.get_j()
print('Ixx + Iyy = {0:.3f}'.format(ixx_c + iyy_c))
print('J = {0:.3f}'.format(j))
