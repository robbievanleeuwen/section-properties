"""
.. _ref_ex_nastran:

Creating a Nastran Section
--------------------------

Calculate section properties of Nastran HAT1 section.

The following example demonstrates how to create a cross-section defined in
a Nastran-based finite element analysis program. The following creates a
HAT1 cross-section and calculates the geometric, warping and plastic properties.
The HAT1 cross-section is meshed with a maximum elemental area of 0.005.

The geometry and mesh are plotted, and the mesh information printed to the terminal
before the analysis is carried out. Detailed time information is printed to the
terminal during the cross-section analysis stage. Once the analysis is complete,
the cross-section properties are printed to the terminal. The centroidal
axis second moments of area and torsion constant are saved to variables and it
is shown that, for non-circular sections, the torsion constant is not equal to the
sum of the second moments of area.
"""

# sphinx_gallery_thumbnail_number = 1

from typing import get_origin
import sectionproperties.pre.library.nastran_sections as nsections
from sectionproperties.analysis.section import Section

# %%
# Create a HAT1 section
geometry = nsections.nastran_hat1(DIM1=4.0, DIM2=2.0, DIM3=1.5, DIM4=0.1875, DIM5=0.375)
geometry.plot_geometry()  # plot the geometry
print(geometry.geom)

# %%
# Create a mesh with a maximum elemental area of 0.005
geometry.create_mesh(mesh_sizes=[0.005])

section = Section(geometry, time_info=True)  # create a Section object
section.display_mesh_info()  # display the mesh information
section.plot_mesh()  # plot the generated mesh`


# %%
# Perform a geometric, warping and plastic analysis, displaying the time info
section.calculate_geometric_properties()
section.calculate_warping_properties()
section.calculate_plastic_properties()

section.display_results()

# %%
# Get the second moments of area and the torsion constant
(ixx_c, iyy_c, ixy_c) = section.get_ic()
j = section.get_j()
print("Ixx + Iyy = {0:.3f}".format(ixx_c + iyy_c))
print("J = {0:.3f}".format(j))
