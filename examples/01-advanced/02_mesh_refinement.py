r"""
.. _ref_ex_advanced2:

Mesh Refinement
---------------
Perform a mesh refinement study.

In this example the convergence of the torsion constant is investigated through an analysis of an
I Section. The mesh is refined both by modifying the mesh size and by specifying the number of
points making up the root radius. The figure below the example code shows that mesh refinement
adjacent to the root radius is a far more efficient method in obtaining fast convergence when
compared to reducing the mesh area size for the entire section.
"""

# sphinx_gallery_thumbnail_number = 1

import numpy as np
import matplotlib.pyplot as plt
import sectionproperties.pre.library.steel_sections as steel_sections
from sectionproperties.analysis.section import Section

# %%
# Define mesh sizes
mesh_size_list = [50, 20, 10, 5, 3, 2, 1]
nr_list = [4, 8, 12, 16, 20, 24, 32, 64]

# %%
# Initialise result lists
mesh_results = []
mesh_elements = []
nr_results = []
nr_elements = []

# %%
# Calculate reference solution
geometry = steel_sections.i_section(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=64)
geometry.create_mesh(mesh_sizes=[0.5])  # create mesh
section = Section(geometry)  # create a Section object
section.calculate_geometric_properties()
section.calculate_warping_properties()
j_reference = section.get_j()  # get the torsion constant

# %%
# Run through mesh_sizes with n_r = 16
for mesh_size in mesh_size_list:
    geometry = steel_sections.i_section(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=16)
    geometry.create_mesh(mesh_sizes=[mesh_size])  # create mesh
    section = Section(geometry)  # create a Section object
    section.calculate_geometric_properties()
    section.calculate_warping_properties()

    mesh_elements.append(len(section.elements))
    mesh_results.append(section.get_j())

# %%
# Run through n_r with mesh_size = 3
for n_r in nr_list:
    geometry = steel_sections.i_section(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=n_r)
    geometry.create_mesh(mesh_sizes=[3])  # create mesh
    section = Section(geometry)  # create a Section object
    section.calculate_geometric_properties()
    section.calculate_warping_properties()

    nr_elements.append(len(section.elements))
    nr_results.append(section.get_j())

# %%
# Convert results to a numpy array and compute the error
mesh_results = np.array(mesh_results)
nr_results = np.array(nr_results)
mesh_error_vals = (mesh_results - j_reference) / mesh_results * 100
nr_error_vals = (nr_results - j_reference) / nr_results * 100

# %%
# Plot the results
(fig, ax) = plt.subplots()
ax.loglog(mesh_elements, mesh_error_vals, "kx-", label="Mesh Size Refinement")
ax.loglog(nr_elements, nr_error_vals, "rx-", label="Root Radius Refinement")
plt.xlabel("Number of Elements")
plt.ylabel("Torsion Constant Error [%]")
plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
plt.tight_layout()
plt.show()
