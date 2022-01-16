r"""
.. _ref_ex_advanced1:

Advanced Plotting
-----------------
Harness some of the plotting features provided by :func:`~sectionproperties.post.post.plotting_context`.
"""

# sphinx_gallery_thumbnail_number = 1

import sectionproperties.pre.library.steel_sections as steel_sections
from sectionproperties.analysis.section import Section
import matplotlib.pyplot as plt

# %%
# The below example creates a 100x6 SHS and plots the geometry, mesh, centroids and torsion vectors
# in a 2x2 subplot arrangement.
geometry = steel_sections.rectangular_hollow_section(d=100, b=100, t=6, r_out=15, n_r=8)

# Plot the geometry
ax = geometry.plot_geometry(nrows=2, ncols=2, figsize=(12, 7), render=False, labels=[])
fig = ax.get_figure()  # get the figure

# Create a mesh and section object, for the mesh, use a maximum area of 2
geometry.create_mesh(mesh_sizes=[2])
section = Section(geometry)
section.plot_mesh(ax=fig.axes[1], materials=False)  # plot the mesh

# Perform a geometry and warping analysis
section.calculate_geometric_properties()
section.calculate_warping_properties()
section.plot_centroids(ax=fig.axes[2])  # plot the cetnroids

# Perform a stress analysis with Mzz = 10 kN.m
stress = section.calculate_stress(Mzz=10e6)
stress.plot_vector_mzz_zxy(
    ax=fig.axes[3], title="Torsion Vectors"
)  # plot the torsion vectors
plt.show()  # show the plot
