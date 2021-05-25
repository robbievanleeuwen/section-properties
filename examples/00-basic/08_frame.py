"""
.. _ref_ex_frame:

Frame Analysis Example
----------------------

Analyse a frame.

The following example demonstrates how *sectionproperties* can be used to
calculate the cross-section properties required for a frame analysis. Using this
method is preferred over executing a geometric and warping analysis as only variables
required for a frame analysis are computed. In this example the torsion constant of
a rectangular section is calculated for a number of different mesh sizes and the
accuracy of the result compared with the time taken to obtain the solution
"""

# sphinx_gallery_thumbnail_number = 1

import time

import matplotlib.pyplot as plt
import numpy as np

import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection


# %%
# create a rectangular section
geometry = sections.RectangularSection(d=100, b=50)

# %%
# create a list of mesh sizes to analyse
mesh_sizes = [1.5, 2, 2.5, 3, 4, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100]
j_calc = []  # list to store torsion constants
t_calc = []  # list to store computation times

# %%
# loop through mesh sizes
for mesh_size in mesh_sizes:
    mesh = geometry.create_mesh(mesh_sizes=[mesh_size])
    section = CrossSection(geometry, mesh)

    start_time = time.time()
    (_, _, _, _, j, _) = section.calculate_frame_properties()
    t = time.time() - start_time
    t_calc.append(t)  # save the time
    j_calc.append(j)  # save the torsion constant

    # print the result
    msg = 'Mesh Size: {0}; '.format(mesh_size)
    msg += 'Solution Time {0:.5f} s; '.format(t)
    msg += 'Torsion Constant: {0:.12e}'.format(j)
    print(msg)

# %%
# Compute the error, assuming that the finest mesh (index 0) gives the 'correct' value
correct_val = j_calc[0]
j_np = np.array(j_calc)
error_vals = (j_calc - correct_val) / j_calc * 100

# %%
# produce a plot of the accuracy of the torsion constant with computation time
plt.loglog(t_calc[1:], error_vals[1:], 'kx-')
plt.xlabel('Solver Time [s]')
plt.ylabel('Torsion Constant Error [%]')
plt.show()
