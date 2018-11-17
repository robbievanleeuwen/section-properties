import time
import numpy as np
import matplotlib.pyplot as plt
import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection

# create a rectangular section
geometry = sections.RectangularSection(d=100, b=50)

# create a list of mesh sizes to analyse
mesh_sizes = [1.5, 2, 2.5, 3, 4, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100]
j_calc = []  # list to store torsion constants
t_calc = []  # list to store computation times

# loop through mesh sizes
for mesh_size in mesh_sizes:
    mesh = geometry.create_mesh(mesh_sizes=[mesh_size])  # create mesh
    section = CrossSection(geometry, mesh)  # create a CrossSection object
    start_time = time.time()  # start timing
    # calculate the frame properties
    (_, _, _, _, j, _) = section.calculate_frame_properties()
    t = time.time() - start_time  # stop timing
    t_calc.append(t)  # save the time
    j_calc.append(j)  # save the torsion constant
    # print the result
    str = "Mesh Size: {0}; ".format(mesh_size)
    str += "Solution Time {0:.5f} s; ".format(t)
    str += "Torsion Constant: {0:.12e}".format(j)
    print(str)

correct_val = j_calc[0]  # assume the finest mesh gives the 'correct' value
j_np = np.array(j_calc)  # convert results to a numpy array
error_vals = (j_calc - correct_val) / j_calc * 100  # compute the error

# produce a plot of the accuracy of the torsion constant with computation time
plt.loglog(t_calc[1:], error_vals[1:], 'kx-')
plt.xlabel("Solver Time [s]")
plt.ylabel("Torsion Constant Error [%]")
plt.show()
