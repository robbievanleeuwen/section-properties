import numpy as np
import matplotlib.pyplot as plt
import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection

# rectangle dimensions
d_list = []
b_list = np.linspace(0.2, 1, 20)
j_list = []  # list holding torsion constant results

# number of elements for each analysis
n = 500

# loop through all the widths
for b in b_list:
    # calculate d assuming area = 1
    d = 1 / b
    d_list.append(d)

    # compute mesh size
    ms = d * b / n

    # perform a warping analysis on rectangle
    geometry = sections.RectangularSection(d=d, b=b)
    mesh = geometry.create_mesh(mesh_sizes=[ms])
    section = CrossSection(geometry, mesh)
    section.calculate_geometric_properties()
    section.calculate_warping_properties()

    # get the torsion constant
    j = section.get_j()
    print("d/b = {0:.3f}; J = {1:.5e}".format(d / b, j))
    j_list.append(j)

# plot the torsion constant as a function of the aspect ratio
(fig, ax) = plt.subplots()
ax.plot(np.array(d_list) / b_list, j_list, 'kx-')
ax.set_xlabel("Aspect Ratio [d/b]")
ax.set_ylabel("Torsion Constant [J]")
ax.set_title("Rectangular Section Torsion Constant")
plt.show()
