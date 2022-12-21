r"""
.. _ref-ex-trapz-approx:

Approximation of Torsion Constant for Trapezoidal Sections
----------------------------------------------------------

Trapezoidal elements or components of a cross-section are
quite common in bridge structures, either concrete or steel
composite construction. However, it's common to determine 
the torsion constant of the trapezoidal section by using a
rectangular approximation. For example, this is done in
the Autodesk Structural Bridge Design software when there
is a haunch in a `Steel Composite Beam
<https://knowledge.autodesk.com/support/structural-bridge-design/learn-explore/caas/CloudHelp/cloudhelp/ENU/ASBD-InProdAU/files/structure/tech-info/ti-torsion/ASBD-InProdAU-structure-tech-info-ti-torsion-Torsion-property-for-SAM-composite-beams-html-html.html>`_

The question then arises, when is it appropriate to make
the rectangular approximation to a trapezoidal section,
and what might the expected error be?
"""

#%%
# Define the Imports
# ==================
# Here we bring in the primative section shapes and also the more generic
# Shapely `Polygon` object.
import numpy as np
import matplotlib.pyplot as plt
from shapely import Polygon
import sectionproperties.pre.geometry as geometry
import sectionproperties.pre.pre as pre
import sectionproperties.pre.library.primitive_sections as sections
from sectionproperties.analysis.section import Section

#%%
# Define the Calculation Engine
# =============================
# It's better to collect the relevant section property calculation in a
# single function. We are only interested in the torsion constant, so this
# is straightforward enough. `geom` is the Section Property Geometry object;
# `ms` is the mesh size.
def get_section_j(geom, ms, plot_geom=False):
    geom.create_mesh(mesh_sizes=[ms])
    section = Section(geom)
    if plot_geom:
        section.plot_mesh()
    section.calculate_geometric_properties()
    section.calculate_warping_properties()
    return section.get_j()


#%%
# Define the Mesh Density
# =======================
# The number of elements per unit area is an important input to the calculations
# even though we are only examining ratios of the results. A nominal value of 100
# is reasonable.
n = 100  # mesh density

#%%
# Create and Analyse the Section
# ==============================
# This function accepts the width `b` and a slope `S` to create the trapezoid.
# Since we are only interested in relative results, the nominal dimensions are
# immaterial. There are a few ways to parametrize the problem, but it has been
# found that setting the middle height of trapezoid (i.e. the average height)
# to a unit value works fine.
def do_section(b, S, d_mid=1, plot_geom=False):
    delta = S * d_mid
    d1 = d_mid - delta
    d2 = d_mid + delta

    # compute mesh size
    ms = d_mid * b / n

    points = []
    points.append([0, 0])
    points.append([0, d1])
    points.append([b, d2])
    points.append([b, 0])
    if S < 1.0:
        trap_geom = geometry.Geometry(Polygon(points))
    else:
        trap_geom = sections.triangular_section(h=d2, b=b)
    jt = get_section_j(trap_geom, ms, plot_geom)

    rect_geom = sections.rectangular_section(d=(d1 + d2) / 2, b=b)
    jr = get_section_j(rect_geom, ms, plot_geom)
    return jt, jr, d1, d2


#%%
# Example Section
# ===============
# The analysis for a particular section looks as follows:
b, S = 4.0, 0.3
jt, jr, d1, d2 = do_section(b, S, plot_geom=True)
print(f"{b=}; {S=}; {jr=}; {jt=}; {jr/jt}")

#%%
# Create Loop Variables
# =====================
# The slope `S` is 0 for a rectangle, and 1 for a triangle and is
# defined per the plot below. A range of `S`, between 0.0 and 1.0 and
# a range of `b` are considered, between 1 and 10 here (but can be extended)
b_list = np.logspace(0, np.log10(10), 10)
S_list = np.linspace(0.0, 1.0, 10)
j_rect = np.zeros((len(b_list), len(S_list)))
j_trap = np.zeros((len(b_list), len(S_list)))

#%%
# The Main Loop
# =============
# Execute the double loop to get the ratios for combinations of `S` and `b`.
#
# An optional deugging line is left in for development but commented out.
for i, b in enumerate(b_list):
    for j, S in enumerate(S_list):
        jt, jr, d1, d2 = do_section(b, S)
        j_trap[i][j] = jt
        j_rect[i][j] = jr

        # print(f"{b=:.3}; {S=:.3}; {d1=:.5}; {d2=:.5}; J rect = {j_rect[i][j]:.5e}; J trap = {j_trap[i][j]:.5e}; J ratio = {j_rect[i][j]/j_trap[i][j]:.5}")
#%%
# Calculate the Ratios
# ====================
# Courtesy of numpy, this is easy:
j_ratio = j_rect / j_trap
#%%
# Plot the Results
# ================
# Here we highlight a few of the contours to illustrate the accuracy and behaviour
# of the approximation.
#
# As expected, when the section is rectangular, the error is small, but as it increases
# towards a triangle the accuracy generally reduces. However, there is an interesting
# line at an aspect ratio of about 2.7 where the rectangular approximation is always
# equal to the trapezoid's torsion constant.
levels = np.arange(start=0.5, stop=1.5, step=0.05)
plt.figure(figsize=(12, 6))
cs = plt.contour(
    S_list,
    b_list,
    j_ratio,
    levels=[0.95, 0.99, 1.00, 1.01, 1.05],
    colors=("k",),
    linestyles=(":",),
    linewidths=(1.2,),
)
plt.clabel(cs, colors="k", fontsize=10)
plt.contourf(S_list, b_list, j_ratio, 25, cmap="Wistia", levels=levels)
# plt.yscale('log')
minor_x_ticks = np.linspace(min(S_list), max(S_list), 10)
# plt.xticks(minor_x_ticks,minor=True)
plt.minorticks_on()
plt.grid(which="both", ls=":")
plt.xlabel(r"Slope $S = (d_2-d_1)/(d_2+d_1); d_2\geq d_1, d_1\geq 0$")
plt.ylabel("Aspect $b/d_{ave}; d_{ave} = (d_1 + d_2)/2$")
plt.colorbar()
plt.title(
    r"Accuracy of rectangular approximation to trapezoid torsion constant $J_{rect}\, /\, J_{trapz}$",
    multialignment="center",
)
plt.show()
