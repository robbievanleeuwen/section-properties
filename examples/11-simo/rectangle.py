r"""

Rectangle
--------------

Calculate section properties of rectangle.
Mesh is refined until relative change of torsion and warping constants
is not more than rtol
"""
import math
import argparse
import sectionproperties.pre.library.primitive_sections as sections
from sectionproperties.analysis.section import Section

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-W","--width", help="width of square",
                    default=0.0032024,type=float)
parser.add_argument("-H","--height", help="height of square",
                    default=0.012377,type=float)
parser.add_argument("-R","--rtol", help="relative tolerance",
                    default=1e-3,type=float)
parser.add_argument("-M","--plot_mesh", help="Plot each mesh",
                    action="store_true")
parser.add_argument("-G","--plot_geometry", help="Plot geometry",
                    action="store_true")
args = parser.parse_args()
print("Rectangle: width = {0:.5g} and height = {1:.5g}, rtol={2:g}".
      format(args.width, args.height,args.rtol))
rtol=args.rtol
geometry = sections.rectangular_section(args.width, args.height)
if args.plot_geometry:
    geometry.plot_geometry()
a=geometry.calculate_area()
j0=a
iw0=a
ms=min(args.width,args.height)
vertices0=0 # sometimes requesting smaller mesh size generates same mesh
while True:
    ms=0.5*ms
    geometry.create_mesh(mesh_sizes=[ms])
    vertices=geometry.mesh.get('vertices').size
    if vertices0==vertices:
        ms=0.5*ms
        continue
    vertices0=vertices
    section = Section(geometry)
    if args.plot_mesh:
        section.plot_mesh()
    section.calculate_geometric_properties()
    section.calculate_warping_properties()
    j = section.get_j()
    if math.isnan(j):
        continue
    iw = section.get_gamma()
    jDiff=abs((j-j0)/j0)
    iwDiff=abs((iw-iw0)/iw0)
    print(("J = {0:.3g}, Iw = {1:.3g}, "+
          "meshSize = {2:.3g}, {5} nodes, {6} elements, "+
          "jDiff = {3:.3g}, iwDiff = {4:.3g}")
          .format(j,iw,ms,jDiff,iwDiff,
                  section.num_nodes,len(section.elements)))
    if(jDiff<rtol and iwDiff<rtol ):
        break
    else:
        j0=j
        iw0=iw
