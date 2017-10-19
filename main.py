import meshpy.triangle as triangle
import numpy as np
import mesh2D
import matplotlib.pyplot as plt

# rectangular hollow section
points = [(0,0), (50,0), (50,100), (0,100), (6,6), (44, 6), (44, 94), (6, 94)]
facets = [(0,1), (1,2), (2,3), (3,0), (4,5), (5,6), (6,7), (7,4)]
info = triangle.MeshInfo()
info.set_points(points)
info.set_holes([(25, 50)])
info.set_facets(facets)

triangularMesh = triangle.build(info, max_volume = 2.5, min_angle = 30) # generate triangular mesh

mesh = mesh2D.triMesh(triangularMesh, 'tri3', 0.2) # create mesh object

mesh.contourPlot(principalAxis = True, z = mesh.tau_torsion, nodes = True)
# mesh.quiverPlot(mesh.tau_zx_torsion, mesh.tau_zy_torsion)
# mesh.contourPlot(False, mesh.Psi)
# mesh.contourPlot(False, mesh.tau_zy_shear)
# mesh.contourPlot(False, mesh.tau_shear)

print "-------------------------"
print "Global xy Axis Properties"
print "-------------------------"
print "Area = {}".format(mesh.area)
print "Qx = {}".format(mesh.Qx)
print "Qy = {}".format(mesh.Qy)
print "cx = {}".format(mesh.cx)
print "cy = {}".format(mesh.cy)
print "Ixx_g = {}".format(mesh.ixx_g)
print "Iyy_g = {}".format(mesh.iyy_g)
print "Ixy_g = {}".format(mesh.ixy_g)
print ""
print "-----------------------------"
print "Centroidal xy Axis Properties"
print "-----------------------------"
print "Ixx_c = {}".format(mesh.ixx_c)
print "Iyy_c = {}".format(mesh.iyy_c)
print "Ixy_c = {}".format(mesh.ixy_c)
print "Zxx_plus = {}".format(mesh.zxx_plus)
print "Zxx_minus = {}".format(mesh.zxx_minus)
print "Zyy_plus = {}".format(mesh.zyy_plus)
print "Zyy_minus = {}".format(mesh.zyy_minus)
print "rx_c = {}".format(mesh.rx_c)
print "ry_c = {}".format(mesh.ry_c)
print ""
print "-----------------------------"
print "Principal Axis Properties"
print "-----------------------------"
print "phi = {}".format(mesh.phi)
print "I11_c = {}".format(mesh.i11_c)
print "I22_c = {}".format(mesh.i22_c)
print "Z11_plus = {}".format(mesh.z11_plus)
print "Z11_minus = {}".format(mesh.z11_minus)
print "Z22_plus = {}".format(mesh.z22_plus)
print "Z22_minus = {}".format(mesh.z22_minus)
print "r1_c = {}".format(mesh.r1_c)
print "r2_c = {}".format(mesh.r2_c)
print ""
print "-----------------------------"
print "Torsional Properties"
print "-----------------------------"
print "J = {}".format(mesh.J)
