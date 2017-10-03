import triangle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import mesh2D
import time

crossSection = {}
crossSection['vertices'] = np.array([[0,0], [50,0], [50,100], [0,100], [6,6], [44, 6], [44, 94], [6, 94]])
crossSection['segments'] = np.array([[0,1], [1,2], [2,3], [3,0], [4,5], [5,6], [6,7], [7,4]])
crossSection['holes'] = np.array([[25,50]])

mesh1 = triangle.triangulate(crossSection, 'pa1q30') # generate triangular mesh
mesh = mesh2D.triMesh(mesh1, 'tri3') # create mesh object

# plot results
plt.figure()
plt.gca().set_aspect('equal')
plt.triplot(mesh1['vertices'][:,0], mesh1['vertices'][:,1], mesh1['triangles'], lw=0.5, color='black')
cmap = cm.get_cmap(name='jet')
v = np.linspace(-300, 300, 11, endpoint=True)
trictr = plt.tricontourf(mesh1['vertices'][:,0], mesh1['vertices'][:,1], mesh1['triangles'], mesh.w, cmap=cmap)
cbar = plt.colorbar(trictr)
d = 50
plt.plot([mesh.cx - d * np.cos(mesh.phi * np.pi / 180), mesh.cx + d * np.cos(mesh.phi * np.pi / 180)], [mesh.cy - d * np.sin(mesh.phi * np.pi / 180), mesh.cy + d * np.sin(mesh.phi * np.pi / 180)])
plt.plot([mesh.cx - d * np.cos(mesh.phi * np.pi / 180 + np.pi / 2), mesh.cx + d * np.cos(mesh.phi * np.pi / 180 + np.pi / 2)], [mesh.cy - d * np.sin(mesh.phi * np.pi / 180 + np.pi / 2), mesh.cy + d * np.sin(mesh.phi * np.pi / 180 + np.pi / 2)])
plt.show()

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
