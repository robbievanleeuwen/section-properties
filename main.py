import triangle
import triangle.plot as plot
import matplotlib.pyplot as plt
import numpy as np
import mesh2D

crossSection = {}
crossSection['vertices'] = np.array([[-10,0], [110,0], [100,10], [55,10], [55,90], [100,90], [110,100], [110,110], [-10,110], [-10,100], [0, 90], [45, 90], [45,10], [-10,10]])
crossSection['segments'] = np.array([[0,1], [1,2], [2,3], [3,4], [4,5], [5,6], [6,7], [7,8], [8,9], [9,10], [10,11], [11,12], [12,13], [13,0]])

ax1 = plt.subplot(121, aspect='equal')
plot.plot(ax1, **crossSection)

mesh1 = triangle.triangulate(crossSection, 'pa20')

ax2 = plt.subplot(122, sharex=ax1, sharey=ax1)
plot.plot(ax2, **mesh1)

# plt.show()
mesh = mesh2D.triMesh(mesh1)

print "Area = {}".format(mesh.area)
print "Qx = {}".format(mesh.Qx)
print "Qy = {}".format(mesh.Qy)
print "cx = {}".format(mesh.cx)
print "cy = {}".format(mesh.cy)
print "Ixx_g = {}".format(mesh.ixx_g)
print "Ixx_c = {}".format(mesh.ixx_c)
print "Iyy_g = {}".format(mesh.iyy_g)
print "Iyy_c = {}".format(mesh.iyy_c)
print "Ixy_g = {}".format(mesh.ixy_g)
print "Ixy_c = {}".format(mesh.ixy_c)
print "rx = {}".format(mesh.rx)
print "ry = {}".format(mesh.ry)
