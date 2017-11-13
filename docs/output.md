# Cross-section property output

The properties printed by the program are defined below:

```
-------------------------
Global xy Axis Properties
-------------------------
Area = Cros-section area
Qx = First moment of area about the global x-axis
Qy = First moment of area about the global x-axis
cx = x-position of the elastic centroid
cy = y-position of the elastic centroid
Ixx_g = Second moment of area about the global x-axis
Iyy_g = Second moment of area about the global y-axis
Ixy_g = Second moment of area about the global xy-axis

-----------------------------
Centroidal xy Axis Properties
-----------------------------
Ixx_c = Second moment of area about the centroidal x-axis
Iyy_c = Second moment of area about the centroidal y-axis
Ixy_c = Second moment of area about the centroidal xy-axis
Zxx+ = Elastic section modulus about the centroidal x-axis with respect to the top fibre
Zxx- = Elastic section modulus about the centroidal x-axis with respect to the bottom fibre
Zyy+ = Elastic section modulus about the centroidal y-axis with respect to the top fibre
Zyy- = Elastic section modulus about the centroidal y-axis with respect to the bottom fibre
rx_c = Radius of gyration about the centroidal x-axis
ry_c = Radius of gyration about the centroidal y-axis

-----------------------------
Principal Axis Properties
-----------------------------
phi = Principal bending axis rotation angle
I11_c = Second moment of area about the 11-axis
I22_c = Second moment of area about the 22-axis
Z11+ = Elastic section modulus about the 11-axis with respect to the top fibre
Z11- = Elastic section modulus about the 11-axis with respect to the bottom fibre
Z22+ = Elastic section modulus about the 22-axis with respect to the top fibre
Z22- = Elastic section modulus about the 22-axis with respect to the bottom fibre
r1_c = Radius of gyration about the centroidal 11-axis
r2_c = Radius of gyration about the centroidal 22-axis

-----------------------------
Plastic Properties
-----------------------------
x_pc = x-position of the plastic centroid with respect to the global axis (global axis bending)
y_pc = y-position of the plastic centroid with respect to the global axis (global axis bending)
Sxx = Plastic section modulus about the centroidal x-axis
Syy = Plastic section modulus about the centroidal y-axis
SF_xx+ = Shape factor for bending about the x-axis with respect to the top fibre
SF_xx- = Shape factor for bending about the x-axis with respect to the bottom fibre
SF_yy+ = Shape factor for bending about the y-axis with respect to the top fibre
SF_yy- = Shape factor for bending about the y-axis with respect to the bottom fibre
x_1_pc = x-position of the plastic centroid with respect to the global axis (principal axis bending)
y_2_pc = y-position of the plastic centroid with respect to the global axis (principal axis bending)
S11 = Plastic section modulus about the 11-axis
S22 = Plastic section modulus about the 22-axis
SF_11+ = Shape factor for bending about the 11-axis with respect to the top fibre
SF_11- = Shape factor for bending about the 11-axis with respect to the bottom fibre
SF_22+ = Shape factor for bending about the 22-axis with respect to the top fibre
SF_22- = Shape factor for bending about the 22-axis with respect to the bottom fibre

-----------------------------
Torsional Properties
-----------------------------
J = St. Venant torsion constant
Iw = Warping constant

-----------------------------
Shear Properties
-----------------------------
x_s,e = x-position of the shear centre (elastic calculation)
y_s,e = y-position of the shear centre (elastic calculation)
x_s,t = x-position of the shear centre (Trefftz's calculation)
y_s,t = y-position of the shear centre (Trefftz's calculation)
x_1_s,e = 11-position of the shear centre (elastic calculation)
y_2_s,e = 22-position of the shear centre (elastic calculation)
A_s,x = Shear area for transverse loading in the x-direction
A_s,y = Shear area for transverse loading in the y-direction
A_s,11 = Shear area for transverse loading in the 11-direction
A_s,22 = Shear area for transverse loading in the 22-direction
```
