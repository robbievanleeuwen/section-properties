import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection

# create a 150x100x6 RHS on its side
geometry = sections.Rhs(d=100, b=150, t=6, r_out=15, n_r=8)

# create a mesh with a maximum area of 2
mesh = geometry.create_mesh(mesh_sizes=[2])

# create a CrossSection object
section = CrossSection(geometry, mesh)

# perform a geometry and warping analysis
section.calculate_geometric_properties()
section.calculate_warping_properties()

# perform a stress analysis with Mx = 5 kN.m; Vx = 10 kN and Mzz = 3 kN.m
case1 = section.calculate_stress(Mxx=5e6, Vx=10e3, Mzz=3e6)

# perform a stress analysis with My = 15 kN.m; Vy = 30 kN and Mzz = 1.5 kN.m
case2 = section.calculate_stress(Myy=15e6, Vy=30e3, Mzz=1.5e6)

case1.plot_stress_m_zz(pause=False)  # plot the bending stress for case1
case1.plot_vector_mzz_zxy(pause=False)  # plot the torsion vectors for case1
case2.plot_stress_v_zxy(pause=False)  # plot the shear stress for case1
case1.plot_stress_vm(pause=False)  # plot the von mises stress for case1
case2.plot_stress_vm()  # plot the von mises stress for case2
