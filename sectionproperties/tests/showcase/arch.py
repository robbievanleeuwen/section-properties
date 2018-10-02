import numpy as np
import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection
from feastruct.fea.frame2d import Frame2D
from feastruct.solvers.linstatic import LinearStatic
from feastruct.solvers.linbuckling import LinearBuckling
from feastruct.post.post import PostProcessor


depths = [180, 190, 200, 210, 220]  # tee-section depths
n = 20  # number of elements
L = 10000  # length of arch
h = 3000  # height of arch
P = 1e5  # load on arch

for d in depths:
    # analyse cross-section
    geometry = sections.TeeSection(d=d, b=100, t_f=12, t_w=8, r=12, n_r=8)
    geometry.plot_geometry()
    mesh = geometry.create_mesh(mesh_sizes=[100])
    section = CrossSection(geometry, mesh)
    section.calculate_geometric_properties()
    area = section.get_area()
    (ixx, iyy, ixy) = section.get_ic()

    # new 2D frame analysis
    analysis = Frame2D()

    # create nodes
    for i in range(n + 1):
        x = i * L / n
        y = h * np.sin(np.pi * x / L)
        analysis.add_node(id=i+1, coord=[x, y])

    # create elements
    for i in range(n):
        analysis.add_element(id=i+1, node_ids=[i+1, i+2], el_type='EB2',
                             E=200e3, A=area, ixx=ixx, rho=7.85e-9)

    # apply boundary conditions
    fc = analysis.add_freedom_case(id=1)  # apply pinned supports on both ends
    fc.add_nodal_support(node_id=1, val=0, dir=1)
    fc.add_nodal_support(node_id=1, val=0, dir=2)
    fc.add_nodal_support(node_id=n+1, val=0, dir=1)
    fc.add_nodal_support(node_id=n+1, val=0, dir=2)
    lc = analysis.add_load_case(id=1)  # apply load at centre of arch
    lc.add_nodal_load(node_id=n/2+1, val=-P, dir=2)
    analysis.add_analysis_case(id=1, fc_id=1, lc_id=1)

    # perform linear static analysis
    LinearStatic(analysis, case_ids=[1]).solve()

    # perform buckling analysis
    LinearBuckling(analysis, case_id=1).solve()

    # display results
    post = PostProcessor(analysis, n_subdiv=5)
    post.plot_buckling_eigenvector(case_id=1, buckling_mode=1)

    # print results
    lf = analysis.nodes[0].get_buckling_results(1, 1)[0]
    print("Depth: {0}; Buckling Load Factor: {1:0.4f}".format(d, lf))
