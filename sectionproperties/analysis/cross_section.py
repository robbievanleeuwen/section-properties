import numpy as np
from scipy.sparse import csc_matrix, coo_matrix, linalg
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sectionproperties.analysis.fea as fea
import sectionproperties.analysis.solver as solver
import sectionproperties.post.post as post

# TODO: fix torsion stress - or related issues???


class CrossSection:
    """Class for structural cross-sections.

    Stores the finite element geometry and mesh and provides methods to compute
    the cross-section properties. The element type used in this program is the
    six-noded quadratic triangular element.

    The constructor extracts information from the provided mesh object and
    creates and stores corresponding tri-6 finite element objects.

    :param geometry: Cross-section geometry object used to generate the mesh
    :type geometry: :class:`~sectionproperties.pre.sections.Geometry`
    :param mesh: Mesh object returned by meshpy
    :type mesh: :class:`meshpy.triangle.MeshInfo`

    The following example creates a
    :class:`~sectionproperties.analysis.cross_section.CrossSection` object of a
    100D x 50W rectangle using a mesh size of 5::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.RectangularSection(d=100, b=50)
            mesh = geometry.create_mesh(mesh_sizes=[5])
            section = CrossSection(geometry, mesh)

    :cvar elements: List of finite element objects describing the cross-section
        mesh
    :vartype elements: list[:class:`~sectionproperties.analysis.fea.Tri6`]
    :cvar int num_nodes: Number of nodes in the finite element mesh
    :cvar geometry: Cross-section geometry object used to generate the mesh
    :vartype geometry: :class:`~sectionproperties.pre.sections.Geometry`
    :cvar mesh: Mesh object returned by meshpy
    :vartype mesh: :class:`meshpy.triangle.MeshInfo`
    :cvar mesh_nodes: Array of node coordinates from the mesh
    :vartype mesh_nodes: :class:`numpy.ndarray`
    :cvar mesh_elements: Array of connectivities from the mesh
    :vartype mesh_elements: :class:`numpy.ndarray`
    :cvar mesh_attributes: Array of attributes from the mesh
    :vartype mesh_attributes: :class:`numpy.ndarray`
    :cvar section_props: Class to store calculated section properties
    :vartype section_props:
        :class:`~sectionproperties.analysis.cross_section.SectionProperties`
    """

    def __init__(self, geometry, mesh):
        """Inits the CrossSection class."""

        self.geometry = geometry  # save geometry data

        # extract mesh data
        nodes = np.array(mesh.points, dtype=np.dtype(float))
        elements = np.array(mesh.elements, dtype=np.dtype(int))
        attributes = np.array(mesh.element_attributes, dtype=np.dtype(int))

        # swap mid-node order to retain node ordering consistency
        elements[:, [3, 4, 5]] = elements[:, [5, 3, 4]]

        self.elements = []  # initialise list holding all element objects

        # build the mesh one element at a time
        for (i, node_ids) in enumerate(elements):
            x1 = nodes[node_ids[0]][0]
            y1 = nodes[node_ids[0]][1]
            x2 = nodes[node_ids[1]][0]
            y2 = nodes[node_ids[1]][1]
            x3 = nodes[node_ids[2]][0]
            y3 = nodes[node_ids[2]][1]
            x4 = nodes[node_ids[3]][0]
            y4 = nodes[node_ids[3]][1]
            x5 = nodes[node_ids[4]][0]
            y5 = nodes[node_ids[4]][1]
            x6 = nodes[node_ids[5]][0]
            y6 = nodes[node_ids[5]][1]

            # create a list containing the vertex and mid-node coordinates
            coords = np.array(
                [[x1, x2, x3, x4, x5, x6], [y1, y2, y3, y4, y5, y6]])

            # add a tri6 element to the mesh
            self.elements.append(fea.Tri6(coords, node_ids, attributes[i]))

        # save total number of nodes in mesh
        self.num_nodes = len(nodes)

        # save mesh input
        self.mesh = mesh
        self.mesh_nodes = nodes
        self.mesh_elements = elements
        self.mesh_attributes = attributes

        # initialise class storing section properties
        self.section_props = SectionProperties()

    def calculate_geometric_properties(self, time_info=False):
        """Calculates all the geometric properties of the cross-section and
        stores them in the
        :class:`~sectionproperties.analysis.cross_section.SectionProperties`
        object contained in section_props.

        :param bool time_info: If set to True, a detailed description of the
            computation and the time cost is printed to the terminal.

        The following section properties are calculated:

        * Area
        * First moments of area
        * Second moments of area about the global axis
        * Second moments of area about the centroidal axis
        * Centroidal section moduli
        * Radii of gyration
        * Principal axis properties

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
        """

        def calculate_geom():
            self.calculate_area()
            self.calculate_qx()
            self.calculate_qy()
            self.calculate_ixx_g()
            self.calculate_iyy_g()
            self.calculate_ixy_g()
            self.section_props.calculate_elastic_centroid()
            self.section_props.calculate_centroidal_properties(self.mesh)

        if time_info:
            text = "--Calculating geometric section properties..."
            solver.function_timer(text, calculate_geom)
            print("")
        else:
            calculate_geom()

    def calculate_warping_properties(self, time_info=False,
                                     solver_type='direct'):
        """Calculates all the warping properties of the cross-section and
        stores them in the
        :class:`~sectionproperties.analysis.cross_section.SectionProperties`
        object contained in section_props.

        :param bool time_info: If set to True, a detailed description of the
            computation and the time cost is printed to the terminal.
        :param string solver_type: Solver used for solving systems of linear
            equations, either using the *'direct'* method or *'cgs'* iterative
            method

        The following section properties are calculated:

        * Torsion constant
        * Shear centre
        * Shear area
        * Warping constant

        Note that the geometric properties must be calculated first for the
        calculation of the warping properties to be correct::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()

        :raises RuntimeError: If the geometric properties have not been
            calculated prior to calling this method
        """

        # check that a geometric analysis has been performed
        if None in [self.section_props.area, self.section_props.ixx_c,
                    self.section_props.cx]:
            err = "Cacluate geometric properties before "
            err += "performing warping analysis."
            raise RuntimeError(err)

        # create a new CrossSection with the origin shifted to the centroid for
        # calculation of the warping properties such that the Lagrangian
        # multiplier approach can be utilised
        warping_section = CrossSection(self.geometry, self.mesh)

        # shift the coordinates of each element
        # N.B. the mesh class attribute remains unshifted!
        for el in warping_section.elements:
            el.coords[0, :] -= self.section_props.cx
            el.coords[1, :] -= self.section_props.cy

        # shift the mesh_nodes
        warping_section.mesh_nodes[:, 0] -= self.section_props.cx
        warping_section.mesh_nodes[:, 1] -= self.section_props.cy

        # assemble stiffness matrix and load vector for warping function
        if time_info:
            text = "--Assembing {0}x{0} stiffness".format(self.num_nodes)
            text += " matrix and load vector..."
            (k, k_lg, f_torsion) = solver.function_timer(
                text, warping_section.assemble_torsion)
        else:
            (k, k_lg, f_torsion) = warping_section.assemble_torsion()

        # ILU decomposition of stiffness matrices
        def ilu_decomp():
            # ILU decomposition on regular stiffness matrix
            k_precond = linalg.LinearOperator(
                (self.num_nodes, self.num_nodes), linalg.spilu(k).solve)

            # ILU decomposition on Lagrangian stiffness matrix
            k_lg_precond = linalg.LinearOperator(
                (self.num_nodes + 1, self.num_nodes + 1),
                linalg.spilu(k_lg).solve)

            return (k_precond, k_lg_precond)

        if solver_type == 'cgs':
            if time_info:
                text = "--Performing ILU decomposition on stiffness"
                text += "matrices..."
                (k_precond, k_lg_precond) = solver.function_timer(
                    text, ilu_decomp)
            else:
                (k_precond, k_lg_precond) = ilu_decomp()

        # solve for warping constant
        def solve_warping():
            if solver_type == 'cgs':
                omega = solver.solve_cgs(k, f_torsion, k_precond)
            elif solver_type == 'direct':
                omega = solver.solve_direct(k, f_torsion)

            return omega

        if time_info:
            text = "--Solving for the warping constant using the "
            text += "{0} solver...".format(solver_type)
            omega = solver.function_timer(text, solve_warping)
        else:
            omega = solve_warping()

        self.section_props.omega = omega

        # determine the torsion constant
        def j_func():
            return (self.section_props.ixx_c + self.section_props.iyy_c -
                    omega.dot(k.dot(np.transpose(omega))))

        if time_info:
            text = "--Computing the torsion constant..."
            self.section_props.j = solver.function_timer(text, j_func)
        else:
            self.section_props.j = j_func()

        # assemble shear function load vectors
        def assemble_shear_load():
            f_psi = np.zeros(self.num_nodes)
            f_phi = np.zeros(self.num_nodes)

            for el in warping_section.elements:
                (f_psi_el, f_phi_el) = el.shear_load_vectors(
                    self.section_props.ixx_c, self.section_props.iyy_c,
                    self.section_props.ixy_c)
                f_psi[el.node_ids] += f_psi_el
                f_phi[el.node_ids] += f_phi_el

            return (f_psi, f_phi)

        if time_info:
            text = "--Assembling shear function load vectors..."
            (f_psi, f_phi) = solver.function_timer(text, assemble_shear_load)
        else:
            (f_psi, f_phi) = assemble_shear_load()

        # solve for shear functions psi and phi
        def solve_shear_functions():
            if solver_type == 'cgs':
                psi_shear = solver.solve_cgs_lagrange(k_lg, f_psi,
                                                      m=k_lg_precond)
                phi_shear = solver.solve_cgs_lagrange(k_lg, f_phi,
                                                      m=k_lg_precond)
            elif solver_type == 'direct':
                psi_shear = solver.solve_direct_lagrange(k_lg, f_psi)
                phi_shear = solver.solve_direct_lagrange(k_lg, f_phi)

            return (psi_shear, phi_shear)

        if time_info:
            text = "--Solving for the shear functions using the "
            text += "{0} solver...".format(solver_type)
            (psi_shear, phi_shear) = solver.function_timer(
                text, solve_shear_functions)
        else:
            (psi_shear, phi_shear) = solve_shear_functions()

        self.section_props.psi_shear = psi_shear
        self.section_props.phi_shear = psi_shear

        # assemble shear centre and warping moment integrals
        def assemle_sc_warping_integrals():
            sc_xint = 0
            sc_yint = 0
            q_omega = 0
            i_omega = 0
            i_xomega = 0
            i_yomega = 0

            for el in warping_section.elements:
                (sc_xint_el, sc_yint_el, q_omega_el, i_omega_el, i_xomega_el,
                 i_yomega_el) = el.shear_warping_integrals(
                    self.section_props.ixx_c, self.section_props.iyy_c,
                    self.section_props.ixy_c, omega[el.node_ids])

                sc_xint += sc_xint_el
                sc_yint += sc_yint_el
                q_omega += q_omega_el
                i_omega += i_omega_el
                i_xomega += i_xomega_el
                i_yomega += i_yomega_el

            return (sc_xint, sc_yint, q_omega, i_omega, i_xomega, i_yomega)

        if time_info:
            text = "--Assembling shear centre and warping moment integrals..."
            (sc_xint, sc_yint, q_omega, i_omega, i_xomega, i_yomega) = (
                solver.function_timer(text, assemle_sc_warping_integrals))
        else:
            (sc_xint, sc_yint, q_omega, i_omega, i_xomega, i_yomega) = (
                assemle_sc_warping_integrals())

        # calculate effective Poissons ratio
        nu = 0  # TODO: implement

        # calculate shear centres
        def shear_centres():
            # calculate shear centres (elasticity approach)
            Delta_s = 2 * (1 + nu) * (
                self.section_props.ixx_c * self.section_props.iyy_c -
                self.section_props.ixy_c ** 2)
            x_se = (1 / Delta_s) * ((nu / 2 * sc_xint) - f_torsion.dot(
                phi_shear))
            y_se = (1 / Delta_s) * ((nu / 2 * sc_yint) + f_torsion.dot(
                psi_shear))
            (x1_se, y2_se) = fea.principal_coordinate(self.section_props.phi,
                                                      x_se, y_se)

            # calculate shear centres (Trefftz's approach)
            x_st = (self.section_props.ixy_c *
                    i_xomega - self.section_props.iyy_c * i_yomega) / (
                self.section_props.ixx_c * self.section_props.iyy_c -
                self.section_props.ixy_c ** 2)
            y_st = (self.section_props.ixx_c *
                    i_xomega - self.section_props.ixy_c * i_yomega) / (
                self.section_props.ixx_c * self.section_props.iyy_c -
                self.section_props.ixy_c ** 2)

            return (Delta_s, x_se, y_se, x1_se, y2_se, x_st, y_st)

        if time_info:
            text = "--Calculating shear centres..."
            (Delta_s, x_se, y_se, x1_se, y2_se, x_st, y_st) = (
                solver.function_timer(text, shear_centres))
        else:
            (Delta_s, x_se, y_se, x1_se, y2_se, x_st, y_st) = shear_centres()

        # save shear centres
        self.section_props.Delta_s = Delta_s
        self.section_props.x_se = x_se
        self.section_props.y_se = y_se
        self.section_props.x1_se = x1_se
        self.section_props.y2_se = y2_se
        self.section_props.x_st = x_st
        self.section_props.y_st = y_st

        # calculate warping constant
        self.section_props.gamma = (
            i_omega - q_omega ** 2 / self.section_props.area -
            y_se * i_xomega + x_se * i_yomega)

        def assemble_shear_deformation():
            # assemble shear deformation coefficients
            kappa_x = 0
            kappa_y = 0
            kappa_xy = 0

            for el in warping_section.elements:
                (kappa_x_el, kappa_y_el, kappa_xy_el) = el.shear_coefficients(
                    self.section_props.ixx_c, self.section_props.iyy_c,
                    self.section_props.ixy_c, psi_shear[el.node_ids],
                    phi_shear[el.node_ids])
                kappa_x += kappa_x_el
                kappa_y += kappa_y_el
                kappa_xy += kappa_xy_el

            return (kappa_x, kappa_y, kappa_xy)

        if time_info:
            text = "--Assembling shear deformation coefficients..."
            (kappa_x, kappa_y, kappa_xy) = (
                solver.function_timer(text, assemble_shear_deformation))
            print("")
        else:
            (kappa_x, kappa_y, kappa_xy) = assemble_shear_deformation()

        # calculate shear areas wrt global axis
        self.section_props.A_sx = Delta_s ** 2 / kappa_x
        self.section_props.A_sy = Delta_s ** 2 / kappa_y
        self.section_props.A_sxy = Delta_s ** 2 / kappa_xy

        # calculate shear areas wrt principal bending axis:
        alpha_xx = kappa_x * self.section_props.area / Delta_s ** 2
        alpha_yy = kappa_y * self.section_props.area / Delta_s ** 2
        alpha_xy = kappa_xy * self.section_props.area / Delta_s ** 2

        # rotate the tensor by the principal axis angle
        phi_rad = self.section_props.phi * np.pi / 180
        R = (np.array([[np.cos(phi_rad),  np.sin(phi_rad)],
                       [-np.sin(phi_rad), np.cos(phi_rad)]]))

        rotatedAlpha = R.dot(np.array(
            [[alpha_xx, alpha_xy],
             [alpha_xy, alpha_yy]])).dot(np.transpose(R))

        # recalculate the shear area based on the rotated alpha value
        self.section_props.A_s11 = self.section_props.area / rotatedAlpha[0, 0]
        self.section_props.A_s22 = self.section_props.area / rotatedAlpha[1, 1]

    def calculate_plastic_properties(self, time_info=False):
        """s"""

        # TODO: implement

    def calculate_stress(self, N=0, Vx=0, Vy=0, Mxx=0, Myy=0, M11=0, M22=0,
                         Mzz=0, time_info=False):
        """Calculates the cross-section stress resulting from design actions
        and returns a
        :class:`~sectionproperties.analysis.cross_section.StressResult` object
        containing the cross-section stresses.

        :param float N: Axial force
        :param float Vx: Shear force acting in the x-direction
        :param float Vy: Shear force acting in the y-direction
        :param float Mxx: Bending moment about the centroidal xx-axis
        :param float Myy: Bending moment about the centroidal yy-axis
        :param float M11: Bending moment about the centroidal 11-axis
        :param float M22: Bending moment about the centroidal 22-axis
        :param float Mzz: Torsion moment about the centroidal zz-axis
        :param bool time_info: If set to True, a detailed description of the
            computation and the time cost is printed to the terminal.
        :return: Object containing the cross-section stresses
        :rtype: :class:`~sectionproperties.analysis.cross_section.StressResult`

        Note that a geometric and warping analysis must be performed before a
        stress analysis is carried out::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_result = section.calculate_stress(N=1e3, Vy=3e3, Mxx=1e6)

        :raises RuntimeError: If a geometric and warping analysis have not been
            performed prior to calling this method
        """

        # check that a geometric and warping analysis has been performed
        if None in [self.section_props.area, self.section_props.ixx_c,
                    self.section_props.cx, self.section_props.j]:
            err = "Perform a geometric and warping analysis before "
            err += "carrying out a stress analysis."
            raise RuntimeError(err)

        def calc_stress():
            # create stress result object
            stress_result = StressResult(self, self.num_nodes)

            # allocate nodal weights vector for nodal averaging
            nodal_weights = np.zeros(self.num_nodes)

            # loop through all elements in the mesh
            for el in self.elements:
                # evaluate stresses at nodes of the current element
                (sig_zz_n_el, sig_zz_mxx_el, sig_zz_myy_el, sig_zz_m11_el,
                    sig_zz_m22_el, sig_zx_mzz_el, sig_zy_mzz_el, sig_zx_vx_el,
                    sig_zy_vx_el, sig_zx_vy_el, sig_zy_vy_el,
                    weights) = el.element_stress(
                    self.section_props.area, self.section_props.cx,
                    self.section_props.cy, self.section_props.ixx_c,
                    self.section_props.iyy_c, self.section_props.ixy_c,
                    self.section_props.i11_c, self.section_props.i22_c,
                    self.section_props.phi, self.section_props.j,
                    self.section_props.omega[el.node_ids],
                    self.section_props.psi_shear[el.node_ids],
                    self.section_props.phi_shear[el.node_ids],
                    self.section_props.Delta_s)

                # add stresses to global vectors
                stress_result.sig_zz_n[el.node_ids] += (
                    N * sig_zz_n_el[:, 0] * weights)
                stress_result.sig_zz_mxx[el.node_ids] += (
                    Mxx * sig_zz_mxx_el * weights)
                stress_result.sig_zz_myy[el.node_ids] += (
                    Myy * sig_zz_myy_el * weights)
                stress_result.sig_zz_m11[el.node_ids] += (
                    M11 * sig_zz_m11_el * weights)
                stress_result.sig_zz_m22[el.node_ids] += (
                    M22 * sig_zz_m22_el * weights)
                stress_result.sig_zx_mzz[el.node_ids] += (
                    Mzz * sig_zx_mzz_el * weights)
                stress_result.sig_zy_mzz[el.node_ids] += (
                    Mzz * sig_zy_mzz_el * weights)
                stress_result.sig_zx_vx[el.node_ids] += (
                    Vx * sig_zx_vx_el * weights)
                stress_result.sig_zy_vx[el.node_ids] += (
                    Vx * sig_zy_vx_el * weights)
                stress_result.sig_zx_vy[el.node_ids] += (
                    Vy * sig_zx_vy_el * weights)
                stress_result.sig_zy_vy[el.node_ids] += (
                    Vy * sig_zy_vy_el * weights)

                # add nodal weights
                nodal_weights[el.node_ids] += weights

            # nodal averaging
            stress_result.sig_zz_n *= 1 / nodal_weights
            stress_result.sig_zz_mxx *= 1 / nodal_weights
            stress_result.sig_zz_myy *= 1 / nodal_weights
            stress_result.sig_zz_m11 *= 1 / nodal_weights
            stress_result.sig_zz_m22 *= 1 / nodal_weights
            stress_result.sig_zx_mzz *= 1 / nodal_weights
            stress_result.sig_zy_mzz *= 1 / nodal_weights
            stress_result.sig_zx_vx *= 1 / nodal_weights
            stress_result.sig_zy_vx *= 1 / nodal_weights
            stress_result.sig_zx_vy *= 1 / nodal_weights
            stress_result.sig_zy_vy *= 1 / nodal_weights

            # calculate combined stresses
            stress_result.calculate_combined_stresses()

            return stress_result

        if time_info:
            text = "--Calculating cross-section stresses..."
            stress_result = solver.function_timer(text, calc_stress)
            print("")
        else:
            stress_result = calc_stress()

        # return the stress result object
        return stress_result

    def calculate_area(self):
        """Calculates the area of the cross-section."""

        self.section_props.area = 0

        for el in self.elements:
            self.section_props.area += el.area

    def calculate_qx(self):
        """Calculates the first moment area of the cross-section about the
        x-axis.
        """

        self.section_props.qx = 0

        for el in self.elements:
            self.section_props.qx += el.qx

    def calculate_qy(self):
        """Calculates the first moment area of the cross-section about the
        y-axis.
        """

        self.section_props.qy = 0

        for el in self.elements:
            self.section_props.qy += el.qy

    def calculate_ixx_g(self):
        """Calculates the second moment area of the cross-section about the
        global x-axis.
        """

        self.section_props.ixx_g = 0

        for el in self.elements:
            self.section_props.ixx_g += el.ixx

    def calculate_iyy_g(self):
        """Calculates the second moment area of the cross-section about the
        global y-axis.
        """

        self.section_props.iyy_g = 0

        for el in self.elements:
            self.section_props.iyy_g += el.iyy

    def calculate_ixy_g(self):
        """Calculates the second moment area of the cross-section about the
        global xy-axis.
        """

        self.section_props.ixy_g = 0

        for el in self.elements:
            self.section_props.ixy_g += el.ixy

    def calculate_elastic_centroid(self):
        """Calculates the elastic centroid of the cross-section.

        :return: Tuple containing the coordinates of the elastic centroid.
        :rtype: tuple(float, float)
        """

        if self.section_props.area is None:
            self.calculate_area()

        if self.section_props.qx is None:
            self.calculate_qx()

        if self.section_props.qy is None:
            self.calculate_qy()

        return(self.section_props.calculate_elastic_centroid())

    def assemble_torsion(self):
        """Assembles stiffness matrices to be used for the computation of
        warping properties and the torsion load vector (f_torsion). Both a
        regular (k) and Lagrangian multiplier (k_lg) stiffness matrix are
        returned. The stiffness matrices are assembled using the sparse COO
        format and returned in the sparse CSC format.

        :return: Regular stiffness matrix, Lagrangian multiplier stiffness
            matrix and torsion load vector *(k, k_lg, f_torsion)*
        :rtype: tuple(:class:`scipy.sparse.csc_matrix`,
            :class:`scipy.sparse.csc_matrix`, :class:`numpy.ndarray`)
        """

        # initialise variables
        N = self.num_nodes  # size of matrix
        row = []  # list holding row indices
        col = []  # list holding column indices
        data = []  # list holding stiffness matrix entries
        f_torsion = np.zeros(N)  # force vector array

        # loop through all elements in the mesh
        for el in self.elements:
            # determine number of nodes in the current element
            n = len(el.node_ids)

            # calculate the element stiffness matrix and torsion load vector
            (k_el, f_el) = el.torsion_properties()

            # assemble the torsion load vector
            f_torsion[el.node_ids] += f_el

            # create row index vector
            r = np.repeat(el.node_ids, n)

            # create column index vector
            c = np.tile(el.node_ids, n)

            # flatten element stiffness matrix
            k = k_el.flatten()

            # add to global arrays
            row = np.hstack((row, r))
            col = np.hstack((col, c))
            data = np.hstack((data, k))

        k = coo_matrix((data, (row, col)), shape=(N, N))

        # construct Lagrangian multiplier matrix:
        # column vector of ones
        row = np.hstack((row, range(N)))
        col = np.hstack((col, np.repeat(N, N)))
        data = np.hstack((data, np.repeat(1, N)))

        # row vector of ones
        row = np.hstack((row, np.repeat(N, N)))
        col = np.hstack((col, range(N)))
        data = np.hstack((data, np.repeat(1, N)))

        # zero in bottom right corner
        row = np.hstack((row, N))
        col = np.hstack((col, N))
        data = np.hstack((data, 0))

        k_lg = coo_matrix((data, (row, col)), shape=(N+1, N+1))

        return (csc_matrix(k), csc_matrix(k_lg), f_torsion)

    def plot_mesh(self, ax=None, pause=True, alpha=1):
        """Plots the finite element mesh. If no axes object is supplied a new
        figure and axis is created.

        :param ax: Axes object on which the mesh is plotted
        :type ax: :class:`matplotlib.axes.Axes`
        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        :param float alpha: Transparency of the mesh: 0 <= alpha <= 1

        The following example plots the mesh generated for a 50D x 100W
        rectangle using a mesh size of 5::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.RectangularSection(d=50, b=100)
            mesh = geometry.create_mesh(mesh_sizes=[5])
            section = CrossSection(geometry, mesh)
            section.plot_mesh()

        ..  figure:: ../images/rectangle_mesh.png
            :align: center

            Finite element mesh generated by the above example.
        """

        # if no axes object is supplied, create and setup the plot
        if ax is None:
            ax_supplied = False
            (fig, ax) = plt.subplots()
            post.setup_plot(ax, pause)
        else:
            ax_supplied = True

        # plot all the elements in the mesh
        ax.triplot(self.mesh_nodes[:, 0], self.mesh_nodes[:, 1],
                   self.mesh_elements[:, 0:3], lw=0.5, color='black',
                   alpha=alpha)

        # if no axes object is supplied, finish the plot
        if not ax_supplied:
            post.finish_plot(ax, pause, title='Finite Element Mesh')

    def plot_centroids(self, pause=True):
        """Plots all the centroids and the shear centre, if they have been
        calculated, on top of the finite element mesh.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.

        The following example analyses a 200 PFC section and displays a plot of
        the centroids::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.PfcSection(d=200, b=75, t_f=12, t_w=6, r=12, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()

            section.plot_centroids()

        ..  figure:: ../images/pfc_centroids.png
            :align: center
            :scale: 75 %

            Plot of the elastic centroid and shear centre for the above example
            generated by
            :func:`~sectionproperties.analysis.cross_section.CrossSection.plot_centroids`

            # TODO: update example with plastic centroids
        """

        # TODO: add plot principal axes

        # create plot and setup the plot
        (fig, ax) = plt.subplots()
        post.setup_plot(ax, pause)

        # plot the finite element mesh
        self.plot_mesh(ax, pause, alpha=0.5)

        # if the elastic centroid has been calculated
        if self.section_props.cx is not None:
            ax.scatter(self.section_props.cx, self.section_props.cy,
                       edgecolors='r', facecolors='none', marker='o', s=100,
                       label='Elastic centroid')

        # if the shear centre has been calculated
        if self.section_props.x_se is not None:
            (x_s, y_s) = self.get_sc()
            ax.scatter(x_s, y_s, c='r', marker='+', s=100,
                       label='Shear centre')

        # if the global plastic centroid has been calculated
        if self.section_props.x_pc is not None:
            ax.scatter(self.section_props.x_pc, self.section_props.y_pc,
                       c='r', marker='x', s=100,
                       label='Global plastic centroid')

        # if the principal plastic centroid has been calculated
        if self.section_props.x1_pc is not None:
            ax.scatter(self.section_props.x1_pc, self.section_props.y2_pc,
                       edgecolors='r', marker='s', s=100,
                       label='Principal plastic centroid')

        # display the legend
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        # finish the plot
        post.finish_plot(ax, pause, title='Centroids')

    def display_mesh_info(self):
        """Prints mesh statistics (number of nodes, elements and regions) to
        the command window.

        The following example displays the mesh statistics for a 100D x 50W
        rectangle using a mesh size of 5::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.RectangularSection(d=100, b=50)
            mesh = geometry.create_mesh(mesh_sizes=[5])
            section = CrossSection(geometry, mesh)
            section.display_mesh_info()

            >>>Mesh Statistics:
            >>>--3282 nodes
            >>>--1591 elements
            >>>--1 region
        """

        print("Mesh Statistics:")
        print("--{0} nodes".format(self.num_nodes))
        print("--{0} elements".format(len(self.elements)))

        regions = max(self.mesh_attributes) + 1
        text = "--{0} region".format(regions)

        if regions == 1:
            text += "\n"
        else:
            text += "s\n"

        print(text)

    def display_results(self, fmt='8.6e'):
        """Prints the results that have been calculated to the terminal.

        :param string fmt: String containing number formatting information

        The following example displays the geometric section properties for a
        100D x 50W rectangle with three digits after the decimal point::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.RectangularSection(d=100, b=50)
            mesh = geometry.create_mesh(mesh_sizes=[5])

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()

            section.display_results(fmt='.3f')
        """

        post.print_results(self, fmt)

    def get_area(self):
        """
        :return: Cross-section area
        :rtype: float

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            area = section.get_area()
        """

        return self.section_props.area

    def get_q(self):
        """
        :return: First moments of area of the cross-section about the global
            axis *(qx, qy)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            (qx, qy) = section.get_q()
        """

        return (self.section_props.qx, self.section_props.qy)

    def get_ig(self):
        """
        :return: Second moments of area of the cross-section about the global
            axis *(ixx_g, iyy_g, ixy_g)*
        :rtype: tuple(float, float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            (ixx_g, iyy_g, ixy_g) = section.get_ig()
        """

        return (self.section_props.ixx_g, self.section_props.iyy_g,
                self.section_props.ixy_g)

    def get_c(self):
        """
        :return: Elastic centroid of the cross-section *(cx, cy)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            (cx, cy) = section.get_elastic_centroid()
        """

        return (self.section_props.cx, self.section_props.cy)

    def get_ic(self):
        """
        :return: Second moments of area of the cross-section about the
            centroidal axis *(ixx_c, iyy_c, ixy_c)*
        :rtype: tuple(float, float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            (ixx_c, iyy_c, ixy_c) = section.get_ic()
        """

        return (self.section_props.ixx_c, self.section_props.iyy_c,
                self.section_props.ixy_c)

    def get_z(self):
        """
        :return: Section moduli about the centroidal axis *(zxx_plus,
            zxx_minus, zyy_plus, zyy_minus)*
        :rtype: tuple(float, float, float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            (zxx_plus, zxx_minus, zyy_plus, zyy_minus) = section.get_z()
        """

        return (self.section_props.zxx_plus, self.section_props.zxx_minus,
                self.section_props.zyy_plus, self.section_props.zyy_minus)

    def get_rc(self):
        """
        :return: Radii of gyration about the centroidal x-axis *(rx, ry)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            (rx, ry) = section.get_rc()
        """

        return (self.section_props.rx_c, self.section_props.ry_c)

    def get_ip(self):
        """
        :return: Second moments of area of the cross-section about the
            principal axis *(i11_c, i22_c)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            (i11_c, i22_c) = section.get_ip()
        """

        return (self.section_props.i11_c, self.section_props.i22_c)

    def get_phi(self):
        """
        :return: Principal bending axis angle
        :rtype: float

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            phi = section.get_phi()
        """

        return self.section_props.phi

    def get_zp(self):
        """
        :return: Section moduli about the principal axis *(z11_plus,
            z11_minus, z22_plus, z22_minus)*
        :rtype: tuple(float, float, float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            (z11_plus, z11_minus, z22_plus, z22_minus) = section.get_zp()
        """

        return (self.section_props.z11_plus, self.section_props.z11_minus,
                self.section_props.z22_plus, self.section_props.z22_minus)

    def get_rp(self):
        """
        :return: Radii of gyration about the principal axis *(r11, r22)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            (r11, r22) = section.get_rp()
        """

        return (self.section_props.r11_c, self.section_props.r22_c)

    def get_j(self):
        """
        :return: Torsion constant
        :rtype: float

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            j = section.get_j()
        """

        return self.section_props.j

    def get_sc(self):
        """
        :return: Shear centre (elasticity approach) *(x_se, y_se)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            (x_se, y_se) = section.get_sc()
        """

        if self.section_props.x_se is None:
            return (None, None)
        else:
            # add centroid location to move section back to original location
            x_se = self.section_props.x_se + self.section_props.cx
            y_se = self.section_props.y_se + self.section_props.cy

        return (x_se, y_se)

    def get_sc_p(self):
        """
        :return: Principal axis shear centre (elasticity approach)
            *(x1_se, y2_se)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            (x1_se, y2_se) = section.get_sc_p()
        """

        if self.section_props.x1_se is None:
            return (None, None)
        else:
            # add centroid location to move section back to original location
            x1_se = self.section_props.x1_se + self.section_props.cx
            y2_se = self.section_props.y2_se + self.section_props.cy

        return (x1_se, y2_se)

    def get_sc_t(self):
        """
        :return: Shear centre (Trefftz's approach) *(x_st, y_st)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            (x_st, y_st) = section.get_sc_t()
        """

        if self.section_props.x_st is None:
            return (None, None)
        else:
            # add centroid location to move section back to original location
            x_st = self.section_props.x_st + self.section_props.cx
            y_st = self.section_props.y_st + self.section_props.cy

        return (x_st, y_st)

    def get_gamma(self):
        """
        :return: Warping constant
        :rtype: float

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            gamma = section.get_gamma()
        """

        return self.section_props.gamma

    def get_As(self):
        """
        :return: Shear areas *(A_sx, A_sy)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            (A_sx, A_sy) = section.get_As()
        """

        return (self.section_props.A_sx, self.section_props.A_sy)

    def get_As_p(self):
        """
        :return: Shear areas about the principal bending axis *(A_s11, A_s22)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            (A_s11, A_s22) = section.get_As_p()
        """

        return (self.section_props.A_s11, self.section_props.A_s22)


class SectionProperties:
    """Class for storing section properties.

    Stores calculated section properties. Also provides methods to calculate
    section properties entirely derived from other section properties.

    :cvar float area: Cross-sectional area
    :cvar float qx: First moment of area about the x-axis
    :cvar float qy: First moment of area about the y-axis
    :cvar float ixx_g: Second moment of area about the global x-axis
    :cvar float iyy_g: Second moment of area about the global y-axis
    :cvar float ixy_g: Second moment of area about the global xy-axis
    :cvar float cx: X coordinate of the elastic centroid
    :cvar float cy: Y coordinate of the elastic centroid
    :cvar float ixx_c: Second moment of area about the centroidal x-axis
    :cvar float iyy_c: Second moment of area about the centroidal y-axis
    :cvar float ixy_c: Second moment of area about the centroidal xy-axis
    :cvar float zxx_plus: Section modulus about the centroidal x-axis for
        stresses at the positive extreme value of y
    :cvar float zxx_minus: Section modulus about the centroidal x-axis for
        stresses at the negative extreme value of y
    :cvar float zyy_plus: Section modulus about the centroidal y-axis for
        stresses at the positive extreme value of x
    :cvar float zyy_minus: Section modulus about the centroidal y-axis for
        stresses at the negative extreme value of x
    :cvar float rx_c: Radius of gyration about the centroidal x-axis.
    :cvar float ry_c: Radius of gyration about the centroidal y-axis.
    :cvar float i11_c: Second moment of area about the centroidal 11-axis
    :cvar float i22_c: Second moment of area about the centroidal 22-axis
    :cvar float phi: Principal axis angle
    :cvar float z11_plus: Section modulus about the principal 11-axis for
        stresses at the positive extreme value of the 22-axis
    :cvar float z11_minus: Section modulus about the principal 11-axis for
        stresses at the negative extreme value of the 22-axis
    :cvar float z22_plus: Section modulus about the principal 22-axis for
        stresses at the positive extreme value of the 11-axis
    :cvar float z22_minus: Section modulus about the principal 22-axis for
        stresses at the negative extreme value of the 11-axis
    :cvar float r11_c: Radius of gyration about the principal 11-axis.
    :cvar float r22_c: Radius of gyration about the principal 22-axis.
    :cvar float xmax: Maximum x-axis coordinate
    :cvar float xmin: Minimum x-axis coordinate
    :cvar float ymax: Maximum y-axis coordinate
    :cvar float ymin: Minimum y-axis coordinate
    :cvar float x1max: Maximum 11-axis coordinate
    :cvar float x1min: Minimum 11-axis coordinate
    :cvar float y2max: Maximum 22-axis coordinate
    :cvar float y2min: Minimum 22-axis coordinate
    :cvar float j: Torsion constant
    :cvar omega: Warping function
    :vartype omega: :class:`numpy.ndarray`
    :cvar psi_shear: Psi shear function
    :vartype psi_shear: :class:`numpy.ndarray`
    :cvar phi_shear: Phi shear function
    :vartype phi_shear: :class:`numpy.ndarray`
    :cvar float Delta_s: Shear factor
    :cvar float x_se: X coordinate of the shear centre (elasticity approach)
    :cvar float y_se: Y coordinate of the shear centre (elasticity approach)
    :cvar float x1_se: 11 coordinate of the shear centre (elasticity approach)
    :cvar float y2_se: 22 coordinate of the shear centre (elasticity approach)
    :cvar float x_st: X coordinate of the shear centre (Trefftz's approach)
    :cvar float y_st: Y coordinate of the shear centre (Trefftz's approach)
    :cvar float gamma: Warping constant
    :cvar float A_sx: Shear area about the x-axis
    :cvar float A_sy: Shear area about the y-axis
    :cvar float A_sxy: Shear area about the xy-axis
    :cvar float A_s11: Shear area about the 11 bending axis
    :cvar float A_s22: Shear area about the 22 bending axis
    :cvar float x_pc: X coordinate of the global plastic centroid
    :cvar float y_pc: Y coordinate of the global plastic centroid
    :cvar float x1_pc: 11 coordinate of the principal plastic centroid
    :cvar float y2_pc: 22 coordinate of the principal plastic centroid
    :cvar float Sxx: Plastic section modulus about the centroidal x-axis
    :cvar float Syy: Plastic section modulus about the centroidal y-axis
    :cvar float SF_xx_plus: Shape factor for bending about the x-axis with
        respect to the top fibre
    :cvar float SF_xx_minus: Shape factor for bending about the x-axis with
        respect to the bottom fibre
    :cvar float SF_yy_plus: Shape factor for bending about the y-axis with
        respect to the top fibre
    :cvar float SF_yy_minus: Shape factor for bending about the y-axis with
        respect to the bottom fibre
    :cvar float S11: Plastic section modulus about the 11-axis
    :cvar float S22: Plastic section modulus about the 22-axis
    :cvar float SF_11_plus: Shape factor for bending about the 11-axis with
        respect to the top fibre
    :cvar float SF_11_minus: Shape factor for bending about the 11-axis with
        respect to the bottom fibre
    :cvar float SF_22_plus: Shape factor for bending about the 22-axis with
        respect to the top fibre
    :cvar float SF_22_minus: Shape factor for bending about the 22-axis with
        respect to the bottom fibre
    """

    def __init__(self):
        """Inits the SectionProperties class."""

        self.area = None
        self.qx = None
        self.qy = None
        self.ixx_g = None
        self.iyy_g = None
        self.ixy_g = None
        self.cx = None
        self.cy = None
        self.ixx_c = None
        self.iyy_c = None
        self.ixy_c = None
        self.zxx_plus = None
        self.zxx_minus = None
        self.zyy_plus = None
        self.zyy_minus = None
        self.rx_c = None
        self.ry_c = None
        self.i11_c = None
        self.i22_c = None
        self.phi = None
        self.z11_plus = None
        self.z11_minus = None
        self.z22_plus = None
        self.z22_minus = None
        self.r11_c = None
        self.r22_c = None
        self.j = None
        self.omega = None
        self.psi_shear = None
        self.phi_shear = None
        self.Delta_s = None
        self.x_se = None
        self.y_se = None
        self.x1_se = None
        self.y2_se = None
        self.x_st = None
        self.y_st = None
        self.gamma = None
        self.A_sx = None
        self.A_sy = None
        self.A_sxy = None
        self.A_s11 = None
        self.A_s22 = None
        self.x_pc = None
        self.y_pc = None
        self.x1_pc = None
        self.y2_pc = None
        self.Sxx = None
        self.Syy = None
        self.SF_xx_plus = None
        self.SF_xx_minus = None
        self.SF_yy_plus = None
        self.SF_yy_minus = None
        self.S11 = None
        self.S22 = None
        self.SF_11_plus = None
        self.SF_11_minus = None
        self.SF_22_plus = None
        self.SF_22_minus = None

    def calculate_elastic_centroid(self):
        """Calculates the elastic centroid based on the cross-section area and
        first moments of area.

        :return: Tuple containing the coordinates of the elastic centroid.
        :rtype: tuple(float, float)
        """

        self.cx = self.qy / self.area
        self.cy = self.qx / self.area

        return (self.cx, self.cy)

    def calculate_centroidal_properties(self, mesh):
        """Calculates the geometric section properties about the centroidal and
        principal axes based on the results about the global axis.
        """

        # calculate second moments of area about the centroidal xy axis
        self.ixx_c = self.ixx_g - self.qx ** 2 / self.area
        self.iyy_c = self.iyy_g - self.qy ** 2 / self.area
        self.ixy_c = self.ixy_g - self.qx * self.qy / self.area

        # calculate section moduli about the centroidal xy axis
        nodes = np.array(mesh.points)
        self.xmax = nodes[:, 0].max()
        self.xmin = nodes[:, 0].min()
        self.ymax = nodes[:, 1].max()
        self.ymin = nodes[:, 1].min()
        self.zxx_plus = self.ixx_c / abs(self.ymax - self.cy)
        self.zxx_minus = self.ixx_c / abs(self.ymin - self.cy)
        self.zyy_plus = self.iyy_c / abs(self.xmax - self.cx)
        self.zyy_minus = self.iyy_c / abs(self.xmin - self.cx)

        # calculate radii of gyration about centroidal xy axis
        self.rx_c = (self.ixx_c / self.area) ** 0.5
        self.ry_c = (self.iyy_c / self.area) ** 0.5

        # calculate prinicpal 2nd moments of area about the centroidal xy axis
        Delta = (((self.ixx_c - self.iyy_c) / 2) ** 2 + self.ixy_c ** 2) ** 0.5
        self.i11_c = (self.ixx_c + self.iyy_c) / 2 + Delta
        self.i22_c = (self.ixx_c + self.iyy_c) / 2 - Delta

        # calculate initial principal axis angle
        if abs(self.ixx_c - self.i11_c) < 1e-12 * self.i11_c:
            self.phi = 0
        else:
            self.phi = np.arctan2(
                self.ixx_c - self.i11_c, self.ixy_c) * 180 / np.pi

        # calculate section moduli about the principal axis
        for (i, pt) in enumerate(nodes):
            x = pt[0] - self.cx
            y = pt[1] - self.cy
            # determine the coordinate of the point wrt the principal axis
            (x1, y2) = fea.principal_coordinate(self.phi, x, y)

            # initialise min, max variables
            if i == 0:
                self.x1max = x1
                self.x1min = x1
                self.y2max = y2
                self.y2min = y2

            # update the mins and maxs where necessary
            self.x1max = max(self.x1max, x1)
            self.x1min = min(self.x1min, x1)
            self.y2max = max(self.y2max, y2)
            self.y2min = min(self.y2min, y2)

        # evaluate principal section moduli
        self.z11_plus = self.i11_c / abs(self.y2max)
        self.z11_minus = self.i11_c / abs(self.y2min)
        self.z22_plus = self.i22_c / abs(self.x1max)
        self.z22_minus = self.i22_c / abs(self.x1min)

        # calculate radii of gyration about centroidal principal axis
        self.r11_c = (self.i11_c / self.area) ** 0.5
        self.r22_c = (self.i22_c / self.area) ** 0.5


class StressResult:
    """Class for storing a stress result.

    Provides a way to store the results from a cross-section stress analysis.
    Also provides methods to calculate combined stresses and post-process the
    stress results.

    :param cross_section: CrossSection object that the stress analysis was
        performed on
    :type cross_section:
        :class:`~sectionproperties.analysis.cross_section.CrossSection`
    :param int num_nodes: Number of nodes in the finite element mesh

    :cvar cross_section: CrossSection object that the stress analysis was
        performed on
    :vartype cross_section:
        :class:`~sectionproperties.analysis.cross_section.CrossSection`
    :cvar sig_zz_n: Normal stress (:math:`\sigma_{zz,N}`) resulting from
        an axial force
    :vartype sig_zz_n: :class:`numpy.ndarray`
    :cvar sig_zz_mxx: Normal stress (:math:`\sigma_{zz,Mxx}`) resulting from
        a bending moment about the xx-axis
    :vartype sig_zz_mxx: :class:`numpy.ndarray`
    :cvar sig_zz_myy: Normal stress (:math:`\sigma_{zz,Myy}`) resulting from
        a bending moment about the yy-axis
    :vartype sig_zz_myy: :class:`numpy.ndarray`
    :cvar sig_zz_m11: Normal stress (:math:`\sigma_{zz,M11}`) resulting from
        a bending moment about the 11-axis
    :vartype sig_zz_m11: :class:`numpy.ndarray`
    :cvar sig_zz_m22: Normal stress (:math:`\sigma_{zz,M22}`) resulting from
        a bending moment about the 22-axis
    :vartype sig_zz_m22: :class:`numpy.ndarray`
    :cvar sig_zx_mzz: Shear stress (:math:`\sigma_{zx,Mzz}`) resulting from
        a torsion moment about the zz-axis
    :vartype sig_zx_mzz: :class:`numpy.ndarray`
    :cvar sig_zy_mzz: Shear stress (:math:`\sigma_{zy,Mzz}`) resulting from
        a torsion moment about the zz-axis
    :vartype sig_zy_mzz: :class:`numpy.ndarray`
    :cvar sig_zx_vx: Shear stress (:math:`\sigma_{zx,Vx}`) resulting from
        a shear force in the x-direction
    :vartype sig_zx_vx: :class:`numpy.ndarray`
    :cvar sig_zy_vx: Shear stress (:math:`\sigma_{zy,Vx}`) resulting from
        a shear force in the x-direction
    :vartype sig_zy_vx: :class:`numpy.ndarray`
    :cvar sig_zx_vy: Shear stress (:math:`\sigma_{zx,Vy}`) resulting from
        a shear force in the y-direction
    :vartype sig_zx_vy: :class:`numpy.ndarray`
    :cvar sig_zy_vy: Shear stress (:math:`\sigma_{zy,Vy}`) resulting from
        a shear force in the y-direction
    :vartype sig_zy_vy: :class:`numpy.ndarray`
    :cvar sig_zz_m: Normal stress (:math:`\sigma_{zz,\Sigma M}`) resulting from
        all applied bending moments
    :vartype sig_zz_m: :class:`numpy.ndarray`
    :cvar sig_zxy_mzz: Resultant shear stress (:math:`\sigma_{zxy,Mzz}`)
        resulting from a torsion moment in the zz-direction
    :vartype sig_zxy_mzz: :class:`numpy.ndarray`
    :cvar sig_zxy_vx: Resultant shear stress (:math:`\sigma_{zxy,Vx}`)
        resulting from a a shear force in the x-direction
    :vartype sig_zxy_vx: :class:`numpy.ndarray`
    :cvar sig_zxy_vy: Resultant shear stress (:math:`\sigma_{zxy,Vy}`)
        resulting from a a shear force in the y-direction
    :vartype sig_zxy_vy: :class:`numpy.ndarray`
    :cvar sig_zx_v: Shear stress (:math:`\sigma_{zx,\Sigma V}`) resulting from
        all applied shear forces
    :vartype sig_zx_v: :class:`numpy.ndarray`
    :cvar sig_zy_v: Shear stress (:math:`\sigma_{zy,\Sigma V}`) resulting from
        all applied shear forces
    :vartype sig_zy_v: :class:`numpy.ndarray`
    :cvar sig_zxy_v: Resultant shear stress (:math:`\sigma_{zxy,\Sigma V}`)
        resulting from all applied shear forces
    :vartype sig_zxy_v: :class:`numpy.ndarray`
    :cvar sig_zz: Combined normal force (:math:`\sigma_{zz}`) resulting from
        all applied actions
    :vartype sig_zz: :class:`numpy.ndarray`
    :cvar sig_zx: Combined shear stress (:math:`\sigma_{zx}`) resulting from
        all applied actions
    :vartype sig_zx: :class:`numpy.ndarray`
    :cvar sig_zy: Combined shear stress (:math:`\sigma_{zy}`) resulting from
        all applied actions
    :vartype sig_zy: :class:`numpy.ndarray`
    :cvar sig_zxy: Combined resultant shear stress (:math:`\sigma_{zxy}`)
        resulting from all applied actions
    :vartype sig_zxy: :class:`numpy.ndarray`
    :cvar sig_vm: von Mises stress (:math:`\sigma_{VM}`) resulting from
        all applied actions
    :vartype sig_vm: :class:`numpy.ndarray`
    """

    def __init__(self, cross_section, num_nodes):
        """Inits the StressResult class."""

        self.cross_section = cross_section

        # allocate stresses arising directly from actions
        self.sig_zz_n = np.zeros(num_nodes)
        self.sig_zz_mxx = np.zeros(num_nodes)
        self.sig_zz_myy = np.zeros(num_nodes)
        self.sig_zz_m11 = np.zeros(num_nodes)
        self.sig_zz_m22 = np.zeros(num_nodes)
        self.sig_zx_mzz = np.zeros(num_nodes)
        self.sig_zy_mzz = np.zeros(num_nodes)
        self.sig_zx_vx = np.zeros(num_nodes)
        self.sig_zy_vx = np.zeros(num_nodes)
        self.sig_zx_vy = np.zeros(num_nodes)
        self.sig_zy_vy = np.zeros(num_nodes)

        # allocate combined stresses
        self.sig_zz_m = np.zeros(num_nodes)
        self.sig_zxy_mzz = np.zeros(num_nodes)
        self.sig_zxy_vx = np.zeros(num_nodes)
        self.sig_zxy_vy = np.zeros(num_nodes)
        self.sig_zx_v = np.zeros(num_nodes)
        self.sig_zy_v = np.zeros(num_nodes)
        self.sig_zxy_v = np.zeros(num_nodes)
        self.sig_zz = np.zeros(num_nodes)
        self.sig_zx = np.zeros(num_nodes)
        self.sig_zy = np.zeros(num_nodes)
        self.sig_zxy = np.zeros(num_nodes)
        self.sig_vm = np.zeros(num_nodes)

    def calculate_combined_stresses(self):
        """Calculates the combined cross-section stresses."""

        self.sig_zz_m = (self.sig_zz_mxx + self.sig_zz_myy + self.sig_zz_m11 +
                         self.sig_zz_m22)
        self.sig_zxy_mzz = (self.sig_zx_mzz ** 2 + self.sig_zy_mzz ** 2) ** 0.5
        self.sig_zxy_vx = (self.sig_zx_vx ** 2 + self.sig_zy_vx ** 2) ** 0.5
        self.sig_zxy_vy = (self.sig_zx_vy ** 2 + self.sig_zy_vy ** 2) ** 0.5
        self.sig_zx_v = self.sig_zx_vx + self.sig_zx_vy
        self.sig_zy_v = self.sig_zy_vx + self.sig_zy_vy
        self.sig_zxy_v = (self.sig_zx_v ** 2 + self.sig_zy_v ** 2) ** 0.5
        self.sig_zz = self.sig_zz_n + self.sig_zz_m
        self.sig_zx = self.sig_zx_mzz + self.sig_zx_v
        self.sig_zy = self.sig_zy_mzz + self.sig_zy_v
        self.sig_zxy = (self.sig_zx ** 2 + self.sig_zy ** 2) ** 0.5
        self.sig_vm = (self.sig_zz ** 2 + 3 * self.sig_zxy ** 2) ** 0.5

    def plot_stress_contour(self, sig, title, pause):
        """Plots filled stress contours over the finite element mesh.

        :param sig: Nodal stress values
        :type sig: :class:`numpy.ndarray`
        :param string title: Plot title
        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        # create plot and setup the plot
        (fig, ax) = plt.subplots()
        post.setup_plot(ax, pause)

        # plot the finite element mesh
        self.cross_section.plot_mesh(ax, pause, alpha=0.5)

        # set up the colormap and contours
        cmap = cm.get_cmap(name='jet')

        # if contour values are not all constant
        if np.amax(sig) - np.amin(sig) > 1e-6:
            v = np.linspace(np.amin(sig), np.amax(sig), 10, endpoint=True)
        else:
            # ten contours
            v = 10

        # plot the filled contour
        trictr = ax.tricontourf(
            self.cross_section.mesh_nodes[:, 0],
            self.cross_section.mesh_nodes[:, 1],
            self.cross_section.mesh_elements[:, 0:3], sig, v, cmap=cmap)
        fig.colorbar(trictr, label='Stress', format='%.4e')

        # finish the plot
        post.finish_plot(ax, pause, title=title)

    def plot_stress_vector(self, sigx, sigy, title, pause):
        """a"""

        pass
        # TODO: implement

    # TODO: implement examples for stress plot methods

    def plot_stress_n_zz(self, pause=True):
        """Produces a contour plot of the normal stress :math:`\sigma_{zz,N}`
        resulting from the applied axial load :math:`N`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zz,N}$'
        self.plot_stress_contour(self.sig_zz_n, title, pause)

    def plot_stress_mxx_zz(self, pause=True):
        """Produces a contour plot of the normal stress :math:`\sigma_{zz,Mxx}`
        resulting from the applied bending moment :math:`M_{xx}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zz,Mxx}$'
        self.plot_stress_contour(self.sig_zz_mxx, title, pause)

    def plot_stress_myy_zz(self, pause=True):
        """Produces a contour plot of the normal stress :math:`\sigma_{zz,Myy}`
        resulting from the applied bending moment :math:`M_{yy}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zz,Myy}$'
        self.plot_stress_contour(self.sig_zz_myy, title, pause)

    def plot_stress_m11_zz(self, pause=True):
        """Produces a contour plot of the normal stress :math:`\sigma_{zz,M11}`
        resulting from the applied bending moment :math:`M_{11}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zz,M11}$'
        self.plot_stress_contour(self.sig_zz_m11, title, pause)

    def plot_stress_m22_zz(self, pause=True):
        """Produces a contour plot of the normal stress :math:`\sigma_{zz,M22}`
        resulting from the applied bending moment :math:`M_{22}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zz,M22}$'
        self.plot_stress_contour(self.sig_zz_m22, title, pause)

    def plot_stress_m_zz(self, pause=True):
        """Produces a contour plot of the normal stress
        :math:`\sigma_{zz,\Sigma M}` resulting from all applied bending moments
        :math:`M_{xx} + M_{yy} + M_{11} + M_{22}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zz,\Sigma M}$'
        self.plot_stress_contour(self.sig_zz_m, title, pause)

    def plot_stress_mzz_zx(self, pause=True):
        """Produces a contour plot of the *x*-component of the shear stress
        :math:`\sigma_{zx,Mzz}` resulting from the applied torsion moment
        :math:`M_{zz}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zx,Mzz}$'
        self.plot_stress_contour(self.sig_zx_mzz, title, pause)

    def plot_stress_mzz_zy(self, pause=True):
        """Produces a contour plot of the *y*-component of the shear stress
        :math:`\sigma_{zy,Mzz}` resulting from the applied torsion moment
        :math:`M_{zz}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zy,Mzz}$'
        self.plot_stress_contour(self.sig_zy_mzz, title, pause)

    def plot_stress_mzz_zxy(self, pause=True):
        """Produces a contour plot of the resultant shear stress
        :math:`\sigma_{zxy,Mzz}` resulting from the applied torsion moment
        :math:`M_{zz}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zxy,Mzz}$'
        self.plot_stress_contour(self.sig_zxy_mzz, title, pause)

    def plot_vector_mzz_zxy(self, pause=True):
        """Produces a vector plot of the resultant shear stress
        :math:`\sigma_{zxy,Mzz}` resulting from the applied torsion moment
        :math:`M_{zz}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Vector Plot - $\sigma_{zxy,Mzz}$'
        self.plot_stress_vector(self.sig_zx_mzz, self.sig_zy_mzz, title, pause)

    def plot_stress_vx_zx(self, pause=True):
        """Produces a contour plot of the *x*-component of the shear stress
        :math:`\sigma_{zx,Vx}` resulting from the applied shear force
        :math:`V_{x}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zx,Vx}$'
        self.plot_stress_contour(self.sig_zx_vx, title, pause)

    def plot_stress_vx_zy(self, pause=True):
        """Produces a contour plot of the *y*-component of the shear stress
        :math:`\sigma_{zy,Vx}` resulting from the applied shear force
        :math:`V_{x}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zy,Vx}$'
        self.plot_stress_contour(self.sig_zy_vx, title, pause)

    def plot_stress_vx_zxy(self, pause=True):
        """Produces a contour plot of the resultant shear stress
        :math:`\sigma_{zxy,Vx}` resulting from the applied shear force
        :math:`V_{x}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zxy,Vx}$'
        self.plot_stress_contour(self.sig_zxy_vx, title, pause)

    def plot_vector_vx_zxy(self, pause=True):
        """Produces a vector plot of the resultant shear stress
        :math:`\sigma_{zxy,Vx}` resulting from the applied shear force
        :math:`V_{x}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Vector Plot - $\sigma_{zxy,Vx}$'
        self.plot_stress_vector(self.sig_zx_vx, self.sig_zx_vy, title, pause)

    def plot_stress_vy_zx(self, pause=True):
        """Produces a contour plot of the *x*-component of the shear stress
        :math:`\sigma_{zx,Vy}` resulting from the applied shear force
        :math:`V_{y}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zx,Vy}$'
        self.plot_stress_contour(self.sig_zx_vy, title, pause)

    def plot_stress_vy_zy(self, pause=True):
        """Produces a contour plot of the *y*-component of the shear stress
        :math:`\sigma_{zy,Vy}` resulting from the applied shear force
        :math:`V_{y}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zy,Vy}$'
        self.plot_stress_contour(self.sig_zy_vy, title, pause)

    def plot_stress_vy_zxy(self, pause=True):
        """Produces a contour plot of the resultant shear stress
        :math:`\sigma_{zxy,Vy}` resulting from the applied shear force
        :math:`V_{y}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zxy,Vy}$'
        self.plot_stress_contour(self.sig_zxy_vy, title, pause)

    def plot_vector_vy_zxy(self, pause=True):
        """Produces a vector plot of the resultant shear stress
        :math:`\sigma_{zxy,Vy}` resulting from the applied shear force
        :math:`V_{y}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Vector Plot - $\sigma_{zxy,Vy}$'
        self.plot_stress_vector(self.sig_zx_vy, self.sig_zy_vy, title, pause)

    def plot_stress_v_zx(self, pause=True):
        """Produces a contour plot of the *x*-component of the shear stress
        :math:`\sigma_{zx,\Sigma V}` resulting from the sum of the applied
        shear forces :math:`V_{x} + V_{y}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zx,\Sigma V}$'
        self.plot_stress_contour(self.sig_zx_v, title, pause)

    def plot_stress_v_zy(self, pause=True):
        """Produces a contour plot of the *y*-component of the shear stress
        :math:`\sigma_{zy,\Sigma V}` resulting from the sum of the applied
        shear forces :math:`V_{x} + V_{y}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zy,\Sigma V}$'
        self.plot_stress_contour(self.sig_zy_v, title, pause)

    def plot_stress_v_zxy(self, pause=True):
        """Produces a contour plot of the resultant shear stress
        :math:`\sigma_{zxy,\Sigma V}` resulting from the sum of the applied
        shear forces :math:`V_{x} + V_{y}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zxy,\Sigma V}$'
        self.plot_stress_contour(self.sig_zxy_v, title, pause)

    def plot_vector_v_zxy(self, pause=True):
        """Produces a vector plot of the resultant shear stress
        :math:`\sigma_{zxy,\Sigma V}` resulting from the sum of the  applied
        shear forces :math:`V_{x} + V_{y}`.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Vector Plot - $\sigma_{zxy,\Sigma V}$'
        self.plot_stress_vector(self.sig_zx_v, self.sig_zy_v, title, pause)

    def plot_stress_zz(self, pause=True):
        """Produces a contour plot of the combined normal stress
        :math:`\sigma_{zz}` resulting from all applied actions.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zz}$'
        self.plot_stress_contour(self.sig_zz, title, pause)

    def plot_stress_zx(self, pause=True):
        """Produces a contour plot of the *x*-component of the shear stress
        :math:`\sigma_{zx}` resulting from all applied actions.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zx}$'
        self.plot_stress_contour(self.sig_zx, title, pause)

    def plot_stress_zy(self, pause=True):
        """Produces a contour plot of the *y*-component of the shear stress
        :math:`\sigma_{zy}` resulting from all applied actions.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zy}$'
        self.plot_stress_contour(self.sig_zy, title, pause)

    def plot_stress_zxy(self, pause=True):
        """Produces a contour plot of the resultant shear stress
        :math:`\sigma_{zxy}` resulting from all applied actions.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{zxy}$'
        self.plot_stress_contour(self.sig_zxy, title, pause)

    def plot_vector_zxy(self, pause=True):
        """Produces a vector plot of the resultant shear stress
        :math:`\sigma_{zxy}` resulting from all applied actions.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Vector Plot - $\sigma_{zxy}$'
        self.plot_stress_vector(self.sig_zx, self.sig_zy, title, pause)

    def plot_stress_vm(self, pause=True):
        """Produces a contour plot of the von Mises :math:`\sigma_{vM}`
        resulting from all applied actions.

        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        title = 'Stress Contour Plot - $\sigma_{vM}$'
        self.plot_stress_contour(self.sig_vm, title, pause)

# def computeSectionProperties(self, points, facets, holes, controlPoints,
#                              materials):
#     """
#     This function computes the the cross section properties for the mesh.
#     """
#
#     if (self.settings.plasticAnalysis):
#         # calculate global plastic properties
#         if (self.settings.outputLog):
#             processText = "-- Calculating global plastic properties..."
#         otherUtilities.functionTimer(
#             processText, self.computeGlobalPlasticProperties, points,
#             facets, holes, controlPoints, materials)
#
#         # calculate principal plastic properties
#         if (self.settings.outputLog):
#             processText = "-- Calculating principal plastic properties..."
#         otherUtilities.functionTimer(
#             processText, self.computePrincipalPlasticProperties, points,
#             facets, holes, controlPoints, materials)
#
#
# def computeGlobalPlasticProperties(self, points, facets, holes,
#                                    controlPoints, materials):
#     """
#     This method computes the plastic centroid and plastic section moduli
#     for global axis bending.
#     """
#
#     # unit vectors in the x & y directions
#     ux = [1, 0]
#     uy = [0, 1]
#
#     # compute plastic centroids and plastic section moduli:
#     # compute x-location of plastic centroid
#     (self.x_pc, topA, botA, topCenX, botCenX) = otherUtilities.pcAlgorithm(
#         self.settings.tol, 1000, uy, self.xmin, self.xmax,
#         points, facets, holes, controlPoints, self.nodes, self.elements,
#         materials, 1)
#
#     # compute y-location of the plastic centroid
#     (self.y_pc, topA, botA, topCenY, botCenY) = otherUtilities.pcAlgorithm(
#         self.settings.tol, 1000, ux, self.ymin, self.ymax,
#         points, facets, holes, controlPoints, self.nodes, self.elements,
#         materials, 2)
#
#     self.Sxx = self.area / 2 * abs(topCenY[1] - botCenY[1])
#     self.Syy = self.area / 2 * abs(topCenX[0] - botCenX[0])
#
#     # compute shape factors
#     self.SF_xx_plus = self.Sxx / self.zxx_plus
#     self.SF_xx_minus = self.Sxx / self.zxx_minus
#     self.SF_yy_plus = self.Syy / self.zyy_plus
#     self.SF_yy_minus = self.Syy / self.zyy_minus
#
# def computePrincipalPlasticProperties(self, points, facets, holes,
#                                       controlPoints, materials):
#     """
#     This method computes the plastic centroid and plastic section moduli
#     for principal axis bending.
#     """
#
#     # unit vectors in the 1 & 2 directions
#     u1 = np.array(
#         [np.cos(self.phi * np.pi / 180), np.sin(self.phi * np.pi / 180)])
#     u2 = np.array(
#         [-np.sin(self.phi * np.pi / 180), np.cos(self.phi * np.pi / 180)])
#
#     # compute plastic centroids and plastic section moduli
#     (x1_pc, topA, botA, topCen1, botCen1) = otherUtilities.pcAlgorithm(
#         self.settings.tol, 1000, -u1, self.y2min, self.y2max,
#         points, facets, holes, controlPoints, self.nodes, self.elements,
#         materials, 1)
#     (y2_pc, topA, botA, topCen2, botCen2) = otherUtilities.pcAlgorithm(
#         self.settings.tol, 1000, -u2, self.x1min, self.x1max,
#         points, facets, holes, controlPoints, self.nodes, self.elements,
#         materials, 2)
#
#     # calculate the area centroids in the principal coordinate system
#     (tc1_1, tc1_2) = otherUtilities.principalCoordinate(
#         self.phi, topCen1[0], topCen1[1])
#     (bc1_1, bc1_2) = otherUtilities.principalCoordinate(
#         self.phi, botCen1[0], botCen1[1])
#     (tc2_1, tc2_2) = otherUtilities.principalCoordinate(
#         self.phi, topCen2[0], topCen2[1])
#     (bc2_1, bc2_2) = otherUtilities.principalCoordinate(
#         self.phi, botCen2[0], botCen2[1])
#
#     self.x1_pc = x1_pc * u2[0] + y2_pc * u1[0]
#     self.y2_pc = x1_pc * u2[1] + y2_pc * u1[1]
#     self.S11 = self.area / 2 * abs(tc1_2 - bc1_2)
#     self.S22 = self.area / 2 * abs(tc2_1 - bc2_1)
#
#     # compute shape factors
#     self.SF_11_plus = self.S11 / self.z11_plus
#     self.SF_11_minus = self.S11 / self.z11_minus
#     self.SF_22_plus = self.S22 / self.z22_plus
#     self.SF_22_minus = self.S22 / self.z22_minus

# def computeAreaSegments(self, u, px, py):
#     """
#     This method computes the area above and below a line defined by unit
#     vector u and point(px, py)
#     """
#
#     # allocate area variables
#     topA = 0
#     botA = 0
#     topQx = 0
#     topQy = 0
#     botQx = 0
#     botQy = 0
#     topCen = [0, 0]
#     botCen = [0, 0]
#
#     # loop through all elements in the mesh
#     for el in self.triElements:
#         # calculate area of element and its first moments of area
#         (elA, Qx, Qy) = el.areaProperties()
#
#         # if the element is not infinitessimally small (meshing artefacts)
#         if elA != 0:
#             # calculate the element centroid
#             elCen = [Qy / elA, Qx / elA]
#         else:
#             elCen = [0, 0]
#
#         # determine location of element and allocate element areas and
#         # first moments of area accordingly
#         if (otherUtilities.pointAboveLine(u, px, py, elCen[0], elCen[1])):
#             topA += elA
#             topQx += Qx
#             topQy += Qy
#         else:
#             botA += elA
#             botQx += Qx
#             botQy += Qy
#
#     # if the element is not infinitessimally small
#     if (topA != 0 and botA != 0):
#         # calculate the centroid of the top and bottom areas
#         topCen = np.array([topQy / topA, topQx / topA])
#         botCen = np.array([botQy / botA, botQx / botA])
#
#     return (topA, botA, topCen, botCen)
#
# def plotResults(self, plots):
#     """
#     This method generates all the plots in the 'plots' list.
#     """
#
#     for plot in plots:
#         # initialise plot variables
#         x = []  # x-component of vector
#         y = []  # y-component of vector
#         z = []  # contour values
#         globalAxis = False
#         principalAxis = False
#         nodes = False
#         plotTitle = ""
#         centroids = False
#         plotType = ""
#
#         if (plot.lower() == "axial"):
#             z = self.axialStress
#             plotTitle = "Axial Stress"
#             plotType = "contour"
#         elif (plot.lower() == "bending"):
#             z = self.bendingStress
#             plotTitle = "Bending Stress"
#             plotType = "contour"
#         elif (plot.lower() == "torsion"):
#             z = self.torsionStress
#             plotTitle = "Torsion Stress"
#             plotType = "contour"
#         elif (plot.lower() == "torsion-vector"):
#             x = self.torsionStress_zx
#             y = self.torsionStress_zy
#             plotTitle = "Torsion Stress Vectors"
#             plotType = "vector"
#         elif (plot.lower() == "shear"):
#             z = self.shearStress
#             plotTitle = "Transverse Shear Stress"
#             plotType = "contour"
#         elif (plot.lower() == "shear-zx"):
#             z = self.shearStress_zx
#             plotTitle = "Transverse Shear (zx) Stress"
#             plotType = "contour"
#         elif (plot.lower() == "shear-zy"):
#             z = self.shearStress_zy
#             plotTitle = "Transverse Shear (zy) Stress"
#             plotType = "contour"
#         elif (plot.lower() == "shear-vector"):
#             x = self.shearStress_zx
#             y = self.shearStress_zy
#             plotTitle = "Transverse Shear Stress Vectors"
#             plotType = "vector"
#         elif (plot.lower() == "combined-normal"):
#             z = self.sigma_zz
#             plotTitle = "Combined Normal Stress"
#             plotType = "contour"
#         elif (plot.lower() == "combined-shear"):
#             z = self.tau
#             plotTitle = "Combined Shear Stress"
#         elif (plot.lower() == "combined-shear-vector"):
#             x = self.tau_zx
#             y = self.tau_zy
#             plotTitle = "Combined Shear Stress Vectors"
#             plotType = "vector"
#         elif (plot.lower() == "von-mises"):
#             z = self.vonMises
#             plotTitle = "von Mises Stress"
#             plotType = "contour"
#         elif (plot.lower() == "centroids"):
#             z = None
#             globalAxis = True
#             principalAxis = True
#             centroids = True
#             plotTitle = "Centroids"
#             plotType = "contour"
#         elif (plot.lower() == "mesh"):
#             z = None
#             globalAxis = True
#             plotTitle = "Mesh"
#             plotType = "contour"
#
#         # if we are displaying a contour plot
#         if (plotType == "contour"):
#             self.contourPlot(
#                 globalAxis=globalAxis, principalAxis=principalAxis, z=z,
#                 nodes=nodes, plotTitle=plotTitle, centroids=centroids)
#         elif (plotType == "vector"):
#             self.quiverPlot(u=x, v=y, plotTitle=plotTitle)
#
# def contourPlot(self, globalAxis=False, principalAxis=False, z=None,
#                 nodes=False, plotTitle="", centroids=False):
#     """
#     This method generates a plot of the mesh with an optional contour
#     plot of results(z). Additional options include displaying the nodes, a
#     plot title, and the principal axis and centroids.
#     """
#
#     fig, ax = plt.subplots()
#     plt.ion()  # interactive mode enabled
#     plt.show()  # show the plot
#     ax.set_aspect("equal")  # set the scale on the x and y axes equal
#
#     # plot the title and axis labels
#     ax.set_title(plotTitle)
#     ax.set_xlabel("x")
#     ax.set_ylabel("y")
#
#     # plot the mesh
#     ax.triplot(self.nodes[:, 0], self.nodes[:, 1], self.elements[:, 0:3],
#                lw=0.5, color='black')
#
#     # plot the global axis as lines
#     if (globalAxis):
#         # determine min and max values of the nodes
#         (xmin, ymin) = np.amin(self.nodes, axis=0)
#         (xmax, ymax) = np.amax(self.nodes, axis=0)
#         xLim = xmax - xmin
#         yLim = ymax - ymin
#
#         # plot x axis
#         ax.plot([xmin - 0.1 * xLim, xmax + 0.1 * xLim],
#                 [-self.cy, -self.cy], label="Global x-axis")
#         # plot y axis
#         ax.plot([-self.cx, -self.cx],
#                 [ymin - 0.1 * yLim, ymax + 0.1 * yLim],
#                 label="Global y-axis")
#
#     # plot the principal axis as lines
#     if principalAxis:
#         start_11 = otherUtilities.globalCoordinate(
#             self.phi, self.x1min, 0)
#         end_11 = otherUtilities.globalCoordinate(
#             self.phi, self.x1max, 0)
#         start_22 = otherUtilities.globalCoordinate(
#             self.phi, 0, self.y2min)
#         end_22 = otherUtilities.globalCoordinate(
#             self.phi, 0, self.y2max)
#
#         lim11_x = end_11[0] - start_11[0]
#         lim11_y = end_11[1] - start_11[1]
#         lim22_x = end_22[0] - start_22[0]
#         lim22_y = end_22[1] - start_22[1]
#
#         ax.plot([start_11[0] - 0.1 * lim11_x, end_11[0] + 0.1 * lim11_x],
#                 [start_11[1] - 0.1 * lim11_y, end_11[1] + 0.1 * lim11_y],
#                 label='Principal 11-axis')
#         ax.plot([start_22[0] - 0.1 * lim22_x, end_22[0] + 0.1 * lim22_x],
#                 [start_22[1] - 0.1 * lim22_y, end_22[1] + 0.1 * lim22_y],
#                 label='Principal 22-axis')
#
#     # plot the locations of the various centroids
#     if centroids:
#         ax.scatter(0, 0, facecolors='None', edgecolors='k', marker='o',
#                    s=100, label='Elastic Centroid')
#         ax.scatter(self.x_se, self.y_se, c='k', marker='+', s=100,
#                    label='Shear Centre')
#
#         if (self.settings.plasticAnalysis):
#             ax.scatter(self.x_pc, self.y_pc, c='k', marker='x', s=100,
#                        label='Global Plastic Centroid')
#             ax.scatter(self.x1_pc, self.y2_pc, facecolors='None',
#                        edgecolors='k', marker='s', s=100,
#                        label='Principal Plastic Centroid')
#
#     # plot a contour of results defined by z
#     if z is not None:
#         cmap = cm.get_cmap(name='jet')
#
#         # if values are not all constant
#         if np.amax(z) - np.amin(z) > 1e-6:
#             v = np.linspace(np.amin(z), np.amax(z), 10, endpoint=True)
#         else:
#             # ten contours
#             v = 10
#
#         trictr = ax.tricontourf(
#             self.nodes[:, 0], self.nodes[:, 1], self.elements[:, 0:3],
#             z, v, cmap=cmap)
#         fig.colorbar(trictr, label='Stress')
#
#     # show the nodes
#     if nodes:
#         ax.plot(self.nodes[:, 0], self.nodes[:, 1], 'ko', markersize=1)
#
#     # show the legend
#     if (globalAxis or principalAxis or centroids):
#         ax.legend()
#
#     ax.grid(True)
#     plt.draw()  # render the figure
#     plt.pause(0.001)
#
#     return fig
#
# def quiverPlot(self, u, v, plotTitle=''):
#     """
#     This method produces a quiver plot of a vector with components u and
#     v, overlaid with the mesh.
#     """
#
#     fig, ax = plt.subplots()
#     plt.ion()  # interactive mode enabled
#     plt.show()  # show the plot
#     ax.set_aspect("equal")  # set the scale on the x and y axes equal
#
#     # plot the title and axis labels
#     ax.set_title(plotTitle)
#     ax.set_xlabel("x")
#     ax.set_ylabel("y")
#
#     # plot the mesh
#     plt.triplot(self.nodes[:, 0], self.nodes[:, 1], self.elements[:, 0:3],
#                 lw=0.5, color='black')
#
#     # scale the colour with respect to the magnitude of the vector
#     c = np.hypot(u, v)
#     cmap = cm.get_cmap(name='jet')
#
#     # generate the quiver plot and apply the colourbar
#     if np.amin(c) != np.amax(c):
#         # only show the quiver plot if there are results
#         quiv = ax.quiver(self.nodes[:, 0], self.nodes[:, 1], u, v, c,
#                          cmap=cmap)
#         v1 = np.linspace(np.amin(c), np.amax(c), 10, endpoint=True)
#         fig.colorbar(quiv, label='Stress', ticks=v1)
#
#     ax.grid(True)
#     plt.draw()  # render the figure
#     plt.pause(0.001)
#
#     return fig
