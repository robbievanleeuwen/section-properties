import numpy as np
from scipy.sparse import csc_matrix, coo_matrix, linalg
import matplotlib.pyplot as plt
import sectionproperties.analysis.fea as fea
import sectionproperties.analysis.solver as solver
import sectionproperties.post.post as post


class CrossSection:
    """Class for structural cross-sections.

    Stores the finite element geometry and mesh and provides methods to compute
    the cross-section properties. The element type used in this program is the
    six-noded quadratic triangular element.

    The constructor extracts information from the provided mesh object and
    creates and stores corresponding tri-6 finite element objects.

    :param geometry: Cross-section geometry object used to generate the mesh
    :type geometry: :class:`sectionproperties.pre.sections.SectionInput`
    :param mesh: Mesh object returned by meshpy
    :type mesh: :class:`meshpy.triangle.MeshInfo`

    The following example creates a CrossSection object of a 100D x 50W
    rectangle using a mesh size of 5::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.RectangularSection(d=100, b=50)
            mesh = geometry.create_mesh(mesh_sizes=[5])
            section = CrossSection(geometry, mesh)

    :cvar elements: List of finite element objects describing the cross-section
        mesh
    :vartype elements: list[:class:`sectionproperties.fea.Tri6`]
    :cvar int num_nodes: Number of nodes in the finite element mesh
    :cvar geometry: Cross-section geometry object used to generate the mesh
    :vartype geometry: :class:`sectionproperties.pre.sections.SectionInput`
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
        :class:`sectionproperties.cross_section.SectionProperties`
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

    def calculate_geometric_properties(self):
        """Calculates all the geometric properties of the cross-section and
        stores them in the
        :class:`sectionproperties.cross_section.SectionProperties` object
        contained in section_props.

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

        self.calculate_area()
        self.calculate_qx()
        self.calculate_qy()
        self.calculate_ixx_g()
        self.calculate_iyy_g()
        self.calculate_ixy_g()
        self.section_props.calculate_elastic_centroid()
        self.section_props.calculate_centroidal_properties(self.mesh)

    def calculate_warping_properties(self):
        """Calculates all the warping properties of the cross-section and
        stores them in the
        :class:`sectionproperties.cross_section.SectionProperties` object
        contained in section_props.

        * Torsion constant
        * Shear centre
        * Shear area
        * Warping constant

        Note that the geometric properties must be calculated first for the
        calculation of the warping properties to be correct::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
        """

        # TODO: if geometric properties are None, raise error

        # create a new CrossSection with the origin shifted to the centroid for
        # calculation of the warping properties such that the Lagrangian
        # multiplier approach can be utilised
        warping_section = CrossSection(self.get_j, self.mesh)

        # shift the coordinates of each element
        # N.B. the mesh class attribute remains unshifted!
        for el in warping_section.elements:
            el.coords[0, :] -= self.section_props.cx
            el.coords[1, :] -= self.section_props.cy

        # shift the mesh_nodes
        warping_section.mesh_nodes[:, 0] -= self.section_props.cx
        warping_section.mesh_nodes[:, 1] -= self.section_props.cy

        # assemble stiffness matrix and load vector for warping function
        (k, k_lg, f_torsion) = warping_section.assemble_torsion()

        # ILU decomposition on regular stiffness matrix
        k_precond = linalg.LinearOperator(
            (self.num_nodes, self.num_nodes), linalg.spilu(k).solve)

        # ILU decomposition on Lagrangian stiffness matrix
        k_lg_precond = linalg.LinearOperator(
            (self.num_nodes + 1, self.num_nodes + 1), linalg.spilu(k_lg).solve)

        # solve for warping constant
        omega = solver.solve_cgs(k, f_torsion, m=k_precond)
        self.section_props.omega = omega

        # determine the torsion constant
        self.section_props.j = (
            self.section_props.ixx_c + self.section_props.iyy_c - omega.dot(
                k.dot(np.transpose(omega))))

        # assemble shear function load vectors
        f_psi = np.zeros(self.num_nodes)
        f_phi = np.zeros(self.num_nodes)

        for el in warping_section.elements:
            (f_psi_el, f_phi_el) = el.shear_load_vectors(
                self.section_props.ixx_c, self.section_props.iyy_c,
                self.section_props.ixy_c)
            f_psi[el.node_ids] += f_psi_el
            f_phi[el.node_ids] += f_phi_el

        # solve for shear functions psi and phi
        psi_shear = solver.solve_cgs_lagrange(k_lg, f_psi, m=k_lg_precond)
        phi_shear = solver.solve_cgs_lagrange(k_lg, f_phi, m=k_lg_precond)
        self.section_props.psi_shear = psi_shear
        self.section_props.phi_shear = psi_shear

        # assemble shear centre and warping moment integrals
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

        # calculate effective Poissons ratio
        nu = 0  # TODO: implement

        # calculate shear centres (elasticity approach)
        Delta_s = 2 * (1 + nu) * (
            self.section_props.ixx_c * self.section_props.iyy_c -
            self.section_props.ixy_c ** 2)
        x_se = (1 / Delta_s) * ((nu / 2 * sc_xint) - f_torsion.dot(phi_shear))
        y_se = (1 / Delta_s) * ((nu / 2 * sc_yint) + f_torsion.dot(psi_shear))
        (x1_se, y2_se) = fea.principal_coordinate(self.section_props.phi, x_se,
                                                  y_se)

        # calculate shear centres (Trefftz's approach)
        x_st = (
            self.section_props.ixy_c * i_xomega - self.section_props.iyy_c *
            i_yomega) / (self.section_props.ixx_c * self.section_props.iyy_c -
                         self.section_props.ixy_c ** 2)
        y_st = (
            self.section_props.ixx_c * i_xomega - self.section_props.ixy_c *
            i_yomega) / (self.section_props.ixx_c * self.section_props.iyy_c -
                         self.section_props.ixy_c ** 2)

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

    def calculate_plastic_properties(self):
        """s"""

        # TODO: implement

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

    def plot_mesh(self, ax=None, pause=True):
        """Plots the finite element mesh. If no axes object is supplied a new
        figure and axis is created.

        :param ax: Axes object on which the mesh is plotted
        :type ax: :class:`matplotlib.axes.Axes`
        :param bool pause: If set to true, the figure pauses the script until
            the window is closed. If set to false, the script continues
            immediately after the window is rendered.
        """

        if ax is None:
            (fig, ax) = plt.subplots()

        post.setup_plot(ax, pause)

        ax.triplot(self.mesh_nodes[:, 0], self.mesh_nodes[:, 1],
                   self.mesh_elements[:, 0:3], lw=0.5, color='black')

        post.finish_plot(pause)

    def display_mesh_info(self):
        """Prints mesh statistics to the command window."""

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
            (i11_c, i22_c) = section.get_ipc()
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
            z22_minus = section.get_J()
        """

        return self.section_props.j

    def get_sc_e(self):
        """
        :return: Shear centre (elasticity approach) *(x_se, y_se)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            (x_se, y_se) = section.get_sc_e()
        """

        # add centroid location to move section back to original location
        x_se = self.section_props.x_se + self.section_props.cx
        y_se = self.section_props.y_se + self.section_props.cy

        return (x_se, y_se)

    def get_sc_p_e(self):
        """
        :return: Principal axis shear centre (elasticity approach)
            *(x1_se, y2_se)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            (x1_se, y2_se) = section.get_sc_p_e()
        """

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
        :return: Shear areas *(A_sx, A_sy, A_sxy)*
        :rtype: tuple(float, float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            (A_sx, A_sy, A_sxy) = section.get_As()
        """

        return (self.section_props.A_sx, self.section_props.A_sy,
                self.section_props.A_sxy)

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
    :cvar float qx: First moment of area of the element about the x-axis
    :cvar float qy: First moment of area of the element about the y-axis
    :cvar float ixx_g: Second moment of area of the element about the global
        x-axis
    :cvar float iyy_g: Second moment of area of the element about the global
        y-axis
    :cvar float ixy_g: Second moment of area of the element about the global
        xy-axis
    :cvar float ixx_c: Second moment of area of the element about the
        centroidal x-axis
    :cvar float iyy_c: Second moment of area of the element about the
        centroidal y-axis
    :cvar float ixy_c: Second moment of area of the element about the
        centroidal xy-axis
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
    :cvar float i11_c: Second moment of area of the element about the
        centroidal 11-axis
    :cvar float i22_c: Second moment of area of the element about the
        centroidal 22-axis
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
    :cvar float x_se: X shear centre (elasticity approach)
    :cvar float y_se: Y shear centre (elasticity approach)
    :cvar float x1_se: 11 shear centre (elasticity approach)
    :cvar float y2_se: 22 shear centre (elasticity approach)
    :cvar float x_st: X shear centre (Trefftz's approach)
    :cvar float y_st: Y shear centre (Trefftz's approach)
    :cvar float gamma: Warping constant
    :cvar float A_sx: Shear area about the x-axis
    :cvar float A_sy: Shear area about the y-axis
    :cvar float A_sxy: Shear area about the xy-axis
    :cvar float A_s11: Shear area about the 11 bending axis
    :cvar float A_s22: Shear area about the 22 bending axis
    """

    def __init__(self):
        """Inits the SectionProperties class."""

        self.area = None
        self.qx = None
        self.qy = None
        self.ixx_g = None
        self.iyy_g = None
        self.ixy_g = None
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
#     # calculate section stresses
#     if (self.settings.outputLog):
#         processText = "-- Calculating cross-section stresses..."
#     otherUtilities.functionTimer(
#         processText, self.calculateStress)
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
# def calculateStress(self):
#     """
#     This method calculates the cross-section stresses as a result of unit
#     loading(axial, moment, shear, torsion).
#     """
#
#     # allocate stress vectors
#     self.sigma_zz_axial = np.zeros(self.noNodes)
#     self.sigma_zz_bending_xx = np.zeros(self.noNodes)
#     self.sigma_zz_bending_yy = np.zeros(self.noNodes)
#     self.sigma_zz_bending_11 = np.zeros(self.noNodes)
#     self.sigma_zz_bending_22 = np.zeros(self.noNodes)
#     self.tau_zx_torsion = np.zeros(self.noNodes)
#     self.tau_zy_torsion = np.zeros(self.noNodes)
#     self.tau_zx_shear_x = np.zeros(self.noNodes)
#     self.tau_zy_shear_x = np.zeros(self.noNodes)
#     self.tau_zx_shear_y = np.zeros(self.noNodes)
#     self.tau_zy_shear_y = np.zeros(self.noNodes)
#
#     # allocate nodal weights vector for nodal averaging
#     nodal_weights = np.zeros(self.noNodes)
#
#     # loop through all elements in the mesh
#     for el in self.triElements:
#         # evaluate stresses at nodes
#         (elSigma_zz_axial, elSigma_zz_bending_xx, elSigma_zz_bending_yy,
#             elSigma_zz_bending_11, elSigma_zz_bending_22, elTau_zx_torsion,
#             elTau_zy_torsion, elTau_shear_zx_x, elTau_shear_zy_x,
#             elTau_shear_zx_y, elTau_shear_zy_y,
#             weights) = el.calculateStress(
#             self.area, self.ixx_c, self.iyy_c, self.ixy_c, self.i11_c,
#             self.i22_c, self.phi, self.omega[el.nodes], self.J,
#             self.Psi[el.nodes], self.Phi[el.nodes], self.Delta_s)
#
#         # add stresses to global vectors
#         self.sigma_zz_axial[el.nodes] += elSigma_zz_axial[:, 0] * weights
#         self.sigma_zz_bending_xx[el.nodes] += (
#             elSigma_zz_bending_xx * weights)
#         self.sigma_zz_bending_yy[el.nodes] += (
#             elSigma_zz_bending_yy * weights)
#         self.sigma_zz_bending_11[el.nodes] += (
#             elSigma_zz_bending_11 * weights)
#         self.sigma_zz_bending_22[el.nodes] += (
#             elSigma_zz_bending_22 * weights)
#         self.tau_zx_torsion[el.nodes] += elTau_zx_torsion * weights
#         self.tau_zy_torsion[el.nodes] += elTau_zy_torsion * weights
#         self.tau_zx_shear_x[el.nodes] += elTau_shear_zx_x * weights
#         self.tau_zy_shear_x[el.nodes] += elTau_shear_zy_x * weights
#         self.tau_zx_shear_y[el.nodes] += elTau_shear_zx_y * weights
#         self.tau_zy_shear_y[el.nodes] += elTau_shear_zy_y * weights
#
#         # increment the nodal count vector
#         nodal_weights[el.nodes] += weights
#
#     # nodal averaging
#     self.sigma_zz_axial *= 1 / nodal_weights
#     self.sigma_zz_bending_xx *= 1 / nodal_weights
#     self.sigma_zz_bending_yy *= 1 / nodal_weights
#     self.sigma_zz_bending_11 *= 1 / nodal_weights
#     self.sigma_zz_bending_22 *= 1 / nodal_weights
#     self.tau_zx_torsion *= 1 / nodal_weights
#     self.tau_zy_torsion *= 1 / nodal_weights
#     self.tau_zx_shear_x *= 1 / nodal_weights
#     self.tau_zy_shear_x *= 1 / nodal_weights
#     self.tau_zx_shear_y *= 1 / nodal_weights
#     self.tau_zy_shear_y *= 1 / nodal_weights
#
# def evaluateSectionStress(self, loadData):
#     """
#     This method scales the results obtained from the calculateStress
#     function by the design actions specified as input.
#     """
#
#     # scale unit stresses by design actions
#     self.axialStress = self.sigma_zz_axial * loadData.Nzz
#     self.bendingStress = (self.sigma_zz_bending_xx * loadData.Mxx +
#                           self.sigma_zz_bending_yy * loadData.Myy +
#                           self.sigma_zz_bending_11 * loadData.M11 +
#                           self.sigma_zz_bending_22 * loadData.M22)
#     self.torsionStress_zx = self.tau_zx_torsion * loadData.Mzz
#     self.torsionStress_zy = self.tau_zy_torsion * loadData.Mzz
#     self.torsionStress = ((self.torsionStress_zx ** 2 +
#                            self.torsionStress_zy ** 2) ** 0.5)
#     self.shearStress_zx = (self.tau_zx_shear_x * loadData.Vx +
#                            self.tau_zx_shear_y * loadData.Vy)
#     self.shearStress_zy = (self.tau_zy_shear_x * loadData.Vx +
#                            self.tau_zy_shear_y * loadData.Vy)
#     self.shearStress = ((self.shearStress_zx ** 2 +
#                          self.shearStress_zy ** 2) ** 0.5)
#
#     # compute combined stresses
#     self.sigma_zz = self.axialStress + self.bendingStress
#     self.tau_zx = self.torsionStress_zx + self.shearStress_zx
#     self.tau_zy = self.torsionStress_zy + self.shearStress_zy
#     self.tau = (self.tau_zx ** 2 + self.tau_zy ** 2) ** 0.5
#     self.vonMises = ((self.sigma_zz ** 2 + 3 * (self.tau ** 2)) ** 0.5)
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
#
# def printResults(self, fmt):
#     '''
#     This function prints the results of the geometric cross-sectional
#     analaysis to the console.
#     '''
#     print("\n-----------------------------")
#     print("Global xy Axis Properties")
#     print("-----------------------------")
#     print("Area\t = {:>{fmt}}".format(self.area, fmt=fmt))
#     print("Qx\t = {:>{fmt}}".format(self.qx, fmt=fmt))
#     print("Qy\t = {:>{fmt}}".format(self.qy, fmt=fmt))
#     print("cx\t = {:>{fmt}}".format(self.cx, fmt=fmt))
#     print("cy\t = {:>{fmt}}".format(self.cy, fmt=fmt))
#     print("Ixx_g\t = {:>{fmt}}".format(self.ixx_g, fmt=fmt))
#     print("Iyy_g\t = {:>{fmt}}".format(self.iyy_g, fmt=fmt))
#     print("Ixy_g\t = {:>{fmt}}\n".format(self.ixy_g, fmt=fmt))
#     print("\n-----------------------------")
#     print("Centroidal xy Axis Properties")
#     print("-----------------------------")
#     print("Ixx_c\t = {:>{fmt}}".format(self.ixx_c, fmt=fmt))
#     print("Iyy_c\t = {:>{fmt}}".format(self.iyy_c, fmt=fmt))
#     print("Ixy_c\t = {:>{fmt}}".format(self.ixy_c, fmt=fmt))
#     print("Zxx+\t = {:>{fmt}}".format(self.zxx_plus, fmt=fmt))
#     print("Zxx-\t = {:>{fmt}}".format(self.zxx_minus, fmt=fmt))
#     print("Zyy+\t = {:>{fmt}}".format(self.zyy_plus, fmt=fmt))
#     print("Zyy-\t = {:>{fmt}}".format(self.zyy_minus, fmt=fmt))
#     print("rx_c\t = {:>{fmt}}".format(self.rx_c, fmt=fmt))
#     print("ry_c\t = {:>{fmt}}\n".format(self.ry_c, fmt=fmt))
#     print("\n-----------------------------")
#     print("Principal Axis Properties")
#     print("-----------------------------")
#     print("phi\t = {:>{fmt}}".format(self.phi, fmt=fmt))
#     print("I11_c\t = {:>{fmt}}".format(self.i11_c, fmt=fmt))
#     print("I22_c\t = {:>{fmt}}".format(self.i22_c, fmt=fmt))
#     print("Z11+\t = {:>{fmt}}".format(self.z11_plus, fmt=fmt))
#     print("Z11-\t = {:>{fmt}}".format(self.z11_minus, fmt=fmt))
#     print("Z22+\t = {:>{fmt}}".format(self.z22_plus, fmt=fmt))
#     print("Z22-\t = {:>{fmt}}".format(self.z22_minus, fmt=fmt))
#     print("r1_c\t = {:>{fmt}}".format(self.r1_c, fmt=fmt))
#     print("r2_c\t = {:>{fmt}}\n".format(self.r2_c, fmt=fmt))
#     print("\n-----------------------------")
#     print("Torsional Properties")
#     print("-----------------------------")
#     print("J\t = {:>{fmt}}".format(self.J, fmt=fmt))
#     print("Iw\t = {:>{fmt}}\n".format(self.Gamma, fmt=fmt))
#     print("\n-----------------------------")
#     print("Shear Properties")
#     print("-----------------------------")
#     print("x_s,e\t = {:>{fmt}}".format(self.x_se, fmt=fmt))
#     print("y_s,e\t = {:>{fmt}}".format(self.y_se, fmt=fmt))
#     print("x_s,t\t = {:>{fmt}}".format(self.x_st, fmt=fmt))
#     print("y_s,t\t = {:>{fmt}}".format(self.y_st, fmt=fmt))
#     print("x1_s,e\t = {:>{fmt}}".format(self.x1_se, fmt=fmt))
#     print("y2_s,e\t = {:>{fmt}}".format(self.y2_se, fmt=fmt))
#     print("A_s,x\t = {:>{fmt}}".format(self.A_sx, fmt=fmt))
#     print("A_s,y\t = {:>{fmt}}".format(self.A_sy, fmt=fmt))
#     print("A_s,11\t = {:>{fmt}}".format(self.A_s11, fmt=fmt))
#     print("A_s,22\t = {:>{fmt}}\n".format(self.A_s22, fmt=fmt))
#
#     if (self.settings.plasticAnalysis):
#         print("\n-----------------------------")
#         print("Plastic Properties")
#         print("-----------------------------")
#         print("x_pc\t = {:>{fmt}}".format(self.x_pc, fmt=fmt))
#         print("y_pc\t = {:>{fmt}}".format(self.y_pc, fmt=fmt))
#         print("Sxx\t = {:>{fmt}}".format(self.Sxx, fmt=fmt))
#         print("Syy\t = {:>{fmt}}".format(self.Syy, fmt=fmt))
#         print("SF_xx+\t = {:>{fmt}}".format(self.SF_xx_plus, fmt=fmt))
#         print("SF_xx-\t = {:>{fmt}}".format(self.SF_xx_minus, fmt=fmt))
#         print("SF_yy+\t = {:>{fmt}}".format(self.SF_yy_plus, fmt=fmt))
#         print("SF_yy-\t = {:>{fmt}}".format(self.SF_yy_minus, fmt=fmt))
#         print("x1_pc\t = {:>{fmt}}".format(self.x1_pc, fmt=fmt))
#         print("y2_pc\t = {:>{fmt}}".format(self.y2_pc, fmt=fmt))
#         print("S11\t = {:>{fmt}}".format(self.S11, fmt=fmt))
#         print("S22\t = {:>{fmt}}".format(self.S22, fmt=fmt))
#         print("SF_11+\t = {:>{fmt}}".format(self.SF_11_plus, fmt=fmt))
#         print("SF_11-\t = {:>{fmt}}".format(self.SF_11_minus, fmt=fmt))
#         print("SF_22+\t = {:>{fmt}}".format(self.SF_22_plus, fmt=fmt))
#         print("SF_22-\t = {:>{fmt}}\n".format(self.SF_22_minus, fmt=fmt))
