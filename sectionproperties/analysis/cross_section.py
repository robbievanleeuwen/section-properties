import copy
import numpy as np
from scipy.sparse import csc_matrix, coo_matrix, linalg
from scipy.optimize import brentq
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.cm as cm
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
import meshpy.triangle as triangle
import sectionproperties.pre.pre as pre
import sectionproperties.analysis.fea as fea
import sectionproperties.analysis.solver as solver
import sectionproperties.post.post as post


class CrossSection:
    """Class for structural cross-sections.

    Stores the finite element geometry, mesh and material information and provides methods to
    compute the cross-section properties. The element type used in this program is the six-noded
    quadratic triangular element.

    The constructor extracts information from the provided mesh object and creates and stores the
    corresponding Tri6 finite element objects.

    :param geometry: Cross-section geometry object used to generate the mesh
    :type geometry: :class:`~sectionproperties.pre.sections.Geometry`
    :param mesh: Mesh object returned by meshpy
    :type mesh: :class:`meshpy.triangle.MeshInfo`
    :param materials: A list of material properties corresponding to various regions in the
        geometry and mesh. Note that if materials are specified, the number of material objects
        ust equal the number of regions in the geometry. If no materials are specified, only a
        purely geometric analysis can take place, and all regions will be assigned a default
        material with an elastic modulus and yield strength equal to 1, and a Poisson's ratio
        equal to 0.
    :type materials: list[:class:`~sectionproperties.pre.pre.Material`]
    :param bool time_info: If set to True, a detailed description of the computation and the time
        cost is printed to the terminal.

    The following example creates a :class:`~sectionproperties.analysis.cross_section.CrossSection`
    object of a 100D x 50W rectangle using a mesh size of 5::

        import sectionproperties.pre.sections as sections
        from sectionproperties.analysis.cross_section import CrossSection

        geometry = sections.RectangularSection(d=100, b=50)
        mesh = geometry.create_mesh(mesh_sizes=[5])
        section = CrossSection(geometry, mesh)

    The following example creates a 100D x 50W rectangle, with the top half of the section
    comprised of timber and the bottom half steel. The timber section is meshed with a maximum area
    of 10 and the steel section mesh with a maximum area of 5::

        import sectionproperties.pre.sections as sections
        from sectionproperties.pre.pre import Material
        from sectionproperties.analysis.cross_section import CrossSection

        geom_steel = sections.RectangularSection(d=50, b=50)
        geom_timber = sections.RectangularSection(d=50, b=50, shift=[0, 50])
        geometry = sections.MergedSection([geom_steel, geom_timber])
        geometry.clean_geometry()

        mesh = geometry.create_mesh(mesh_sizes=[5, 10])

        steel = Material(
            name='Steel', elastic_modulus=200e3, poissons_ratio=0.3, yield_strength=250,
            color='grey'
        )
        timber = Material(
            name='Timber', elastic_modulus=8e3, poissons_ratio=0.35, yield_strength=20,
            color='burlywood'
        )

        section = CrossSection(geometry, mesh, [steel, timber])
        section.plot_mesh(materials=True, alpha=0.5)

    :cvar elements: List of finite element objects describing the cross-section mesh
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
    :cvar materials: List of materials
    :type materials: list[:class:`~sectionproperties.pre.pre.Material`]
    :cvar material_groups: List of objects containing the elements in each defined material
    :type material_groups: list[:class:`~sectionproperties.pre.pre.MaterialGroup`]
    :cvar section_props: Class to store calculated section properties
    :vartype section_props: :class:`~sectionproperties.analysis.cross_section.SectionProperties`

    :raises AssertionError: If the number of materials does not equal the number of regions
    """

    def __init__(self, geometry, mesh, materials=None, time_info=False):
        """Inits the CrossSection class."""

        def init():
            self.geometry = geometry  # save geometry data

            # extract mesh data
            nodes = np.array(mesh.points, dtype=np.dtype(float))
            elements = np.array(mesh.elements, dtype=np.dtype(int))
            attributes = np.array(mesh.element_attributes, dtype=np.dtype(int))

            # swap mid-node order to retain node ordering consistency
            elements[:, [3, 4, 5]] = elements[:, [5, 3, 4]]

            # save total number of nodes in mesh
            self.num_nodes = len(nodes)

            # initialise material_sections variable
            self.material_groups = []

            # if materials are specified, check that the right number of material properties are
            # specified and then populate material_groups list
            if materials is not None:
                msg = "Number of materials ({0}), ".format(len(materials))
                msg += "should match the number of regions ({0}).".format(max(attributes) + 1)
                assert(len(materials) == max(attributes) + 1), msg

                # add a MaterialGroup object to the material_groups list for each uniquely
                # encountered material
                for (i, material) in enumerate(materials):
                    # add the first material to the list
                    if i == 0:
                        self.material_groups.append(MaterialGroup(material, self.num_nodes))
                    else:
                        # if the material hasn't been encountered
                        if material not in materials[:i]:
                            self.material_groups.append(MaterialGroup(material, self.num_nodes))
            # if there are no materials defined, add only the default material
            else:
                default_material = pre.Material('default', 1, 0, 1)
                self.material_groups.append(MaterialGroup(default_material, self.num_nodes))

            self.materials = materials  # save the input materials list

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
                coords = np.array([[x1, x2, x3, x4, x5, x6], [y1, y2, y3, y4, y5, y6]])

                # if materials are specified, get the material
                if materials is not None:
                    # get attribute index of current element
                    att_el = attributes[i]

                    # fetch the material
                    material = materials[att_el]
                # if there are no materials specified, use a default material
                else:
                    material = default_material

                # add tri6 elements to the mesh
                new_element = fea.Tri6(i, coords, node_ids, material)
                self.elements.append(new_element)

                # add element to relevant MaterialGroup
                for group in self.material_groups:
                    if material is group.material:
                        group.add_element(new_element)
                        break

            # save mesh input
            self.mesh = mesh
            self.mesh_nodes = nodes
            self.mesh_elements = elements
            self.mesh_attributes = attributes

            # initialise class storing section properties
            self.section_props = SectionProperties()

        if time_info:
            text = "--Initialising the CrossSection class..."
            solver.function_timer(text, init)
            print("")
        else:
            init()

    def calculate_geometric_properties(self, time_info=False):
        """Calculates the geometric properties of the cross-section and stores them in the
        :class:`~sectionproperties.analysis.cross_section.SectionProperties` object contained in
        the ``section_props`` class variable.

        :param bool time_info: If set to True, a detailed description of the computation and the
            time cost is printed to the terminal.

        The following geometric section properties are calculated:

        * Cross-sectional area
        * Cross-sectional perimeter
        * Modulus weighted area (axial rigidity)
        * First moments of area
        * Second moments of area about the global axis
        * Second moments of area about the centroidal axis
        * Elastic centroid
        * Centroidal section moduli
        * Radii of gyration
        * Principal axis properties

        If materials are specified for the cross-section, the moments of area and section moduli
        are elastic modulus weighted.

        The following example demonstrates the use of this method::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
        """

        def calculate_geom():
            # initialise properties
            self.section_props.area = 0
            self.section_props.perimeter = 0
            self.section_props.ea = 0
            self.section_props.ga = 0
            self.section_props.qx = 0
            self.section_props.qy = 0
            self.section_props.ixx_g = 0
            self.section_props.iyy_g = 0
            self.section_props.ixy_g = 0

            # calculate perimeter
            self.section_props.perimeter = self.geometry.calculate_perimeter()

            # calculate global geometric properties
            for el in self.elements:
                (area, qx, qy, ixx_g, iyy_g, ixy_g, e, g) = el.geometric_properties()

                self.section_props.area += area
                self.section_props.ea += area * e
                self.section_props.ga += area * g
                self.section_props.qx += qx * e
                self.section_props.qy += qy * e
                self.section_props.ixx_g += ixx_g * e
                self.section_props.iyy_g += iyy_g * e
                self.section_props.ixy_g += ixy_g * e

            self.section_props.nu_eff = self.section_props.ea / (2 * self.section_props.ga) - 1
            self.section_props.calculate_elastic_centroid()
            self.section_props.calculate_centroidal_properties(self.mesh)

        if time_info:
            text = "--Calculating geometric section properties..."
            solver.function_timer(text, calculate_geom)
            print("")
        else:
            calculate_geom()

    def calculate_warping_properties(self, time_info=False, solver_type='direct'):
        """Calculates all the warping properties of the cross-section and stores them in the
        :class:`~sectionproperties.analysis.cross_section.SectionProperties` object contained in
        the ``section_props`` class variable.

        :param bool time_info: If set to True, a detailed description of the computation and the
            time cost is printed to the terminal.
        :param string solver_type: Solver used for solving systems of linear equations, either
            using the *'direct'* method or *'cgs'* iterative method

        The following warping section properties are calculated:

        * Torsion constant
        * Shear centre
        * Shear area
        * Warping constant
        * Monosymmetry constant

        If materials are specified, the values calculated for the torsion constant, warping
        constant and shear area are elastic modulus weighted.

        Note that the geometric properties must be calculated first for the calculation of the
        warping properties to be correct::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()

        :raises RuntimeError: If the geometric properties have not been
            calculated prior to calling this method
        """

        # check that a geometric analysis has been performed
        if None in [self.section_props.area, self.section_props.ixx_c, self.section_props.cx]:
            err = "Calculate geometric properties before performing a warping analysis."
            raise RuntimeError(err)

        # create a new CrossSection with the origin shifted to the centroid for calculation of the
        # warping properties such that the Lagrangian multiplier approach can be utilised
        warping_section = CrossSection(self.geometry, self.mesh, self.materials)

        # shift the coordinates of each element N.B. the mesh class attribute remains unshifted!
        for el in warping_section.elements:
            el.coords[0, :] -= self.section_props.cx
            el.coords[1, :] -= self.section_props.cy

        # assemble stiffness matrix and load vector for warping function
        if time_info:
            text = "--Assembling {0}x{0} stiffness matrix and load vector...".format(self.num_nodes)
            (k, k_lg, f_torsion) = solver.function_timer(text, warping_section.assemble_torsion)
        else:
            (k, k_lg, f_torsion) = warping_section.assemble_torsion()

        # ILU decomposition of stiffness matrices
        def ilu_decomp():
            # ILU decomposition on regular stiffness matrix
            k_precond = linalg.LinearOperator(
                (self.num_nodes, self.num_nodes), linalg.spilu(k).solve
            )

            # ILU decomposition on Lagrangian stiffness matrix
            k_lg_precond = linalg.LinearOperator(
                (self.num_nodes + 1, self.num_nodes + 1), linalg.spilu(k_lg).solve
            )

            return (k_precond, k_lg_precond)

        # if the cgs method is used, perform ILU decomposition
        if solver_type == 'cgs':
            if time_info:
                text = "--Performing ILU decomposition on the stiffness matrices..."
                (k_precond, k_lg_precond) = solver.function_timer(text, ilu_decomp)
            else:
                (k_precond, k_lg_precond) = ilu_decomp()

        # solve for warping function
        def solve_warping():
            if solver_type == 'cgs':
                omega = solver.solve_cgs(k, f_torsion, k_precond)
            elif solver_type == 'direct':
                omega = solver.solve_direct(k, f_torsion)

            return omega

        if time_info:
            text = "--Solving for the warping function using the {0} solver...".format(solver_type)
            omega = solver.function_timer(text, solve_warping)
        else:
            omega = solve_warping()

        # save the warping function
        self.section_props.omega = omega

        # determine the torsion constant
        def j_func():
            return (
                self.section_props.ixx_c + self.section_props.iyy_c
                - omega.dot(k.dot(np.transpose(omega)))
            )

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
                    self.section_props.ixy_c, self.section_props.nu_eff)
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
                psi_shear = solver.solve_cgs_lagrange(k_lg, f_psi, m=k_lg_precond)
                phi_shear = solver.solve_cgs_lagrange(k_lg, f_phi, m=k_lg_precond)
            elif solver_type == 'direct':
                psi_shear = solver.solve_direct_lagrange(k_lg, f_psi)
                phi_shear = solver.solve_direct_lagrange(k_lg, f_phi)

            return (psi_shear, phi_shear)

        if time_info:
            text = "--Solving for the shear functions using the {0} solver...".format(solver_type)
            (psi_shear, phi_shear) = solver.function_timer(text, solve_shear_functions)
        else:
            (psi_shear, phi_shear) = solve_shear_functions()

        # save the shear functions
        self.section_props.psi_shear = psi_shear
        self.section_props.phi_shear = phi_shear

        # assemble shear centre and warping moment integrals
        def assemble_sc_warping_integrals():
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
                    self.section_props.ixy_c, omega[el.node_ids]
                )

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
                solver.function_timer(text, assemble_sc_warping_integrals))
        else:
            (sc_xint, sc_yint, q_omega, i_omega, i_xomega, i_yomega) = (
                assemble_sc_warping_integrals())

        # calculate shear centres
        def shear_centres():
            # calculate shear centres (elasticity approach)
            Delta_s = (
                2 * (1 + self.section_props.nu_eff) * (
                    self.section_props.ixx_c * self.section_props.iyy_c
                    - self.section_props.ixy_c ** 2)
            )
            x_se = (
                (1 / Delta_s) * ((self.section_props.nu_eff / 2
                                  * sc_xint) - f_torsion.dot(phi_shear))
            )
            y_se = (
                (1 / Delta_s) * ((self.section_props.nu_eff / 2
                                  * sc_yint) + f_torsion.dot(psi_shear))
            )
            (x11_se, y22_se) = fea.principal_coordinate(self.section_props.phi, x_se, y_se)

            # calculate shear centres (Trefftz's approach)
            x_st = (
                (self.section_props.ixy_c * i_xomega - self.section_props.iyy_c * i_yomega) / (
                    self.section_props.ixx_c * self.section_props.iyy_c
                    - self.section_props.ixy_c ** 2)
            )
            y_st = (
                (self.section_props.ixx_c * i_xomega - self.section_props.ixy_c * i_yomega) / (
                    self.section_props.ixx_c * self.section_props.iyy_c
                    - self.section_props.ixy_c ** 2)
            )

            return (Delta_s, x_se, y_se, x11_se, y22_se, x_st, y_st)

        if time_info:
            text = "--Calculating shear centres..."
            (Delta_s, x_se, y_se, x11_se, y22_se, x_st, y_st) = solver.function_timer(
                text, shear_centres)
        else:
            (Delta_s, x_se, y_se, x11_se, y22_se, x_st, y_st) = shear_centres()

        # save shear centres
        self.section_props.Delta_s = Delta_s
        self.section_props.x_se = x_se
        self.section_props.y_se = y_se
        self.section_props.x11_se = x11_se
        self.section_props.y22_se = y22_se
        self.section_props.x_st = x_st
        self.section_props.y_st = y_st

        # calculate warping constant
        self.section_props.gamma = (
            i_omega - q_omega ** 2 / self.section_props.ea - y_se * i_xomega + x_se * i_yomega
        )

        def assemble_shear_deformation():
            # assemble shear deformation coefficients
            kappa_x = 0
            kappa_y = 0
            kappa_xy = 0

            for el in warping_section.elements:
                (kappa_x_el, kappa_y_el, kappa_xy_el) = el.shear_coefficients(
                    self.section_props.ixx_c, self.section_props.iyy_c,
                    self.section_props.ixy_c, psi_shear[el.node_ids], phi_shear[el.node_ids],
                    self.section_props.nu_eff
                )

                kappa_x += kappa_x_el
                kappa_y += kappa_y_el
                kappa_xy += kappa_xy_el

            return (kappa_x, kappa_y, kappa_xy)

        if time_info:
            text = "--Assembling shear deformation coefficients..."
            (kappa_x, kappa_y, kappa_xy) = solver.function_timer(text, assemble_shear_deformation)
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
        R = np.array([
            [np.cos(phi_rad), np.sin(phi_rad)],
            [-np.sin(phi_rad), np.cos(phi_rad)]
        ])

        rotatedAlpha = R.dot(np.array([
            [alpha_xx, alpha_xy],
            [alpha_xy, alpha_yy]
        ])).dot(np.transpose(R))

        # recalculate the shear area based on the rotated alpha value
        self.section_props.A_s11 = self.section_props.area / rotatedAlpha[0, 0]
        self.section_props.A_s22 = self.section_props.area / rotatedAlpha[1, 1]

        # calculate the monosymmetry consants
        def calculate_monosymmetry_integrals():
            int_x = 0
            int_y = 0
            int_11 = 0
            int_22 = 0

            for el in warping_section.elements:
                (int_x_el, int_y_el, int_11_el, int_22_el) = el.monosymmetry_integrals(
                    self.section_props.phi
                )

                int_x += int_x_el
                int_y += int_y_el
                int_11 += int_11_el
                int_22 += int_22_el

            return (int_x, int_y, int_11, int_22)

        if time_info:
            text = "--Assembling monosymmetry integrals..."
            (int_x, int_y, int_11, int_22) = solver.function_timer(
                text, calculate_monosymmetry_integrals
            )
            print("")
        else:
            (int_x, int_y, int_11, int_22) = calculate_monosymmetry_integrals()

        # calculate the monosymmetry constants
        self.section_props.beta_x_plus = (
            -int_x / self.section_props.ixx_c + 2 * self.section_props.y_se
        )
        self.section_props.beta_x_minus = (
            int_x / self.section_props.ixx_c - 2 * self.section_props.y_se
        )
        self.section_props.beta_y_plus = (
            -int_y / self.section_props.iyy_c + 2 * self.section_props.x_se
        )
        self.section_props.beta_y_minus = (
            int_y / self.section_props.iyy_c - 2 * self.section_props.x_se
        )
        self.section_props.beta_11_plus = (
            -int_11 / self.section_props.i11_c + 2 * self.section_props.y22_se
        )
        self.section_props.beta_11_minus = (
            int_11 / self.section_props.i11_c - 2 * self.section_props.y22_se
        )
        self.section_props.beta_22_plus = (
            -int_22 / self.section_props.i22_c + 2 * self.section_props.x11_se
        )
        self.section_props.beta_22_minus = (
            int_22 / self.section_props.i22_c - 2 * self.section_props.x11_se
        )

    def calculate_frame_properties(self, time_info=False, solver_type='direct'):
        """Calculates and returns the properties required for a frame analysis. The properties are
        also stored in the :class:`~sectionproperties.analysis.cross_section.SectionProperties`
        object contained in the ``section_props`` class variable.

        :param bool time_info: If set to True, a detailed description of the computation and the
            time cost is printed to the terminal.
        :param string solver_type: Solver used for solving systems of linear equations, either
            using the *'direct'* method or *'cgs'* iterative method

        :return: Cross-section properties to be used for a frame analysis *(area, ixx, iyy, ixy, j,
            phi)*
        :rtype: tuple(float, float, float, float, float, float)

        The following section properties are calculated:

        * Cross-sectional area *(area)*
        * Second moments of area about the centroidal axis *(ixx, iyy, ixy)*
        * Torsion constant *(j)*
        * Principal axis angle *(phi)*

        If materials are specified for the cross-section, the area, second moments of area and
        torsion constant are elastic modulus weighted.

        The following example demonstrates the use of this method::

            section = CrossSection(geometry, mesh)
            (area, ixx, iyy, ixy, j, phi) = section.calculate_frame_properties()
        """

        def calculate_frame():
            # initialise geometric properties
            self.section_props.area = 0
            self.section_props.ea = 0
            self.section_props.qx = 0
            self.section_props.qy = 0
            self.section_props.ixx_g = 0
            self.section_props.iyy_g = 0
            self.section_props.ixy_g = 0
            self.section_props.ixx_c = 0
            self.section_props.iyy_c = 0
            self.section_props.ixy_c = 0
            self.section_props.j = 0
            self.section_props.phi = 0

            # calculate global geometric properties
            for el in self.elements:
                (area, qx, qy, ixx_g,
                 iyy_g, ixy_g, e, _) = el.geometric_properties()

                self.section_props.area += area
                self.section_props.ea += area * e
                self.section_props.qx += qx * e
                self.section_props.qy += qy * e
                self.section_props.ixx_g += ixx_g * e
                self.section_props.iyy_g += iyy_g * e
                self.section_props.ixy_g += ixy_g * e

            # calculate elastic centroid location
            self.section_props.calculate_elastic_centroid()

            # calculate second moments of area about the centroidal xy axis
            self.section_props.ixx_c = (
                self.section_props.ixx_g - self.section_props.qx ** 2 / self.section_props.ea
            )
            self.section_props.iyy_c = (
                self.section_props.iyy_g - self.section_props.qy ** 2 / self.section_props.ea
            )
            self.section_props.ixy_c = (
                self.section_props.ixy_g - self.section_props.qx * self.section_props.qy
                / self.section_props.ea
            )

            # calculate the principal axis angle
            Delta = (
                ((self.section_props.ixx_c - self.section_props.iyy_c) / 2) ** 2
                + self.section_props.ixy_c ** 2
            ) ** 0.5

            i11_c = (
                (self.section_props.ixx_c + self.section_props.iyy_c) / 2 + Delta
            )

            # calculate initial principal axis angle
            if abs(self.section_props.ixx_c - i11_c) < 1e-12 * i11_c:
                self.section_props.phi = 0
            else:
                self.section_props.phi = np.arctan2(
                    self.section_props.ixx_c - i11_c,
                    self.section_props.ixy_c
                ) * 180 / np.pi

            # create a new CrossSection with the origin shifted to the centroid for calculation of
            # the warping properties
            warping_section = CrossSection(self.geometry, self.mesh, self.materials)

            # shift the coordinates of each element N.B. the mesh class attribute remains unshifted
            for el in warping_section.elements:
                el.coords[0, :] -= self.section_props.cx
                el.coords[1, :] -= self.section_props.cy

            (k, _, f) = warping_section.assemble_torsion(lg=False)

            # if the cgs method is used, perform ILU decomposition
            if solver_type == 'cgs':
                k_precond = linalg.LinearOperator(
                    (self.num_nodes, self.num_nodes), linalg.spilu(k).solve
                )

            # solve for warping function
            if solver_type == 'cgs':
                omega = solver.solve_cgs(k, f, k_precond)
            elif solver_type == 'direct':
                omega = solver.solve_direct(k, f)

            # calculate the torsion constant
            self.section_props.j = (
                self.section_props.ixx_c + self.section_props.iyy_c - omega.dot(k.dot(
                    np.transpose(omega)))
            )

        if time_info:
            text = "--Calculating frame section properties..."
            solver.function_timer(text, calculate_frame)
            print("")
        else:
            calculate_frame()

        return (
            self.section_props.ea, self.section_props.ixx_c,
            self.section_props.iyy_c, self.section_props.ixy_c,
            self.section_props.j, self.section_props.phi)

    def calculate_plastic_properties(self, time_info=False, verbose=False, debug=False):
        """Calculates the plastic properties of the cross-section and stores the, in the
        :class:`~sectionproperties.analysis.cross_section.SectionProperties` object contained in
        the ``section_props`` class variable.

        :param bool time_info: If set to True, a detailed description of the computation and the
            time cost is printed to the terminal.
        :param bool verbose: If set to True, the number of iterations required for each plastic
            axis is printed to the terminal.
        :param bool debug: If set to True, the geometry is plotted each time a new mesh is
            generated by the plastic centroid algorithm.

        The following warping section properties are calculated:

        * Plastic centroid for bending about the centroidal and principal axes
        * Plastic section moduli for bending about the centroidal and principal axes
        * Shape factors for bending about the centroidal and principal axes

        If materials are specified for the cross-section, the plastic section moduli are displayed
        as plastic moments (i.e :math:`M_p = f_y S`) and the shape factors are not calculated.

        Note that the geometric properties must be calculated before the plastic properties are
        calculated::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_plastic_properties()

        :raises RuntimeError: If the geometric properties have not been calculated prior to calling
            this method
        """

        # check that a geometric analysis has been performed
        if self.section_props.cx is None:
            err = "Calculate geometric properties before performing a plastic analysis."
            raise RuntimeError(err)

        def calc_plastic():
            plastic_section = PlasticSection(self.geometry, self.materials, debug)

            # calculate plastic properties
            try:
                plastic_section.calculate_plastic_properties(self, verbose)
            except ValueError as exp:
                msg = "Plastic section properties calculation failed. Contact "
                msg += "robbie.vanleeuwen@gmail.com with your analysis parameters."
                raise RuntimeError(msg) from exp

        if time_info:
            text = "--Calculating plastic properties..."
            solver.function_timer(text, calc_plastic)
            print("")
        else:
            calc_plastic()

    def calculate_stress(self, N=0, Vx=0, Vy=0, Mxx=0, Myy=0, M11=0, M22=0,
                         Mzz=0, time_info=False):
        """Calculates the cross-section stress resulting from design actions and returns a
        :class:`~sectionproperties.analysis.cross_section.StressPost` object allowing
        post-processing of the stress results.

        :param float N: Axial force
        :param float Vx: Shear force acting in the x-direction
        :param float Vy: Shear force acting in the y-direction
        :param float Mxx: Bending moment about the centroidal xx-axis
        :param float Myy: Bending moment about the centroidal yy-axis
        :param float M11: Bending moment about the centroidal 11-axis
        :param float M22: Bending moment about the centroidal 22-axis
        :param float Mzz: Torsion moment about the centroidal zz-axis
        :param bool time_info: If set to True, a detailed description of the computation and the
            time cost is printed to the terminal.
        :return: Object for post-processing cross-section stresses
        :rtype: :class:`~sectionproperties.analysis.cross_section.StressPost`

        Note that a geometric and warping analysis must be performed before a stress analysis is
        carried out::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(N=1e3, Vy=3e3, Mxx=1e6)

        :raises RuntimeError: If a geometric and warping analysis have not been performed prior to
            calling this method
        """

        # check that a geometric and warping analysis has been performed
        if None in [
            self.section_props.area, self.section_props.ixx_c, self.section_props.cx,
            self.section_props.j
        ]:
            err = "Perform a geometric and warping analysis before carrying out a stress analysis."
            raise RuntimeError(err)

        def calc_stress():
            # create stress post object
            stress_post = StressPost(self)

            # get relevant section properties
            ea = self.section_props.ea
            cx = self.section_props.cx
            cy = self.section_props.cy
            ixx = self.section_props.ixx_c
            iyy = self.section_props.iyy_c
            ixy = self.section_props.ixy_c
            i11 = self.section_props.i11_c
            i22 = self.section_props.i22_c
            phi = self.section_props.phi
            j = self.section_props.j
            Delta_s = self.section_props.Delta_s
            nu = self.section_props.nu_eff

            # loop through all material groups
            for group in stress_post.material_groups:
                # allocate nodal weights vector for nodal averaging
                nodal_weights = np.zeros(self.num_nodes)

                # loop through all elements in the material group
                for el in group.elements:
                    (
                        sig_zz_n_el, sig_zz_mxx_el, sig_zz_myy_el, sig_zz_m11_el, sig_zz_m22_el,
                        sig_zx_mzz_el, sig_zy_mzz_el, sig_zx_vx_el, sig_zy_vx_el, sig_zx_vy_el,
                        sig_zy_vy_el, weights
                    ) = el.element_stress(
                        N, Mxx, Myy, M11, M22, Mzz, Vx, Vy, ea, cx, cy, ixx, iyy, ixy, i11, i22,
                        phi, j, nu, self.section_props.omega[el.node_ids],
                        self.section_props.psi_shear[el.node_ids],
                        self.section_props.phi_shear[el.node_ids], Delta_s
                    )

                    # add stresses to global vectors
                    group.stress_result.sig_zz_n[el.node_ids] += sig_zz_n_el * weights
                    group.stress_result.sig_zz_mxx[el.node_ids] += sig_zz_mxx_el * weights
                    group.stress_result.sig_zz_myy[el.node_ids] += sig_zz_myy_el * weights
                    group.stress_result.sig_zz_m11[el.node_ids] += sig_zz_m11_el * weights
                    group.stress_result.sig_zz_m22[el.node_ids] += sig_zz_m22_el * weights
                    group.stress_result.sig_zx_mzz[el.node_ids] += sig_zx_mzz_el * weights
                    group.stress_result.sig_zy_mzz[el.node_ids] += sig_zy_mzz_el * weights
                    group.stress_result.sig_zx_vx[el.node_ids] += sig_zx_vx_el * weights
                    group.stress_result.sig_zy_vx[el.node_ids] += sig_zy_vx_el * weights
                    group.stress_result.sig_zx_vy[el.node_ids] += sig_zx_vy_el * weights
                    group.stress_result.sig_zy_vy[el.node_ids] += sig_zy_vy_el * weights

                    # add nodal weights
                    nodal_weights[el.node_ids] += weights

                # nodal averaging
                for (i, weight) in enumerate(nodal_weights):
                    if weight != 0:
                        group.stress_result.sig_zz_n[i] *= 1 / weight
                        group.stress_result.sig_zz_mxx[i] *= 1 / weight
                        group.stress_result.sig_zz_myy[i] *= 1 / weight
                        group.stress_result.sig_zz_m11[i] *= 1 / weight
                        group.stress_result.sig_zz_m22[i] *= 1 / weight
                        group.stress_result.sig_zx_mzz[i] *= 1 / weight
                        group.stress_result.sig_zy_mzz[i] *= 1 / weight
                        group.stress_result.sig_zx_vx[i] *= 1 / weight
                        group.stress_result.sig_zy_vx[i] *= 1 / weight
                        group.stress_result.sig_zx_vy[i] *= 1 / weight
                        group.stress_result.sig_zy_vy[i] *= 1 / weight

                # calculate combined stresses
                group.stress_result.calculate_combined_stresses()

            return stress_post

        if time_info:
            text = "--Calculating cross-section stresses..."
            stress_post = solver.function_timer(text, calc_stress)
            print("")
        else:
            stress_post = calc_stress()

        # return the stress_post object
        return stress_post

    def assemble_torsion(self, lg=True):
        """Assembles stiffness matrices to be used for the computation of warping properties and
        the torsion load vector (f_torsion). Both a regular (k) and Lagrangian multiplier (k_lg)
        stiffness matrix are returned. The stiffness matrices are assembled using the sparse COO
        format and returned in the sparse CSC format.

        :param bool lg: Whether or not to calculate the Lagrangian multiplier stiffness matrix

        :return: Regular stiffness matrix, Lagrangian multiplier stiffness matrix and torsion load
            vector *(k, k_lg, f_torsion)*
        :rtype: tuple(:class:`scipy.sparse.csc_matrix`, :class:`scipy.sparse.csc_matrix`,
            :class:`numpy.ndarray`)
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

        if not lg:
            return (csc_matrix(k), None, f_torsion)

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

        k_lg = coo_matrix((data, (row, col)), shape=(N + 1, N + 1))

        return (csc_matrix(k), csc_matrix(k_lg), f_torsion)

    def plot_mesh(self, ax=None, pause=True, alpha=1, materials=False, mask=None):
        """Plots the finite element mesh. If no axes object is supplied a new figure and axis is
        created.

        :param ax: Axes object on which the mesh is plotted
        :type ax: :class:`matplotlib.axes.Axes`
        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.
        :param float alpha: Transparency of the mesh outlines: :math:`0 \leq \\alpha \leq 1`
        :param bool materials: If set to true and material properties have been provided to the
            :class:`~sectionproperties.analysis.cross_section.CrossSection` object, shades the
            elements with the specified material colours
        :param mask: Mask array, of length ``num_nodes``, to mask out triangles
        :type mask: list[bool]

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots the mesh generated for the second example
        listed under the :class:`~sectionproperties.analysis.cross_section.CrossSection` object
        definition::

            import sectionproperties.pre.sections as sections
            from sectionproperties.pre.pre import Material
            from sectionproperties.analysis.cross_section import CrossSection

            geom_steel = sections.RectangularSection(d=50, b=50)
            geom_timber = sections.RectangularSection(d=50, b=50, shift=[50, 0])
            geometry = sections.MergedSection([geom_steel, geom_timber])
            geometry.clean_geometry()

            mesh = geometry.create_mesh(mesh_sizes=[5, 10])

            steel = Material(
                name='Steel', elastic_modulus=200e3, poissons_ratio=0.3, yield_strength=250,
                color='grey'
            )
            timber = Material(
                name='Timber', elastic_modulus=8e3, poissons_ratio=0.35, yield_strength=20,
                color='burlywood'
            )

            section = CrossSection(geometry, mesh, [steel, timber])
            section.plot_mesh(materials=True, alpha=0.5)

        ..  figure:: ../images/composite_mesh.png
            :align: center
            :scale: 75 %

            Finite element mesh generated by the above example.
        """

        # if no axes object is supplied, create and setup the plot
        if ax is None:
            ax_supplied = False
            (fig, ax) = plt.subplots()
            post.setup_plot(pause)
        else:
            ax_supplied = True

        # plot the mesh
        ax.triplot(
            self.mesh_nodes[:, 0], self.mesh_nodes[:, 1], self.mesh_elements[:, 0:3], lw=0.5,
            color='black', alpha=alpha, mask=mask
        )

        # if the material colours are to be displayed
        if materials and self.materials is not None:
            color_array = []
            legend_list = []

            # create an array of finite element colours
            for element in self.elements:
                color_array.append(element.material.color)

            # create a list of unique material legend entries
            for (i, material) in enumerate(self.materials):
                # if the material has not be entered yet
                if i == 0 or material not in self.materials[0:i]:
                    # add the material colour and name to the legend list
                    legend_list.append(mpatches.Patch(color=material.color, label=material.name))

            cmap = ListedColormap(color_array)  # custom colormap
            c = np.arange(len(color_array))  # indices of elements

            # plot the mesh colours
            ax.tripcolor(
                self.mesh_nodes[:, 0], self.mesh_nodes[:, 1], self.mesh_elements[:, 0:3], c,
                cmap=cmap
            )

            # display the legend
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), handles=legend_list)

        if ax_supplied:
            # if an axis is supplied, return None for figure and axes to indicate that it is not
            # yet finished
            return None, None

        # if no axes object is supplied, finish the plot and return the figure and axes
        post.finish_plot(ax, pause, title='Finite Element Mesh')
        return (fig, ax)

    def plot_centroids(self, pause=True):
        """Plots the elastic centroid, the shear centre, the plastic centroids and the principal
        axis, if they have been calculated, on top of the finite element mesh.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example analyses a 200 PFC section and displays a plot of
        the centroids::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.PfcSection(d=200, b=75, t_f=12, t_w=6, r=12, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            section.calculate_plastic_properties()

            section.plot_centroids()

        ..  figure:: ../images/pfc_centroids.png
            :align: center
            :scale: 75 %

            Plot of the centroids generated by the above example.

        The following example analyses a 150x90x12 UA section and displays a plot of the
        centroids::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            section.calculate_plastic_properties()

            section.plot_centroids()

        ..  figure:: ../images/angle_centroids.png
            :align: center
            :scale: 75 %

            Plot of the centroids generated by the above example.
        """

        # create plot and setup the plot
        (fig, ax) = plt.subplots()
        post.setup_plot(pause)

        # plot the finite element mesh
        self.plot_mesh(ax, pause, alpha=0.5)

        # if the elastic centroid has been calculated
        if self.section_props.cx is not None:
            ax.scatter(
                self.section_props.cx, self.section_props.cy, edgecolors='r', facecolors='none',
                marker='o', s=100, label='Elastic centroid'
            )

        # if the shear centre has been calculated
        if self.section_props.x_se is not None:
            (x_s, y_s) = self.get_sc()
            ax.scatter(x_s, y_s, c='r', marker='+', s=100, label='Shear centre')

        # if the global plastic centroid has been calculated
        if self.section_props.x_pc is not None:
            (x_pc, y_pc) = self.get_pc()
            ax.scatter(x_pc, y_pc, c='r', marker='x', s=100, label='Global plastic centroid')

        # if the principal plastic centroid has been calculated
        if self.section_props.x11_pc is not None:
            (x11_pc, y22_pc) = self.get_pc_p()
            ax.scatter(
                x11_pc, y22_pc, edgecolors='r', facecolors='none', marker='s', s=100,
                label='Principal plastic centroid'
            )

        # if the principal axis has been calculated
        if self.section_props.phi is not None:
            post.draw_principal_axis(
                ax, self.section_props.phi * np.pi / 180, self.section_props.cx,
                self.section_props.cy
            )

        # display the legend
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        # finish the plot
        post.finish_plot(ax, pause, title='Centroids')

        return (fig, ax)

    def display_mesh_info(self):
        """Prints mesh statistics (number of nodes, elements and regions) to the command window.

        The following example displays the mesh statistics for a Tee section merged from two
        rectangles::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            rec1 = sections.RectangularSection(d=100, b=25, shift=[-12.5, 0])
            rec2 = sections.RectangularSection(d=25, b=100, shift=[-50, 100])
            geometry = sections.MergedSection([rec1, rec2])
            mesh = geometry.create_mesh(mesh_sizes=[5, 2.5])
            section = CrossSection(geometry, mesh)
            section.display_mesh_info()

            >>>Mesh Statistics:
            >>>--4920 nodes
            >>>--2365 elements
            >>>--2 regions
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

        :param string fmt: Number formatting string

        The following example displays the geometric section properties for a 100D x 50W rectangle
        with three digits after the decimal point::

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

    def get_perimeter(self):
        """
        :return: Cross-section perimeter
        :rtype: float

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            perimeter = section.get_perimeter()
        """

        return self.section_props.perimeter

    def get_ea(self):
        """
        :return: Modulus weighted area (axial rigidity)
        :rtype: float

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            ea = section.get_ea()
        """

        return self.section_props.ea

    def get_q(self):
        """
        :return: First moments of area about the global axis *(qx, qy)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            (qx, qy) = section.get_q()
        """

        return (self.section_props.qx, self.section_props.qy)

    def get_ig(self):
        """
        :return: Second moments of area about the global axis *(ixx_g, iyy_g, ixy_g)*
        :rtype: tuple(float, float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            (ixx_g, iyy_g, ixy_g) = section.get_ig()
        """

        return (self.section_props.ixx_g, self.section_props.iyy_g, self.section_props.ixy_g)

    def get_c(self):
        """
        :return: Elastic centroid *(cx, cy)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            (cx, cy) = section.get_c()
        """

        return (self.section_props.cx, self.section_props.cy)

    def get_ic(self):
        """
        :return: Second moments of area centroidal axis *(ixx_c, iyy_c, ixy_c)*
        :rtype: tuple(float, float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            (ixx_c, iyy_c, ixy_c) = section.get_ic()
        """

        return (self.section_props.ixx_c, self.section_props.iyy_c, self.section_props.ixy_c)

    def get_z(self):
        """
        :return: Elastic section moduli about the centroidal axis with respect to the top and
            bottom fibres *(zxx_plus, zxx_minus, zyy_plus, zyy_minus)*
        :rtype: tuple(float, float, float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            (zxx_plus, zxx_minus, zyy_plus, zyy_minus) = section.get_z()
        """

        return (
            self.section_props.zxx_plus, self.section_props.zxx_minus, self.section_props.zyy_plus,
            self.section_props.zyy_minus
        )

    def get_rc(self):
        """
        :return: Radii of gyration about the centroidal axis *(rx, ry)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            (rx, ry) = section.get_rc()
        """

        return (self.section_props.rx_c, self.section_props.ry_c)

    def get_ip(self):
        """
        :return: Second moments of area about the principal axis *(i11_c, i22_c)*
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
        :return: Elastic section moduli about the principal axis with respect to the top and bottom
            fibres *(z11_plus, z11_minus, z22_plus, z22_minus)*
        :rtype: tuple(float, float, float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            (z11_plus, z11_minus, z22_plus, z22_minus) = section.get_zp()
        """

        return (
            self.section_props.z11_plus, self.section_props.z11_minus, self.section_props.z22_plus,
            self.section_props.z22_minus
        )

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
        :return: St. Venant torsion constant
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
        :return: Centroidal axis shear centre (elasticity approach) *(x_se, y_se)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            (x_se, y_se) = section.get_sc()
        """

        if self.section_props.x_se is None:
            return (None, None)

        # add centroid location to move section back to original location
        x_se = self.section_props.x_se + self.section_props.cx
        y_se = self.section_props.y_se + self.section_props.cy

        return (x_se, y_se)

    def get_sc_p(self):
        """
        :return: Principal axis shear centre (elasticity approach) *(x11_se, y22_se)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            (x11_se, y22_se) = section.get_sc_p()
        """

        if self.section_props.x11_se is None:
            return (None, None)

        x11_se = self.section_props.x11_se
        y22_se = self.section_props.y22_se

        return (x11_se, y22_se)

    def get_sc_t(self):
        """
        :return: Centroidal axis shear centre (Trefftz's approach) *(x_st, y_st)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            (x_st, y_st) = section.get_sc_t()
        """

        if self.section_props.x_st is None:
            return (None, None)

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
        :return: Shear area for loading about the centroidal axis *(A_sx, A_sy)*
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
        :return: Shear area for loading about the principal bending axis *(A_s11, A_s22)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            (A_s11, A_s22) = section.get_As_p()
        """

        return (self.section_props.A_s11, self.section_props.A_s22)

    def get_beta(self):
        """
        :return: Monosymmetry constant for bending about both global axes *(beta_x_plus,
            beta_x_minus, beta_y_plus, beta_y_minus)*. The *plus* value relates to the top flange
            in compression and the *minus* value relates to the bottom flange in compression.
        :rtype: tuple(float, float, float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            (beta_x_plus, beta_x_minus, beta_y_plus, beta_y_minus) = section.get_beta()
        """

        return (
            self.section_props.beta_x_plus, self.section_props.beta_x_minus,
            self.section_props.beta_y_plus, self.section_props.beta_y_minus
        )

    def get_beta_p(self):
        """
        :return: Monosymmetry constant for bending about both principal axes *(beta_11_plus,
            beta_11_minus, beta_22_plus, beta_22_minus)*. The *plus* value relates to the top
            flange in compression and the *minus* value relates to the bottom flange in
            compression.
        :rtype: tuple(float, float, float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            (beta_11_plus, beta_11_minus, beta_22_plus, beta_22_minus) = section.get_beta_p()
        """

        return (
            self.section_props.beta_11_plus, self.section_props.beta_11_minus,
            self.section_props.beta_22_plus, self.section_props.beta_22_minus
        )

    def get_pc(self):
        """
        :return: Centroidal axis plastic centroid *(x_pc, y_pc)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_plastic_properties()
            (x_pc, y_pc) = section.get_pc()
        """

        if self.section_props.x_pc is None:
            return (None, None)

        # add centroid location to move section back to original location
        x_pc = self.section_props.x_pc + self.section_props.cx
        y_pc = self.section_props.y_pc + self.section_props.cy

        return (x_pc, y_pc)

    def get_pc_p(self):
        """
        :return: Principal bending axis plastic centroid *(x11_pc, y22_pc)*
        :rtype: tuple(float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_plastic_properties()
            (x11_pc, y22_pc) = section.get_pc_p()
        """

        if self.section_props.x11_pc is None:
            return (None, None)

        # determine the position of the plastic centroid in the global axis
        (x_pc, y_pc) = fea.global_coordinate(
            self.section_props.phi, self.section_props.x11_pc, self.section_props.y22_pc
        )

        # add centroid location to move section back to original location
        return (x_pc + self.section_props.cx, y_pc + self.section_props.cy)

    def get_s(self):
        """
        :return: Plastic section moduli about the centroidal axis *(sxx, syy)*
        :rtype: tuple(float, float)

        If material properties have been specified, returns the plastic moment :math:`M_p = f_y S`.

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_plastic_properties()
            (sxx, syy) = section.get_s()
        """

        return (self.section_props.sxx, self.section_props.syy)

    def get_sp(self):
        """
        :return: Plastic section moduli about the principal bending axis *(s11, s22)*
        :rtype: tuple(float, float)

        If material properties have been specified, returns the plastic moment
        :math:`M_p = f_y S`.

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_plastic_properties()
            (s11, s22) = section.get_sp()
        """

        return (self.section_props.s11, self.section_props.s22)

    def get_sf(self):
        """
        :return: Centroidal axis shape factors with respect to the top and bottom fibres
            *(sf_xx_plus, sf_xx_minus, sf_yy_plus, sf_yy_minus)*
        :rtype: tuple(float, float, float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_plastic_properties()
            (sf_xx_plus, sf_xx_minus, sf_yy_plus, sf_yy_minus) = section.get_sf()
        """

        return (
            self.section_props.sf_xx_plus, self.section_props.sf_xx_minus,
            self.section_props.sf_yy_plus, self.section_props.sf_yy_minus
        )

    def get_sf_p(self):
        """
        :return: Principal bending axis shape factors with respect to the top and bottom fibres
            *(sf_11_plus, sf_11_minus, sf_22_plus, sf_22_minus)*
        :rtype: tuple(float, float, float, float)

        ::

            section = CrossSection(geometry, mesh)
            section.calculate_geometric_properties()
            section.calculate_plastic_properties()
            (sf_11_plus, sf_11_minus, sf_22_plus, sf_22_minus) = section.get_sf_p()
        """

        return (
            self.section_props.sf_11_plus, self.section_props.sf_11_minus,
            self.section_props.sf_22_plus, self.section_props.sf_22_minus
        )


class PlasticSection:
    """Class for the plastic analysis of cross-sections.

    Stores the finite element geometry and material information and provides methods to compute the
    plastic section properties.

    :param geometry: Cross-section geometry object
    :type geometry: :class:`~sectionproperties.pre.sections.Geometry`
    :param materials: A list of material properties corresponding to various regions in the
        geometry and mesh.
    :type materials: list[:class:`~sectionproperties.pre.pre.Material`]
    :param bool debug: If set to True, the geometry is plotted each time a new mesh is generated by
        the plastic centroid algorithm.

    :cvar geometry: Deep copy of the cross-section geometry object provided to the constructor
    :vartype geometry: :class:`~sectionproperties.pre.sections.Geometry`
    :cvar materials: A list of material properties corresponding to various regions in the geometry
        and mesh.
    :vartype materials: list[:class:`~sectionproperties.pre.pre.Material`]
    :cvar bool debug: If set to True, the geometry is plotted each time a new mesh is generated by
        the plastic centroid algorithm.
    :cvar mesh: Mesh object returned by meshpy
    :vartype mesh: :class:`meshpy.triangle.MeshInfo`
    :cvar mesh_nodes: Array of node coordinates from the mesh
    :vartype mesh_nodes: :class:`numpy.ndarray`
    :cvar mesh_elements: Array of connectivities from the mesh
    :vartype mesh_elements: :class:`numpy.ndarray`
    :cvar elements: List of finite element objects describing the cross-section mesh
    :vartype elements: list[:class:`~sectionproperties.analysis.fea.Tri6`]
    :cvar float f_top: Current force in the top region
    :cvar c_top: Centroid of the force in the top region *(c_top_x, c_top_y)*
    :type c_top: list[float, float]
    :cvar c_bot: Centroid of the force in the bottom region *(c_bot_x, c_bot_y)*
    :type c_bot: list[float, float]
    """

    def __init__(self, geometry, materials, debug):
        """Inits the PlasticSection class."""

        # make a deepcopy of the geometry & materials so that we can modify it
        self.geometry = copy.deepcopy(geometry)
        self.materials = copy.deepcopy(materials)
        self.debug = debug

        # initialize variables to be defined later within calculate_plastic_force
        self.c_top = [0.0, 0.0]
        self.c_bot = [0.0, 0.0]
        self.f_top = 0.0
        self.f_bot = 0.0

        if self.materials is not None:
            # create dummy control point at the start of the list
            (x_min, _, y_min, _) = geometry.calculate_extents()
            self.geometry.control_points.insert(0, [x_min - 1, y_min - 1])

            # create matching dummy material
            self.materials.insert(0, pre.Material('default', 1, 0, 1))

        # create simple mesh of the geometry
        mesh = self.create_plastic_mesh()

        # get the elements of the mesh
        (_, _, elements) = self.get_elements(mesh)

        # calculate centroid of the mesh
        (cx, cy) = self.calculate_centroid(elements)

        # shift geometry such that the origin is at the centroid
        self.geometry.shift = [-cx, -cy]
        self.geometry.shift_section()

        # remesh the geometry and store the mesh
        self.mesh = self.create_plastic_mesh()

        # store the nodes, elements and list of elements in the mesh
        (self.mesh_nodes, self.mesh_elements, self.elements) = self.get_elements(self.mesh)

    def get_elements(self, mesh):
        """Extracts finite elements from the provided mesh and returns Tri6 finite elements with
        their associated material properties.

        :param mesh: Mesh object returned by meshpy
        :type mesh: :class:`meshpy.triangle.MeshInfo`
        :return: A tuple containing an array of the nodes locations, element indices and a list of
            the finite elements.
        :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`,
            list[:class:`~sectionproperties.analysis.fea.Tri6`])
        """

        # extract mesh data
        nodes = np.array(mesh.points, dtype=np.dtype(float))
        elements = np.array(mesh.elements, dtype=np.dtype(int))
        attributes = np.array(mesh.element_attributes, dtype=np.dtype(int))

        # swap mid-node order to retain node ordering consistency
        elements[:, [3, 4, 5]] = elements[:, [5, 3, 4]]

        # initialise list of Tri6 elements
        element_list = []

        # build the element list one element at a time
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
            coords = np.array([[x1, x2, x3, x4, x5, x6], [y1, y2, y3, y4, y5, y6]])

            # if materials are specified, get the material
            if self.materials is not None:
                # get attribute index of current element
                att_el = attributes[i]

                # if the current element is assigned the default attribute
                if att_el == 0:
                    # determine point within current element (centroid)
                    pt = [(x1 + x2 + x3) / 3, (y1 + y2 + y3) / 3]

                    # search within original elements - find coinciding element
                    for el in self.elements:
                        # if the point lies within the current element
                        if el.point_within_element(pt):
                            material = el.material
                            break
                else:
                    # fetch the material
                    material = self.materials[att_el]
            # if there are no materials specified, use a default material
            else:
                material = pre.Material('default', 1, 0, 1)

            # add tri6 elements to the element list
            element_list.append(fea.Tri6(i, coords, node_ids, material))

        return (nodes, elements, element_list)

    @staticmethod
    def calculate_centroid(elements):
        """Calculates the elastic centroid from a list of finite elements.

        :param elements: A list of Tri6 finite elements.
        :type elements: list[:class:`~sectionproperties.analysis.fea.Tri6`]
        :return: A tuple containing the x and y location of the elastic centroid.
        :rtype: tuple(float, float)
        """

        ea = 0
        qx = 0
        qy = 0

        # loop through all the elements
        for el in elements:
            (area, qx_el, qy_el, _, _, _, e, _) = el.geometric_properties()
            ea += area * e
            qx += qx_el * e
            qy += qy_el * e

        return (qy / ea, qx / ea)

    def calculate_plastic_properties(self, cross_section, verbose):
        """Calculates the location of the plastic centroid with respect to the centroidal and
        principal bending axes, the plastic section moduli and shape factors and stores the results
        to the supplied :class:`~sectionproperties.analysis.cross_section.CrossSection` object.

        :param cross_section: Cross section object that uses the same geometry and materials
            specified in the class constructor
        :type cross_section: :class:`~sectionproperties.analysis.cross_section.CrossSection`
        :param bool verbose: If set to True, the number of iterations required for each plastic
            axis is printed to the terminal.
        """

        # 1) Calculate plastic properties for centroidal axis
        # calculate distances to the extreme fibres
        fibres = self.calculate_extreme_fibres(0)

        # 1a) Calculate x-axis plastic centroid
        (y_pc, r, f, c_top, c_bot) = self.pc_algorithm(np.array([1, 0]), fibres[2:], 1, verbose)

        self.check_convergence(r, 'x-axis')
        cross_section.section_props.y_pc = y_pc
        cross_section.section_props.sxx = f * abs(c_top[1] - c_bot[1])

        if verbose:
            self.print_verbose(y_pc, r, 'x-axis')

        # 1b) Calculate y-axis plastic centroid
        (x_pc, r, f, c_top, c_bot) = self.pc_algorithm(np.array([0, 1]), fibres[0:2], 2, verbose)

        self.check_convergence(r, 'y-axis')
        cross_section.section_props.x_pc = x_pc
        cross_section.section_props.syy = f * abs(c_top[0] - c_bot[0])

        if verbose:
            self.print_verbose(x_pc, r, 'y-axis')

        # 2) Calculate plastic properties for principal axis
        # convert principal axis angle to radians
        angle = cross_section.section_props.phi * np.pi / 180

        # unit vectors in the axis directions
        ux = np.array([np.cos(angle), np.sin(angle)])
        uy = np.array([-np.sin(angle), np.cos(angle)])

        # calculate distances to the extreme fibres in the principal axis
        fibres = self.calculate_extreme_fibres(cross_section.section_props.phi)

        # 2a) Calculate 11-axis plastic centroid
        (y22_pc, r, f, c_top, c_bot) = self.pc_algorithm(ux, fibres[2:], 1, verbose)

        # calculate the centroids in the principal coordinate system
        c_top_p = fea.principal_coordinate(cross_section.section_props.phi, c_top[0], c_top[1])
        c_bot_p = fea.principal_coordinate(cross_section.section_props.phi, c_bot[0], c_bot[1])

        self.check_convergence(r, '11-axis')
        cross_section.section_props.y22_pc = y22_pc
        cross_section.section_props.s11 = f * abs(c_top_p[1] - c_bot_p[1])

        if verbose:
            self.print_verbose(y22_pc, r, '11-axis')

        # 2b) Calculate 22-axis plastic centroid
        (x11_pc, r, f, c_top, c_bot) = self.pc_algorithm(uy, fibres[0:2], 2, verbose)

        # calculate the centroids in the principal coordinate system
        c_top_p = fea.principal_coordinate(cross_section.section_props.phi, c_top[0], c_top[1])
        c_bot_p = fea.principal_coordinate(cross_section.section_props.phi, c_bot[0], c_bot[1])

        self.check_convergence(r, '22-axis')
        cross_section.section_props.x11_pc = x11_pc
        cross_section.section_props.s22 = f * abs(c_top_p[0] - c_bot_p[0])

        if verbose:
            self.print_verbose(x11_pc, r, '22-axis')

        # if there are no materials specified, calculate shape factors
        if cross_section.materials is None:
            cross_section.section_props.sf_xx_plus = (
                cross_section.section_props.sxx / cross_section.section_props.zxx_plus
            )
            cross_section.section_props.sf_xx_minus = (
                cross_section.section_props.sxx / cross_section.section_props.zxx_minus
            )
            cross_section.section_props.sf_yy_plus = (
                cross_section.section_props.syy / cross_section.section_props.zyy_plus
            )
            cross_section.section_props.sf_yy_minus = (
                cross_section.section_props.syy / cross_section.section_props.zyy_minus
            )

            cross_section.section_props.sf_11_plus = (
                cross_section.section_props.s11 / cross_section.section_props.z11_plus
            )
            cross_section.section_props.sf_11_minus = (
                cross_section.section_props.s11 / cross_section.section_props.z11_minus
            )
            cross_section.section_props.sf_22_plus = (
                cross_section.section_props.s22 / cross_section.section_props.z22_plus
            )
            cross_section.section_props.sf_22_minus = (
                cross_section.section_props.s22 / cross_section.section_props.z22_minus
            )

    @staticmethod
    def check_convergence(root_result, axis):
        """Checks that the function solver converged and if not, raises a helpful error.

        :param root_result: Result object from the root finder
        :type root_result: :class:`scipy.optimize.RootResults`
        :param string axis: Axis being considered by the function solver
        :raises RuntimeError: If the function solver did not converge
        """

        if not root_result.converged:
            msg = "Plastic centroid calculation about the {0}".format(axis)
            msg += " failed. Contact robbie.vanleeuwen@gmail.com with your"
            msg += " analysis parameters. Termination flag: {0}".format(root_result.flag)

            raise RuntimeError(msg)

    @staticmethod
    def print_verbose(d, root_result, axis):
        """Prints information related to the function solver convergence to the terminal.

        :param float d: Location of the plastic centroid axis
        :param root_result: Result object from the root finder
        :type root_result: :class:`scipy.optimize.RootResults`
        :param string axis: Axis being considered by the function solver
        """

        msg = "---{0} plastic centroid calculation converged at ".format(axis)
        msg += "{0:.5e} in {1} iterations.".format(d, root_result.iterations)
        print(msg)

    def calculate_extreme_fibres(self, angle):
        """Calculates the locations of the extreme fibres along and perpendicular to the axis
        specified by 'angle' using the elements stored in `self.elements`.

        :param float angle: Angle (in radians) along which to calculate the extreme fibre locations
        :return: The location of the extreme fibres parallel (u) and perpendicular (v) to the axis
            *(u_min, u_max, v_min, v_max)*
        :rtype: tuple(float, float, float, float)
        """

        # loop through all nodes in the mesh
        for (i, pt) in enumerate(self.mesh_nodes):
            # determine the coordinate of the point wrt the axis
            (u, v) = fea.principal_coordinate(angle, pt[0], pt[1])

            # initialise min, max variables
            if i == 0:
                u_min = u
                u_max = u
                v_min = v
                v_max = v

            # update the mins and maxes where necessary
            u_min = min(u_min, u)
            u_max = max(u_max, u)
            v_min = min(v_min, v)
            v_max = max(v_max, v)

        return (u_min, u_max, v_min, v_max)

    def evaluate_force_eq(self, d, u, u_p, verbose):
        """Given a distance *d* from the centroid to an axis (defined by unit vector *u*), creates
        a mesh including the new and axis and calculates the force equilibrium. The resultant
        force, as a ratio of the total force, is returned.

        :param float d: Distance from the centroid to current axis
        :param u: Unit vector defining the direction of the axis
        :type u: :class:`numpy.ndarray`
        :param u_p: Unit vector perpendicular to the direction of the axis
        :type u_p: :class:`numpy.ndarray`
        :param bool verbose: If set to True, the number of iterations required for each plastic
            axis is printed to the terminal.
        :return: The force equilibrium norm
        :rtype: float
        """

        p = np.array([d * u_p[0], d * u_p[1]])

        # create a mesh with the axis included
        mesh = self.create_plastic_mesh([p, u])
        (nodes, elements, element_list) = self.get_elements(mesh)

        if self.debug:
            self.plot_mesh(nodes, elements, element_list, self.materials)

        # calculate force equilibrium
        (f_top, f_bot) = self.calculate_plastic_force(element_list, u, p)

        # calculate the force norm
        f_norm = (f_top - f_bot) / (f_top + f_bot)

        # print verbose results
        if verbose:
            print("d = {0}; f_norm = {1}".format(d, f_norm))

        # return the force norm
        return f_norm

    def pc_algorithm(self, u, dlim, axis, verbose):
        """An algorithm used for solving for the location of the plastic centroid. The algorithm
        searches for the location of the axis, defined by unit vector *u* and within the section
        depth, that satisfies force equilibrium.

        :param u: Unit vector defining the direction of the axis
        :type u: :class:`numpy.ndarray`
        :param dlim: List [dmax, dmin] containing the distances from the centroid to the extreme
            fibres perpendicular to the axis
        :type dlim: list[float, float]
        :param int axis: The current axis direction: 1 (e.g. x or 11) or 2 (e.g. y or 22)
        :param bool verbose: If set to True, the number of iterations required for each plastic
            axis is printed to the terminal.
        :return: The distance to the plastic centroid axis *d*, the result object *r*, the force in
            the top of the section *f_top* and the location of the centroids of the top and bottom
            areas *c_top* and *c_bottom*
        :rtype: tuple(float, :class:`scipy.optimize.RootResults`, float, list[float, float],
            list[float, float])
        """

        # calculate vector perpendicular to u
        if axis == 1:
            u_p = np.array([-u[1], u[0]])
        else:
            u_p = np.array([u[1], -u[0]])

        a = dlim[0]
        b = dlim[1]

        (d, r) = brentq(
            self.evaluate_force_eq, a, b, args=(u, u_p, verbose), full_output=True, disp=False,
            xtol=1e-6, rtol=1e-6
        )

        return (d, r, self.f_top, self.c_top, self.c_bot)

    def calculate_plastic_force(self, elements, u, p):
        """Sums the forces above and below the axis defined by unit vector *u* and point *p*. Also
        returns the force centroid of the forces above and below the axis.

        :param elements: A list of Tri6 finite elements.
        :type elements: list[:class:`~sectionproperties.analysis.fea.Tri6`]
        :param u: Unit vector defining the direction of the axis
        :type u: :class:`numpy.ndarray`
        :param p: Point on the axis
        :type p: :class:`numpy.ndarray`
        :return: Force in the top and bottom areas *(f_top, f_bot)*
        :rtype: tuple(float, float)
        """

        # initialise variables
        (f_top, f_bot) = (0, 0)
        (ea_top, ea_bot) = (0, 0)
        (qx_top, qx_bot) = (0, 0)
        (qy_top, qy_bot) = (0, 0)

        # loop through all elements in the mesh
        for el in elements:
            # calculate element force and area properties
            (f_el, ea_el, qx_el, qy_el, is_above) = el.plastic_properties(u, p)

            # assign force and area properties to the top or bottom segments
            if is_above:
                f_top += f_el
                ea_top += ea_el
                qx_top += qx_el
                qy_top += qy_el
            else:
                f_bot += f_el
                ea_bot += ea_el
                qx_bot += qx_el
                qy_bot += qy_el

        # if there are no elements in the top/bottom prevent division by zero N.B. the algorithm
        # will never converge at this point, this is purely done to ensure a 100% search range
        if ea_top == 0:
            ea_top = 1
        if ea_bot == 0:
            ea_bot = 1

        # calculate the centroid of the top and bottom segments and save
        self.c_top = [qy_top / ea_top, qx_top / ea_top]
        self.c_bot = [qy_bot / ea_bot, qx_bot / ea_bot]
        self.f_top = f_top
        self.f_bot = f_bot

        return (f_top, f_bot)

    def create_plastic_mesh(self, new_line=None):
        """Generates a triangular mesh of a deep copy of the geometry stored in `self.geometry`.
        Optionally, a line can be added to the copied geometry, which is defined by a point *p* and
        a unit vector *u*.

        :param new_line: A point p and a unit vector u defining a line to add to the mesh
            (new_line: p -> p + u) [*p*, *u*]
        :type new_line: list[:class:`numpy.ndarray`, :class:`numpy.ndarray`]
        :param mesh: Mesh object returned by meshpy
        :type mesh: :class:`meshpy.triangle.MeshInfo`
        """

        # start with the initial geometry
        geom = copy.deepcopy(self.geometry)

        # add line at new_line
        if new_line is not None:
            self.add_line(geom, new_line)

            # fast clean the geometry after adding the line
            geom.zip_points()
            geom.remove_zero_length_facets()
            geom.remove_unused_points()

        if self.debug:
            if new_line is not None:
                geom.plot_geometry(labels=True)

        # build mesh object
        mesh = triangle.MeshInfo()  # create mesh info object
        mesh.set_points(geom.points)  # set points
        mesh.set_facets(geom.facets)  # set facets
        mesh.set_holes(geom.holes)  # set holes

        # set regions
        mesh.regions.resize(len(geom.control_points))
        region_id = 0  # initialise region ID variable

        for (i, cp) in enumerate(geom.control_points):
            mesh.regions[i] = [cp[0], cp[1], region_id, 1]
            region_id += 1

        mesh = triangle.build(mesh, mesh_order=2, quality_meshing=False, attributes=True)

        return mesh

    def add_line(self, geometry, line):
        """Adds a line a geometry object. Finds the intersection points of the line with the
        current facets and splits the existing facets to accommodate the new line.

        :param geometry: Cross-section geometry object used to generate the mesh
        :type geometry: :class:`~sectionproperties.pre.sections.Geometry`
        :param line: A point p and a unit vector u defining a line to add to the mesh
            (line: p -> p + u)
        :type line: list[:class:`numpy.ndarray`, :class:`numpy.ndarray`]
        """

        # initialise intersection points and facet index list
        int_pts = []
        fct_idx = []

        # get current number of points in the geometry object
        num_pts = len(geometry.points)

        # line: p -> p + r
        p = line[0]
        r = line[1]

        # loop through all the facets in the geometry to find intersection pts
        for (idx, fct) in enumerate(geometry.facets):
            # facet: q -> q + s
            q = np.array(geometry.points[fct[0]])
            s = geometry.points[fct[1]] - q

            # calculate intersection point between p -> p + r and q -> q + s N.B. make line
            # p -> p + r infinitely long to find all intersects if the lines are not parallel
            if np.cross(r, s) != 0:
                # calculate t and u
                t = np.cross(q - p, s) / np.cross(r, s)
                u = np.cross(p - q, r) / np.cross(s, r)

                new_pt = p + t * r

                # if the line lies within q -> q + s and the point hasn't already been added
                # (ignore t as it is infinitely long)
                if (0 <= u <= 1 and list(new_pt) not in [list(item) for item in int_pts]):
                    int_pts.append(new_pt)
                    fct_idx.append(idx)

        # if less than 2 intersection points are found, we are at the edge of the section,
        # therefore no line to add
        if len(int_pts) < 2:
            return

        # sort intersection points and facet list first by x, then by y
        int_pts = np.array(int_pts)
        idx_sort = np.lexsort((int_pts[:, 0], int_pts[:, 1]))
        int_pts = int_pts[idx_sort]
        fct_idx = list(np.array(fct_idx)[idx_sort])

        # add points to the geometry object
        for pt in int_pts:
            geometry.points.append([pt[0], pt[1]])

        # add new facets by looping from the second facet index to the end
        for (i, idx) in enumerate(fct_idx[1:]):
            # get mid-point of proposed new facet
            mid_pt = 0.5 * (int_pts[i] + int_pts[i + 1])

            # check to see if the mid-point is not in a hole
            # add the facet
            if self.point_within_element(mid_pt):
                geometry.facets.append([num_pts + i, num_pts + i + 1])

            # rebuild the intersected facet
            self.rebuild_parent_facet(geometry, idx, num_pts + i + 1)

            # rebuild the first facet the looped skipped
            if i == 0:
                self.rebuild_parent_facet(geometry, fct_idx[0], num_pts + i)

        # sort list of facet indices (to be removed) in reverse order so as not to compromise the
        # indices during deletion
        idx_to_remove = sorted(fct_idx, reverse=True)

        for idx in idx_to_remove:
            geometry.facets.pop(idx)

    @staticmethod
    def rebuild_parent_facet(geometry, fct_idx, pt_idx):
        """Splits and rebuilds a facet at a given point.

        :param geometry: Cross-section geometry object used to generate the mesh
        :type geometry: :class:`~sectionproperties.pre.sections.Geometry`
        :param int fct_idx: Index of the facet to be split
        :param int pt_idx: Index of the point to insert into the facet
        """

        # get current facet
        fct = geometry.facets[fct_idx]

        # rebuild facet
        geometry.facets.append([fct[0], pt_idx])
        geometry.facets.append([pt_idx, fct[1]])

    def point_within_element(self, pt):
        """Determines whether a point lies within an element in the mesh stored in
        `self.mesh_elements`.

        :param pt: Point to check
        :type pt: :class:`numpy.ndarray`
        :return: Whether the point lies within an element
        :rtype: bool
        """

        px = pt[0]
        py = pt[1]

        # loop through elements in the mesh
        for el in self.mesh_elements:
            # get coordinates of corner points
            x1 = self.mesh_nodes[el[0]][0]
            y1 = self.mesh_nodes[el[0]][1]
            x2 = self.mesh_nodes[el[1]][0]
            y2 = self.mesh_nodes[el[1]][1]
            x3 = self.mesh_nodes[el[2]][0]
            y3 = self.mesh_nodes[el[2]][1]

            # compute variables alpha, beta and gamma
            alpha = (
                ((y2 - y3) * (px - x3) + (x3 - x2) * (py - y3))
                / ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3))
            )
            beta = (
                ((y3 - y1) * (px - x3) + (x1 - x3) * (py - y3))
                / ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3))
            )
            gamma = 1.0 - alpha - beta

            # if the point lies within an element
            if alpha >= 0 and beta >= 0 and gamma >= 0:
                return True

        return False

    @staticmethod
    def plot_mesh(nodes, elements, element_list, materials):
        """Watered down implementation of the CrossSection method to plot the finite element mesh,
        showing material properties."""

        (_, ax) = plt.subplots()
        post.setup_plot(True)

        # plot the mesh
        ax.triplot(nodes[:, 0], nodes[:, 1], elements[:, 0:3], lw=0.5,
                   color='black')

        color_array = []
        legend_list = []

        if materials is not None:
            # create an array of finite element colours
            for el in element_list:
                color_array.append(el.material.color)

            # create a list of unique material legend entries
            for (i, mat) in enumerate(materials):
                # if the material has not be entered yet
                if i == 0 or mat not in materials[0:i]:
                    # add the material colour and name to the legend list
                    legend_list.append(mpatches.Patch(color=mat.color, label=mat.name))

            cmap = ListedColormap(color_array)  # custom colormap
            c = np.arange(len(color_array))  # indices of elements

            # plot the mesh colours
            ax.tripcolor(nodes[:, 0], nodes[:, 1], elements[:, 0:3], c, cmap=cmap)

            # display the legend
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), handles=legend_list)

        # finish the plot
        post.finish_plot(ax, True, title='Finite Element Mesh')


class StressPost:
    """Class for post-processing finite element stress results.

    A StressPost object is created when a stress analysis is carried out and is returned as an
    object to allow post-processing of the results. The StressPost object creates a deep copy of
    the MaterialGroups within the cross-section to allow the calculation of stresses for each
    material. Methods for post-processing the calculated stresses are provided.

    :param cross_section: Cross section object for stress calculation
    :type cross_section: :class:`~sectionproperties.analysis.cross_section.CrossSection`

    :cvar cross_section: Cross section object for stress calculation
    :vartype cross_section: :class:`~sectionproperties.analysis.cross_section.CrossSection`
    :cvar material_groups: A deep copy of the `cross_section` material groups to allow a new stress
        analysis
    :vartype material_groups: list[:class:`~sectionproperties.pre.pre.MaterialGroup`]
    """

    def __init__(self, cross_section):
        """Inits the StressPost class."""

        self.cross_section = cross_section

        # make a deep copy of the material groups to the StressPost object such that stress results
        # can be saved to a new material group
        self.material_groups = copy.deepcopy(cross_section.material_groups)

    def plot_stress_contour(self, sigs, title, pause):
        """Plots filled stress contours over the finite element mesh.

        :param sigs: List of nodal stress values for each material
        :type sigs: list[:class:`numpy.ndarray`]
        :param string title: Plot title
        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)
        """

        # create plot and setup the plot
        (fig, ax) = plt.subplots()
        post.setup_plot(pause)

        # plot the finite element mesh
        self.cross_section.plot_mesh(ax, pause, alpha=0.5)

        # set up the colormap
        cmap = cm.get_cmap(name='jet')

        # create triangulation
        triang = tri.Triangulation(
            self.cross_section.mesh_nodes[:, 0], self.cross_section.mesh_nodes[:, 1],
            self.cross_section.mesh_elements[:, 0:3]
        )

        # determine minimum and maximum stress values for the contour list
        sig_min = min([min(x) for x in sigs])
        sig_max = max([max(x) for x in sigs])
        v = np.linspace(sig_min, sig_max, 15, endpoint=True)

        if np.isclose(v[0], v[-1], atol=1e-12):
            v = 15
            ticks = None
        else:
            ticks = v

        # plot the filled contour, looping through the materials
        for (i, sig) in enumerate(sigs):
            # create and set the mask for the current material
            mask_array = np.ones(len(self.cross_section.elements), dtype=bool)
            mask_array[self.material_groups[i].el_ids] = False
            triang.set_mask(mask_array)

            # plot the filled contour
            trictr = ax.tricontourf(triang, sig, v, cmap=cmap)

        # display the colourbar
        fig.colorbar(trictr, label='Stress', format='%.4e', ticks=ticks)

        # TODO: display stress values in the toolbar (format_coord)

        # finish the plot
        post.finish_plot(ax, pause, title)

        return (fig, ax)

    def plot_stress_vector(self, sigxs, sigys, title, pause):
        """Plots stress vectors over the finite element mesh.

        :param sigxs: List of x-components of the nodal stress values for each material
        :type sigxs: list[:class:`numpy.ndarray`]
        :param sigys: List of y-components of the nodal stress values for each material
        :type sigys: list[:class:`numpy.ndarray`]
        :param string title: Plot title
        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)
        """

        # create plot and setup the plot
        (fig, ax) = plt.subplots()
        post.setup_plot(pause)

        # plot the finite element mesh
        self.cross_section.plot_mesh(ax, pause, alpha=0.5)

        # set up the colormap
        cmap = cm.get_cmap(name='jet')

        # initialise quiver plot list max scale
        quiv_list = []
        max_scale = 0

        # plot the vectors
        for (i, sigx) in enumerate(sigxs):
            sigy = sigys[i]

            # scale the colour with respect to the magnitude of the vector
            c = np.hypot(sigx, sigy)

            quiv = ax.quiver(
                self.cross_section.mesh_nodes[:, 0], self.cross_section.mesh_nodes[:, 1], sigx,
                sigy, c, cmap=cmap
            )

            # get the scale and store the max value
            quiv._init()
            max_scale = max(max_scale, quiv.scale)
            quiv_list.append(quiv)

            # update the colormap values
            if i == 0:
                c_min = min(c)
                c_max = max(c)
            else:
                c_min = min(c_min, min(c))
                c_max = max(c_max, max(c))

        # apply the scale
        for quiv_plot in quiv_list:
            quiv_plot.scale = max_scale

        # apply the colourbar
        v1 = np.linspace(c_min, c_max, 15, endpoint=True)
        fig.colorbar(quiv, label='Stress', ticks=v1, format='%.4e')

        # finish the plot
        post.finish_plot(ax, pause, title=title)

        return (fig, ax)

    def get_stress(self):
        """Returns the stresses within each material belonging to the current
        :class:`~sectionproperties.analysis.cross_section.StressPost` object.

        :return: A list of dictionaries containing the cross-section stresses for each material.
        :rtype: list[dict]

        A dictionary is returned for each material in the cross-section, containing the following
        keys and values:

        * *'Material'*: Material name
        * *'sig_zz_n'*: Normal stress :math:`\sigma_{zz,N}` resulting from the axial load :math:`N`
        * *'sig_zz_mxx'*: Normal stress :math:`\sigma_{zz,Mxx}` resulting from the bending moment
          :math:`M_{xx}`
        * *'sig_zz_myy'*: Normal stress :math:`\sigma_{zz,Myy}` resulting from the bending moment
          :math:`M_{yy}`
        * *'sig_zz_m11'*: Normal stress :math:`\sigma_{zz,M11}` resulting from the bending moment
          :math:`M_{11}`
        * *'sig_zz_m22'*: Normal stress :math:`\sigma_{zz,M22}` resulting from the bending moment
          :math:`M_{22}`
        * *'sig_zz_m'*: Normal stress :math:`\sigma_{zz,\Sigma M}` resulting from all bending
          moments
        * *'sig_zx_mzz'*: *x*-component of the shear stress :math:`\sigma_{zx,Mzz}` resulting from
          the torsion moment
        * *'sig_zy_mzz'*: *y*-component of the shear stress :math:`\sigma_{zy,Mzz}` resulting from
          the torsion moment
        * *'sig_zxy_mzz'*: Resultant shear stress :math:`\sigma_{zxy,Mzz}` resulting from the
          torsion moment
        * *'sig_zx_vx'*: *x*-component of the shear stress :math:`\sigma_{zx,Vx}` resulting from
          the shear force :math:`V_{x}`
        * *'sig_zy_vx'*: *y*-component of the shear stress :math:`\sigma_{zy,Vx}` resulting from
          the shear force :math:`V_{x}`
        * *'sig_zxy_vx'*: Resultant shear stress :math:`\sigma_{zxy,Vx}` resulting from the shear
          force :math:`V_{x}`
        * *'sig_zx_vy'*: *x*-component of the shear stress :math:`\sigma_{zx,Vy}` resulting from
          the shear force :math:`V_{y}`
        * *'sig_zy_vy'*: *y*-component of the shear stress :math:`\sigma_{zy,Vy}` resulting from
          the shear force :math:`V_{y}`
        * *'sig_zxy_vy'*: Resultant shear stress :math:`\sigma_{zxy,Vy}` resulting from the shear
          force :math:`V_{y}`
        * *'sig_zx_v'*: *x*-component of the shear stress :math:`\sigma_{zx,\Sigma V}` resulting
          from all shear forces
        * *'sig_zy_v'*: *y*-component of the shear stress :math:`\sigma_{zy,\Sigma V}` resulting
          from all shear forces
        * *'sig_zxy_v'*: Resultant shear stress :math:`\sigma_{zxy,\Sigma V}` resulting from all
          shear forces
        * *'sig_zz'*: Combined normal stress :math:`\sigma_{zz}` resulting from all actions
        * *'sig_zx'*: *x*-component of the shear stress :math:`\sigma_{zx}` resulting from all
          actions
        * *'sig_zy'*: *y*-component of the shear stress :math:`\sigma_{zy}` resulting from all
          actions
        * *'sig_zxy'*: Resultant shear stress :math:`\sigma_{zxy}` resulting from all actions
        * *'sig_vm'*: von Mises stress :math:`\sigma_{vM}` resulting from all actions

        The following example returns the normal stress within a 150x90x12 UA section resulting
        from an axial force of 10 kN::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(N=10e3)

            stresses = stress_post.get_stress()
            print('Material: {0}'.format(stresses[0]['Material']))
            print('Axial Stresses: {0}'.format(stresses[0]['sig_zz_n']))

            $ Material: default
            $ Axial Stresses: [3.6402569 3.6402569 3.6402569 ... 3.6402569 3.6402569 3.6402569]
        """

        stress = []

        for group in self.material_groups:
            stress.append({
                'Material': group.material.name,
                'sig_zz_n': group.stress_result.sig_zz_n,
                'sig_zz_mxx': group.stress_result.sig_zz_mxx,
                'sig_zz_myy': group.stress_result.sig_zz_myy,
                'sig_zz_m11': group.stress_result.sig_zz_m11,
                'sig_zz_m22': group.stress_result.sig_zz_m22,
                'sig_zz_m': group.stress_result.sig_zz_m,
                'sig_zx_mzz': group.stress_result.sig_zx_mzz,
                'sig_zy_mzz': group.stress_result.sig_zy_mzz,
                'sig_zxy_mzz': group.stress_result.sig_zxy_mzz,
                'sig_zx_vx': group.stress_result.sig_zx_vx,
                'sig_zy_vx': group.stress_result.sig_zy_vx,
                'sig_zxy_vx': group.stress_result.sig_zxy_vx,
                'sig_zx_vy': group.stress_result.sig_zx_vy,
                'sig_zy_vy': group.stress_result.sig_zy_vy,
                'sig_zxy_vy': group.stress_result.sig_zxy_vy,
                'sig_zx_v': group.stress_result.sig_zx_v,
                'sig_zy_v': group.stress_result.sig_zy_v,
                'sig_zxy_v': group.stress_result.sig_zxy_v,
                'sig_zz': group.stress_result.sig_zz,
                'sig_zx': group.stress_result.sig_zx,
                'sig_zy': group.stress_result.sig_zy,
                'sig_zxy': group.stress_result.sig_zxy,
                'sig_vm': group.stress_result.sig_vm
            })

        return stress

    def plot_stress_n_zz(self, pause=True):
        """Produces a contour plot of the normal stress :math:`\sigma_{zz,N}` resulting from the
        axial load :math:`N`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots the normal stress within a 150x90x12 UA section resulting from
        an axial force of 10 kN::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(N=10e3)

            stress_post.plot_stress_n_zz()

        ..  figure:: ../images/stress/stress_n_zz.png
            :align: center
            :scale: 75 %

            Contour plot of the axial stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zz,N}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zz_n)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_stress_mxx_zz(self, pause=True):
        """Produces a contour plot of the normal stress :math:`\sigma_{zz,Mxx}` resulting from the
        bending moment :math:`M_{xx}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots the normal stress within a 150x90x12 UA section resulting from
        a bending moment about the x-axis of 5 kN.m::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mxx=5e6)

            stress_post.plot_stress_mxx_zz()

        ..  figure:: ../images/stress/stress_mxx_zz.png
            :align: center
            :scale: 75 %

            Contour plot of the bending stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zz,Mxx}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zz_mxx)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_stress_myy_zz(self, pause=True):
        """Produces a contour plot of the normal stress :math:`\sigma_{zz,Myy}` resulting from the
        bending moment :math:`M_{yy}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots the normal stress within a 150x90x12 UA section resulting from
        a bending moment about the y-axis of 2 kN.m::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Myy=2e6)

            stress_post.plot_stress_myy_zz()

        ..  figure:: ../images/stress/stress_myy_zz.png
            :align: center
            :scale: 75 %

            Contour plot of the bending stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zz,Myy}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zz_myy)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_stress_m11_zz(self, pause=True):
        """Produces a contour plot of the normal stress :math:`\sigma_{zz,M11}` resulting from the
        bending moment :math:`M_{11}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots the normal stress within a 150x90x12 UA section resulting from
        a bending moment about the 11-axis of 5 kN.m::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(M11=5e6)

            stress_post.plot_stress_m11_zz()

        ..  figure:: ../images/stress/stress_m11_zz.png
            :align: center
            :scale: 75 %

            Contour plot of the bending stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zz,M11}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zz_m11)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_stress_m22_zz(self, pause=True):
        """Produces a contour plot of the normal stress :math:`\sigma_{zz,M22}` resulting from the
        bending moment :math:`M_{22}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots the normal stress within a 150x90x12 UA section resulting from
        a bending moment about the 22-axis of 2 kN.m::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(M22=5e6)

            stress_post.plot_stress_m22_zz()

        ..  figure:: ../images/stress/stress_m22_zz.png
            :align: center
            :scale: 75 %

            Contour plot of the bending stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zz,M22}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zz_m22)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_stress_m_zz(self, pause=True):
        """Produces a contour plot of the normal stress :math:`\sigma_{zz,\Sigma M}` resulting from
        all bending moments :math:`M_{xx} + M_{yy} + M_{11} + M_{22}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots the normal stress within a 150x90x12 UA section resulting from
        a bending moment about the x-axis of 5 kN.m, a bending moment about the y-axis of 2 kN.m
        and a bending moment of 3 kN.m about the 11-axis::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mxx=5e6, Myy=2e6, M11=3e6)

            stress_post.plot_stress_m_zz()

        ..  figure:: ../images/stress/stress_m_zz.png
            :align: center
            :scale: 75 %

            Contour plot of the bending stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zz,\Sigma M}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zz_m)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_stress_mzz_zx(self, pause=True):
        """Produces a contour plot of the *x*-component of the shear stress :math:`\sigma_{zx,Mzz}`
        resulting from the torsion moment :math:`M_{zz}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots the x-component of the shear stress within a 150x90x12 UA
        section resulting from a torsion moment of 1 kN.m::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mzz=1e6)

            stress_post.plot_stress_mzz_zx()

        ..  figure:: ../images/stress/stress_mzz_zx.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zx,Mzz}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zx_mzz)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_stress_mzz_zy(self, pause=True):
        """Produces a contour plot of the *y*-component of the shear stress :math:`\sigma_{zy,Mzz}`
        resulting from the torsion moment :math:`M_{zz}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots the y-component of the shear stress within a 150x90x12 UA
        section resulting from a torsion moment of 1 kN.m::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mzz=1e6)

            stress_post.plot_stress_mzz_zy()

        ..  figure:: ../images/stress/stress_mzz_zy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zy,Mzz}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zy_mzz)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_stress_mzz_zxy(self, pause=True):
        """Produces a contour plot of the resultant shear stress :math:`\sigma_{zxy,Mzz}` resulting
        from the torsion moment :math:`M_{zz}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots a contour of the resultant shear stress within a 150x90x12 UA
        section resulting from a torsion moment of 1 kN.m::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mzz=1e6)

            stress_post.plot_stress_mzz_zxy()

        ..  figure:: ../images/stress/stress_mzz_zxy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zxy,Mzz}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zxy_mzz)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_vector_mzz_zxy(self, pause=True):
        """Produces a vector plot of the resultant shear stress :math:`\sigma_{zxy,Mzz}` resulting
        from the torsion moment :math:`M_{zz}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example generates a vector plot of the shear stress within a 150x90x12 UA
        section resulting from a torsion moment of 1 kN.m::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mzz=1e6)

            stress_post.plot_vector_mzz_zxy()

        ..  figure:: ../images/stress/vector_mzz_zxy.png
            :align: center
            :scale: 75 %

            Vector plot of the shear stress.
        """

        title = 'Stress Vector Plot - $\sigma_{zxy,Mzz}$'
        sigxs = []
        sigys = []

        for group in self.material_groups:
            sigxs.append(group.stress_result.sig_zx_mzz)
            sigys.append(group.stress_result.sig_zy_mzz)

        return self.plot_stress_vector(sigxs, sigys, title, pause)

    def plot_stress_vx_zx(self, pause=True):
        """Produces a contour plot of the *x*-component of the shear stress :math:`\sigma_{zx,Vx}`
        resulting from the shear force :math:`V_{x}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots the x-component of the shear stress within a 150x90x12 UA
        section resulting from a shear force in the x-direction of 15 kN::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vx=15e3)

            stress_post.plot_stress_vx_zx()

        ..  figure:: ../images/stress/stress_vx_zx.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zx,Vx}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zx_vx)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_stress_vx_zy(self, pause=True):
        """Produces a contour plot of the *y*-component of the shear stress :math:`\sigma_{zy,Vx}`
        resulting from the shear force :math:`V_{x}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots the y-component of the shear stress within a 150x90x12 UA
        section resulting from a shear force in the x-direction of 15 kN::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vx=15e3)

            stress_post.plot_stress_vx_zy()

        ..  figure:: ../images/stress/stress_vx_zy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zy,Vx}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zy_vx)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_stress_vx_zxy(self, pause=True):
        """Produces a contour plot of the resultant shear stress :math:`\sigma_{zxy,Vx}` resulting
        from the shear force :math:`V_{x}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots a contour of the resultant shear stress within a 150x90x12 UA
        section resulting from a shear force in the x-direction of 15 kN::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vx=15e3)

            stress_post.plot_stress_vx_zxy()

        ..  figure:: ../images/stress/stress_vx_zxy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zxy,Vx}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zxy_vx)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_vector_vx_zxy(self, pause=True):
        """Produces a vector plot of the resultant shear stress :math:`\sigma_{zxy,Vx}` resulting
        from the shear force :math:`V_{x}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example generates a vector plot of the shear stress within a 150x90x12 UA
        section resulting from a shear force in the x-direction of 15 kN::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vx=15e3)

            stress_post.plot_vector_vx_zxy()

        ..  figure:: ../images/stress/vector_vx_zxy.png
            :align: center
            :scale: 75 %

            Vector plot of the shear stress.
        """

        title = 'Stress Vector Plot - $\sigma_{zxy,Vx}$'
        sigxs = []
        sigys = []

        for group in self.material_groups:
            sigxs.append(group.stress_result.sig_zx_vx)
            sigys.append(group.stress_result.sig_zy_vx)

        return self.plot_stress_vector(sigxs, sigys, title, pause)

    def plot_stress_vy_zx(self, pause=True):
        """Produces a contour plot of the *x*-component of the shear stress :math:`\sigma_{zx,Vy}`
        resulting from the shear force :math:`V_{y}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots the x-component of the shear stress within a 150x90x12 UA
        section resulting from a shear force in the y-direction of 30 kN::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vy=30e3)

            stress_post.plot_stress_vy_zx()

        ..  figure:: ../images/stress/stress_vy_zx.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zx,Vy}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zx_vy)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_stress_vy_zy(self, pause=True):
        """Produces a contour plot of the *y*-component of the shear stress :math:`\sigma_{zy,Vy}`
        resulting from the shear force :math:`V_{y}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots the y-component of the shear stress within a 150x90x12 UA
        section resulting from a shear force in the y-direction of 30 kN::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vy=30e3)

            stress_post.plot_stress_vy_zy()

        ..  figure:: ../images/stress/stress_vy_zy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zy,Vy}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zy_vy)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_stress_vy_zxy(self, pause=True):
        """Produces a contour plot of the resultant shear stress :math:`\sigma_{zxy,Vy}` resulting
        from the shear force :math:`V_{y}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots a contour of the resultant shear stress within a 150x90x12 UA
        section resulting from a shear force in the y-direction of 30 kN::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vy=30e3)

            stress_post.plot_stress_vy_zxy()

        ..  figure:: ../images/stress/stress_vy_zxy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zxy,Vy}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zxy_vy)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_vector_vy_zxy(self, pause=True):
        """Produces a vector plot of the resultant shear stress :math:`\sigma_{zxy,Vy}` resulting
        from the shear force :math:`V_{y}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example generates a vector plot of the shear stress within a 150x90x12 UA
        section resulting from a shear force in the y-direction of 30 kN::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vy=30e3)

            stress_post.plot_vector_vy_zxy()

        ..  figure:: ../images/stress/vector_vy_zxy.png
            :align: center
            :scale: 75 %

            Vector plot of the shear stress.
        """

        title = 'Stress Vector Plot - $\sigma_{zxy,Vy}$'
        sigxs = []
        sigys = []

        for group in self.material_groups:
            sigxs.append(group.stress_result.sig_zx_vy)
            sigys.append(group.stress_result.sig_zy_vy)

        return self.plot_stress_vector(sigxs, sigys, title, pause)

    def plot_stress_v_zx(self, pause=True):
        """Produces a contour plot of the *x*-component of the shear stress
        :math:`\sigma_{zx,\Sigma V}` resulting from the sum of the applied shear forces
        :math:`V_{x} + V_{y}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots the x-component of the shear stress within a 150x90x12 UA
        section resulting from a shear force of 15 kN in the x-direction and 30 kN in the
        y-direction::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vx=15e3, Vy=30e3)

            stress_post.plot_stress_v_zx()

        ..  figure:: ../images/stress/stress_v_zx.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zx,\Sigma V}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zx_v)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_stress_v_zy(self, pause=True):
        """Produces a contour plot of the *y*-component of the shear stress
        :math:`\sigma_{zy,\Sigma V}` resulting from the sum of the applied shear forces
        :math:`V_{x} + V_{y}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots the y-component of the shear stress within a 150x90x12 UA
        section resulting from a shear force of 15 kN in the x-direction and 30 kN in the
        y-direction::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vx=15e3, Vy=30e3)

            stress_post.plot_stress_v_zy()

        ..  figure:: ../images/stress/stress_v_zy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zy,\Sigma V}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zy_v)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_stress_v_zxy(self, pause=True):
        """Produces a contour plot of the resultant shear stress
        :math:`\sigma_{zxy,\Sigma V}` resulting from the sum of the applied shear forces
        :math:`V_{x} + V_{y}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots a contour of the resultant shear stress within a 150x90x12 UA
        section resulting from a shear force of 15 kN in the x-direction and 30 kN in the
        y-direction::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vx=15e3, Vy=30e3)

            stress_post.plot_stress_v_zxy()

        ..  figure:: ../images/stress/stress_v_zxy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zxy,\Sigma V}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zxy_v)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_vector_v_zxy(self, pause=True):
        """Produces a vector plot of the resultant shear stress
        :math:`\sigma_{zxy,\Sigma V}` resulting from the sum of the  applied shear forces
        :math:`V_{x} + V_{y}`.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example generates a vector plot of the shear stress within a 150x90x12 UA
        section resulting from a shear force of 15 kN in the x-direction and 30 kN in the
        y-direction::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vx=15e3, Vy=30e3)

            stress_post.plot_vector_v_zxy()

        ..  figure:: ../images/stress/vector_v_zxy.png
            :align: center
            :scale: 75 %

            Vector plot of the shear stress.
        """

        title = 'Stress Vector Plot - $\sigma_{zxy,\Sigma V}$'
        sigxs = []
        sigys = []

        for group in self.material_groups:
            sigxs.append(group.stress_result.sig_zx_v)
            sigys.append(group.stress_result.sig_zy_v)

        return self.plot_stress_vector(sigxs, sigys, title, pause)

    def plot_stress_zz(self, pause=True):
        """Produces a contour plot of the combined normal stress :math:`\sigma_{zz}` resulting from
        all actions.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots the normal stress within a 150x90x12 UA section resulting from
        an axial force of 100 kN, a bending moment about the x-axis of 5 kN.m and a bending moment
        about the y-axis of 2 kN.m::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(N=100e3, Mxx=5e6, Myy=2e6)

            stress_post.plot_stress_zz()

        ..  figure:: ../images/stress/stress_zz.png
            :align: center
            :scale: 75 %

            Contour plot of the normal stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zz}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zz)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_stress_zx(self, pause=True):
        """Produces a contour plot of the *x*-component of the shear stress :math:`\sigma_{zx}`
        resulting from all actions.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots the x-component of the shear stress within a 150x90x12 UA
        section resulting from a torsion moment of 1 kN.m and a shear force of 30 kN in the
        y-direction::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mzz=1e6, Vy=30e3)

            stress_post.plot_stress_zx()

        ..  figure:: ../images/stress/stress_zx.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zx}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zx)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_stress_zy(self, pause=True):
        """Produces a contour plot of the *y*-component of the shear stress :math:`\sigma_{zy}`
        resulting from all actions.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots the y-component of the shear stress within a 150x90x12 UA
        section resulting from a torsion moment of 1 kN.m and a shear force of 30 kN in the
        y-direction::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mzz=1e6, Vy=30e3)

            stress_post.plot_stress_zy()

        ..  figure:: ../images/stress/stress_zy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zy}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zy)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_stress_zxy(self, pause=True):
        """Produces a contour plot of the resultant shear stress :math:`\sigma_{zxy}` resulting
        from all actions.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots a contour of the resultant shear stress within a 150x90x12 UA
        section resulting from a torsion moment of 1 kN.m and a shear force of 30 kN in the
        y-direction::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mzz=1e6, Vy=30e3)

            stress_post.plot_stress_zxy()

        ..  figure:: ../images/stress/stress_zxy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        title = 'Stress Contour Plot - $\sigma_{zxy}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zxy)

        return self.plot_stress_contour(sigs, title, pause)

    def plot_vector_zxy(self, pause=True):
        """Produces a vector plot of the resultant shear stress :math:`\sigma_{zxy}` resulting
        from all actions.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example generates a vector plot of the shear stress within a 150x90x12 UA
        section resulting from a torsion moment of 1 kN.m and a shear force of 30 kN in the
        y-direction::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mzz=1e6, Vy=30e3)

            stress_post.plot_vector_zxy()

        ..  figure:: ../images/stress/vector_zxy.png
            :align: center
            :scale: 75 %

            Vector plot of the shear stress.
        """

        title = 'Stress Vector Plot - $\sigma_{zxy}$'
        sigxs = []
        sigys = []

        for group in self.material_groups:
            sigxs.append(group.stress_result.sig_zx)
            sigys.append(group.stress_result.sig_zy)

        return self.plot_stress_vector(sigxs, sigys, title, pause)

    def plot_stress_vm(self, pause=True):
        """Produces a contour plot of the von Mises stress :math:`\sigma_{vM}` resulting from all
        actions.

        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example plots a contour of the von Mises stress within a 150x90x12 UA section
        resulting from the following actions:

        * :math:`N = 50` kN
        * :math:`M_{xx} = -5` kN.m
        * :math:`M_{22} = 2.5` kN.m
        * :math:`M_{zz} = 1.5` kN.m
        * :math:`V_{x} = 10` kN
        * :math:`V_{y} = 5` kN

        ::

            import sectionproperties.pre.sections as sections
            from sectionproperties.analysis.cross_section import CrossSection

            geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(
                N=50e3, Mxx=-5e6, M22=2.5e6, Mzz=0.5e6, Vx=10e3, Vy=5e3
            )

            stress_post.plot_stress_vm()

        ..  figure:: ../images/stress/stress_vm.png
            :align: center
            :scale: 75 %

            Contour plot of the von Mises stress.
        """

        title = 'Stress Contour Plot - $\sigma_{vM}$'
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_vm)

        return self.plot_stress_contour(sigs, title, pause)


class MaterialGroup:
    """Class for storing elements of different materials.

    A MaterialGroup object contains the finite element objects for a specified `material`. The
    `stress_result` variable provides storage for stresses related each material.

    :param material: Material object for the current MaterialGroup
    :type material: :class:`~sectionproperties.pre.pre.Material`
    :param int num_nods: Number of nodes for the entire cross-section

    :cvar material: Material object for the current MaterialGroup
    :vartype material: :class:`~sectionproperties.pre.pre.Material`
    :cvar stress_result: A StressResult object for saving the stresses of the current material
    :vartype stress_result: :class:`~sectionproperties.analysis.cross_section.StressResult`
    :cvar elements: A list of finite element objects that are of the current material type
    :vartype elements: list[:class:`~sectionproperties.analysis.fea.Tri6`]
    :cvar el_ids: A list of the element IDs of the elements that are of the current material type
    :vartype el_ids: list[int]
    """

    def __init__(self, material, num_nodes):
        """Inits the MaterialGroup class."""

        self.material = material
        self.stress_result = StressResult(num_nodes)
        self.elements = []
        self.el_ids = []

    def add_element(self, element):
        """Adds an element and its element ID to the MaterialGroup.

        :param element: Element to add to the MaterialGroup
        :type element: :class:`~sectionproperties.analysis.fea.Tri6`
        """

        # add Tri6 element to the list of elements
        self.elements.append(element)
        self.el_ids.append(element.el_id)


class StressResult:
    """Class for storing a stress result.

    Provides variables to store the results from a cross-section stress analysis. Also provides a
    method to calculate combined stresses.

    :param int num_nodes: Number of nodes in the finite element mesh

    :cvar sig_zz_n: Normal stress (:math:`\sigma_{zz,N}`) resulting from an axial force
    :vartype sig_zz_n: :class:`numpy.ndarray`
    :cvar sig_zz_mxx: Normal stress (:math:`\sigma_{zz,Mxx}`) resulting from a bending moment about
        the xx-axis
    :vartype sig_zz_mxx: :class:`numpy.ndarray`
    :cvar sig_zz_myy: Normal stress (:math:`\sigma_{zz,Myy}`) resulting from a bending moment about
        the yy-axis
    :vartype sig_zz_myy: :class:`numpy.ndarray`
    :cvar sig_zz_m11: Normal stress (:math:`\sigma_{zz,M11}`) resulting from a bending moment about
        the 11-axis
    :vartype sig_zz_m11: :class:`numpy.ndarray`
    :cvar sig_zz_m22: Normal stress (:math:`\sigma_{zz,M22}`) resulting from a bending moment about
        the 22-axis
    :vartype sig_zz_m22: :class:`numpy.ndarray`
    :cvar sig_zx_mzz: Shear stress (:math:`\sigma_{zx,Mzz}`) resulting from a torsion moment about
        the zz-axis
    :vartype sig_zx_mzz: :class:`numpy.ndarray`
    :cvar sig_zy_mzz: Shear stress (:math:`\sigma_{zy,Mzz}`) resulting from a torsion moment about
        the zz-axis
    :vartype sig_zy_mzz: :class:`numpy.ndarray`
    :cvar sig_zx_vx: Shear stress (:math:`\sigma_{zx,Vx}`) resulting from a shear force in the
        x-direction
    :vartype sig_zx_vx: :class:`numpy.ndarray`
    :cvar sig_zy_vx: Shear stress (:math:`\sigma_{zy,Vx}`) resulting from a shear force in the
        x-direction
    :vartype sig_zy_vx: :class:`numpy.ndarray`
    :cvar sig_zx_vy: Shear stress (:math:`\sigma_{zx,Vy}`) resulting from a shear force in the
        y-direction
    :vartype sig_zx_vy: :class:`numpy.ndarray`
    :cvar sig_zy_vy: Shear stress (:math:`\sigma_{zy,Vy}`) resulting from a shear force in the
        y-direction
    :vartype sig_zy_vy: :class:`numpy.ndarray`
    :cvar sig_zz_m: Normal stress (:math:`\sigma_{zz,\Sigma M}`) resulting from all bending moments
    :vartype sig_zz_m: :class:`numpy.ndarray`
    :cvar sig_zxy_mzz: Resultant shear stress (:math:`\sigma_{zxy,Mzz}`) resulting from a torsion
        moment in the zz-direction
    :vartype sig_zxy_mzz: :class:`numpy.ndarray`
    :cvar sig_zxy_vx: Resultant shear stress (:math:`\sigma_{zxy,Vx}`) resulting from a a shear
        force in the x-direction
    :vartype sig_zxy_vx: :class:`numpy.ndarray`
    :cvar sig_zxy_vy: Resultant shear stress (:math:`\sigma_{zxy,Vy}`) resulting from a a shear
        force in the y-direction
    :vartype sig_zxy_vy: :class:`numpy.ndarray`
    :cvar sig_zx_v: Shear stress (:math:`\sigma_{zx,\Sigma V}`) resulting from all shear forces
    :vartype sig_zx_v: :class:`numpy.ndarray`
    :cvar sig_zy_v: Shear stress (:math:`\sigma_{zy,\Sigma V}`) resulting from all shear forces
    :vartype sig_zy_v: :class:`numpy.ndarray`
    :cvar sig_zxy_v: Resultant shear stress (:math:`\sigma_{zxy,\Sigma V}`) resulting from all
        shear forces
    :vartype sig_zxy_v: :class:`numpy.ndarray`
    :cvar sig_zz: Combined normal force (:math:`\sigma_{zz}`) resulting from all actions
    :vartype sig_zz: :class:`numpy.ndarray`
    :cvar sig_zx: Combined shear stress (:math:`\sigma_{zx}`) resulting from all actions
    :vartype sig_zx: :class:`numpy.ndarray`
    :cvar sig_zy: Combined shear stress (:math:`\sigma_{zy}`) resulting from all actions
    :vartype sig_zy: :class:`numpy.ndarray`
    :cvar sig_zxy: Combined resultant shear stress (:math:`\sigma_{zxy}`) resulting from all
        actions
    :vartype sig_zxy: :class:`numpy.ndarray`
    :cvar sig_vm: von Mises stress (:math:`\sigma_{VM}`) resulting from all actions
    :vartype sig_vm: :class:`numpy.ndarray`
    """

    def __init__(self, num_nodes):
        """Inits the StressResult class."""

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

        self.sig_zz_m = self.sig_zz_mxx + self.sig_zz_myy + self.sig_zz_m11 + self.sig_zz_m22
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


class SectionProperties:
    """Class for storing section properties.

    Stores calculated section properties. Also provides methods to calculate section properties
    entirely derived from other section properties.

    :cvar float area: Cross-sectional area
    :cvar float perimeter: Cross-sectional perimeter
    :cvar float ea: Modulus weighted area (axial rigidity)
    :cvar float ga: Modulus weighted product of shear modulus and area
    :cvar float nu_eff: Effective Poisson's ratio
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
    :cvar float zxx_plus: Section modulus about the centroidal x-axis for stresses at the positive
        extreme value of y
    :cvar float zxx_minus: Section modulus about the centroidal x-axis for stresses at the negative
        extreme value of y
    :cvar float zyy_plus: Section modulus about the centroidal y-axis for stresses at the positive
        extreme value of x
    :cvar float zyy_minus: Section modulus about the centroidal y-axis for stresses at the negative
        extreme value of x
    :cvar float rx_c: Radius of gyration about the centroidal x-axis.
    :cvar float ry_c: Radius of gyration about the centroidal y-axis.
    :cvar float i11_c: Second moment of area about the centroidal 11-axis
    :cvar float i22_c: Second moment of area about the centroidal 22-axis
    :cvar float phi: Principal axis angle
    :cvar float z11_plus: Section modulus about the principal 11-axis for stresses at the positive
        extreme value of the 22-axis
    :cvar float z11_minus: Section modulus about the principal 11-axis for stresses at the negative
        extreme value of the 22-axis
    :cvar float z22_plus: Section modulus about the principal 22-axis for stresses at the positive
        extreme value of the 11-axis
    :cvar float z22_minus: Section modulus about the principal 22-axis for stresses at the negative
        extreme value of the 11-axis
    :cvar float r11_c: Radius of gyration about the principal 11-axis.
    :cvar float r22_c: Radius of gyration about the principal 22-axis.
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
    :cvar float x11_se: 11 coordinate of the shear centre (elasticity approach)
    :cvar float y22_se: 22 coordinate of the shear centre (elasticity approach)
    :cvar float x_st: X coordinate of the shear centre (Trefftz's approach)
    :cvar float y_st: Y coordinate of the shear centre (Trefftz's approach)
    :cvar float gamma: Warping constant
    :cvar float A_sx: Shear area about the x-axis
    :cvar float A_sy: Shear area about the y-axis
    :cvar float A_sxy: Shear area about the xy-axis
    :cvar float A_s11: Shear area about the 11 bending axis
    :cvar float A_s22: Shear area about the 22 bending axis
    :cvar float beta_x_plus: Monosymmetry constant for bending about the x-axis with the top flange
        in compression
    :cvar float beta_x_minus: Monosymmetry constant for bending about the x-axis with the bottom
        flange in compression
    :cvar float beta_y_plus: Monosymmetry constant for bending about the y-axis with the top flange
        in compression
    :cvar float beta_y_minus: Monosymmetry constant for bending about the y-axis with the bottom
        flange in compression
    :cvar float beta_11_plus: Monosymmetry constant for bending about the 11-axis with the top
        flange in compression
    :cvar float beta_11_minus: Monosymmetry constant for bending about the 11-axis with the bottom
        flange in compression
    :cvar float beta_22_plus: Monosymmetry constant for bending about the 22-axis with the top
        flange in compression
    :cvar float beta_22_minus: Monosymmetry constant for bending about the 22-axis with the bottom
        flange in compression
    :cvar float x_pc: X coordinate of the global plastic centroid
    :cvar float y_pc: Y coordinate of the global plastic centroid
    :cvar float x11_pc: 11 coordinate of the principal plastic centroid
    :cvar float y22_pc: 22 coordinate of the principal plastic centroid
    :cvar float sxx: Plastic section modulus about the centroidal x-axis
    :cvar float syy: Plastic section modulus about the centroidal y-axis
    :cvar float sf_xx_plus: Shape factor for bending about the x-axis with respect to the top fibre
    :cvar float sf_xx_minus: Shape factor for bending about the x-axis with respect to the bottom
        fibre
    :cvar float sf_yy_plus: Shape factor for bending about the y-axis with respect to the top fibre
    :cvar float sf_yy_minus: Shape factor for bending about the y-axis with respect to the bottom
        fibre
    :cvar float s11: Plastic section modulus about the 11-axis
    :cvar float s22: Plastic section modulus about the 22-axis
    :cvar float sf_11_plus: Shape factor for bending about the 11-axis with respect to the top
        fibre
    :cvar float sf_11_minus: Shape factor for bending about the 11-axis with respect to the bottom
        fibre
    :cvar float sf_22_plus: Shape factor for bending about the 22-axis with respect to the top
        fibre
    :cvar float sf_22_minus: Shape factor for bending about the 22-axis with respect to the bottom
        fibre
    """

    def __init__(self):
        """Inits the SectionProperties class."""

        self.area = None
        self.perimeter = None
        self.ea = None
        self.ga = None
        self.nu_eff = None
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
        self.x11_se = None
        self.y22_se = None
        self.x_st = None
        self.y_st = None
        self.gamma = None
        self.A_sx = None
        self.A_sy = None
        self.A_sxy = None
        self.A_s11 = None
        self.A_s22 = None
        self.beta_x_plus = None
        self.beta_x_minus = None
        self.beta_y_plus = None
        self.beta_y_minus = None
        self.beta_11_plus = None
        self.beta_11_minus = None
        self.beta_22_plus = None
        self.beta_22_minus = None
        self.x_pc = None
        self.y_pc = None
        self.x11_pc = None
        self.y22_pc = None
        self.sxx = None
        self.syy = None
        self.sf_xx_plus = None
        self.sf_xx_minus = None
        self.sf_yy_plus = None
        self.sf_yy_minus = None
        self.s11 = None
        self.s22 = None
        self.sf_11_plus = None
        self.sf_11_minus = None
        self.sf_22_plus = None
        self.sf_22_minus = None

    def calculate_elastic_centroid(self):
        """Calculates the elastic centroid based on the cross-section area and first moments of
        area.
        """

        self.cx = self.qy / self.ea
        self.cy = self.qx / self.ea

    def calculate_centroidal_properties(self, mesh):
        """Calculates the geometric section properties about the centroidal and principal axes
        based on the results about the global axis.
        """

        # calculate second moments of area about the centroidal xy axis
        self.ixx_c = self.ixx_g - self.qx ** 2 / self.ea
        self.iyy_c = self.iyy_g - self.qy ** 2 / self.ea
        self.ixy_c = self.ixy_g - self.qx * self.qy / self.ea

        # calculate section moduli about the centroidal xy axis
        nodes = np.array(mesh.points)
        xmax = nodes[:, 0].max()
        xmin = nodes[:, 0].min()
        ymax = nodes[:, 1].max()
        ymin = nodes[:, 1].min()
        self.zxx_plus = self.ixx_c / abs(ymax - self.cy)
        self.zxx_minus = self.ixx_c / abs(ymin - self.cy)
        self.zyy_plus = self.iyy_c / abs(xmax - self.cx)
        self.zyy_minus = self.iyy_c / abs(xmin - self.cx)

        # calculate radii of gyration about centroidal xy axis
        self.rx_c = (self.ixx_c / self.ea) ** 0.5
        self.ry_c = (self.iyy_c / self.ea) ** 0.5

        # calculate principal 2nd moments of area about the centroidal xy axis
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
                x1max = x1
                x1min = x1
                y2max = y2
                y2min = y2

            # update the mins and maxes where necessary
            x1max = max(x1max, x1)
            x1min = min(x1min, x1)
            y2max = max(y2max, y2)
            y2min = min(y2min, y2)

        # evaluate principal section moduli
        self.z11_plus = self.i11_c / abs(y2max)
        self.z11_minus = self.i11_c / abs(y2min)
        self.z22_plus = self.i22_c / abs(x1max)
        self.z22_minus = self.i22_c / abs(x1min)

        # calculate radii of gyration about centroidal principal axis
        self.r11_c = (self.i11_c / self.ea) ** 0.5
        self.r22_c = (self.i22_c / self.ea) ** 0.5
