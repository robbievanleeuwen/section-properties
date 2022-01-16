from typing import Union, Optional, Tuple

import copy
from dataclasses import dataclass, asdict

import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.cm as cm
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap, CenteredNorm
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

import numpy as np
from scipy.sparse import csc_matrix, coo_matrix, linalg
from scipy.optimize import brentq
import sectionproperties.pre.pre as pre
import sectionproperties.pre.geometry as section_geometry
import sectionproperties.analysis.fea as fea
import sectionproperties.analysis.solver as solver
import sectionproperties.post.post as post


class Section:
    """Class for structural cross-sections.

    Stores the finite element geometry, mesh and material information and provides methods to
    compute the cross-section properties. The element type used in this program is the six-noded
    quadratic triangular element.

    The constructor extracts information from the provided mesh object and creates and stores the
    corresponding Tri6 finite element objects.

    :param geometry: Cross-section geometry object used to generate the mesh
    :type geometry: :class:`~sectionproperties.pre.geometry.Geometry`
    :param bool time_info: If set to True, a detailed description of the computation and the time
        cost is printed to the terminal for every computation performed.

    The following example creates a :class:`~sectionproperties.analysis.section.Section`
    object of a 100D x 50W rectangle using a mesh size of 5::

        import sectionproperties.pre.library.primitive_sections as primitive_sections
        from sectionproperties.analysis.section import Section

        geometry = primitive_sections.rectangular_section(d=100, b=50)
        geometry.create_mesh(mesh_sizes=[5])
        section = Section(geometry)

    :cvar elements: List of finite element objects describing the cross-section mesh
    :vartype elements: list[:class:`~sectionproperties.analysis.fea.Tri6`]
    :cvar int num_nodes: Number of nodes in the finite element mesh
    :cvar geometry: Cross-section geometry object used to generate the mesh
    :vartype geometry: :class:`~sectionproperties.pre.geometry.Geometry`
    :cvar mesh: Mesh dict returned by triangle
    :vartype mesh: dict(mesh)
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
    :vartype section_props: :class:`~sectionproperties.analysis.section.SectionProperties`

    :raises AssertionError: If the number of materials does not equal the number of regions
    :raises ValueError: If geometry does not contain a mesh
    """

    def __init__(
        self,
        geometry: Union[section_geometry.Geometry, section_geometry.CompoundGeometry],
        time_info: bool = False,
    ):
        """Inits the Section class."""
        if not hasattr(geometry, "mesh") or not geometry.mesh:
            raise ValueError(
                "Selected Geometry or CompoundGeometry "
                "object does not contain a mesh.\n"
                "Try running {geometry}.create_mesh() before adding to "
                "a Section object for analysis."
            )
        self.geometry = geometry
        self.time_info = time_info
        self.mesh = geometry.mesh
        self.materials = []
        mesh = self.mesh

        def init():
            if isinstance(self.geometry, section_geometry.CompoundGeometry):
                self.materials = [geom.material for geom in self.geometry.geoms]
            else:
                self.materials = [self.geometry.material]

            # extract mesh data
            nodes = np.array(mesh["vertices"], dtype=np.dtype(float))
            elements = np.array(mesh["triangles"], dtype=np.dtype(int))
            attributes = np.array(mesh["triangle_attributes"].T[0], dtype=np.dtype(int))

            # swap mid-node order to retain node ordering consistency
            elements[:, [3, 4, 5]] = elements[:, [5, 3, 4]]

            # save total number of nodes in mesh
            self.num_nodes = len(nodes)

            # initialise material_sections variable
            self.material_groups = []

            # if materials are specified, check that the right number of material properties are
            # specified and then populate material_groups list
            if self.materials:
                msg = "Number of materials ({0}), ".format(len(self.materials))
                msg += "should match the number of regions ({0}).".format(
                    max(attributes) + 1
                )
                assert len(self.materials) == max(attributes) + 1, msg

                # add a MaterialGroup object to the material_groups list for each uniquely
                # encountered material
                for (i, material) in enumerate(self.materials):
                    # add the first material to the list
                    if i == 0:
                        self.material_groups.append(
                            MaterialGroup(material, self.num_nodes)
                        )
                    else:
                        # if the material hasn't been encountered
                        if material not in self.materials[:i]:
                            self.material_groups.append(
                                MaterialGroup(material, self.num_nodes)
                            )

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
                if self.materials:
                    # get attribute index of current element
                    att_el = attributes[i]

                    # fetch the material
                    material = self.materials[att_el]
                # if there are no materials specified, use a default material
                else:  # Should not happen but included as failsafe
                    material = pre.DEFAULT_MATERIAL

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

        if self.time_info:
            text = "--Initialising the Section class..."
            solver.function_timer(text, init)
        else:
            init()

    def calculate_geometric_properties(self):
        """Calculates the geometric properties of the cross-section and stores them in the
        :class:`~sectionproperties.analysis.section.SectionProperties` object contained in
        the ``section_props`` class variable.

        The following geometric section properties are calculated:

        * Cross-sectional area
        * Cross-sectional perimeter
        * Cross-sectional mass
        * Area weighted material properties, composite only :math:`E_{eff}`, :math:`G_{eff}`, :math:`{nu}_{eff}`
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

            section = Section(geometry)
            section.calculate_geometric_properties()
        """

        def calculate_geom():
            # initialise properties
            self.section_props.area = 0
            self.section_props.perimeter = 0
            self.section_props.mass = 0
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
                (
                    area,
                    qx,
                    qy,
                    ixx_g,
                    iyy_g,
                    ixy_g,
                    e,
                    g,
                    rho,
                ) = el.geometric_properties()
                self.section_props.area += area
                self.section_props.mass += area * rho
                self.section_props.ea += area * e
                self.section_props.ga += area * g
                self.section_props.qx += qx * e
                self.section_props.qy += qy * e
                self.section_props.ixx_g += ixx_g * e
                self.section_props.iyy_g += iyy_g * e
                self.section_props.ixy_g += ixy_g * e

            self.section_props.nu_eff = (
                self.section_props.ea / (2 * self.section_props.ga) - 1
            )
            self.section_props.e_eff = self.section_props.ea / self.section_props.area
            self.section_props.g_eff = self.section_props.ga / self.section_props.area
            self.section_props.calculate_elastic_centroid()
            self.section_props.calculate_centroidal_properties(self.mesh)

        if self.time_info:
            text = "--Calculating geometric section properties..."
            solver.function_timer(text, calculate_geom)
        else:
            calculate_geom()

    def calculate_warping_properties(self, solver_type="direct"):
        """Calculates all the warping properties of the cross-section and stores them in the
        :class:`~sectionproperties.analysis.section.SectionProperties` object contained in
        the ``section_props`` class variable.

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

        Note that the geometric properties must be calculated prior to the calculation of the
        warping properties::

            section = Section(geometry)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()

        :raises RuntimeError: If the geometric properties have not been
            calculated prior to calling this method
        """

        # check that a geometric analysis has been performed
        if None in [
            self.section_props.area,
            self.section_props.ixx_c,
            self.section_props.cx,
        ]:
            err = "Calculate geometric properties before performing a warping analysis."
            raise RuntimeError(err)

        # create a new Section with the origin shifted to the centroid for calculation of the
        # warping properties such that the Lagrangian multiplier approach can be utilised
        warping_section = Section(self.geometry)

        # shift the coordinates of each element N.B. the mesh class attribute remains unshifted!
        for el in warping_section.elements:
            el.coords[0, :] -= self.section_props.cx
            el.coords[1, :] -= self.section_props.cy

        # assemble stiffness matrix and load vector for warping function
        if self.time_info:
            text = "--Assembling {0}x{0} stiffness matrix and load vector...".format(
                self.num_nodes
            )
            (k, k_lg, f_torsion) = solver.function_timer(
                text, warping_section.assemble_torsion
            )
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
        if solver_type == "cgs":
            if self.time_info:
                text = "--Performing ILU decomposition on the stiffness matrices..."
                (k_precond, k_lg_precond) = solver.function_timer(text, ilu_decomp)
            else:
                (k_precond, k_lg_precond) = ilu_decomp()

        # solve for warping function
        def solve_warping():
            if solver_type == "cgs":
                omega = solver.solve_cgs(k, f_torsion, k_precond)
            elif solver_type == "direct":
                omega = solver.solve_direct(k, f_torsion)

            return omega

        if self.time_info:
            text = "--Solving for the warping function using the {0} solver...".format(
                solver_type
            )
            omega = solver.function_timer(text, solve_warping)
        else:
            omega = solve_warping()

        # save the warping function
        self.section_props.omega = omega

        # determine the torsion constant
        def j_func():
            return (
                self.section_props.ixx_c
                + self.section_props.iyy_c
                - omega.dot(k.dot(np.transpose(omega)))
            )

        if self.time_info:
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
                    self.section_props.ixx_c,
                    self.section_props.iyy_c,
                    self.section_props.ixy_c,
                    self.section_props.nu_eff,
                )
                f_psi[el.node_ids] += f_psi_el
                f_phi[el.node_ids] += f_phi_el

            return (f_psi, f_phi)

        if self.time_info:
            text = "--Assembling shear function load vectors..."
            (f_psi, f_phi) = solver.function_timer(text, assemble_shear_load)
        else:
            (f_psi, f_phi) = assemble_shear_load()

        # solve for shear functions psi and phi
        def solve_shear_functions():
            if solver_type == "cgs":
                psi_shear = solver.solve_cgs_lagrange(k_lg, f_psi, m=k_lg_precond)
                phi_shear = solver.solve_cgs_lagrange(k_lg, f_phi, m=k_lg_precond)
            elif solver_type == "direct":
                psi_shear = solver.solve_direct_lagrange(k_lg, f_psi)
                phi_shear = solver.solve_direct_lagrange(k_lg, f_phi)

            return (psi_shear, phi_shear)

        if self.time_info:
            text = "--Solving for the shear functions using the {0} solver...".format(
                solver_type
            )
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
                (
                    sc_xint_el,
                    sc_yint_el,
                    q_omega_el,
                    i_omega_el,
                    i_xomega_el,
                    i_yomega_el,
                ) = el.shear_warping_integrals(
                    self.section_props.ixx_c,
                    self.section_props.iyy_c,
                    self.section_props.ixy_c,
                    omega[el.node_ids],
                )

                sc_xint += sc_xint_el
                sc_yint += sc_yint_el
                q_omega += q_omega_el
                i_omega += i_omega_el
                i_xomega += i_xomega_el
                i_yomega += i_yomega_el

            return (sc_xint, sc_yint, q_omega, i_omega, i_xomega, i_yomega)

        if self.time_info:
            text = "--Assembling shear centre and warping moment integrals..."
            (
                sc_xint,
                sc_yint,
                q_omega,
                i_omega,
                i_xomega,
                i_yomega,
            ) = solver.function_timer(text, assemble_sc_warping_integrals)
        else:
            (
                sc_xint,
                sc_yint,
                q_omega,
                i_omega,
                i_xomega,
                i_yomega,
            ) = assemble_sc_warping_integrals()

        # calculate shear centres
        def shear_centres():
            # calculate shear centres (elasticity approach)
            Delta_s = (
                2
                * (1 + self.section_props.nu_eff)
                * (
                    self.section_props.ixx_c * self.section_props.iyy_c
                    - self.section_props.ixy_c ** 2
                )
            )
            x_se = (1 / Delta_s) * (
                (self.section_props.nu_eff / 2 * sc_xint) - f_torsion.dot(phi_shear)
            )
            y_se = (1 / Delta_s) * (
                (self.section_props.nu_eff / 2 * sc_yint) + f_torsion.dot(psi_shear)
            )
            (x11_se, y22_se) = fea.principal_coordinate(
                self.section_props.phi, x_se, y_se
            )

            # calculate shear centres (Trefftz's approach)
            x_st = (
                self.section_props.ixy_c * i_xomega
                - self.section_props.iyy_c * i_yomega
            ) / (
                self.section_props.ixx_c * self.section_props.iyy_c
                - self.section_props.ixy_c ** 2
            )
            y_st = (
                self.section_props.ixx_c * i_xomega
                - self.section_props.ixy_c * i_yomega
            ) / (
                self.section_props.ixx_c * self.section_props.iyy_c
                - self.section_props.ixy_c ** 2
            )

            return (Delta_s, x_se, y_se, x11_se, y22_se, x_st, y_st)

        if self.time_info:
            text = "--Calculating shear centres..."
            (Delta_s, x_se, y_se, x11_se, y22_se, x_st, y_st) = solver.function_timer(
                text, shear_centres
            )
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
            i_omega
            - q_omega ** 2 / self.section_props.ea
            - y_se * i_xomega
            + x_se * i_yomega
        )

        def assemble_shear_deformation():
            # assemble shear deformation coefficients
            kappa_x = 0
            kappa_y = 0
            kappa_xy = 0

            for el in warping_section.elements:
                (kappa_x_el, kappa_y_el, kappa_xy_el) = el.shear_coefficients(
                    self.section_props.ixx_c,
                    self.section_props.iyy_c,
                    self.section_props.ixy_c,
                    psi_shear[el.node_ids],
                    phi_shear[el.node_ids],
                    self.section_props.nu_eff,
                )

                kappa_x += kappa_x_el
                kappa_y += kappa_y_el
                kappa_xy += kappa_xy_el

            return (kappa_x, kappa_y, kappa_xy)

        if self.time_info:
            text = "--Assembling shear deformation coefficients..."
            (kappa_x, kappa_y, kappa_xy) = solver.function_timer(
                text, assemble_shear_deformation
            )
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
        R = np.array(
            [[np.cos(phi_rad), np.sin(phi_rad)], [-np.sin(phi_rad), np.cos(phi_rad)]]
        )

        rotatedAlpha = R.dot(
            np.array([[alpha_xx, alpha_xy], [alpha_xy, alpha_yy]])
        ).dot(np.transpose(R))

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

        if self.time_info:
            text = "--Assembling monosymmetry integrals..."
            (int_x, int_y, int_11, int_22) = solver.function_timer(
                text, calculate_monosymmetry_integrals
            )
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

    def calculate_frame_properties(self, solver_type="direct"):
        """Calculates and returns the properties required for a frame analysis. The properties are
        also stored in the :class:`~sectionproperties.analysis.section.SectionProperties`
        object contained in the ``section_props`` class variable.

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

            section = Section(geometry)
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
                (area, qx, qy, ixx_g, iyy_g, ixy_g, e, _, _) = el.geometric_properties()

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
                self.section_props.ixx_g
                - self.section_props.qx ** 2 / self.section_props.ea
            )
            self.section_props.iyy_c = (
                self.section_props.iyy_g
                - self.section_props.qy ** 2 / self.section_props.ea
            )
            self.section_props.ixy_c = (
                self.section_props.ixy_g
                - self.section_props.qx * self.section_props.qy / self.section_props.ea
            )

            # calculate the principal axis angle
            Delta = (
                ((self.section_props.ixx_c - self.section_props.iyy_c) / 2) ** 2
                + self.section_props.ixy_c ** 2
            ) ** 0.5

            i11_c = (self.section_props.ixx_c + self.section_props.iyy_c) / 2 + Delta

            # calculate initial principal axis angle
            if abs(self.section_props.ixx_c - i11_c) < 1e-12 * i11_c:
                self.section_props.phi = 0
            else:
                self.section_props.phi = (
                    np.arctan2(
                        self.section_props.ixx_c - i11_c, self.section_props.ixy_c
                    )
                    * 180
                    / np.pi
                )

            # create a new Section with the origin shifted to the centroid for calculation of
            # the warping properties
            warping_section = Section(self.geometry, self.time_info)

            # shift the coordinates of each element N.B. the mesh class attribute remains unshifted
            for el in warping_section.elements:
                el.coords[0, :] -= self.section_props.cx
                el.coords[1, :] -= self.section_props.cy

            (k, _, f) = warping_section.assemble_torsion(lg=False)

            # if the cgs method is used, perform ILU decomposition
            if solver_type == "cgs":
                k_precond = linalg.LinearOperator(
                    (self.num_nodes, self.num_nodes), linalg.spilu(k).solve
                )

            # solve for warping function
            if solver_type == "cgs":
                omega = solver.solve_cgs(k, f, k_precond)
            elif solver_type == "direct":
                omega = solver.solve_direct(k, f)

            # calculate the torsion constant
            self.section_props.j = (
                self.section_props.ixx_c
                + self.section_props.iyy_c
                - omega.dot(k.dot(np.transpose(omega)))
            )

        if self.time_info:
            text = "--Calculating frame section properties..."
            solver.function_timer(text, calculate_frame)
        else:
            calculate_frame()

        return (
            self.section_props.ea,
            self.section_props.ixx_c,
            self.section_props.iyy_c,
            self.section_props.ixy_c,
            self.section_props.j,
            self.section_props.phi,
        )

    def calculate_plastic_properties(self, verbose=False, debug=False):
        """Calculates the plastic properties of the cross-section and stores the, in the
        :class:`~sectionproperties.analysis.section.SectionProperties` object contained in
        the ``section_props`` class variable.

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

            section = Section(geometry)
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
            plastic_section = PlasticSection(self.geometry)
            plastic_section.calculate_plastic_properties(self, verbose)

        if self.time_info:
            text = "--Calculating plastic properties..."
            solver.function_timer(text, calc_plastic)
        else:
            calc_plastic()

    def calculate_stress(self, N=0, Vx=0, Vy=0, Mxx=0, Myy=0, M11=0, M22=0, Mzz=0):
        """Calculates the cross-section stress resulting from design actions and returns a
        :class:`~sectionproperties.analysis.section.StressPost` object allowing
        post-processing of the stress results.

        :param float N: Axial force
        :param float Vx: Shear force acting in the x-direction
        :param float Vy: Shear force acting in the y-direction
        :param float Mxx: Bending moment about the centroidal xx-axis
        :param float Myy: Bending moment about the centroidal yy-axis
        :param float M11: Bending moment about the centroidal 11-axis
        :param float M22: Bending moment about the centroidal 22-axis
        :param float Mzz: Torsion moment about the centroidal zz-axis
        :return: Object for post-processing cross-section stresses
        :rtype: :class:`~sectionproperties.analysis.section.StressPost`

        Note that a geometric and warping analysis must be performed before a stress analysis is
        carried out::

            section = Section(geometry)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(N=1e3, Vy=3e3, Mxx=1e6)

        :raises RuntimeError: If a geometric and warping analysis have not been performed prior to
            calling this method
        """

        # check that a geometric and warping analysis has been performed
        if (
            None
            in [
                self.section_props.area,
                self.section_props.ixx_c,
                self.section_props.cx,
                self.section_props.j,
            ]
            and self.section_props.omega is None
        ):
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
                        sig_zz_n_el,
                        sig_zz_mxx_el,
                        sig_zz_myy_el,
                        sig_zz_m11_el,
                        sig_zz_m22_el,
                        sig_zx_mzz_el,
                        sig_zy_mzz_el,
                        sig_zx_vx_el,
                        sig_zy_vx_el,
                        sig_zx_vy_el,
                        sig_zy_vy_el,
                        weights,
                    ) = el.element_stress(
                        N,
                        Mxx,
                        Myy,
                        M11,
                        M22,
                        Mzz,
                        Vx,
                        Vy,
                        ea,
                        cx,
                        cy,
                        ixx,
                        iyy,
                        ixy,
                        i11,
                        i22,
                        phi,
                        j,
                        nu,
                        self.section_props.omega[el.node_ids],
                        self.section_props.psi_shear[el.node_ids],
                        self.section_props.phi_shear[el.node_ids],
                        Delta_s,
                    )

                    # add stresses to global vectors
                    group.stress_result.sig_zz_n[el.node_ids] += sig_zz_n_el * weights
                    group.stress_result.sig_zz_mxx[el.node_ids] += (
                        sig_zz_mxx_el * weights
                    )
                    group.stress_result.sig_zz_myy[el.node_ids] += (
                        sig_zz_myy_el * weights
                    )
                    group.stress_result.sig_zz_m11[el.node_ids] += (
                        sig_zz_m11_el * weights
                    )
                    group.stress_result.sig_zz_m22[el.node_ids] += (
                        sig_zz_m22_el * weights
                    )
                    group.stress_result.sig_zx_mzz[el.node_ids] += (
                        sig_zx_mzz_el * weights
                    )
                    group.stress_result.sig_zy_mzz[el.node_ids] += (
                        sig_zy_mzz_el * weights
                    )
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

        if self.time_info:
            text = "--Calculating cross-section stresses..."
            stress_post = solver.function_timer(text, calc_stress)
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

    def plot_mesh(
        self,
        alpha=0.5,
        materials=True,
        mask=None,
        title="Finite Element Mesh",
        **kwargs,
    ):
        r"""Plots the finite element mesh.

        :param float alpha: Transparency of the mesh outlines: :math:`0 \leq \alpha \leq 1`
        :param bool materials: If set to true shades the elements with the specified material
            colours
        :param mask: Mask array, of length ``num_nodes``, to mask out triangles
        :type mask: list[bool]
        :param string title: Plot title
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots the mesh generated for the second example
        listed under the :class:`~sectionproperties.analysis.section.Section` object
        definition::

            import sectionproperties.pre.library.primitive_sections as primitive_sections
            from sectionproperties.pre.pre import Material
            from sectionproperties.analysis.section import Section

            steel = Material(
                name='Steel', elastic_modulus=200e3, poissons_ratio=0.3, density=7.85e-6,
                yield_strength=250, color='grey'
            )
            timber = Material(
                name='Timber', elastic_modulus=8e3, poissons_ratio=0.35, density=6.5e-7,
                yield_strength=20, color='burlywood'
            )

            geom_steel = primitive_sections.rectangular_section(d=50, b=50, material=steel)
            geom_timber = primitive_sections.rectangular_section(d=50, b=50, material=timber)
            geometry = geom_timber.align_to(geom_steel, on="right") + geom_steel

            geometry.create_mesh(mesh_sizes=[10, 5])

            section = Section(geometry)
            section.plot_mesh(materials=True, alpha=0.5)

        ..  figure:: ../images/composite_mesh.png
            :align: center
            :scale: 75 %

            Finite element mesh generated by the above example.
        """

        # create plot and setup the plot
        with post.plotting_context(title=title, **kwargs) as (fig, ax):
            # if the material colours are to be displayed
            if materials and self.materials is not None:
                color_array = []
                legend_labels = []
                c = []  # Indices of elements for mapping colors

                # create an array of finite element colours
                for idx, element in enumerate(self.elements):
                    color_array.append(element.material.color)
                    c.append(idx)

                # create a list of unique material legend entries
                for (i, material) in enumerate(self.materials):
                    # if the material has not be entered yet
                    if i == 0 or material not in self.materials[0:i]:
                        # add the material colour and name to the legend list
                        patch = mpatches.Patch(
                            color=material.color, label=material.name
                        )
                        legend_labels.append(patch)

                cmap = ListedColormap(color_array)  # custom colormap

                # plot the mesh colours
                ax.tripcolor(
                    self.mesh_nodes[:, 0],
                    self.mesh_nodes[:, 1],
                    self.mesh_elements[:, 0:3],
                    c,
                    cmap=cmap,
                )

                # display the legend
                ax.legend(
                    loc="center left", bbox_to_anchor=(1, 0.5), handles=legend_labels
                )

            # plot the mesh
            ax.triplot(
                self.mesh_nodes[:, 0],
                self.mesh_nodes[:, 1],
                self.mesh_elements[:, 0:3],
                lw=0.5,
                color="black",
                alpha=alpha,
                mask=mask,
            )

        return ax

    def plot_centroids(self, title="Centroids", **kwargs):
        """Plots the elastic centroid, the shear centre, the plastic centroids and the principal
        axis, if they have been calculated, on top of the finite element mesh.

        :param string title: Plot title
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example analyses a 200 PFC section and displays a plot of
        the centroids::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.channel_section(d=200, b=75, t_f=12, t_w=6, r=12, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])

            section = Section(geometry)
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

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])

            section = Section(geometry)
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
        with post.plotting_context(title=title, **kwargs) as (fig, ax):
            # plot the finite element mesh
            self.plot_mesh(**dict(kwargs, ax=ax))

            # if the elastic centroid has been calculated
            if self.section_props.cx is not None:
                ax.scatter(
                    self.section_props.cx,
                    self.section_props.cy,
                    edgecolors="r",
                    facecolors="none",
                    marker="o",
                    s=100,
                    label="Elastic centroid",
                )

            # if the shear centre has been calculated
            if self.section_props.x_se is not None:
                (x_s, y_s) = self.get_sc()
                ax.scatter(x_s, y_s, c="r", marker="+", s=100, label="Shear centre")

            # if the global plastic centroid has been calculated
            if self.section_props.x_pc is not None:
                (x_pc, y_pc) = self.get_pc()
                ax.scatter(
                    x_pc,
                    y_pc,
                    c="r",
                    marker="x",
                    s=100,
                    label="Global plastic centroid",
                )

            # if the principal plastic centroid has been calculated
            if self.section_props.x11_pc is not None:
                (x11_pc, y22_pc) = self.get_pc_p()
                ax.scatter(
                    x11_pc,
                    y22_pc,
                    edgecolors="r",
                    facecolors="none",
                    marker="s",
                    s=100,
                    label="Principal plastic centroid",
                )

            # if the principal axis has been calculated
            if self.section_props.phi is not None:
                post.draw_principal_axis(
                    ax,
                    self.section_props.phi * np.pi / 180,
                    self.section_props.cx,
                    self.section_props.cy,
                )

            # display the legend
            ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

        return ax

    def display_mesh_info(self):
        """Prints mesh statistics (number of nodes, elements and regions) to the command window.

        The following example displays the mesh statistics for a Tee section merged from two
        rectangles::

            import sectionproperties.pre.library.primitive_sections as primitive_sections
            from sectionproperties.analysis.section import Section

            rec1 = primitive_sections.rectangular_section(d=100, b=25)
            rec2 = primitive_sections.rectangular_section(d=25, b=100)
            rec1 = rec1.shift_section(x_offset=-12.5)
            rec2 = rec2.shift_section(x_offset=-50, y_offset=100)

            geometry = rec1 + rec2
            geometry.create_mesh(mesh_sizes=[5, 2.5])
            section = Section(geometry)
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

    def display_results(self, fmt="8.6e"):
        """Prints the results that have been calculated to the terminal.

        :param string fmt: Number formatting string

        The following example displays the geometric section properties for a 100D x 50W rectangle
        with three digits after the decimal point::

            import sectionproperties.pre.library.primitive_sections as primitive_sections
            from sectionproperties.analysis.section import Section

            geometry = primitive_sections.rectangular_section(d=100, b=50)
            geometry.create_mesh(mesh_sizes=[5])

            section = Section(geometry)
            section.calculate_geometric_properties()

            section.display_results(fmt='.3f')
        """

        post.print_results(self, fmt)

    def get_area(self):
        """
        :return: Cross-section area
        :rtype: float

        ::

            section = Section(geometry)
            section.calculate_geometric_properties()
            area = section.get_area()
        """

        return self.section_props.area

    def get_perimeter(self):
        """
        :return: Cross-section perimeter
        :rtype: float

        ::

            section = Section(geometry)
            section.calculate_geometric_properties()
            perimeter = section.get_perimeter()
        """

        return self.section_props.perimeter

    def get_mass(self):
        """
        :return: Cross-section mass
        :rtype: float

        ::

            section = Section(geometry)
            section.calculate_geometric_properties()
            perimeter = section.get_mass()
        """

        return self.section_props.mass

    def get_ea(self):
        """
        :return: Modulus weighted area (axial rigidity)
        :rtype: float

        ::

            section = Section(geometry)
            section.calculate_geometric_properties()
            ea = section.get_ea()
        """

        return self.section_props.ea

    def get_q(self):
        """
        :return: First moments of area about the global axis *(qx, qy)*
        :rtype: tuple(float, float)

        ::

            section = Section(geometry)
            section.calculate_geometric_properties()
            (qx, qy) = section.get_q()
        """

        return (self.section_props.qx, self.section_props.qy)

    def get_ig(self):
        """
        :return: Second moments of area about the global axis *(ixx_g, iyy_g, ixy_g)*
        :rtype: tuple(float, float, float)

        ::

            section = Section(geometry)
            section.calculate_geometric_properties()
            (ixx_g, iyy_g, ixy_g) = section.get_ig()
        """

        return (
            self.section_props.ixx_g,
            self.section_props.iyy_g,
            self.section_props.ixy_g,
        )

    def get_c(self):
        """
        :return: Elastic centroid *(cx, cy)*
        :rtype: tuple(float, float)

        ::

            section = Section(geometry)
            section.calculate_geometric_properties()
            (cx, cy) = section.get_c()
        """

        return (self.section_props.cx, self.section_props.cy)

    def get_ic(self):
        """
        :return: Second moments of area centroidal axis *(ixx_c, iyy_c, ixy_c)*
        :rtype: tuple(float, float, float)

        ::

            section = Section(geometry)
            section.calculate_geometric_properties()
            (ixx_c, iyy_c, ixy_c) = section.get_ic()
        """

        return (
            self.section_props.ixx_c,
            self.section_props.iyy_c,
            self.section_props.ixy_c,
        )

    def get_z(self):
        """
        :return: Elastic section moduli about the centroidal axis with respect to the top and
            bottom fibres *(zxx_plus, zxx_minus, zyy_plus, zyy_minus)*
        :rtype: tuple(float, float, float, float)

        ::

            section = Section(geometry)
            section.calculate_geometric_properties()
            (zxx_plus, zxx_minus, zyy_plus, zyy_minus) = section.get_z()
        """

        return (
            self.section_props.zxx_plus,
            self.section_props.zxx_minus,
            self.section_props.zyy_plus,
            self.section_props.zyy_minus,
        )

    def get_rc(self):
        """
        :return: Radii of gyration about the centroidal axis *(rx, ry)*
        :rtype: tuple(float, float)

        ::

            section = Section(geometry)
            section.calculate_geometric_properties()
            (rx, ry) = section.get_rc()
        """

        return (self.section_props.rx_c, self.section_props.ry_c)

    def get_ip(self):
        """
        :return: Second moments of area about the principal axis *(i11_c, i22_c)*
        :rtype: tuple(float, float)

        ::

            section = Section(geometry)
            section.calculate_geometric_properties()
            (i11_c, i22_c) = section.get_ip()
        """

        return (self.section_props.i11_c, self.section_props.i22_c)

    def get_phi(self):
        """
        :return: Principal bending axis angle
        :rtype: float

        ::

            section = Section(geometry)
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

            section = Section(geometry)
            section.calculate_geometric_properties()
            (z11_plus, z11_minus, z22_plus, z22_minus) = section.get_zp()
        """

        return (
            self.section_props.z11_plus,
            self.section_props.z11_minus,
            self.section_props.z22_plus,
            self.section_props.z22_minus,
        )

    def get_rp(self):
        """
        :return: Radii of gyration about the principal axis *(r11, r22)*
        :rtype: tuple(float, float)

        ::

            section = Section(geometry)
            section.calculate_geometric_properties()
            (r11, r22) = section.get_rp()
        """

        return (self.section_props.r11_c, self.section_props.r22_c)

    def get_nu_eff(self):
        """
        :return: Effective Poisson's ratio
        :rtype: float

        ::

            section = Section(geometry)
            section.calculate_geometric_properties()
            nu_eff = section.get_nu_eff()
        """

        return self.section_props.nu_eff

    def get_e_eff(self):
        """
        :return: Effective elastic modulus based on area
        :rtype: float

        ::

            section = Section(geometry)
            section.calculate_warping_properties()
            e_eff = section.get_e_eff()
        """

        return self.section_props.e_eff

    def get_g_eff(self):
        """
        :return: Effective shear modulus based on area
        :rtype: float

        ::

            section = Section(geometry)
            section.calculate_geometric_properties()
            g_eff = section.get_g_eff()
        """

        return self.section_props.g_eff

    def get_j(self):
        """
        :return: St. Venant torsion constant
        :rtype: float

        ::

            section = Section(geometry)
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

            section = Section(geometry)
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
        :return: Principal axis shear centre (elasticity approach) *(x11_se, y22_se)*
        :rtype: tuple(float, float)

        ::

            section = Section(geometry)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            (x11_se, y22_se) = section.get_sc_p()
        """

        if self.section_props.x11_se is None:
            return (None, None)
        else:
            x11_se = self.section_props.x11_se
            y22_se = self.section_props.y22_se

        return (x11_se, y22_se)

    def get_sc_t(self):
        """
        :return: Centroidal axis shear centre (Trefftz's approach) *(x_st, y_st)*
        :rtype: tuple(float, float)

        ::

            section = Section(geometry)
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

            section = Section(geometry)
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

            section = Section(geometry)
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

            section = Section(geometry)
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

            section = Section(geometry)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            (beta_x_plus, beta_x_minus, beta_y_plus, beta_y_minus) = section.get_beta()
        """

        return (
            self.section_props.beta_x_plus,
            self.section_props.beta_x_minus,
            self.section_props.beta_y_plus,
            self.section_props.beta_y_minus,
        )

    def get_beta_p(self):
        """
        :return: Monosymmetry constant for bending about both principal axes *(beta_11_plus,
            beta_11_minus, beta_22_plus, beta_22_minus)*. The *plus* value relates to the top
            flange in compression and the *minus* value relates to the bottom flange in
            compression.
        :rtype: tuple(float, float, float, float)

        ::

            section = Section(geometry)
            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            (beta_11_plus, beta_11_minus, beta_22_plus, beta_22_minus) = section.get_beta_p()
        """

        return (
            self.section_props.beta_11_plus,
            self.section_props.beta_11_minus,
            self.section_props.beta_22_plus,
            self.section_props.beta_22_minus,
        )

    def get_pc(self):
        """
        :return: Centroidal axis plastic centroid *(x_pc, y_pc)*
        :rtype: tuple(float, float)

        ::

            section = Section(geometry)
            section.calculate_geometric_properties()
            section.calculate_plastic_properties()
            (x_pc, y_pc) = section.get_pc()
        """

        if self.section_props.x_pc is None:
            return (None, None)
        else:
            # add centroid location to move section back to original location
            x_pc = self.section_props.x_pc + self.section_props.cx
            y_pc = self.section_props.y_pc + self.section_props.cy

        return (x_pc, y_pc)

    def get_pc_p(self):
        """
        :return: Principal bending axis plastic centroid *(x11_pc, y22_pc)*
        :rtype: tuple(float, float)

        ::

            section = Section(geometry)
            section.calculate_geometric_properties()
            section.calculate_plastic_properties()
            (x11_pc, y22_pc) = section.get_pc_p()
        """

        if self.section_props.x11_pc is None:
            return (None, None)
        else:
            # determine the position of the plastic centroid in the global axis
            (x_pc, y_pc) = fea.global_coordinate(
                self.section_props.phi,
                self.section_props.x11_pc,
                self.section_props.y22_pc,
            )

            # add centroid location to move section back to original location
            return (x_pc + self.section_props.cx, y_pc + self.section_props.cy)

    def get_s(self):
        """
        :return: Plastic section moduli about the centroidal axis *(sxx, syy)*
        :rtype: tuple(float, float)

        If material properties have been specified, returns the plastic moment :math:`M_p = f_y S`.

        ::

            section = Section(geometry)
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

            section = Section(geometry)
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

            section = Section(geometry)
            section.calculate_geometric_properties()
            section.calculate_plastic_properties()
            (sf_xx_plus, sf_xx_minus, sf_yy_plus, sf_yy_minus) = section.get_sf()
        """

        return (
            self.section_props.sf_xx_plus,
            self.section_props.sf_xx_minus,
            self.section_props.sf_yy_plus,
            self.section_props.sf_yy_minus,
        )

    def get_sf_p(self):
        """
        :return: Principal bending axis shape factors with respect to the top and bottom fibres
            *(sf_11_plus, sf_11_minus, sf_22_plus, sf_22_minus)*
        :rtype: tuple(float, float, float, float)

        ::

            section = Section(geometry)
            section.calculate_geometric_properties()
            section.calculate_plastic_properties()
            (sf_11_plus, sf_11_minus, sf_22_plus, sf_22_minus) = section.get_sf_p()
        """

        return (
            self.section_props.sf_11_plus,
            self.section_props.sf_11_minus,
            self.section_props.sf_22_plus,
            self.section_props.sf_22_minus,
        )


class PlasticSection:
    """Class for the plastic analysis of cross-sections.

    Stores the finite element geometry and material information and provides methods to compute the
    plastic section properties.

    :param section: Section object
    :type section: :class:`~sectionproperties.analysis.section.Section`

    :param bool debug: If set to True, the geometry is plotted each time a new mesh is generated by
        the plastic centroid algorithm.

    :cvar geometry: Deep copy of the Section geometry object provided to the constructor
    :vartype geometry: :class:`~sectionproperties.pre.geometry.Geometry`
    :cvar materials: A list of material properties corresponding to various regions in the geometry
        and mesh.
    :vartype materials: list[:class:`~sectionproperties.pre.pre.Material`]
    :cvar bool debug: If set to True, the geometry is plotted each time a new mesh is generated by
        the plastic centroid algorithm.
    :cvar mesh: Mesh dict returned by triangle
    :vartype mesh: dict(mesh)
    :cvar mesh_nodes: Array of node coordinates from the mesh
    :vartype mesh_nodes: :class:`numpy.ndarray`
    :cvar mesh_elements: Array of connectivities from the mesh
    :vartype mesh_elements: :class:`numpy.ndarray`
    :cvar elements: List of finite element objects describing the cross-section mesh
    :vartype elements: list[:class:`~sectionproperties.analysis.fea.Tri6`]
    :cvar f_top: Current force in the top region
    :type f_top: float
    :cvar c_top: Centroid of the force in the top region *(c_top_x, c_top_y)*
    :type c_top: list[float, float]
    :cvar c_bot: Centroid of the force in the bottom region *(c_bot_x, c_bot_y)*
    :type c_bot: list[float, float]
    """

    def __init__(
        self, geom: Union[section_geometry.Geometry, section_geometry.CompoundGeometry]
    ):
        """Inits the PlasticSection class."""
        self.geometry = geom.align_center()
        self.geometry.compile_geometry()

        # initialize variables to be defined later within calculate_plastic_force
        self.c_top = [0.0, 0.0]
        self.c_bot = [0.0, 0.0]
        self.f_top = 0.0
        self.f_bot = 0.0

    def calculate_centroid(self, elements):
        """Calculates the elastic centroid from a list of finite elements.

        :param elements: A list of Tri6 finite elements.
        :type elements: list[:class:`~sectionproperties.analysis.fea.Tri6`]
        :return: A tuple containing the x and y location of the elastic centroid.
        :rtype: tuple(float, float)
        """
        ea = 0
        qx = 0
        qy = 0

        # General case
        # loop through all the geometries
        if isinstance(self.geometry, section_geometry.CompoundGeometry):
            for geom in self.geometry.geoms:
                e = geom.material.elastic_modulus
                area = geom.calculate_area()
                cx, cy = geom.calculate_centroid()

                ea += area * e
                qx += cx * e
                qy += cy * e

        # Special case (just one geometry w/ one material)
        else:
            e = self.geometry.material.elastic_modulus
            area = self.geometry.calculate_area()
            cx, cy = self.geometry.calculate_centroid()

            ea += area * e
            qx += cx * area * e
            qy += cy * area * e

        return (qy / ea, qx, ea)

    def calculate_plastic_properties(self, section, verbose):
        """Calculates the location of the plastic centroid with respect to the centroidal and
        principal bending axes, the plastic section moduli and shape factors and stores the results
        to the supplied :class:`~sectionproperties.analysis.section.Section` object.

        :param section: Cross section object that uses the same geometry and materials
            specified in the class constructor
        :type section: :class:`~sectionproperties.analysis.section.Section`
        :param bool verbose: If set to True, the number of iterations required for each plastic
            axis is printed to the terminal.
        """
        # 1) Calculate plastic properties for centroidal axis
        # calculate distances to the extreme fibres
        fibres = self.calculate_extreme_fibres(
            0
        )  # Angle = 0 ; returns bounding box Provides upper/lower bound for brentq algo

        # 1a) Calculate x-axis plastic centroid
        (y_pc, r, f, c_top, c_bot) = self.pc_algorithm(
            np.array([1, 0]), fibres[2:], 1, verbose
        )  # fibres[2:] = ymin, ymax

        self.check_convergence(r, "x-axis")
        section.section_props.y_pc = y_pc
        section.section_props.sxx = f * abs(c_top[1] - c_bot[1])

        if verbose:
            self.print_verbose(y_pc, r, "x-axis")  # Location of axis for each iteration

        # 1b) Calculate y-axis plastic centroid
        (x_pc, r, f, c_top, c_bot) = self.pc_algorithm(
            np.array([0, 1]), fibres[0:2], 2, verbose
        )  # fibres[0:2] = xmin, xmax

        self.check_convergence(r, "y-axis")
        section.section_props.x_pc = x_pc
        section.section_props.syy = f * abs(c_top[0] - c_bot[0])

        if verbose:
            self.print_verbose(x_pc, r, "y-axis")  # Location of axis for each iteration

        # 2) Calculate plastic properties for principal axis
        # convert principal axis angle to radians
        angle = section.section_props.phi * np.pi / 180

        # unit vectors in the axis directions
        ux = np.array([np.cos(angle), np.sin(angle)])  # Unit vectors
        uy = np.array([-np.sin(angle), np.cos(angle)])  # Unit vectors

        # calculate distances to the extreme fibres in the principal axis
        fibres = self.calculate_extreme_fibres(section.section_props.phi)

        # 2a) Calculate 11-axis plastic centroid
        (y22_pc, r, f, c_top, c_bot) = self.pc_algorithm(
            ux, fibres[2:], 1, verbose
        )  # ux

        # calculate the centroids in the principal coordinate system
        c_top_p = fea.principal_coordinate(
            section.section_props.phi, c_top[0], c_top[1]
        )
        c_bot_p = fea.principal_coordinate(
            section.section_props.phi, c_bot[0], c_bot[1]
        )

        self.check_convergence(r, "11-axis")
        section.section_props.y22_pc = y22_pc
        section.section_props.s11 = f * abs(c_top_p[1] - c_bot_p[1])

        if verbose:
            self.print_verbose(y22_pc, r, "11-axis")

        # 2b) Calculate 22-axis plastic centroid
        (x11_pc, r, f, c_top, c_bot) = self.pc_algorithm(
            uy, fibres[0:2], 2, verbose
        )  # uy

        # calculate the centroids in the principal coordinate system
        c_top_p = fea.principal_coordinate(
            section.section_props.phi, c_top[0], c_top[1]
        )
        c_bot_p = fea.principal_coordinate(
            section.section_props.phi, c_bot[0], c_bot[1]
        )

        self.check_convergence(r, "22-axis")
        section.section_props.x11_pc = x11_pc
        section.section_props.s22 = f * abs(c_top_p[0] - c_bot_p[0])

        if verbose:
            self.print_verbose(x11_pc, r, "22-axis")

        # if there are no materials specified, calculate shape factors
        if list(set(section.materials)) == [pre.DEFAULT_MATERIAL]:
            section.section_props.sf_xx_plus = (
                section.section_props.sxx / section.section_props.zxx_plus
            )
            section.section_props.sf_xx_minus = (
                section.section_props.sxx / section.section_props.zxx_minus
            )
            section.section_props.sf_yy_plus = (
                section.section_props.syy / section.section_props.zyy_plus
            )
            section.section_props.sf_yy_minus = (
                section.section_props.syy / section.section_props.zyy_minus
            )

            section.section_props.sf_11_plus = (
                section.section_props.s11 / section.section_props.z11_plus
            )
            section.section_props.sf_11_minus = (
                section.section_props.s11 / section.section_props.z11_minus
            )
            section.section_props.sf_22_plus = (
                section.section_props.s22 / section.section_props.z22_plus
            )
            section.section_props.sf_22_minus = (
                section.section_props.s22 / section.section_props.z22_minus
            )

    def check_convergence(self, root_result, axis):
        """Checks that the function solver converged and if not, raises a helpful error.

        :param root_result: Result object from the root finder
        :type root_result: :class:`scipy.optimize.RootResults`
        :param string axis: Axis being considered by the function solver
        :raises RuntimeError: If the function solver did not converge
        """

        if not root_result.converged:
            msg = "Plastic centroid calculation about the {0}".format(axis)
            msg += " failed. Contact robbie.vanleeuwen@gmail.com with your"
            msg += " analysis parameters. Termination flag: {0}".format(
                root_result.flag
            )

            raise RuntimeError(msg)

    def print_verbose(self, d, root_result, axis):
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
        for (i, pt) in enumerate(self.geometry.points):
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

        p = np.array(
            [d * u_p[0], d * u_p[1]]
        )  # p finding a point on the axis by scaling the perpendicular

        # calculate force equilibrium
        (f_top, f_bot) = self.calculate_plastic_force(u, p)

        # calculate the force norm
        f_norm = (f_top - f_bot) / (f_top + f_bot)  # Going down to zero

        # print verbose results
        if verbose:
            print("d = {0}; f_norm = {1}".format(d, f_norm))

        # return the force norm
        return f_norm  # This is the result that is the target for the root finding algorithm

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

        a = dlim[0]  # Upper and lower bound for algorithm
        b = dlim[1]  # Upper and lower bound for algorithm

        (
            d,
            r,
        ) = brentq(  # d = neutral axis depth, measured from the origin; r is the convegence indicator
            self.evaluate_force_eq,
            a,
            b,
            args=(u, u_p, verbose),
            full_output=True,
            disp=False,  # Unit vector w/ unit vector perpendicular
            xtol=1e-6,
            rtol=1e-6,
        )
        return (d, r, self.f_top, self.c_top, self.c_bot)

    def calculate_plastic_force(
        self, u: np.ndarray, p: np.ndarray
    ) -> Tuple[float, float]:
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

        top_geoms, bot_geoms = self.geometry.split_section(point_i=p, vector=u)

        if top_geoms:
            for top_geom in top_geoms:
                e = top_geom.material.elastic_modulus
                f_y = top_geom.material.yield_strength
                area_top = top_geom.calculate_area()
                ea_top += e * area_top
                cx, cy = top_geom.calculate_centroid()
                qx_top += cy * area_top
                qy_top += cx * area_top
                f_top += f_y * area_top

        if bot_geoms:
            for bot_geom in bot_geoms:
                e = bot_geom.material.elastic_modulus
                f_y = bot_geom.material.yield_strength
                area_bot = bot_geom.calculate_area()
                ea_bot += e * area_bot
                cx, cy = bot_geom.calculate_centroid()
                qx_bot += cy * area_bot
                qy_bot += cx * area_bot
                f_bot += f_y * area_bot

        try:
            self.c_top = [qy_top / ea_top, qx_top / ea_top]
            self.f_top = f_top
        except ZeroDivisionError:
            self.c_top = [0, 0]
            self.f_top = 0

        try:
            self.c_bot = [qy_bot / ea_bot, qx_bot / ea_bot]
        except ZeroDivisionError:
            self.c_bot = [0, 0]

        return (f_top, f_bot)


class StressPost:
    """Class for post-processing finite element stress results.

    A StressPost object is created when a stress analysis is carried out and is returned as an
    object to allow post-processing of the results. The StressPost object creates a deep copy of
    the MaterialGroups within the cross-section to allow the calculation of stresses for each
    material. Methods for post-processing the calculated stresses are provided.

    :param section: Cross section object for stress calculation
    :type section: :class:`~sectionproperties.analysis.section.Section`

    :cvar section: Cross section object for stress calculation
    :vartype section: :class:`~sectionproperties.analysis.section.Section`
    :cvar material_groups: A deep copy of the `section` material groups to allow a new stress
        analysis
    :vartype material_groups: list[:class:`~sectionproperties.pre.pre.MaterialGroup`]
    """

    def __init__(self, section):
        """Inits the StressPost class."""

        self.section = section

        # make a deep copy of the material groups to the StressPost object such that stress results
        # can be saved to a new material group
        self.material_groups = copy.deepcopy(section.material_groups)

    def plot_stress_contour(self, sigs, title, cmap, normalize, **kwargs):
        """Plots filled stress contours over the finite element mesh.

        :param sigs: List of nodal stress values for each material
        :type sigs: list[:class:`numpy.ndarray`]
        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axe object
        :rtype: :class:`matplotlib.axes`
        """

        # create plot and setup the plot
        with post.plotting_context(title=title, **kwargs) as (fig, ax):
            # set up the colormap
            cmap = cm.get_cmap(name=cmap)

            # create triangulation
            triang = tri.Triangulation(
                self.section.mesh_nodes[:, 0],
                self.section.mesh_nodes[:, 1],
                self.section.mesh_elements[:, 0:3],
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

            norm = None
            if normalize:
                norm = CenteredNorm()

            # plot the filled contour, looping through the materials
            for (i, sig) in enumerate(sigs):
                # create and set the mask for the current material
                mask_array = np.ones(len(self.section.elements), dtype=bool)
                mask_array[self.material_groups[i].el_ids] = False
                triang.set_mask(mask_array)

                # plot the filled contour
                trictr = ax.tricontourf(triang, sig, v, cmap=cmap, norm=norm)

            # display the colourbar
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.1)

            fig.colorbar(trictr, label="Stress", format="%.4e", ticks=ticks, cax=cax)

            # plot the finite element mesh
            self.section.plot_mesh(materials=False, **dict(kwargs, ax=ax))

        return ax

    def plot_stress_vector(self, sigxs, sigys, title, cmap, normalize, **kwargs):
        r"""Plots stress vectors over the finite element mesh.

        :param sigxs: List of x-components of the nodal stress values for each material
        :type sigxs: list[:class:`numpy.ndarray`]
        :param sigys: List of y-components of the nodal stress values for each material
        :type sigys: list[:class:`numpy.ndarray`]
        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`
        """

        # create plot and setup the plot
        with post.plotting_context(title=title, **kwargs) as (fig, ax):
            # set up the colormap
            cmap = cm.get_cmap(name=cmap)

            # initialise quiver plot list max scale
            quiv_list = []
            max_scale = 0

            norm = None
            if normalize:
                norm = CenteredNorm()

            # plot the vectors
            for (i, sigx) in enumerate(sigxs):
                sigy = sigys[i]

                # scale the colour with respect to the magnitude of the vector
                c = np.hypot(sigx, sigy)

                quiv = ax.quiver(
                    self.section.mesh_nodes[:, 0],
                    self.section.mesh_nodes[:, 1],
                    sigx,
                    sigy,
                    c,
                    cmap=cmap,
                    norm=norm,
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
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.1)

            fig.colorbar(quiv, label="Stress", format="%.4e", ticks=v1, cax=cax)

            # plot the finite element mesh
            self.section.plot_mesh(materials=False, **dict(kwargs, ax=ax))

        return ax

    def get_stress(self):
        r"""Returns the stresses within each material belonging to the current
        :class:`~sectionproperties.analysis.section.StressPost` object.

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
        * *'sig_1'*: Major principal stress :math:`\sigma_{1}` resulting from all actions
        * *'sig_3'*: Minor principal stress :math:`\sigma_{3}` resulting from all actions
        * *'sig_vm'*: von Mises stress :math:`\sigma_{vM}` resulting from all actions

        The following example returns stresses for each material within a composite section, note
        that a result is generated for each node in the mesh for all materials irrespective of
        whether the materials exists at that point or not.

        ::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(
                N=50e3, Mxx=-5e6, M22=2.5e6, Mzz=0.5e6, Vx=10e3, Vy=5e3
            )
            stresses = stress_post.get_stress()

            print("Number of nodes: {0}".format(section.num_nodes))

            for stress in stresses:
                print('Material: {0}'.format(stress['Material']))
                print('List Size: {0}'.format(len(stress['sig_zz_n'])))
                print('Normal Stresses: {0}'.format(stress['sig_zz_n']))
                print('von Mises Stresses: {0}'.format(stress['sig_vm']))

        ::

            $ Number of nodes: 2465

            $ Material: Timber
            $ List Size: 2465
            $ Normal Stresses: [0.76923077 0.76923077 0.76923077 ... 0.76923077 0.76923077 0.76923077]
            $ von Mises Stresses: [7.6394625  5.38571866 3.84784964 ... 3.09532948 3.66992556 2.81976647]

            $ Material: Steel
            $ List Size: 2465
            $ Normal Stresses: [19.23076923 0. 0. ... 0. 0. 0.]
            $ von Mises Stresses: [134.78886419 0. 0. ... 0. 0. 0.]
        """

        stress = []

        for group in self.material_groups:
            stress.append(
                {
                    "Material": group.material.name,
                    "sig_zz_n": group.stress_result.sig_zz_n,
                    "sig_zz_mxx": group.stress_result.sig_zz_mxx,
                    "sig_zz_myy": group.stress_result.sig_zz_myy,
                    "sig_zz_m11": group.stress_result.sig_zz_m11,
                    "sig_zz_m22": group.stress_result.sig_zz_m22,
                    "sig_zz_m": group.stress_result.sig_zz_m,
                    "sig_zx_mzz": group.stress_result.sig_zx_mzz,
                    "sig_zy_mzz": group.stress_result.sig_zy_mzz,
                    "sig_zxy_mzz": group.stress_result.sig_zxy_mzz,
                    "sig_zx_vx": group.stress_result.sig_zx_vx,
                    "sig_zy_vx": group.stress_result.sig_zy_vx,
                    "sig_zxy_vx": group.stress_result.sig_zxy_vx,
                    "sig_zx_vy": group.stress_result.sig_zx_vy,
                    "sig_zy_vy": group.stress_result.sig_zy_vy,
                    "sig_zxy_vy": group.stress_result.sig_zxy_vy,
                    "sig_zx_v": group.stress_result.sig_zx_v,
                    "sig_zy_v": group.stress_result.sig_zy_v,
                    "sig_zxy_v": group.stress_result.sig_zxy_v,
                    "sig_zz": group.stress_result.sig_zz,
                    "sig_zx": group.stress_result.sig_zx,
                    "sig_zy": group.stress_result.sig_zy,
                    "sig_zxy": group.stress_result.sig_zxy,
                    "sig_1": group.stress_result.sig_1,
                    "sig_3": group.stress_result.sig_3,
                    "sig_vm": group.stress_result.sig_vm,
                }
            )

        return stress

    def plot_stress_n_zz(
        self,
        title=r"Stress Contour Plot - $\sigma_{zz,N}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the normal stress :math:`\sigma_{zz,N}` resulting from the
        axial load :math:`N`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots the normal stress within a 150x90x12 UA section resulting from
        an axial force of 10 kN::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(N=10e3)

            stress_post.plot_stress_n_zz()

        ..  figure:: ../images/stress/stress_n_zz.png
            :align: center
            :scale: 75 %

            Contour plot of the axial stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zz_n)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_stress_mxx_zz(
        self,
        title=r"Stress Contour Plot - $\sigma_{zz,Mxx}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the normal stress :math:`\sigma_{zz,Mxx}` resulting from the
        bending moment :math:`M_{xx}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots the normal stress within a 150x90x12 UA section resulting from
        a bending moment about the x-axis of 5 kN.m::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mxx=5e6)

            stress_post.plot_stress_mxx_zz()

        ..  figure:: ../images/stress/stress_mxx_zz.png
            :align: center
            :scale: 75 %

            Contour plot of the bending stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zz_mxx)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_stress_myy_zz(
        self,
        title=r"Stress Contour Plot - $\sigma_{zz,Myy}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the normal stress :math:`\sigma_{zz,Myy}` resulting from the
        bending moment :math:`M_{yy}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots the normal stress within a 150x90x12 UA section resulting from
        a bending moment about the y-axis of 2 kN.m::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Myy=2e6)

            stress_post.plot_stress_myy_zz()

        ..  figure:: ../images/stress/stress_myy_zz.png
            :align: center
            :scale: 75 %

            Contour plot of the bending stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zz_myy)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_stress_m11_zz(
        self,
        title=r"Stress Contour Plot - $\sigma_{zz,M11}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the normal stress :math:`\sigma_{zz,M11}` resulting from the
        bending moment :math:`M_{11}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots the normal stress within a 150x90x12 UA section resulting from
        a bending moment about the 11-axis of 5 kN.m::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(M11=5e6)

            stress_post.plot_stress_m11_zz()

        ..  figure:: ../images/stress/stress_m11_zz.png
            :align: center
            :scale: 75 %

            Contour plot of the bending stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zz_m11)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_stress_m22_zz(
        self,
        title=r"Stress Contour Plot - $\sigma_{zz,M22}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the normal stress :math:`\sigma_{zz,M22}` resulting from the
        bending moment :math:`M_{22}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots the normal stress within a 150x90x12 UA section resulting from
        a bending moment about the 22-axis of 2 kN.m::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(M22=5e6)

            stress_post.plot_stress_m22_zz()

        ..  figure:: ../images/stress/stress_m22_zz.png
            :align: center
            :scale: 75 %

            Contour plot of the bending stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zz_m22)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_stress_m_zz(
        self,
        title=r"Stress Contour Plot - $\sigma_{zz,\Sigma M}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the normal stress :math:`\sigma_{zz,\Sigma M}` resulting from
        all bending moments :math:`M_{xx} + M_{yy} + M_{11} + M_{22}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots the normal stress within a 150x90x12 UA section resulting from
        a bending moment about the x-axis of 5 kN.m, a bending moment about the y-axis of 2 kN.m
        and a bending moment of 3 kN.m about the 11-axis::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mxx=5e6, Myy=2e6, M11=3e6)

            stress_post.plot_stress_m_zz()

        ..  figure:: ../images/stress/stress_m_zz.png
            :align: center
            :scale: 75 %

            Contour plot of the bending stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zz_m)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_stress_mzz_zx(
        self,
        title=r"Stress Contour Plot - $\sigma_{zx,Mzz}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the *x*-component of the shear stress :math:`\sigma_{zx,Mzz}`
        resulting from the torsion moment :math:`M_{zz}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots the x-component of the shear stress within a 150x90x12 UA
        section resulting from a torsion moment of 1 kN.m::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mzz=1e6)

            stress_post.plot_stress_mzz_zx()

        ..  figure:: ../images/stress/stress_mzz_zx.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zx_mzz)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_stress_mzz_zy(
        self,
        title=r"Stress Contour Plot - $\sigma_{zy,Mzz}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the *y*-component of the shear stress :math:`\sigma_{zy,Mzz}`
        resulting from the torsion moment :math:`M_{zz}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots the y-component of the shear stress within a 150x90x12 UA
        section resulting from a torsion moment of 1 kN.m::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mzz=1e6)

            stress_post.plot_stress_mzz_zy()

        ..  figure:: ../images/stress/stress_mzz_zy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zy_mzz)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_stress_mzz_zxy(
        self,
        title=r"Stress Contour Plot - $\sigma_{zxy,Mzz}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the resultant shear stress :math:`\sigma_{zxy,Mzz}` resulting
        from the torsion moment :math:`M_{zz}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots a contour of the resultant shear stress within a 150x90x12 UA
        section resulting from a torsion moment of 1 kN.m::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mzz=1e6)

            stress_post.plot_stress_mzz_zxy()

        ..  figure:: ../images/stress/stress_mzz_zxy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zxy_mzz)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_vector_mzz_zxy(
        self,
        title=r"Stress Vector Plot - $\sigma_{zxy,Mzz}$",
        cmap="YlOrBr",
        normalize=False,
        **kwargs,
    ):
        r"""Produces a vector plot of the resultant shear stress :math:`\sigma_{zxy,Mzz}` resulting
        from the torsion moment :math:`M_{zz}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example generates a vector plot of the shear stress within a 150x90x12 UA
        section resulting from a torsion moment of 1 kN.m::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mzz=1e6)

            stress_post.plot_vector_mzz_zxy()

        ..  figure:: ../images/stress/vector_mzz_zxy.png
            :align: center
            :scale: 75 %

            Vector plot of the shear stress.
        """

        sigxs = []
        sigys = []

        for group in self.material_groups:
            sigxs.append(group.stress_result.sig_zx_mzz)
            sigys.append(group.stress_result.sig_zy_mzz)

        return self.plot_stress_vector(sigxs, sigys, title, cmap, normalize, **kwargs)

    def plot_stress_vx_zx(
        self,
        title=r"Stress Contour Plot - $\sigma_{zx,Vx}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the *x*-component of the shear stress :math:`\sigma_{zx,Vx}`
        resulting from the shear force :math:`V_{x}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots the x-component of the shear stress within a 150x90x12 UA
        section resulting from a shear force in the x-direction of 15 kN::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vx=15e3)

            stress_post.plot_stress_vx_zx()

        ..  figure:: ../images/stress/stress_vx_zx.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zx_vx)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_stress_vx_zy(
        self,
        title=r"Stress Contour Plot - $\sigma_{zy,Vx}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the *y*-component of the shear stress :math:`\sigma_{zy,Vx}`
        resulting from the shear force :math:`V_{x}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots the y-component of the shear stress within a 150x90x12 UA
        section resulting from a shear force in the x-direction of 15 kN::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vx=15e3)

            stress_post.plot_stress_vx_zy()

        ..  figure:: ../images/stress/stress_vx_zy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zy_vx)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_stress_vx_zxy(
        self,
        title=r"Stress Contour Plot - $\sigma_{zz,Myy}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the resultant shear stress :math:`\sigma_{zxy,Vx}` resulting
        from the shear force :math:`V_{x}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots a contour of the resultant shear stress within a 150x90x12 UA
        section resulting from a shear force in the x-direction of 15 kN::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vx=15e3)

            stress_post.plot_stress_vx_zxy()

        ..  figure:: ../images/stress/stress_vx_zxy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        title = r"Stress Contour Plot - $\sigma_{zxy,Vx}$"
        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zxy_vx)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_vector_vx_zxy(
        self,
        title=r"Stress Vector Plot - $\sigma_{zxy,Vx}$",
        cmap="YlOrBr",
        normalize=False,
        **kwargs,
    ):
        r"""Produces a vector plot of the resultant shear stress :math:`\sigma_{zxy,Vx}` resulting
        from the shear force :math:`V_{x}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example generates a vector plot of the shear stress within a 150x90x12 UA
        section resulting from a shear force in the x-direction of 15 kN::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vx=15e3)

            stress_post.plot_vector_vx_zxy()

        ..  figure:: ../images/stress/vector_vx_zxy.png
            :align: center
            :scale: 75 %

            Vector plot of the shear stress.
        """

        sigxs = []
        sigys = []

        for group in self.material_groups:
            sigxs.append(group.stress_result.sig_zx_vx)
            sigys.append(group.stress_result.sig_zy_vx)

        return self.plot_stress_vector(sigxs, sigys, title, cmap, normalize, **kwargs)

    def plot_stress_vy_zx(
        self,
        title=r"Stress Contour Plot - $\sigma_{zx,Vy}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the *x*-component of the shear stress :math:`\sigma_{zx,Vy}`
        resulting from the shear force :math:`V_{y}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots the x-component of the shear stress within a 150x90x12 UA
        section resulting from a shear force in the y-direction of 30 kN::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vy=30e3)

            stress_post.plot_stress_vy_zx()

        ..  figure:: ../images/stress/stress_vy_zx.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zx_vy)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_stress_vy_zy(
        self,
        title=r"Stress Contour Plot - $\sigma_{zy,Vy}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the *y*-component of the shear stress :math:`\sigma_{zy,Vy}`
        resulting from the shear force :math:`V_{y}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots the y-component of the shear stress within a 150x90x12 UA
        section resulting from a shear force in the y-direction of 30 kN::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vy=30e3)

            stress_post.plot_stress_vy_zy()

        ..  figure:: ../images/stress/stress_vy_zy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zy_vy)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_stress_vy_zxy(
        self,
        title=r"Stress Contour Plot - $\sigma_{zxy,Vy}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the resultant shear stress :math:`\sigma_{zxy,Vy}` resulting
        from the shear force :math:`V_{y}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots a contour of the resultant shear stress within a 150x90x12 UA
        section resulting from a shear force in the y-direction of 30 kN::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vy=30e3)

            stress_post.plot_stress_vy_zxy()

        ..  figure:: ../images/stress/stress_vy_zxy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zxy_vy)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_vector_vy_zxy(
        self,
        title=r"Stress Vector Plot - $\sigma_{zxy,Vy}$",
        cmap="YlOrBr",
        normalize=False,
        **kwargs,
    ):
        r"""Produces a vector plot of the resultant shear stress :math:`\sigma_{zxy,Vy}` resulting
        from the shear force :math:`V_{y}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example generates a vector plot of the shear stress within a 150x90x12 UA
        section resulting from a shear force in the y-direction of 30 kN::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vy=30e3)

            stress_post.plot_vector_vy_zxy()

        ..  figure:: ../images/stress/vector_vy_zxy.png
            :align: center
            :scale: 75 %

            Vector plot of the shear stress.
        """

        sigxs = []
        sigys = []

        for group in self.material_groups:
            sigxs.append(group.stress_result.sig_zx_vy)
            sigys.append(group.stress_result.sig_zy_vy)

        return self.plot_stress_vector(sigxs, sigys, title, cmap, normalize, **kwargs)

    def plot_stress_v_zx(
        self,
        title=r"Stress Contour Plot - $\sigma_{zx,\Sigma V}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the *x*-component of the shear stress
        :math:`\sigma_{zx,\Sigma V}` resulting from the sum of the applied shear forces
        :math:`V_{x} + V_{y}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots the x-component of the shear stress within a 150x90x12 UA
        section resulting from a shear force of 15 kN in the x-direction and 30 kN in the
        y-direction::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vx=15e3, Vy=30e3)

            stress_post.plot_stress_v_zx()

        ..  figure:: ../images/stress/stress_v_zx.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zx_v)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_stress_v_zy(
        self,
        title=r"Stress Contour Plot - $\sigma_{zy,\Sigma V}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the *y*-component of the shear stress
        :math:`\sigma_{zy,\Sigma V}` resulting from the sum of the applied shear forces
        :math:`V_{x} + V_{y}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots the y-component of the shear stress within a 150x90x12 UA
        section resulting from a shear force of 15 kN in the x-direction and 30 kN in the
        y-direction::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vx=15e3, Vy=30e3)

            stress_post.plot_stress_v_zy()

        ..  figure:: ../images/stress/stress_v_zy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zy_v)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_stress_v_zxy(
        self,
        title=r"Stress Contour Plot - $\sigma_{zxy,\Sigma V}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the resultant shear stress
        :math:`\sigma_{zxy,\Sigma V}` resulting from the sum of the applied shear forces
        :math:`V_{x} + V_{y}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots a contour of the resultant shear stress within a 150x90x12 UA
        section resulting from a shear force of 15 kN in the x-direction and 30 kN in the
        y-direction::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vx=15e3, Vy=30e3)

            stress_post.plot_stress_v_zxy()

        ..  figure:: ../images/stress/stress_v_zxy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zxy_v)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_vector_v_zxy(
        self,
        title=r"Stress Vector Plot - $\sigma_{zxy,\Sigma V}$",
        cmap="YlOrBr",
        normalize=False,
        **kwargs,
    ):
        r"""Produces a vector plot of the resultant shear stress
        :math:`\sigma_{zxy,\Sigma V}` resulting from the sum of the  applied shear forces
        :math:`V_{x} + V_{y}`.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example generates a vector plot of the shear stress within a 150x90x12 UA
        section resulting from a shear force of 15 kN in the x-direction and 30 kN in the
        y-direction::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Vx=15e3, Vy=30e3)

            stress_post.plot_vector_v_zxy()

        ..  figure:: ../images/stress/vector_v_zxy.png
            :align: center
            :scale: 75 %

            Vector plot of the shear stress.
        """

        sigxs = []
        sigys = []

        for group in self.material_groups:
            sigxs.append(group.stress_result.sig_zx_v)
            sigys.append(group.stress_result.sig_zy_v)

        return self.plot_stress_vector(sigxs, sigys, title, cmap, normalize, **kwargs)

    def plot_stress_zz(
        self,
        title=r"Stress Contour Plot - $\sigma_{zz}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the combined normal stress :math:`\sigma_{zz}` resulting from
        all actions.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots the normal stress within a 150x90x12 UA section resulting from
        an axial force of 100 kN, a bending moment about the x-axis of 5 kN.m and a bending moment
        about the y-axis of 2 kN.m::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(N=100e3, Mxx=5e6, Myy=2e6)

            stress_post.plot_stress_zz()

        ..  figure:: ../images/stress/stress_zz.png
            :align: center
            :scale: 75 %

            Contour plot of the normal stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zz)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_stress_zx(
        self,
        title=r"Stress Contour Plot - $\sigma_{zx}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the *x*-component of the shear stress :math:`\sigma_{zx}`
        resulting from all actions.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots the x-component of the shear stress within a 150x90x12 UA
        section resulting from a torsion moment of 1 kN.m and a shear force of 30 kN in the
        y-direction::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mzz=1e6, Vy=30e3)

            stress_post.plot_stress_zx()

        ..  figure:: ../images/stress/stress_zx.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zx)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_stress_zy(
        self,
        title=r"Stress Contour Plot - $\sigma_{zy}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the *y*-component of the shear stress :math:`\sigma_{zy}`
        resulting from all actions.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots the y-component of the shear stress within a 150x90x12 UA
        section resulting from a torsion moment of 1 kN.m and a shear force of 30 kN in the
        y-direction::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mzz=1e6, Vy=30e3)

            stress_post.plot_stress_zy()

        ..  figure:: ../images/stress/stress_zy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zy)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_stress_zxy(
        self,
        title=r"Stress Contour Plot - $\sigma_{zxy}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the resultant shear stress :math:`\sigma_{zxy}` resulting
        from all actions.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots a contour of the resultant shear stress within a 150x90x12 UA
        section resulting from a torsion moment of 1 kN.m and a shear force of 30 kN in the
        y-direction::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mzz=1e6, Vy=30e3)

            stress_post.plot_stress_zxy()

        ..  figure:: ../images/stress/stress_zxy.png
            :align: center
            :scale: 75 %

            Contour plot of the shear stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_zxy)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_vector_zxy(
        self,
        title=r"Stress Vector Plot - $\sigma_{zxy}$",
        cmap="YlOrBr",
        normalize=False,
        **kwargs,
    ):
        r"""Produces a vector plot of the resultant shear stress :math:`\sigma_{zxy}` resulting
        from all actions.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example generates a vector plot of the shear stress within a 150x90x12 UA
        section resulting from a torsion moment of 1 kN.m and a shear force of 30 kN in the
        y-direction::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(Mzz=1e6, Vy=30e3)

            stress_post.plot_vector_zxy()

        ..  figure:: ../images/stress/vector_zxy.png
            :align: center
            :scale: 75 %

            Vector plot of the shear stress.
        """

        sigxs = []
        sigys = []

        for group in self.material_groups:
            sigxs.append(group.stress_result.sig_zx)
            sigys.append(group.stress_result.sig_zy)

        return self.plot_stress_vector(sigxs, sigys, title, cmap, normalize, **kwargs)

    def plot_stress_1(
        self,
        title=r"Stress Contour Plot - $\sigma_{1}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the major principal stress :math:`\sigma_{1}` resulting from all
        actions.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots a contour of the major principal stress within a 150x90x12 UA
        section resulting from the following actions:

        * :math:`N = 50` kN
        * :math:`M_{xx} = -5` kN.m
        * :math:`M_{22} = 2.5` kN.m
        * :math:`M_{zz} = 1.5` kN.m
        * :math:`V_{x} = 10` kN
        * :math:`V_{y} = 5` kN

        ::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(
                N=50e3, Mxx=-5e6, M22=2.5e6, Mzz=0.5e6, Vx=10e3, Vy=5e3
            )

            stress_post.plot_stress_1()

        ..  figure:: ../images/stress/stress_1.png
            :align: center
            :scale: 50 %

            Contour plot of the major principal stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_1)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_stress_3(
        self,
        title=r"Stress Contour Plot - $\sigma_{3}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        r"""Produces a contour plot of the Minor principal stress :math:`\sigma_{3}` resulting from all
        actions.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots a contour of the Minor principal stress within a 150x90x12 UA
        section resulting from the following actions:

        * :math:`N = 50` kN
        * :math:`M_{xx} = -5` kN.m
        * :math:`M_{22} = 2.5` kN.m
        * :math:`M_{zz} = 1.5` kN.m
        * :math:`V_{x} = 10` kN
        * :math:`V_{y} = 5` kN

        ::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = CrossSection(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(
                N=50e3, Mxx=-5e6, M22=2.5e6, Mzz=0.5e6, Vx=10e3, Vy=5e3
            )

            stress_post.plot_stress_3()

        ..  figure:: ../images/stress/stress_3.png
            :align: center
            :scale: 50 %

            Contour plot of the minor principal stress.
        """

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_3)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_stress_vm(
        self,
        title=r"Stress Contour Plot - $\sigma_{vM}$",
        cmap="coolwarm",
        normalize=True,
        **kwargs,
    ):
        """Produces a contour plot of the von Mises stress :math:`\sigma_{vM}` resulting from all
        actions.

        :param string title: Plot title
        :param string cmap: Matplotlib color map.
        :param bool normalize: If set to true, the CenteredNorm is used to scale the colormap.
            If set to false, the default linear scaling is used.
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots a contour of the von Mises stress within a 150x90x12 UA section
        resulting from the following actions:

        * :math:`N = 50` kN
        * :math:`M_{xx} = -5` kN.m
        * :math:`M_{22} = 2.5` kN.m
        * :math:`M_{zz} = 1.5` kN.m
        * :math:`V_{x} = 10` kN
        * :math:`V_{y} = 5` kN

        ::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            geometry.create_mesh(mesh_sizes=[20])
            section = Section(geometry)

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

        sigs = []

        for group in self.material_groups:
            sigs.append(group.stress_result.sig_vm)

        return self.plot_stress_contour(sigs, title, cmap, normalize, **kwargs)

    def plot_mohrs_circles(self, x, y, title=None, **kwargs):
        r"""Plots Mohr's Circles of the 3D stress state at position x, y

        :param float x: x-coordinate of the point to draw Mohr's Circle
        :param float y: y-coordinate of the point to draw Mohr's Circle
        :param string title: Plot title
        :param \**kwargs: Passed to :func:`~sectionproperties.post.post.plotting_context`

        :return: Matplotlib axes object
        :rtype: :class:`matplotlib.axes`

        The following example plots the Mohr's Circles for the 3D stress state within a 150x90x12
        UA section resulting from the following actions:

        * :math:`N = 50` kN
        * :math:`M_{xx} = -5` kN.m
        * :math:`M_{22} = 2.5` kN.m
        * :math:`M_{zz} = 1.5` kN.m
        * :math:`V_{x} = 10` kN
        * :math:`V_{y} = 5` kN

        at the point (10, 88.9).

        ::

            import sectionproperties.pre.library.steel_sections as steel_sections
            from sectionproperties.analysis.section import Section

            geometry = steel_sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])
            section = Section(geometry, mesh)

            section.calculate_geometric_properties()
            section.calculate_warping_properties()
            stress_post = section.calculate_stress(
                N=50e3, Mxx=-5e6, M22=2.5e6, Mzz=0.5e6, Vx=10e3, Vy=5e3
            )

            stress_post.plot_mohrs_circles(10, 88.9)

        ..  figure:: ../images/stress/mohrs_circles.png
            :align: center
            :scale: 100 %

            Mohr's Circles of the 3D stress state at (10, 88.9).

        """

        pt = (x, y)
        nodes = self.section.mesh_nodes
        ele = self.section.mesh_elements
        triang = tri.Triangulation(nodes[:, 0], nodes[:, 1], ele[:, 0:3])

        # Find in which material group the point lies
        pt_group = None
        for group in self.material_groups:
            mask_array = np.ones(len(self.section.elements), dtype=bool)
            mask_array[group.el_ids] = False
            triang.set_mask(mask_array)
            trifinder = triang.get_trifinder()
            if trifinder(*pt) != -1:
                pt_group = group
            triang.set_mask(None)

        if pt_group is None:
            raise Exception(f"Point {(*pt,)} is not within mesh")

        # Assesmble the stress results from the relevant material group
        sigma_zz_v = pt_group.stress_result.sig_zz
        tau_xz_v = pt_group.stress_result.sig_zx
        tau_yz_v = pt_group.stress_result.sig_zy

        # Get the interpolators
        sigma_zz_interp = tri.LinearTriInterpolator(triang, sigma_zz_v)
        tau_xz_interp = tri.LinearTriInterpolator(triang, tau_xz_v)
        tau_yz_interp = tri.LinearTriInterpolator(triang, tau_yz_v)

        # Get the stresses at the point
        sigma_zz = sigma_zz_interp(*pt).item()
        tau_xz = tau_xz_interp(*pt).item()
        tau_yz = tau_yz_interp(*pt).item()

        # Assemble the stress tensor
        sigma_xx = 0
        sigma_yy = 0
        tau_xy = 0
        sigma = np.array(
            [
                [sigma_xx, tau_xy, tau_xz],
                [tau_xy, sigma_yy, tau_yz],
                [tau_xz, tau_yz, sigma_zz],
            ]
        )

        # Solve for the principal stresses using the general approach
        s, n = np.linalg.eig(sigma)
        sigma_3, sigma_2, sigma_1 = np.sort(s)
        s_l = s.tolist()
        sigma_1i, sigma_2i, sigma_3i = (
            s_l.index(sigma_1),
            s_l.index(sigma_2),
            s_l.index(sigma_3),
        )
        n1, n2, n3 = n[:, sigma_1i], n[:, sigma_2i], n[:, sigma_3i]

        # The tractions on each plane in cartesian coords wrt principal axes
        n_inv = np.linalg.inv(n)
        tractions = []
        for col in range(3):
            ss = n_inv[:, col].T @ np.diag(s) @ n_inv[:, col]
            ts = np.sqrt(np.linalg.norm(np.diag(s) @ n_inv[:, col]) ** 2 - ss ** 2)
            tractions.append((ss, ts))

        def plot_circle(ax, c, R, col, label=None, fill=None):
            circ = plt.Circle(c, R, fill=fill, ec=col, label=label)
            ax.add_patch(circ)
            ax.set_aspect(1)
            ax.autoscale_view()

        if title is None:
            title = f"Mohr's Circles for 3D Stress State at {(*pt,)}"

        # create plot and setup the plot
        with post.plotting_context(title=title, **kwargs) as (fig, ax):
            plot_circle(
                ax,
                (0.5 * (sigma_2 + sigma_3), 0),
                0.5 * (sigma_2 - sigma_3),
                "r",
                r"C1: ($\sigma_2$, $\sigma_3$)",
            )
            plot_circle(
                ax,
                (0.5 * (sigma_1 + sigma_3), 0),
                0.5 * (sigma_1 - sigma_3),
                "b",
                r"C2: ($\sigma_1$, $\sigma_3$)",
            )
            plot_circle(
                ax,
                (0.5 * (sigma_1 + sigma_2), 0),
                0.5 * (sigma_1 - sigma_2),
                "k",
                r"C3: ($\sigma_1$, $\sigma_2$)",
            )

            for i, plane, col in zip(range(3), ["X", "Y", "Z"], ["r", "b", "k"]):
                ax.plot(*tractions[i], f"{col}.", label=rf"{plane}-face")

            ax.set_axisbelow(True)
            ax.grid(which="both")
            ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))

            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)

            ax.set_ylabel(r"Shear stress $\tau$ (MPa)")
            ax.set_xlabel(r"Direct stress $\sigma$ (MPa)")

            ax.xaxis.set_tick_params(bottom=True, top=False, direction="inout")
            ax.yaxis.set_tick_params(left=True, right=False, direction="inout")

            #
            # The following is just to get the labels positioned outside the axes
            #
            # Store default label positions
            x_lbl_pos = ax.xaxis.label.get_position()
            y_lbl_pos = ax.yaxis.label.get_position()
            # Make spines pass through zero of the other axis
            ax.spines["bottom"].set_position("zero")
            ax.spines["left"].set_position("zero")
            # Now set the coords
            ax.xaxis.set_label_coords(*x_lbl_pos)
            ax.yaxis.set_label_coords(*y_lbl_pos)

        return ax


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
    :vartype stress_result: :class:`~sectionproperties.analysis.section.StressResult`
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


@dataclass
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
    :cvar sig_1: Major principal stress (:math:`\sigma_{1}`) resulting from all actions
    :vartype sig_1: :class:`numpy.ndarray`
    :cvar sig_3: Minor principal stress (:math:`\sigma_{3}`) resulting from all actions
    :vartype sig_3: :class:`numpy.ndarray`
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
        self.sig_1 = np.zeros(num_nodes)
        self.sig_3 = np.zeros(num_nodes)
        self.sig_vm = np.zeros(num_nodes)

    def calculate_combined_stresses(self):
        """Calculates the combined cross-section stresses."""

        self.sig_zz_m = (
            self.sig_zz_mxx + self.sig_zz_myy + self.sig_zz_m11 + self.sig_zz_m22
        )
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
        self.sig_1 = self.sig_zz / 2 + np.sqrt(
            (self.sig_zz / 2) ** 2 + self.sig_zxy ** 2
        )
        self.sig_3 = self.sig_zz / 2 - np.sqrt(
            (self.sig_zz / 2) ** 2 + self.sig_zxy ** 2
        )
        self.sig_vm = (self.sig_zz ** 2 + 3 * self.sig_zxy ** 2) ** 0.5


@dataclass
class SectionProperties:
    """Class for storing section properties.

    Stores calculated section properties. Also provides methods to calculate section properties
    entirely derived from other section properties.

    :cvar float area: Cross-sectional area
    :cvar float perimeter: Cross-sectional perimeter
    :cvar float mass: Cross-sectional mass
    :cvar float ea: Modulus weighted area (axial rigidity)
    :cvar float ga: Modulus weighted product of shear modulus and area
    :cvar float nu_eff: Effective Poisson's ratio
    :cvar float e_eff: Effective elastic modulus
    :cvar float g_eff: Effective shear modulus
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

    area: Optional[float] = None
    perimeter: Optional[float] = None
    mass: Optional[float] = None
    ea: Optional[float] = None
    ga: Optional[float] = None
    nu_eff: Optional[float] = None
    e_eff: Optional[float] = None
    g_eff: Optional[float] = None
    qx: Optional[float] = None
    qy: Optional[float] = None
    ixx_g: Optional[float] = None
    iyy_g: Optional[float] = None
    ixy_g: Optional[float] = None
    cx: Optional[float] = None
    cy: Optional[float] = None
    ixx_c: Optional[float] = None
    iyy_c: Optional[float] = None
    ixy_c: Optional[float] = None
    zxx_plus: Optional[float] = None
    zxx_minus: Optional[float] = None
    zyy_plus: Optional[float] = None
    zyy_minus: Optional[float] = None
    rx_c: Optional[float] = None
    ry_c: Optional[float] = None
    i11_c: Optional[float] = None
    i22_c: Optional[float] = None
    phi: Optional[float] = None
    z11_plus: Optional[float] = None
    z11_minus: Optional[float] = None
    z22_plus: Optional[float] = None
    z22_minus: Optional[float] = None
    r11_c: Optional[float] = None
    r22_c: Optional[float] = None
    j: Optional[float] = None
    omega: Optional[np.ndarray] = None
    psi_shear: Optional[np.ndarray] = None
    phi_shear: Optional[np.ndarray] = None
    Delta_s: Optional[float] = None
    x_se: Optional[float] = None
    y_se: Optional[float] = None
    x11_se: Optional[float] = None
    y22_se: Optional[float] = None
    x_st: Optional[float] = None
    y_st: Optional[float] = None
    gamma: Optional[float] = None
    A_sx: Optional[float] = None
    A_sy: Optional[float] = None
    A_sxy: Optional[float] = None
    A_s11: Optional[float] = None
    A_s22: Optional[float] = None
    beta_x_plus: Optional[float] = None
    beta_x_minus: Optional[float] = None
    beta_y_plus: Optional[float] = None
    beta_y_minus: Optional[float] = None
    beta_11_plus: Optional[float] = None
    beta_11_minus: Optional[float] = None
    beta_22_plus: Optional[float] = None
    beta_22_minus: Optional[float] = None
    x_pc: Optional[float] = None
    y_pc: Optional[float] = None
    x11_pc: Optional[float] = None
    y22_pc: Optional[float] = None
    sxx: Optional[float] = None
    syy: Optional[float] = None
    sf_xx_plus: Optional[float] = None
    sf_xx_minus: Optional[float] = None
    sf_yy_plus: Optional[float] = None
    sf_yy_minus: Optional[float] = None
    s11: Optional[float] = None
    s22: Optional[float] = None
    sf_11_plus: Optional[float] = None
    sf_11_minus: Optional[float] = None
    sf_22_plus: Optional[float] = None
    sf_22_minus: Optional[float] = None

    def asdict(self):
        """
        Returns the SectionProperties dataclass object as a dictionary.
        """
        return asdict(self)

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
        nodes = np.array(mesh["vertices"])
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
            self.phi = np.arctan2(self.ixx_c - self.i11_c, self.ixy_c) * 180 / np.pi

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
