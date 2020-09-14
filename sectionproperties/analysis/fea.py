import numpy as np


class Tri6:
    """Class for a six noded quadratic triangular element.

    Provides methods for the calculation of section properties based on the finite element method.

    :param int el_id: Unique element id
    :param coords: A 2 x 6 array of the coordinates of the tri-6 nodes. The first three columns
        relate to the vertices of the triangle and the last three columns correspond to the
        mid-nodes.
    :type coords: :class:`numpy.ndarray`
    :param node_ids: A list of the global node ids for the current element
    :type node_ids: list[int]
    :param material: Material object for the current finite element.
    :type material: :class:`~sectionproperties.pre.pre.Material`

    :cvar int el_id: Unique element id
    :cvar coords: A 2 x 6 array of the coordinates of the tri-6 nodes. The first three columns
        relate to the vertices of the triangle and the last three columns correspond to the
        mid-nodes.
    :vartype coords: :class:`numpy.ndarray`
    :cvar node_ids: A list of the global node ids for the current element
    :vartype node_ids: list[int]
    :cvar material: Material of the current finite element.
    :vartype material: :class:`~sectionproperties.pre.pre.Material`
    """

    def __init__(self, el_id, coords, node_ids, material):
        """Inits the Tri6 class."""

        self.el_id = el_id
        self.coords = coords
        self.node_ids = node_ids
        self.material = material

    def geometric_properties(self):
        """Calculates the geometric properties for the current finite element.

        :return: Tuple containing the geometric properties and the elastic and shear moduli of the
            element: *(area, qx, qy, ixx, iyy, ixy, e, g)*
        :rtype: tuple(float)
        """

        # initialise geometric properties
        area = 0
        qx = 0
        qy = 0
        ixx = 0
        iyy = 0
        ixy = 0

        # Gauss points for 6 point Gaussian integration
        gps = gauss_points(6)

        # loop through each Gauss point
        for gp in gps:
            # determine shape function, shape function derivative and jacobian
            (N, _, j) = shape_function(self.coords, gp)

            area += gp[0] * j
            qx += gp[0] * np.dot(N, np.transpose(self.coords[1, :])) * j
            qy += gp[0] * np.dot(N, np.transpose(self.coords[0, :])) * j
            ixx += gp[0] * np.dot(N, np.transpose(self.coords[1, :])) ** 2 * j
            iyy += gp[0] * np.dot(N, np.transpose(self.coords[0, :])) ** 2 * j
            ixy += (
                gp[0] * np.dot(N, np.transpose(self.coords[1, :])) * np.dot(
                    N, np.transpose(self.coords[0, :])) * j
            )

        return (
            area, qx, qy, ixx, iyy, ixy, self.material.elastic_modulus, self.material.shear_modulus
        )

    def torsion_properties(self):
        """Calculates the element stiffness matrix used for warping analysis and the torsion load
        vector.

        :return: Element stiffness matrix *(k_el)* and element torsion load vector *(f_el)*
        :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`)
        """

        # initialise stiffness matrix and load vector
        k_el = 0
        f_el = 0

        # Gauss points for 6 point Gaussian integration
        gps = gauss_points(6)

        for gp in gps:
            # determine shape function, shape function derivative and jacobian
            (N, B, j) = shape_function(self.coords, gp)

            # determine x and y position at Gauss point
            Nx = np.dot(N, np.transpose(self.coords[0, :]))
            Ny = np.dot(N, np.transpose(self.coords[1, :]))

            # calculated modulus weighted stiffness matrix and load vector
            k_el += gp[0] * np.dot(np.transpose(B), B) * j * (self.material.elastic_modulus)
            f_el += (
                gp[0] * np.dot(np.transpose(B), np.transpose(np.array([Ny, -Nx])))
                * j * self.material.elastic_modulus
            )

        return (k_el, f_el)

    def shear_load_vectors(self, ixx, iyy, ixy, nu):
        """Calculates the element shear load vectors used to evaluate the shear functions.

        :param float ixx: Second moment of area about the centroidal x-axis
        :param float iyy: Second moment of area about the centroidal y-axis
        :param float ixy: Second moment of area about the centroidal xy-axis
        :param float nu: Effective Poisson's ratio for the cross-section

        :return: Element shear load vector psi *(f_psi)* and phi *(f_phi)*
        :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`)
        """

        # initialise force vectors
        f_psi = 0
        f_phi = 0

        # Gauss points for 6 point Gaussian integration
        gps = gauss_points(6)

        for gp in gps:
            # determine shape function, shape function derivative and jacobian
            (N, B, j) = shape_function(self.coords, gp)

            # determine x and y position at Gauss point
            Nx = np.dot(N, np.transpose(self.coords[0, :]))
            Ny = np.dot(N, np.transpose(self.coords[1, :]))

            # determine shear parameters
            r = Nx ** 2 - Ny ** 2
            q = 2 * Nx * Ny
            d1 = ixx * r - ixy * q
            d2 = ixy * r + ixx * q
            h1 = -ixy * r + iyy * q
            h2 = -iyy * r - ixy * q

            f_psi += (
                gp[0] * (nu / 2 * np.transpose(np.transpose(B).dot(np.array([[d1], [d2]])))[0]
                         + 2 * (1 + nu) * np.transpose(N) * (ixx * Nx - ixy * Ny)) * j
                * self.material.elastic_modulus
            )
            f_phi += (
                gp[0] * (nu / 2 * np.transpose(np.transpose(B).dot(np.array([[h1], [h2]])))[0]
                         + 2 * (1 + nu) * np.transpose(N) * (iyy * Ny - ixy * Nx)) * j
                * self.material.elastic_modulus
            )

        return (f_psi, f_phi)

    def shear_warping_integrals(self, ixx, iyy, ixy, omega):
        """Calculates the element shear centre and warping integrals required for shear analysis of
        the cross-section.

        :param float ixx: Second moment of area about the centroidal x-axis
        :param float iyy: Second moment of area about the centroidal y-axis
        :param float ixy: Second moment of area about the centroidal xy-axis
        :param omega: Values of the warping function at the element nodes
        :type omega: :class:`numpy.ndarray`

        :return: Shear centre integrals about the x and y-axes *(sc_xint, sc_yint)*, warping
            integrals *(q_omega, i_omega, i_xomega, i_yomega)*
        :rtype: tuple(float, float, float, float, float, float)
        """

        # initialise integrals
        sc_xint = 0
        sc_yint = 0
        q_omega = 0
        i_omega = 0
        i_xomega = 0
        i_yomega = 0

        # Gauss points for 6 point Gaussian integration
        gps = gauss_points(6)

        for gp in gps:
            # determine shape function, shape function derivative and jacobian
            (N, B, j) = shape_function(self.coords, gp)

            # determine x and y position at Gauss point
            Nx = np.dot(N, np.transpose(self.coords[0, :]))
            Ny = np.dot(N, np.transpose(self.coords[1, :]))
            Nomega = np.dot(N, np.transpose(omega))

            sc_xint += (
                gp[0] * (iyy * Nx + ixy * Ny) * (Nx ** 2 + Ny ** 2)
                * j * self.material.elastic_modulus
            )
            sc_yint += (
                gp[0] * (ixx * Ny + ixy * Nx) * (Nx ** 2 + Ny ** 2)
                * j * self.material.elastic_modulus
            )
            q_omega += gp[0] * Nomega * j * self.material.elastic_modulus
            i_omega += gp[0] * Nomega ** 2 * j * self.material.elastic_modulus
            i_xomega += gp[0] * Nx * Nomega * j * self.material.elastic_modulus
            i_yomega += gp[0] * Ny * Nomega * j * self.material.elastic_modulus

        return (sc_xint, sc_yint, q_omega, i_omega, i_xomega, i_yomega)

    def shear_coefficients(self, ixx, iyy, ixy, psi_shear, phi_shear, nu):
        """Calculates the variables used to determine the shear deformation coefficients.

        :param float ixx: Second moment of area about the centroidal x-axis
        :param float iyy: Second moment of area about the centroidal y-axis
        :param float ixy: Second moment of area about the centroidal xy-axis
        :param psi_shear: Values of the psi shear function at the element nodes
        :type psi_shear: :class:`numpy.ndarray`
        :param phi_shear: Values of the phi shear function at the element nodes
        :type phi_shear: :class:`numpy.ndarray`
        :param float nu: Effective Poisson's ratio for the cross-section

        :return: Shear deformation variables *(kappa_x, kappa_y, kappa_xy)*
        :rtype: tuple(float, float, float)
        """

        # initialise properties
        kappa_x = 0
        kappa_y = 0
        kappa_xy = 0

        # Gauss points for 6 point Gaussian integration
        gps = gauss_points(6)

        for gp in gps:
            # determine shape function, shape function derivative and jacobian
            (N, B, j) = shape_function(self.coords, gp)

            # determine x and y position at Gauss point
            Nx = np.dot(N, np.transpose(self.coords[0, :]))
            Ny = np.dot(N, np.transpose(self.coords[1, :]))

            # determine shear parameters
            r = Nx ** 2 - Ny ** 2
            q = 2 * Nx * Ny
            d1 = ixx * r - ixy * q
            d2 = ixy * r + ixx * q
            h1 = -ixy * r + iyy * q
            h2 = -iyy * r - ixy * q

            kappa_x += (
                gp[0] * (psi_shear.dot(np.transpose(B)) - nu / 2 * np.array([d1, d2])).dot(
                    B.dot(psi_shear) - nu / 2 * np.array([d1, d2])) * j
                * self.material.elastic_modulus
            )
            kappa_y += (
                gp[0] * (phi_shear.dot(np.transpose(B)) - nu / 2 * np.array([h1, h2])).dot(
                    B.dot(phi_shear) - nu / 2 * np.array([h1, h2])) * j
                * self.material.elastic_modulus
            )
            kappa_xy += (
                gp[0] * (psi_shear.dot(np.transpose(B)) - nu / 2 * np.array([d1, d2])).dot(
                    B.dot(phi_shear) - nu / 2 * np.array([h1, h2])) * j
                * self.material.elastic_modulus
            )

        return (kappa_x, kappa_y, kappa_xy)

    def monosymmetry_integrals(self, phi):
        """Calculates the integrals used to evaluate the monosymmetry constant about both global
        axes and both prinicipal axes.

        :param float phi: Principal bending axis angle

        :return: Integrals used to evaluate the monosymmetry constants *(int_x, int_y, int_11,
            int_22)*
        :rtype: tuple(float, float, float, float)
        """

        # initialise properties
        int_x = 0
        int_y = 0
        int_11 = 0
        int_22 = 0

        # Gauss points for 6 point Gaussian integration
        gps = gauss_points(6)

        for gp in gps:
            # determine shape function and jacobian
            (N, _, j) = shape_function(self.coords, gp)

            # determine x and y position at Gauss point
            Nx = np.dot(N, np.transpose(self.coords[0, :]))
            Ny = np.dot(N, np.transpose(self.coords[1, :]))

            # determine 11 and 22 position at Gauss point
            (Nx_11, Ny_22) = principal_coordinate(phi, Nx, Ny)

            # weight the monosymmetry integrals by the section elastic modulus
            int_x += gp[0] * (Nx * Nx * Ny + Ny * Ny * Ny) * j * self.material.elastic_modulus
            int_y += gp[0] * (Ny * Ny * Nx + Nx * Nx * Nx) * j * self.material.elastic_modulus
            int_11 += (
                gp[0] * (Nx_11 * Nx_11 * Ny_22 + Ny_22 * Ny_22 * Ny_22) * j
                * self.material.elastic_modulus
            )
            int_22 += (
                gp[0] * (Ny_22 * Ny_22 * Nx_11 + Nx_11 * Nx_11 * Nx_11) * j
                * self.material.elastic_modulus
            )

        return (int_x, int_y, int_11, int_22)

    def plastic_properties(self, u, p):
        """Calculates total force resisted by the element when subjected to a stress equal to the
        yield strength. Also returns the modulus weighted area and first moments of area, and
        determines whether or not the element is above or below the line defined by the unit
        vector *u* and point *p*.

        :param u: Unit vector in the direction of the line
        :type u: :class:`numpy.ndarray`
        :param p: Point on the line
        :type p: :class:`numpy.ndarray`

        :return: Element force *(force)*, modulus weighted area properties *(ea, e.qx, e.qy)* and
            whether or not the element is above the line
        :rtype: tuple(float, float, float, float, bool)
        """

        # initialise geometric properties
        e = self.material.elastic_modulus
        area = 0
        qx = 0
        qy = 0
        force = 0

        # Gauss points for 3 point Gaussian integration
        gps = gauss_points(3)

        # loop through each Gauss point
        for gp in gps:
            # determine shape function, shape function derivative and jacobian
            (N, _, j) = shape_function(self.coords, gp)

            area += gp[0] * j
            qx += gp[0] * np.dot(N, np.transpose(self.coords[1, :])) * j
            qy += gp[0] * np.dot(N, np.transpose(self.coords[0, :])) * j
            force += gp[0] * j * self.material.yield_strength

        # calculate element centroid
        (cx, cy) = (qy / area, qx / area)

        # determine if the element is above the line p + u
        is_above = point_above_line(u, p[0], p[1], cx, cy)

        return (force, area * e, qx * e, qy * e, is_above)

    def element_stress(self, N, Mxx, Myy, M11, M22, Mzz, Vx, Vy, ea, cx, cy, ixx, iyy, ixy, i11,
                       i22, phi, j, nu, omega, psi_shear, phi_shear, Delta_s):
        """Calculates the stress within an element resulting from a specified loading. Also returns
        the shape function weights.

        :param float N: Axial force
        :param float Mxx: Bending moment about the centroidal xx-axis
        :param float Myy: Bending moment about the centroidal yy-axis
        :param float M11: Bending moment about the centroidal 11-axis
        :param float M22: Bending moment about the centroidal 22-axis
        :param float Mzz: Torsion moment about the centroidal zz-axis
        :param float Vx: Shear force acting in the x-direction
        :param float Vy: Shear force acting in the y-direction
        :param float ea: Modulus weighted area
        :param float cx: x position of the elastic centroid
        :param float cy: y position of the elastic centroid
        :param float ixx: Second moment of area about the centroidal x-axis
        :param float iyy: Second moment of area about the centroidal y-axis
        :param float ixy: Second moment of area about the centroidal xy-axis
        :param float i11: Second moment of area about the principal 11-axis
        :param float i22: Second moment of area about the principal 22-axis
        :param float phi: Principal bending axis angle
        :param float j: St. Venant torsion constant
        :param float nu: Effective Poisson's ratio for the cross-section
        :param omega: Values of the warping function at the element nodes
        :type omega: :class:`numpy.ndarray`
        :param psi_shear: Values of the psi shear function at the element nodes
        :type psi_shear: :class:`numpy.ndarray`
        :param phi_shear: Values of the phi shear function at the element nodes
        :type phi_shear: :class:`numpy.ndarray`
        :param float Delta_s: Cross-section shear factor
        :return: Tuple containing element stresses and integration weights
            (:math:`\sigma_{zz,n}`, :math:`\sigma_{zz,mxx}`,
            :math:`\sigma_{zz,myy}`, :math:`\sigma_{zz,m11}`,
            :math:`\sigma_{zz,m22}`, :math:`\sigma_{zx,mzz}`,
            :math:`\sigma_{zy,mzz}`, :math:`\sigma_{zx,vx}`,
            :math:`\sigma_{zy,vx}`, :math:`\sigma_{zx,vy}`,
            :math:`\sigma_{zy,vy}`, :math:`w_i`)
        :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`, ...)
        """

        # calculate axial stress
        sig_zz_n = N * np.ones(6) * self.material.elastic_modulus / ea

        # initialise stresses at the gauss points
        sig_zz_mxx_gp = np.zeros((6, 1))
        sig_zz_myy_gp = np.zeros((6, 1))
        sig_zz_m11_gp = np.zeros((6, 1))
        sig_zz_m22_gp = np.zeros((6, 1))
        sig_zxy_mzz_gp = np.zeros((6, 2))
        sig_zxy_vx_gp = np.zeros((6, 2))
        sig_zxy_vy_gp = np.zeros((6, 2))

        # Gauss points for 6 point Gaussian integration
        gps = gauss_points(6)

        for (i, gp) in enumerate(gps):
            # determine x and y positions with respect to the centroidal axis
            coords_c = np.zeros((2, 6))
            coords_c[0, :] = self.coords[0, :] - cx
            coords_c[1, :] = self.coords[1, :] - cy

            # determine shape function, shape function derivative and jacobian
            (N, B, _) = shape_function(coords_c, gp)

            # determine x and y position at Gauss point
            Nx = np.dot(N, np.transpose(coords_c[0, :]))
            Ny = np.dot(N, np.transpose(coords_c[1, :]))

            # determine 11 and 22 position at Gauss point
            (Nx_11, Ny_22) = principal_coordinate(phi, Nx, Ny)

            # determine shear parameters
            r = Nx ** 2 - Ny ** 2
            q = 2 * Nx * Ny
            d1 = ixx * r - ixy * q
            d2 = ixy * r + ixx * q
            h1 = -ixy * r + iyy * q
            h2 = -iyy * r - ixy * q

            # calculate element stresses
            sig_zz_mxx_gp[i, :] = (
                self.material.elastic_modulus * (-(ixy * Mxx) / (ixx * iyy - ixy ** 2) * Nx + (
                    iyy * Mxx) / (ixx * iyy - ixy ** 2) * Ny)
            )
            sig_zz_myy_gp[i, :] = (
                self.material.elastic_modulus * (-(ixx * Myy) / (ixx * iyy - ixy ** 2) * Nx + (
                    ixy * Myy) / (ixx * iyy - ixy ** 2) * Ny)
            )
            sig_zz_m11_gp[i, :] = self.material.elastic_modulus * M11 / i11 * Ny_22
            sig_zz_m22_gp[i, :] = self.material.elastic_modulus * -M22 / i22 * Nx_11

            if Mzz != 0:
                sig_zxy_mzz_gp[i, :] = (
                    self.material.elastic_modulus * Mzz / j * (B.dot(omega) - np.array([Ny, -Nx]))
                )

            if Vx != 0:
                sig_zxy_vx_gp[i, :] = (
                    self.material.elastic_modulus * Vx / Delta_s * (
                        B.dot(psi_shear) - nu / 2 * np.array([d1, d2]))
                )

            if Vy != 0:
                sig_zxy_vy_gp[i, :] = (
                    self.material.elastic_modulus * Vy / Delta_s * (
                        B.dot(phi_shear) - nu / 2 * np.array([h1, h2]))
                )

        # extrapolate results to nodes
        sig_zz_mxx = extrapolate_to_nodes(sig_zz_mxx_gp[:, 0])
        sig_zz_myy = extrapolate_to_nodes(sig_zz_myy_gp[:, 0])
        sig_zz_m11 = extrapolate_to_nodes(sig_zz_m11_gp[:, 0])
        sig_zz_m22 = extrapolate_to_nodes(sig_zz_m22_gp[:, 0])
        sig_zx_mzz = extrapolate_to_nodes(sig_zxy_mzz_gp[:, 0])
        sig_zy_mzz = extrapolate_to_nodes(sig_zxy_mzz_gp[:, 1])
        sig_zx_vx = extrapolate_to_nodes(sig_zxy_vx_gp[:, 0])
        sig_zy_vx = extrapolate_to_nodes(sig_zxy_vx_gp[:, 1])
        sig_zx_vy = extrapolate_to_nodes(sig_zxy_vy_gp[:, 0])
        sig_zy_vy = extrapolate_to_nodes(sig_zxy_vy_gp[:, 1])

        return (sig_zz_n, sig_zz_mxx, sig_zz_myy, sig_zz_m11, sig_zz_m22, sig_zx_mzz, sig_zy_mzz,
                sig_zx_vx, sig_zy_vx, sig_zx_vy, sig_zy_vy, gps[:, 0])

    def point_within_element(self, pt):
        """Determines whether a point lies within the current element.

        :param pt: Point to check *(x, y)*
        :type pt: list[float, float]
        :return: Whether the point lies within an element
        :rtype: bool
        """

        px = pt[0]
        py = pt[1]

        # get coordinates of corner points
        x1 = self.coords[0][0]
        y1 = self.coords[1][0]
        x2 = self.coords[0][1]
        y2 = self.coords[1][1]
        x3 = self.coords[0][2]
        y3 = self.coords[1][2]

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
        else:
            return False


def gauss_points(n):
    """Returns the Gaussian weights and locations for *n* point Gaussian integration of a quadratic
    triangular element.

    :param int n: Number of Gauss points (1, 3 or 6)
    :return: An *n x 4* matrix consisting of the integration weight and the eta, xi and zeta
        locations for *n* Gauss points
    :rtype: :class:`numpy.ndarray`
    """

    if n == 1:
        # one point gaussian integration
        return np.array([[1, 1.0 / 3, 1.0 / 3, 1.0 / 3]])

    elif n == 3:
        # three point gaussian integration
        return np.array([
            [1.0 / 3, 2.0 / 3, 1.0 / 6, 1.0 / 6],
            [1.0 / 3, 1.0 / 6, 2.0 / 3, 1.0 / 6],
            [1.0 / 3, 1.0 / 6, 1.0 / 6, 2.0 / 3]
        ])
    elif n == 6:
        # six point gaussian integration
        g1 = 1.0 / 18 * (8 - np.sqrt(10) + np.sqrt(38 - 44 * np.sqrt(2.0 / 5)))
        g2 = 1.0 / 18 * (8 - np.sqrt(10) - np.sqrt(38 - 44 * np.sqrt(2.0 / 5)))
        w1 = (620 + np.sqrt(213125 - 53320 * np.sqrt(10))) / 3720
        w2 = (620 - np.sqrt(213125 - 53320 * np.sqrt(10))) / 3720

        return np.array([
            [w2, 1 - 2 * g2, g2, g2],
            [w2, g2, 1 - 2 * g2, g2],
            [w2, g2, g2, 1 - 2 * g2],
            [w1, g1, g1, 1 - 2 * g1],
            [w1, 1 - 2 * g1, g1, g1],
            [w1, g1, 1 - 2 * g1, g1]
        ])


def shape_function(coords, gauss_point):
    """Computes shape functions, shape function derivatives and the determinant of the Jacobian
    matrix for a tri 6 element at a given Gauss point.

    :param coords: Global coordinates of the quadratic triangle vertices [2 x 6]
    :type coords: :class:`numpy.ndarray`
    :param gauss_point: Gaussian weight and isoparametric location of the Gauss point
    :type gauss_point: :class:`numpy.ndarray`
    :return: The value of the shape functions *N(i)* at the given Gauss point [1 x 6], the
        derivative of the shape functions in the j-th global direction *B(i,j)* [2 x 6] and the
        determinant of the Jacobian matrix *j*
    :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`, float)
    """

    # location of isoparametric co-ordinates for each Gauss point
    eta = gauss_point[1]
    xi = gauss_point[2]
    zeta = gauss_point[3]

    # value of the shape functions
    N = np.array([
        eta * (2 * eta - 1),
        xi * (2 * xi - 1),
        zeta * (2 * zeta - 1),
        4 * eta * xi,
        4 * xi * zeta,
        4 * eta * zeta
    ])

    # derivatives of the shape functions wrt the isoparametric co-ordinates
    B_iso = np.array([
        [4 * eta - 1, 0, 0, 4 * xi, 0, 4 * zeta],
        [0, 4 * xi - 1, 0, 4 * eta, 4 * zeta, 0],
        [0, 0, 4 * zeta - 1, 0, 4 * xi, 4 * eta]
    ])

    # form Jacobian matrix
    J_upper = np.array([[1, 1, 1]])
    J_lower = np.dot(coords, np.transpose(B_iso))
    J = np.vstack((J_upper, J_lower))

    # calculate the jacobian
    j = 0.5 * np.linalg.det(J)

    # if the area of the element is not zero
    if j != 0:
        # cacluate the P matrix
        P = np.dot(np.linalg.inv(J), np.array([[0, 0], [1, 0], [0, 1]]))

        # calculate the B matrix in terms of cartesian co-ordinates
        B = np.transpose(np.dot(np.transpose(B_iso), P))
    else:
        B = np.zeros((2, 6))  # empty B matrix

    return (N, B, j)


def extrapolate_to_nodes(w):
    """Extrapolates results at six Gauss points to the six noes of a quadratic triangular element.

    :param w: Result at the six Gauss points [1 x 6]
    :type w: :class:`numpy.ndarray`
    :return: Extrapolated nodal values at the six nodes [1 x 6]
    :rtype: :class:`numpy.ndarray`
    """

    H_inv = np.array([
        [1.87365927351160, 0.138559587411935, 0.138559587411935,
         -0.638559587411936, 0.126340726488397, -0.638559587411935],
        [0.138559587411935, 1.87365927351160, 0.138559587411935,
         -0.638559587411935, -0.638559587411935, 0.126340726488397],
        [0.138559587411935, 0.138559587411935, 1.87365927351160,
         0.126340726488396, -0.638559587411935, -0.638559587411935],
        [0.0749010751157440, 0.0749010751157440, 0.180053080734478,
         1.36051633430762, -0.345185782636792, -0.345185782636792],
        [0.180053080734478, 0.0749010751157440, 0.0749010751157440,
         -0.345185782636792, 1.36051633430762, -0.345185782636792],
        [0.0749010751157440, 0.180053080734478, 0.0749010751157440,
         -0.345185782636792, -0.345185782636792, 1.36051633430762]
    ])

    return H_inv.dot(w)


def principal_coordinate(phi, x, y):
    """Determines the coordinates of the cartesian point *(x, y)* in the
    principal axis system given an axis rotation angle phi.

    :param float phi: Prinicpal bending axis angle (degrees)
    :param float x: x coordinate in the global axis
    :param float y: y coordinate in the global axis
    :return: Principal axis coordinates *(x1, y2)*
    :rtype: tuple(float, float)
    """

    # convert principal axis angle to radians
    phi_rad = phi * np.pi / 180

    # form rotation matrix
    R = np.array([
        [np.cos(phi_rad), np.sin(phi_rad)],
        [-np.sin(phi_rad), np.cos(phi_rad)]
    ])

    # calculate rotated x and y coordinates
    x_rotated = R.dot(np.array([x, y]))

    return (x_rotated[0], x_rotated[1])


def global_coordinate(phi, x11, y22):
    """Determines the global coordinates of the principal axis point *(x1, y2)* given principal
    axis rotation angle phi.

    :param float phi: Prinicpal bending axis angle (degrees)
    :param float x11: 11 coordinate in the principal axis
    :param float y22: 22 coordinate in the principal axis
    :return: Global axis coordinates *(x, y)*
    :rtype: tuple(float, float)
    """

    # convert principal axis angle to radians
    phi_rad = phi * np.pi / 180

    # form transposed rotation matrix
    R = np.array([
        [np.cos(phi_rad), -np.sin(phi_rad)],
        [np.sin(phi_rad), np.cos(phi_rad)]
    ])
    # calculate rotated x_1 and y_2 coordinates
    x_rotated = R.dot(np.array([x11, y22]))

    return (x_rotated[0], x_rotated[1])


def point_above_line(u, px, py, x, y):
    """Determines whether a point *(x, y)* is a above or below the line defined by the parallel
    unit vector *u* and the point *(px, py)*.

    :param u: Unit vector parallel to the line [1 x 2]
    :type u: :class:`numpy.ndarray`
    :param float px: x coordinate of a point on the line
    :param float py: y coordinate of a point on the line
    :param float x: x coordinate of the point to be tested
    :param float y: y coordinate of the point to be tested
    :return: This method returns *True* if the point is above the line or *False* if the point is
        below the line
    :rtype: bool
    """

    # vector from point to point on line
    PQ = np.array([px - x, py - y])
    return np.cross(PQ, u) > 0
