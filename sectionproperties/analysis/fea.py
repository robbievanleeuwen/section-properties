import numpy as np


class Tri6:
    """Class for a six noded quadratic triangular element.

    Provides methods for the calculation of section properties based on the
    finite element method.

    :cvar coords: A 2 x 6 array of the coordinates of the tri-6 nodes. The
        first three columns relate to the vertices of the triangle and the last
        three columns correspond to the mid-nodes.
    :vartype coords: :class:`numpy.ndarray`
    :cvar node_ids: A list of the global node ids for the current element
    :vartype node_ids: list[int]
    :cvar int attribute: Attribute number for the element
    :cvar float area: Area of the element
    :cvar float qx: First moment of area of the element about the x-axis
    :cvar float qy: First moment of area of the element about the y-axis
    :cvar float ixx: Second moment of area of the element about the x-axis
    :cvar float iyy: Second moment of area of the element about the y-axis
    :cvar float ixy: Second moment of area of the element about the xy-axis
    """

    def __init__(self, coords, node_ids, attribute):
        """Inits the Tri6 class. Calculates and stores the area properties of
        the element.

        :param coords: A 2 x 6 array of the coordinates of the tri-6 nodes. The
            first three columns relate to the vertices of the triangle and the
            last three columns correspond to the mid-nodes
        :type coords: :class:`numpy.ndarray`
        :param node_ids: A list of the global node ids for the current element
        :type node_ids: list[int]
        :param int attribute: Attribute number for the element
        """

        self.coords = coords
        self.node_ids = node_ids
        self.attribute = attribute
        self.nu = 0  # TODO: implement Poissons ratio

        # initialise area properties
        self.area = 0
        self.qx = 0
        self.qy = 0
        self.ixx = 0
        self.iyy = 0
        self.ixy = 0

        # Gauss points for 6 point Gaussian integration
        gps = gauss_points(6)

        # loop through each Gauss point
        for gp in gps:
            # determine shape function, shape function derivative and jacobian
            (N, B, j) = shape_function(self.coords, gp)

            self.area += gp[0] * j
            self.qx += gp[0] * np.dot(N, np.transpose(self.coords[1, :])) * j
            self.qy += gp[0] * np.dot(N, np.transpose(self.coords[0, :])) * j
            self.ixx += gp[0] * (
                np.dot(N, np.transpose(self.coords[1, :]))) ** 2 * j
            self.iyy += gp[0] * (
                np.dot(N, np.transpose(self.coords[0, :]))) ** 2 * j
            self.ixy += (gp[0] * np.dot(N, np.transpose(self.coords[1, :])) *
                         np.dot(N, np.transpose(self.coords[0, :])) * j)

    def torsion_properties(self):
        """Calculates the element stiffness matrix used for warping analysis
        and the torsion load vector.

        :return: Element stiffness matrix (k_el) and element torsion load
            vector (f_el)
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

            k_el += gp[0] * np.dot(np.transpose(B), B) * j
            f_el += gp[0] * np.dot(np.transpose(B), np.transpose(np.array(
                [Ny, -Nx]))) * j

        return (k_el, f_el)

    def shear_load_vectors(self, ixx, iyy, ixy):
        """Calculates the element shear load vectors used to evaluate the shear
        functions.

        :return: Element shear load vector psi (f_psi) and phi (f_phi)
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

            f_psi += gp[0] * (self.nu / 2 * np.transpose(
                np.transpose(B).dot(np.array([[d1], [d2]])))[0] + 2 *
                (1 + self.nu) * np.transpose(N) * (ixx * Nx - ixy * Ny)) * j
            f_phi += gp[0] * (self.nu / 2 * np.transpose(
                np.transpose(B).dot(np.array([[h1], [h2]])))[0] + 2 *
                (1 + self.nu) * np.transpose(N) * (iyy * Ny - ixy * Nx)) * j

        return (f_psi, f_phi)

    def shear_warping_integrals(self, ixx, iyy, ixy, omega):
        """Calculates the element shear centre and warping integrals required
        for shear analysis of the cross-section.

        :return: Shear centre integrals about the x and y-axes (sc_xint,
            sc_yint), warping integrals (q_omega, i_omega, i_xomega, i_yomega)
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

            sc_xint += gp[0] * (iyy * Nx + ixy * Ny) * (Nx ** 2 + Ny ** 2) * j
            sc_yint += gp[0] * (ixx * Ny + ixy * Nx) * (Nx ** 2 + Ny ** 2) * j
            q_omega += gp[0] * Nomega * j
            i_omega += gp[0] * Nomega ** 2 * j
            i_xomega += gp[0] * Nx * Nomega * j
            i_yomega += gp[0] * Ny * Nomega * j

        return (sc_xint, sc_yint, q_omega, i_omega, i_xomega, i_yomega)

    def shear_coefficients(self, ixx, iyy, ixy, psi_shear, phi_shear):
        """Calculates the variables used to determine the shear
        deformation coefficients.

        :return: Shear deformation variables (kappa_x, kappa_y, kappa_xy)
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

            kappa_x += gp[0] * (
                psi_shear.dot(np.transpose(B)) - self.nu / 2 * np.array(
                    [d1, d2])).dot(B.dot(psi_shear) - self.nu / 2 * np.array(
                        [d1, d2])) * j
            kappa_y += gp[0] * (
                phi_shear.dot(np.transpose(B)) - self.nu / 2 * np.array(
                    [h1, h2])).dot(B.dot(phi_shear) - self.nu / 2 * np.array(
                        [h1, h2])) * j
            kappa_xy += gp[0] * (
                psi_shear.dot(np.transpose(B)) - self.nu / 2 * np.array(
                    [d1, d2])).dot(B.dot(phi_shear) - self.nu / 2 * np.array(
                        [h1, h2])) * j

        return (kappa_x, kappa_y, kappa_xy)

    # def calculateStress(self, area, ixx, iyy, ixy, i11, i22, phi, omega, J,
    #                     Psi, Phi, Delta_s):
    #     """
    #     This method calculates the stresses within an an element resulting
    #     from unit loading.
    #     """
    #
    #     # initialise stress vectors
    #     sigma_zz_axial = np.ones((6, 1)) * 1 / area
    #     sigma_zz_bending_xx_gp = np.zeros((6, 1))
    #     sigma_zz_bending_yy_gp = np.zeros((6, 1))
    #     sigma_zz_bending_11_gp = np.zeros((6, 1))
    #     sigma_zz_bending_22_gp = np.zeros((6, 1))
    #     tau_torsion_gp = np.zeros((6, 2))
    #     tau_shear_x_gp = np.zeros((6, 2))
    #     tau_shear_y_gp = np.zeros((6, 2))
    #
    #     # Gauss points for 6 point Gaussian integration
    #     gps = femUtilities.gaussPoints(6)
    #
    #     for (i, gp) in enumerate(gps):
    #         # determine shape function, shape function derivative and jacobian
    #         (N, B, j) = femUtilities.shapeFunction(self.xy, gp)
    #
    #         # determine x and y position at Gauss point
    #         Nx = np.dot(N, np.transpose(self.xy[0, :]))
    #         Ny = np.dot(N, np.transpose(self.xy[1, :]))
    #
    #         # determine 11 and 22 position at Gauss point
    #         (Nx_1, Ny_2) = principalCoordinate(phi, Nx, Ny)
    #
    #         # determine shear parameters
    #         r = Nx ** 2 - Ny ** 2
    #         q = 2 * Nx * Ny
    #         d1 = ixx * r - ixy * q
    #         d2 = ixy * r + ixx * q
    #         h1 = -ixy * r + iyy * q
    #         h2 = -iyy * r - ixy * q
    #
    #         sigma_zz_bending_xx_gp[i, :] = (
    #             -(ixy * 1) / (ixx * iyy - ixy ** 2) * Nx + (iyy * 1) /
    #             (ixx * iyy - ixy ** 2) * Ny)
    #         sigma_zz_bending_yy_gp[i, :] = (
    #             -(ixx * 1) / (ixx * iyy - ixy ** 2) * Nx + (ixy * 1) /
    #             (ixx * iyy - ixy ** 2) * Ny)
    #         sigma_zz_bending_11_gp[i, :] = 1 / i11 * Ny_2
    #         sigma_zz_bending_22_gp[i, :] = -1 / i22 * Nx_1
    #         tau_torsion_gp[i, :] = (1 / J * (
    #             B.dot(omega) - np.array([Ny, -Nx])))
    #         tau_shear_x_gp[i, :] = (
    #             1 / Delta_s * (B.dot(Psi) - self.nu / 2 * np.array([d1, d2])))
    #         tau_shear_y_gp[i, :] = (
    #             1 / Delta_s * (B.dot(Phi) - self.nu / 2 * np.array([h1, h2])))
    #
    #     # extrapolate results to nodes
    #     sigma_zz_bending_xx = femUtilities.extrapolateToNodes(
    #         sigma_zz_bending_xx_gp[:, 0])
    #     sigma_zz_bending_yy = femUtilities.extrapolateToNodes(
    #         sigma_zz_bending_yy_gp[:, 0])
    #     sigma_zz_bending_11 = femUtilities.extrapolateToNodes(
    #         sigma_zz_bending_11_gp[:, 0])
    #     sigma_zz_bending_22 = femUtilities.extrapolateToNodes(
    #         sigma_zz_bending_22_gp[:, 0])
    #     tau_zx_torsion = femUtilities.extrapolateToNodes(tau_torsion_gp[:, 0])
    #     tau_zy_torsion = femUtilities.extrapolateToNodes(tau_torsion_gp[:, 1])
    #     tau_shear_zx_x = femUtilities.extrapolateToNodes(tau_shear_x_gp[:, 0])
    #     tau_shear_zy_x = femUtilities.extrapolateToNodes(tau_shear_x_gp[:, 1])
    #     tau_shear_zx_y = femUtilities.extrapolateToNodes(tau_shear_y_gp[:, 0])
    #     tau_shear_zy_y = femUtilities.extrapolateToNodes(tau_shear_y_gp[:, 1])
    #
    #     return (sigma_zz_axial, sigma_zz_bending_xx, sigma_zz_bending_yy,
    #             sigma_zz_bending_11, sigma_zz_bending_22, tau_zx_torsion,
    #             tau_zy_torsion, tau_shear_zx_x, tau_shear_zy_x,
    #             tau_shear_zx_y, tau_shear_zy_y, gps[:, 0])


def gauss_points(n):
    """Returns the Gaussian weights and locations for *n* point Gaussian
    integration of a quadratic triangular element.

    :param int n: Number of Gauss points (1, 3 or 6)
    :return: An *n x 4* matrix consisting of the integration weight and the
        eta, xi and zeta locations for *n* Gauss points
    :rtype: :class:`numpy.ndarray`
    """

    if n == 1:
        # one point gaussian integration
        return np.array([[1, 1.0 / 3, 1.0 / 3, 1.0 / 3]])

    elif n == 3:
        # three point gaussian integration
        return np.array([[1.0 / 3, 2.0 / 3, 1.0 / 6, 1.0 / 6],
                         [1.0 / 3, 1.0 / 6, 2.0 / 3, 1.0 / 6],
                         [1.0 / 3, 1.0 / 6, 1.0 / 6, 2.0 / 3]])
    elif n == 6:
        # six point gaussian integration
        g1 = 1.0 / 18 * (8 - np.sqrt(10) + np.sqrt(38 - 44 * np.sqrt(2.0 / 5)))
        g2 = 1.0 / 18 * (8 - np.sqrt(10) - np.sqrt(38 - 44 * np.sqrt(2.0 / 5)))
        w1 = (620 + np.sqrt(213125 - 53320 * np.sqrt(10))) / 3720
        w2 = (620 - np.sqrt(213125 - 53320 * np.sqrt(10))) / 3720

        return np.array([[w2, 1 - 2 * g2, g2, g2],
                         [w2, g2, 1 - 2 * g2, g2],
                         [w2, g2, g2, 1 - 2 * g2],
                         [w1, g1, g1, 1 - 2 * g1],
                         [w1, 1 - 2 * g1, g1, g1],
                         [w1, g1, 1 - 2 * g1, g1]])


def shape_function(coords, gauss_point):
    """Computes shape functions, shape function derivatives and the determinant
    of the Jacobian matrix for a tri 6 element at a given Gauss point.

    :param coords: Global coordinates of the quadratic triangle vertices
        [2 x 6]
    :type coords: :class:`numpy.ndarray`
    :param gauss_point: Gaussian weight and isoparametric location of the Gauss
        point
    :type gauss_point: :class:`numpy.ndarray`
    :return: The value of the shape functions *N(i)* at the given Gauss point
        [1 x 6], the derivative of the shape functions in the j-th global
        direction *B(i,j)* [2 x 6] and the determinant of the Jacobian matrix
        *j*
    :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`, float)
    """

    # location of isoparametric co-ordinates for each Gauss point
    eta = gauss_point[1]
    xi = gauss_point[2]
    zeta = gauss_point[3]

    # value of the shape functions
    N = np.array([eta * (2 * eta - 1),
                  xi * (2 * xi - 1),
                  zeta * (2 * zeta - 1),
                  4 * eta * xi,
                  4 * xi * zeta,
                  4 * eta * zeta])

    # derivatives of the shape functions wrt the isoparametric co-ordinates
    B_iso = np.array([[4 * eta - 1, 0, 0, 4 * xi, 0, 4 * zeta],
                      [0, 4 * xi - 1, 0, 4 * eta, 4 * zeta, 0],
                      [0, 0, 4 * zeta - 1, 0, 4 * xi, 4 * eta]])

    # form Jacobian matrix
    J_upper = np.array([[1, 1, 1]])
    J_lower = np.dot(coords, np.transpose(B_iso))
    J = np.vstack((J_upper, J_lower))

    # calculate the jacobian
    try:
        j = 0.5 * np.linalg.det(J)
    except ValueError:
        # handle warning if area is zero during plastic centroid algorithm
        j = 0

    # cacluate the P matrix
    if j != 0:
        P = np.dot(np.linalg.inv(J), np.array([[0, 0], [1, 0], [0, 1]]))
        # calculate the B matrix in terms of cartesian co-ordinates
        B = np.transpose(np.dot(np.transpose(B_iso), P))

    return (N, B, j)


def extrapolate_to_nodes(w):
    """Extrapolates results at six Gauss points to the six noes of a quadratic
    triangular element.

    :param w: Result at the six Gauss points [1 x 6]
    :type w: :class:`numpy.ndarray`
    :return: Extrapolated nodal values at the six nodes [1 x 6]
    :rtype: :class:`numpy.ndarray`
    """

    H_inv = np.array(
        [[1.87365927351160,	0.138559587411935, 0.138559587411935,
          -0.638559587411936, 0.126340726488397, -0.638559587411935],
         [0.138559587411935, 1.87365927351160, 0.138559587411935,
          -0.638559587411935, -0.638559587411935, 0.126340726488397],
         [0.138559587411935, 0.138559587411935, 1.87365927351160,
          0.126340726488396, -0.638559587411935, -0.638559587411935],
         [0.0749010751157440, 0.0749010751157440, 0.180053080734478,
          1.36051633430762,	-0.345185782636792, -0.345185782636792],
         [0.180053080734478, 0.0749010751157440, 0.0749010751157440,
          -0.345185782636792, 1.36051633430762, -0.345185782636792],
         [0.0749010751157440, 0.180053080734478,  0.0749010751157440,
          -0.345185782636792, -0.345185782636792, 1.36051633430762]])

    return H_inv.dot(w)


def principal_coordinate(phi, x, y):
    """Determines the coordinates of the cartesian point *(x, y)* in the
    principal axis system given an axis rotation angle phi.

    :param float phi: Prinicpal axis rotation angle
    :param float x: x coordinate in the global axis
    :param float y: y coordinate in the global axis
    :return: Principal axis coordinates *(x1, y2)*
    :rtype: tuple(float, float)
    """

    # convert principal axis angle to radians
    phi_rad = phi * np.pi / 180

    # form rotation matrix
    R = np.array([[np.cos(phi_rad), np.sin(phi_rad)],
                  [-np.sin(phi_rad), np.cos(phi_rad)]])

    # calculate rotated x and y coordinates
    x_rotated = R.dot(np.array([x, y]))

    return (x_rotated[0], x_rotated[1])


def global_coordinate(phi, x1, y2):
    """Determines the global coordinates of the principal axis point *(x1, y2)*
    given principal axis rotation angle phi.

    :param float phi: Prinicpal axis rotation angle
    :param float x1: 11 coordinate in the principal axis
    :param float y2: 22 coordinate in the principal axis
    :return: Global axis coordinates *(x, y)*
    :rtype: tuple(float, float)
    """

    # convert principal axis angle to radians
    phi_rad = phi * np.pi / 180

    # form transposed rotation matrix
    R = (np.array([[np.cos(phi_rad), -np.sin(phi_rad)], [np.sin(phi_rad),
                                                         np.cos(phi_rad)]]))
    # calculate rotated x_1 and y_2 coordinates
    x_rotated = R.dot(np.array([x1, y2]))

    return (x_rotated[0], x_rotated[1])


def point_above_line(u, px, py, x, y):
    """Determines whether a point *(x, y)* is a above or below the line defined
    by the parallel unit vector *u* and the point *(px, py)*.

    :param u: Unit vector parallel to the line [1 x 2]
    :type u: :class:`numpy.ndarray`
    :param float px: x coordinate of a point on the line
    :param float py: y coordinate of a point on the line
    :param float x: x coordinate of the point to be tested
    :param float y: y coordinate of the point to be tested
    :return: This method returns *True* if the point is above the line or
        *False* if the point is below the line
    :rtype: bool
    """

    # vector from point to point on line
    PQ = np.array([px - x, py - y])
    return np.cross(PQ, u) > 0
