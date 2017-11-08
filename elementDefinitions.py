'''
This module contains the definition of the tri6 element class and includes the
initialisation of the tri6 element and provides functions used at element level.
'''

import numpy as np
import femFunctions

class tri6:
    '''
    Six noded quadratic triangular element class, input a 2 x 6 array containing
    the coordinates of the vertices, the node numbers of the vertices and the
    poissons ratio.
    '''

    def __init__(self, vertices, nodes, nu):
        self.xy = vertices # triangle vertex co-ordinates [2 x 6]
        self.nodes = nodes # array of node numbers [1 x 6]
        self.nu = nu # poissons ratio of material

    def areaProperties(self):
        '''
        This function calculates the simple area propeties of element.
        '''
        # initialise properties
        area = 0
        Qx = 0
        Qy = 0

        # Gauss points for 3 point Gaussian integration
        gps = femFunctions.gaussPoints(3)

        for gp in gps:
            # determine shape function, shape function derivative and jacobian
            (N, B, j) = femFunctions.shapeFunction(self.xy, gp)

            area += gp[0] * j
            Qx += gp[0] * np.dot(N, np.transpose(self.xy[1,:])) * j
            Qy += gp[0] * np.dot(N, np.transpose(self.xy[0,:])) * j

        return (area, Qx, Qy)

    def geometricProperties(self):
        '''
        This function calculates the geometric propeties of the element.
        '''
        # initialise properties
        area = 0
        Qx = 0
        Qy = 0
        ixx = 0
        iyy = 0
        ixy = 0

        # Gauss points for 6 point Gaussian integration
        gps = femFunctions.gaussPoints(6)

        for gp in gps:
            # determine shape function, shape function derivative and jacobian
            (N, B, j) = femFunctions.shapeFunction(self.xy, gp)

            area += gp[0] * j
            Qx += gp[0] * np.dot(N, np.transpose(self.xy[1,:])) * j
            Qy += gp[0] * np.dot(N, np.transpose(self.xy[0,:])) * j
            ixx += gp[0] * (np.dot(N, np.transpose(self.xy[1,:]))) ** 2 * j
            iyy += gp[0] * (np.dot(N, np.transpose(self.xy[0,:]))) ** 2 * j
            ixy += (gp[0] * np.dot(N, np.transpose(self.xy[1,:])) *
                np.dot(N, np.transpose(self.xy[0,:])) * j)

        return (area, Qx, Qy, ixx, iyy, ixy)

    def torsionProperties(self):
        '''
        The function calculates the element stiffness matrix used for warping
        analysis and the torsion load vector.
        '''
        # initialise properties
        shearKe = 0
        torsionFe = 0

        # Gauss points for 6 point Gaussian integration
        gps = femFunctions.gaussPoints(6)

        for gp in gps:
            # determine shape function, shape function derivative and jacobian
            (N, B, j) = femFunctions.shapeFunction(self.xy, gp)

            # determine x and y position at Gauss point
            Nx = np.dot(N, np.transpose(self.xy[0,:]))
            Ny = np.dot(N, np.transpose(self.xy[1,:]))

            shearKe += gp[0] * np.dot(np.transpose(B), B) * j
            torsionFe += (gp[0] * np.dot(np.transpose(B),
                    np.transpose(np.array([Ny, -Nx]))) * j)

        return (shearKe, torsionFe)

    def shearProperties(self, ixx, iyy, ixy, omega):
        '''
        This function calculates the shear load vectors and also integrals used
        for shear analysis.
        '''
        # initialise properties
        shearFPsi = 0
        shearFPhi = 0
        shearCentreXInt = 0
        shearCentreYInt = 0
        Q_omega = 0
        i_omega = 0
        i_xomega = 0
        i_yomega = 0

        # Gauss points for 6 point Gaussian integration
        gps = femFunctions.gaussPoints(6)

        for gp in gps:
            # determine shape function, shape function derivative and jacobian
            (N, B, j) = femFunctions.shapeFunction(self.xy, gp)

            # determine x and y position at Gauss point
            Nx = np.dot(N, np.transpose(self.xy[0,:]))
            Ny = np.dot(N, np.transpose(self.xy[1,:]))

            # determine shear parameters
            r = Nx ** 2 - Ny ** 2
            q = 2 * Nx * Ny
            d1 = ixx * r - ixy * q
            d2 = ixy * r + ixx * q
            h1 = -ixy * r + iyy * q
            h2 = -iyy * r - ixy * q
            Nomega = np.dot(N, np.transpose(omega))

            shearFPsi += (gp[0] * (self.nu / 2 * np.transpose(
                np.transpose(B).dot(np.array([[d1], [d2]])))[0] + 2 * (1 +
                self.nu) *np.transpose(N) * (ixx * Nx - ixy * Ny)) * j)
            shearFPhi += (gp[0] * (self.nu / 2 * np.transpose(
                np.transpose(B).dot(np.array([[h1], [h2]])))[0] + 2 * (1 +
                self.nu) *np.transpose(N) * (iyy * Ny - ixy * Nx)) * j)
            shearCentreXInt += (gp[0] * (iyy * Nx + ixy * Ny) * (Nx ** 2 +
                Ny ** 2) * j)
            shearCentreYInt += (gp[0] * (ixx * Ny + ixy * Nx) * (Nx ** 2 +
                Ny ** 2) * j)
            Q_omega += gp[0] * Nomega * j
            i_omega += gp[0] * Nomega ** 2 * j
            i_xomega += gp[0] * Nx * Nomega * j
            i_yomega += gp[0] * Ny * Nomega * j

        return ((shearFPsi, shearFPhi, shearCentreXInt, shearCentreYInt,
            Q_omega, i_omega, i_xomega, i_yomega))

    def shearCoefficients(self, ixx, iyy, ixy, Psi, Phi):
        '''
        This functions calculates the variables used to determine the shear
        deformation coefficients.
        '''
        # initialise properties
        kappa_x = 0
        kappa_y = 0
        kappa_xy = 0

        # Gauss points for 6 point Gaussian integration
        gps = femFunctions.gaussPoints(6)

        for gp in gps:
            # determine shape function, shape function derivative and jacobian
            (N, B, j) = femFunctions.shapeFunction(self.xy, gp)

            # determine x and y position at Gauss point
            Nx = np.dot(N, np.transpose(self.xy[0,:]))
            Ny = np.dot(N, np.transpose(self.xy[1,:]))

            # determine shear parameters
            r = Nx ** 2 - Ny ** 2
            q = 2 * Nx * Ny
            d1 = ixx * r - ixy * q
            d2 = ixy * r + ixx * q
            h1 = -ixy * r + iyy * q
            h2 = -iyy * r - ixy * q

            kappa_x += (gp[0] * (Psi.dot(np.transpose(B)) - self.nu / 2 *
                np.array([d1, d2])).dot(B.dot(Psi) - self.nu / 2 *
                np.array([d1, d2])) * j)
            kappa_y += (gp[0] * (Phi.dot(np.transpose(B)) - self.nu / 2 *
                np.array([h1, h2])).dot(B.dot(Phi) - self.nu / 2 *
                np.array([h1, h2])) * j)
            kappa_xy += (gp[0] * (Psi.dot(np.transpose(B)) - self.nu / 2 *
                np.array([d1, d2])).dot(B.dot(Phi) - self.nu / 2 *
                np.array([h1, h2])) * j)

        return (kappa_x, kappa_y, kappa_xy)

    def calculateStress(self, area, ixx, iyy, ixy, i11, i22, phi, omega, J, Psi,
        Phi, Delta_s):
        '''
        This function calculates stresses within ana element as a result of unit
        loading.
        '''
        # initialise stress vectors
        sigma_zz_axial = np.ones((6,1)) * 1 / area
        sigma_zz_bending_xx_gp = np.zeros((6,1))
        sigma_zz_bending_yy_gp = np.zeros((6,1))
        sigma_zz_bending_11_gp = np.zeros((6,1))
        sigma_zz_bending_22_gp = np.zeros((6,1))
        tau_torsion_gp = np.zeros((6,2))
        tau_shear_x_gp = np.zeros((6,2))
        tau_shear_y_gp = np.zeros((6,2))

        # Gauss points for 6 point Gaussian integration
        gps = femFunctions.gaussPoints(6)

        for (i, gp) in enumerate(gps):
            # determine shape function, shape function derivative and jacobian
            (N, B, j) = femFunctions.shapeFunction(self.xy, gp)

            # determine x and y position at Gauss point
            Nx = np.dot(N, np.transpose(self.xy[0,:]))
            Ny = np.dot(N, np.transpose(self.xy[1,:]))

            # determine 11 and 22 position at Gauss point
            (Nx_1, Ny_2) = principalCoordinate(phi, Nx, Ny)

            # determine shear parameters
            r = Nx ** 2 - Ny ** 2
            q = 2 * Nx * Ny
            d1 = ixx * r - ixy * q
            d2 = ixy * r + ixx * q
            h1 = -ixy * r + iyy * q
            h2 = -iyy * r - ixy * q

            sigma_zz_bending_xx_gp[i,:] = (-(ixy * 1) / (ixx * iyy - ixy ** 2) *
                Nx + (iyy * 1) / (ixx * iyy - ixy ** 2) * Ny)
            sigma_zz_bending_yy_gp[i,:] = (-(ixx * 1) / (ixx * iyy - ixy ** 2) *
                Nx + (ixy * 1) / (ixx * iyy - ixy ** 2) * Ny)
            sigma_zz_bending_11_gp[i,:] = 1 / i11 * Ny_2
            sigma_zz_bending_22_gp[i,:] = -1 / i22 * Nx_1
            tau_torsion_gp[i,:] = (1 / J * (B.dot(omega) - np.array([Ny, -Nx])))
            tau_shear_x_gp[i,:] = (1 / Delta_s * (B.dot(Psi) - self.nu / 2 *
                np.array([d1, d2])))
            tau_shear_y_gp[i,:] = (1 / Delta_s * (B.dot(Phi) - self.nu / 2 *
                np.array([h1, h2])))

        # extrapolate results to nodes
        sigma_zz_bending_xx = (femFunctions.extrapolateToNodes(
            sigma_zz_bending_xx_gp[:,0]))
        sigma_zz_bending_yy = (femFunctions.extrapolateToNodes(
            sigma_zz_bending_yy_gp[:,0]))
        sigma_zz_bending_11 = (femFunctions.extrapolateToNodes(
            sigma_zz_bending_11_gp[:,0]))
        sigma_zz_bending_22 = (femFunctions.extrapolateToNodes(
            sigma_zz_bending_22_gp[:,0]))
        tau_zx_torsion = femFunctions.extrapolateToNodes(tau_torsion_gp[:,0])
        tau_zy_torsion = femFunctions.extrapolateToNodes(tau_torsion_gp[:,1])
        tau_shear_zx_x = femFunctions.extrapolateToNodes(tau_shear_x_gp[:,0])
        tau_shear_zy_x = femFunctions.extrapolateToNodes(tau_shear_x_gp[:,1])
        tau_shear_zx_y = femFunctions.extrapolateToNodes(tau_shear_y_gp[:,0])
        tau_shear_zy_y = femFunctions.extrapolateToNodes(tau_shear_y_gp[:,1])

        return ((sigma_zz_axial, sigma_zz_bending_xx, sigma_zz_bending_yy,
            sigma_zz_bending_11, sigma_zz_bending_22, tau_zx_torsion, tau_zy_torsion,
            tau_shear_zx_x, tau_shear_zy_x, tau_shear_zx_y, tau_shear_zy_y))

# ------------------------------------------------------------------------------
# GENERAL METHODS:
# ------------------------------------------------------------------------------

def principalCoordinate(phi, x, y):
    '''
    This function determines the coordinates of the point (x,y) in the principal
    axis system given rotation phi.
    '''
    # convert principal axis angle to radians
    phi_rad = phi * np.pi / 180
    # form rotation matrix
    R = (np.array([[np.cos(phi_rad), np.sin(phi_rad)], [-np.sin(phi_rad),
        np.cos(phi_rad)]]))
    # calculate rotated x and y coordinates
    x_rotated = R.dot(np.array([x,y]))

    return (x_rotated[0], x_rotated[1])
