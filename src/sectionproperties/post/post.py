"""Post processing routines."""

import contextlib

import matplotlib.pyplot as plt
import numpy as np


@contextlib.contextmanager
def plotting_context(
    ax=None, pause=True, title='', filename='', render=True, axis_index=None, **kwargs
):
    r"""Executes code required to set up a matplotlib figure.

    Parameters
    ----------
    ax : :class:`matplotlib.axes.Axes`
        Axes object on which to plot
    pause : bool
        If set to true, the figure pauses the script until the window is closed. If set to false,
        the script continues immediately after the window is rendered.
    title : str
        Plot title
    filename : str
        Pass a non-empty string or path to save the image as. If this option is used, the figure is
        closed after the file is saved.
    render : bool
        If set to False, the image is not popped up. This may be useful if the figure or axes will
        be embedded or further edited before being displayed.
    axis_index : Union[None, int, Tuple(int)]
        If more than 1 axes is created by subplot, then this is the axis to plot on. This may be a
        tuple if a 2D array of plots is returned.  The default value of None will select the top
        left plot.
    \**kwargs
        Passed to :func:`matplotlib.pyplot.subplots`
    """
    if filename:
        render = False
    if ax is None:
        if not render:
            plt.ioff()
        elif pause:
            plt.ioff()
        else:
            plt.ion()
        ax_supplied = False
        (fig, ax) = plt.subplots(**kwargs)
        try:
            if axis_index is None:
                axis_index = (0,) * ax.ndim
            ax = ax[axis_index]
        except (AttributeError, TypeError):
            pass  # only 1 axis, not an array
        except IndexError as exc:
            raise ValueError(
                f'axis_index={axis_index} is not compatible with arguments to subplots: {kwargs}'
            ) from exc
    else:
        fig = ax.get_figure()
        ax_supplied = True
        if not render:
            plt.ioff()

    yield fig, ax
    ax.set_aspect('equal', anchor='C')

    if ax_supplied:
        # if an axis was supplied, don't continue with displaying or configuring the plot
        return

    # if no axes was supplied, finish the plot and return the figure and axes
    ax.set_title(title)
    plt.tight_layout()

    if filename:
        fig.savefig(filename, dpi=fig.dpi)
        plt.close(fig)  # close the figure to free the memory
        return  # if the figure was to be saved, then don't show it also

    if render:
        if pause:
            plt.show()
        else:
            plt.draw()
            plt.pause(0.001)


def draw_principal_axis(ax, phi, cx, cy):
    """Draws the principal axis on a plot.

    Parameters
    ----------
    ax : :class:`matplotlib.axes.Axes`
        Axes object on which to plot
    phi : float
        Principal axis angle in radians
    cx : float
        x-location of the centroid
    cy : float
        y-location of the centroid
    """
    # get current axis limits
    (xmin, xmax) = ax.get_xlim()
    (ymin, ymax) = ax.get_ylim()
    lims = [xmin, xmax, ymin, ymax]

    # form rotation matrix
    R = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])

    # get basis vectors in the directions of the principal axes
    x11_basis = R.dot(np.array([1, 0]))
    y22_basis = R.dot(np.array([0, 1]))

    def add_point(vec, basis, centroid, num, denom):
        """Adds a point to the list *vec* if there is an intersection."""
        if denom != 0:
            point = basis * num / denom + centroid
            vec.append([point[0], point[1]])

    def get_principal_points(basis, lims, centroid):
        """Return points to plot the principal axes.

        Determines the intersections of the principal axis with the four lines defining a
        bounding box around the limits of the cross-section. The middle two intersection points are
        returned for plotting.

        Parameters
        ----------
        basis : :class:`numpy.ndarray`
            Basis (unit) vector in the direction of the principal axis
        lims : tuple(float, float, float, float)
            Tuple containing the axis limits *(xmin, xmax, ymin, ymax)*
        centroid : list[float, float]
            Centroid *(cx, cy)* of the cross-section, through which the principal axis passes
        """
        pts = []  # initialise list containing the intersection points

        # add intersection points to the list
        add_point(pts, basis, centroid, lims[0] - centroid[0], basis[0])
        add_point(pts, basis, centroid, lims[1] - centroid[0], basis[0])
        add_point(pts, basis, centroid, lims[2] - centroid[1], basis[1])
        add_point(pts, basis, centroid, lims[3] - centroid[1], basis[1])

        # sort point vector
        pts = np.array(pts)
        pts = pts[pts[:, 0].argsort()]  # stackoverflow sort numpy array by col

        # if there are four points, take the middle two points
        if len(pts) == 4:
            return pts[1:3, :]

        return pts

    # get intersection points for the 11 and 22 axes
    x11 = get_principal_points(x11_basis, lims, [cx, cy])
    y22 = get_principal_points(y22_basis, lims, [cx, cy])

    # plot the principal axis
    ax.plot(x11[:, 0], x11[:, 1], 'k--', alpha=0.5, label='11-axis')
    ax.plot(y22[:, 0], y22[:, 1], 'k-.', alpha=0.5, label='22-axis')


def print_results(cross_section, fmt):
    """Prints the results that have been calculated to the terminal.

    Parameters
    ----------
    cross_section : :class:`~sectionproperties.analysis.cross_section.CrossSection`
        Structural cross-section object
    fmt : str
        Number format
    """
    if cross_section.materials is not None:
        prefix = "E."
    else:
        prefix = ""

    def sprint(label, value):
        if value >= 0:
            print(f'{label:<7} =  {value:>{fmt}}')
        else:
            print(f'{label:<7} = {value:>{fmt}}')

    area = cross_section.get_area()
    if area is not None:
        print("Section Properties:")
        sprint("A", area)

    perimeter = cross_section.get_perimeter()
    if perimeter is not None:
        sprint("Perim.", perimeter)

    if cross_section.materials is not None:
        ea = cross_section.get_ea()
        if ea is not None:
            sprint("E.A", ea)

    (qx, qy) = cross_section.get_q()
    if qx is not None:
        sprint(prefix + "Qx", qx)
        sprint(prefix + "Qy", qy)

    (cx, cy) = cross_section.get_c()
    if cx is not None:
        sprint("cx", cx)
        sprint("cy", cy)

    (ixx_g, iyy_g, ixy_g) = cross_section.get_ig()
    if ixx_g is not None:
        sprint(prefix + "Ixx_g", ixx_g)
        sprint(prefix + "Iyy_g", iyy_g)
        sprint(prefix + "Ixy_g", ixy_g)

    (ixx_c, iyy_c, ixy_c) = cross_section.get_ic()
    if ixx_c is not None:
        sprint(prefix + "Ixx_c", ixx_c)
        sprint(prefix + "Iyy_c", iyy_c)
        sprint(prefix + "Ixy_c", ixy_c)

    (zxx_plus, zxx_minus, zyy_plus, zyy_minus) = cross_section.get_z()
    if zxx_plus is not None:
        sprint(prefix + "Zxx+", zxx_plus)
        sprint(prefix + "Zxx-", zxx_minus)
        sprint(prefix + "Zyy+", zyy_plus)
        sprint(prefix + "Zyy-", zyy_minus)

    (rx, ry) = cross_section.get_rc()
    if rx is not None:
        sprint("rx", rx)
        sprint("ry", ry)

    phi = cross_section.get_phi()
    (i11_c, i22_c) = cross_section.get_ip()
    if phi is not None:
        sprint("phi", phi)
        sprint(prefix + "I11_c", i11_c)
        sprint(prefix + "I22_c", i22_c)

    (z11_plus, z11_minus, z22_plus, z22_minus) = cross_section.get_zp()
    if z11_plus is not None:
        sprint(prefix + "Z11+", z11_plus)
        sprint(prefix + "Z11-", z11_minus)
        sprint(prefix + "Z22+", z22_plus)
        sprint(prefix + "Z22-", z22_minus)

    (r11, r22) = cross_section.get_rp()
    if r11 is not None:
        sprint("r11", r11)
        sprint("r22", r22)

    j = cross_section.get_j()
    if j is not None:
        if cross_section.materials is not None:
            sprint("G.J", j / (2 * (1 + cross_section.section_props.nu_eff)))
        else:
            sprint("J", j)

    gamma = cross_section.get_gamma()
    if gamma is not None:
        if cross_section.materials is not None:
            sprint("G.Iw", gamma / (2 * (1 + cross_section.section_props.nu_eff)))
        else:
            sprint("Iw", gamma)

    (x_se, y_se) = cross_section.get_sc()
    if x_se is not None:
        sprint("x_se", x_se)
        sprint("y_se", y_se)

    (x_st, y_st) = cross_section.get_sc_t()
    if x_se is not None:
        sprint("x_st", x_st)
        sprint("y_st", y_st)

    (x1_se, y2_se) = cross_section.get_sc_p()
    if x1_se is not None:
        sprint("x1_se", x1_se)
        sprint("y2_se", y2_se)

    (A_sx, A_sy) = cross_section.get_As()
    if A_sx is not None:
        if cross_section.materials is not None:
            sprint("A_sx", A_sx * cross_section.section_props.area / cross_section.section_props.ea)
            sprint("A_sy", A_sy * cross_section.section_props.area / cross_section.section_props.ea)
        else:
            sprint("A_sx", A_sx)
            sprint("A_sy", A_sy)

    (A_s11, A_s22) = cross_section.get_As_p()
    if A_s11 is not None:
        if cross_section.materials is not None:
            sprint(
                "A_s11", A_s11 * cross_section.section_props.area / cross_section.section_props.ea
            )
            sprint(
                "A_s22", A_s22 * cross_section.section_props.area / cross_section.section_props.ea
            )
        else:
            sprint("A_s11", A_s11)
            sprint("A_s22", A_s22)

    (beta_x_plus, beta_x_minus, beta_y_plus, beta_y_minus) = cross_section.get_beta()
    if beta_x_plus is not None:
        sprint("betax+", beta_x_plus)
        sprint("betax-", beta_x_minus)
        sprint("betay+", beta_y_plus)
        sprint("betay-", beta_y_minus)

    (beta_11_plus, beta_11_minus, beta_22_plus, beta_22_minus) = cross_section.get_beta_p()
    if beta_x_plus is not None:
        sprint("beta11+", beta_11_plus)
        sprint("beta11-", beta_11_minus)
        sprint("beta22+", beta_22_plus)
        sprint("beta22-", beta_22_minus)

    (x_pc, y_pc) = cross_section.get_pc()
    if x_pc is not None:
        sprint("x_pc", x_pc)
        sprint("y_pc", y_pc)

    (sxx, syy) = cross_section.get_s()
    if sxx is not None:
        if cross_section.materials is not None:
            sprint("M_p,xx", sxx)
            sprint("M_p,yy", syy)
        else:
            sprint("Sxx", sxx)
            sprint("Syy", syy)

    (sf_xx_plus, sf_xx_minus, sf_yy_plus, sf_yy_minus) = cross_section.get_sf()
    if sf_xx_plus is not None:
        sprint("SF_xx+", sf_xx_plus)
        sprint("SF_xx-", sf_xx_minus)
        sprint("SF_yy+", sf_yy_plus)
        sprint("SF_yy-", sf_yy_minus)

    (x11_pc, y22_pc) = cross_section.get_pc_p()
    if x_pc is not None:
        sprint("x11_pc", x11_pc)
        sprint("y22_pc", y22_pc)

    (s11, s22) = cross_section.get_sp()
    if s11 is not None:
        if cross_section.materials is not None:
            sprint("M_p,11", s11)
            sprint("M_p,22", s22)
        else:
            sprint("S11", s11)
            sprint("S22", s22)

    (sf_11_plus, sf_11_minus, sf_22_plus, sf_22_minus) = cross_section.get_sf_p()
    if sf_11_plus is not None:
        sprint("SF_11+", sf_11_plus)
        sprint("SF_11-", sf_11_minus)
        sprint("SF_22+", sf_22_plus)
        sprint("SF_22-", sf_22_minus)

    print("")
