"""Post-processor methods and classes."""

from __future__ import annotations

import contextlib
from dataclasses import asdict, dataclass
from typing import TYPE_CHECKING

import matplotlib.axes
import matplotlib.pyplot as plt
import numpy as np
from rich.console import Console
from rich.table import Table

from sectionproperties.analysis.fea import principal_coordinate
from sectionproperties.pre.pre import DEFAULT_MATERIAL


if TYPE_CHECKING:
    from sectionproperties.analysis.section import Section


@dataclass
class SectionProperties:
    """Class for storing section properties.

    Stores calculated section properties. Also provides methods to calculate section
    properties entirely derived from other section properties.

    Attributes:
        area: Cross-sectional area
        perimeter: Cross-sectional perimeter
        mass: Cross-sectional mass
        ea: Modulus weighted area (axial rigidity)
        ga: Modulus weighted product of shear modulus and area
        nu_eff: Effective Poisson's ratio
        e_eff: Effective elastic modulus
        g_eff: Effective shear modulus
        qx: First moment of area about the x-axis
        qy: First moment of area about the y-axis
        ixx_g: Second moment of area about the global x-axis
        iyy_g: Second moment of area about the global y-axis
        ixy_g: Second moment of area about the global xy-axis
        cx: X coordinate of the elastic centroid
        cy: Y coordinate of the elastic centroid
        ixx_c: Second moment of area about the centroidal x-axis
        iyy_c: Second moment of area about the centroidal y-axis
        ixy_c: Second moment of area about the centroidal xy-axis
        zxx_plus: Section modulus about the centroidal x-axis for stresses at the
            positive extreme value of y
        zxx_minus: Section modulus about the centroidal x-axis for stresses at the
            negative extreme value of y
        zyy_plus: Section modulus about the centroidal y-axis for stresses at the
            positive extreme value of x
        zyy_minus: Section modulus about the centroidal y-axis for stresses at the
            negative extreme value of x
        rx_c: Radius of gyration about the centroidal x-axis.
        ry_c: Radius of gyration about the centroidal y-axis.
        i11_c: Second moment of area about the centroidal 11-axis
        i22_c: Second moment of area about the centroidal 22-axis
        phi: Principal axis angle
        z11_plus: Section modulus about the principal 11-axis for stresses at the
            positive extreme value of the 22-axis
        z11_minus: Section modulus about the principal 11-axis for stresses at the
            negative extreme value of the 22-axis
        z22_plus: Section modulus about the principal 22-axis for stresses at the
            positive extreme value of the 11-axis
        z22_minus: Section modulus about the principal 22-axis for stresses at the
            negative extreme value of the 11-axis
        r11_c: Radius of gyration about the principal 11-axis.
        r22_c: Radius of gyration about the principal 22-axis.
        j: Torsion constant
        omega: Warping function
        psi_shear: Psi shear function
        phi_shear: Phi shear function
        Delta_s: Shear factor
        x_se: x-coordinate of the shear centre (elasticity approach)
        y_se: y-coordinate of the shear centre (elasticity approach)
        x11_se: 11-coordinate of the shear centre (elasticity approach)
        y22_se: 22-coordinate of the shear centre (elasticity approach)
        x_st: x-coordinate of the shear centre (Trefftz's approach)
        y_st: y-coordinate of the shear centre (Trefftz's approach)
        gamma: Warping constant
        A_sx: Shear area about the x-axis
        A_sy: Shear area about the y-axis
        A_sxy: Shear area about the xy-axis
        A_s11: Shear area about the 11-bending axis
        A_s22: Shear area about the 22-bending axis
        beta_x_plus: Monosymmetry constant for bending about the x-axis with the top
            flange in compression
        beta_x_minus: Monosymmetry constant for bending about the x-axis with the bottom
            flange in compression
        beta_y_plus: Monosymmetry constant for bending about the y-axis with the top
            flange in compression
        beta_y_minus: Monosymmetry constant for bending about the y-axis with the bottom
            flange in compression
        beta_11_plus: Monosymmetry constant for bending about the 11-axis with the top
            flange in compression
        beta_11_minus: Monosymmetry constant for bending about the 11-axis with the
            bottom flange in compression
        beta_22_plus: Monosymmetry constant for bending about the 22-axis with the top
            flange in compression
        beta_22_minus: Monosymmetry constant for bending about the 22-axis with the
            bottom flange in compression
        x_pc: x-coordinate of the global plastic centroid
        y_pc: y-coordinate of the global plastic centroid
        x11_pc: 11-coordinate of the principal plastic centroid
        y22_pc: 22-coordinate of the principal plastic centroid
        sxx: Plastic section modulus about the centroidal x-axis
        syy: Plastic section modulus about the centroidal y-axis
        sf_xx_plus: Shape factor for bending about the x-axis with respect to the top
            fibre
        sf_xx_minus: Shape factor for bending about the x-axis with respect to the
            bottom fibre
        sf_yy_plus: Shape factor for bending about the y-axis with respect to the top
            fibre
        sf_yy_minus: Shape factor for bending about the y-axis with respect to the
            bottom fibre
        s11: Plastic section modulus about the 11-axis
        s22: Plastic section modulus about the 22-axis
        sf_11_plus: Shape factor for bending about the 11-axis with respect to the top
            fibre
        sf_11_minus: Shape factor for bending about the 11-axis with respect to the
            bottom fibre
        sf_22_plus: Shape factor for bending about the 22-axis with respect to the top
            fibre
        sf_22_minus: Shape factor for bending about the 22-axis with respect to the
            bottom fibre
    """

    area: float | None = None
    perimeter: float | None = None
    mass: float | None = None
    ea: float | None = None
    ga: float | None = None
    nu_eff: float | None = None
    e_eff: float | None = None
    g_eff: float | None = None
    qx: float | None = None
    qy: float | None = None
    ixx_g: float | None = None
    iyy_g: float | None = None
    ixy_g: float | None = None
    cx: float | None = None
    cy: float | None = None
    ixx_c: float | None = None
    iyy_c: float | None = None
    ixy_c: float | None = None
    zxx_plus: float | None = None
    zxx_minus: float | None = None
    zyy_plus: float | None = None
    zyy_minus: float | None = None
    rx_c: float | None = None
    ry_c: float | None = None
    i11_c: float | None = None
    i22_c: float | None = None
    phi: float | None = None
    z11_plus: float | None = None
    z11_minus: float | None = None
    z22_plus: float | None = None
    z22_minus: float | None = None
    r11_c: float | None = None
    r22_c: float | None = None
    j: float | None = None
    omega: np.ndarray | None = None
    psi_shear: np.ndarray | None = None
    phi_shear: np.ndarray | None = None
    Delta_s: float | None = None
    x_se: float | None = None
    y_se: float | None = None
    x11_se: float | None = None
    y22_se: float | None = None
    x_st: float | None = None
    y_st: float | None = None
    gamma: float | None = None
    A_sx: float | None = None
    A_sy: float | None = None
    A_sxy: float | None = None
    A_s11: float | None = None
    A_s22: float | None = None
    beta_x_plus: float | None = None
    beta_x_minus: float | None = None
    beta_y_plus: float | None = None
    beta_y_minus: float | None = None
    beta_11_plus: float | None = None
    beta_11_minus: float | None = None
    beta_22_plus: float | None = None
    beta_22_minus: float | None = None
    x_pc: float | None = None
    y_pc: float | None = None
    x11_pc: float | None = None
    y22_pc: float | None = None
    sxx: float | None = None
    syy: float | None = None
    sf_xx_plus: float | None = None
    sf_xx_minus: float | None = None
    sf_yy_plus: float | None = None
    sf_yy_minus: float | None = None
    s11: float | None = None
    s22: float | None = None
    sf_11_plus: float | None = None
    sf_11_minus: float | None = None
    sf_22_plus: float | None = None
    sf_22_minus: float | None = None

    def asdict(self) -> dict:
        """Returns the SectionProperties dataclass object as a dictionary.

        Return:
            Dictionary of the SectionProperties class
        """
        return asdict(self)

    def calculate_elastic_centroid(self) -> None:
        """Calculates and stores the elastic centroid."""
        if self.qx and self.qy and self.ea:
            self.cx = self.qy / self.ea
            self.cy = self.qx / self.ea
        else:
            raise RuntimeError("Calculate geometric properties first.")

    def calculate_centroidal_properties(
        self,
        node_list: list[list[float]],
    ) -> None:
        """Calculates and stores derived geometric properties.

        Args:
            node_list: List of mesh node coordinates
        """
        # calculate second moments of area about the centroidal xy axis
        if self.qx and self.qy and self.ea and self.ixx_g and self.iyy_g and self.ixy_g:
            self.ixx_c = self.ixx_g - self.qx**2 / self.ea
            self.iyy_c = self.iyy_g - self.qy**2 / self.ea
            self.ixy_c = self.ixy_g - self.qx * self.qy / self.ea

            # calculate section moduli about the centroidal xy axis
            nodes = np.array(node_list, dtype=float)
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
            delta = (((self.ixx_c - self.iyy_c) / 2) ** 2 + self.ixy_c**2) ** 0.5
            self.i11_c = (self.ixx_c + self.iyy_c) / 2 + delta
            self.i22_c = (self.ixx_c + self.iyy_c) / 2 - delta

            # calculate initial principal axis angle
            if abs(self.ixx_c - self.i11_c) < 1e-12 * self.i11_c:
                self.phi = 0
            else:
                self.phi = np.arctan2(self.ixx_c - self.i11_c, self.ixy_c) * 180 / np.pi

            # initialise min, max variables
            x1, y2 = principal_coordinate(
                phi=self.phi,
                x=nodes[0][0] - self.cx,
                y=nodes[0][1] - self.cy,
            )
            x1max = x1
            x1min = x1
            y2max = y2
            y2min = y2

            # calculate section moduli about the principal axis
            for pt in nodes[1:]:
                x = pt[0] - self.cx
                y = pt[1] - self.cy

                # determine the coordinate of the point wrt the principal axis
                x1, y2 = principal_coordinate(phi=self.phi, x=x, y=y)

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
        else:
            raise RuntimeError("Calculate geometric properties first.")


@contextlib.contextmanager
def plotting_context(
    ax: matplotlib.axes.Axes | None = None,
    pause: bool = True,
    title: str = "",
    filename: str = "",
    render: bool = True,
    axis_index: int | tuple[int, int] | None = None,
    **kwargs,
):
    """Executes code required to set up a matplotlib figure.

    Args:
        ax: Axes object on which to plot
        pause: If set to True, the figure pauses the script until the window is closed.
            If set to False, the script continues immediately after the window is
            rendered.
        title: Plot title
        filename: Pass a non-empty string or path to save the image as. If this option
            is used, the figure is closed after the file is saved.
        render: If set to False, the image is not displayed. This may be useful if the
            figure or axes will be embedded or further edited before being displayed.
        axis_index: If more than 1 axis is created by subplot, then this is the axis to
            plot on. This may be a tuple if a 2D array of plots is returned.  The
            default value of None will select the top left plot.
        kwargs: Passed to :func:`matplotlib.pyplot.subplots`

    Yields:
        Matplotlib figure and axes
    """
    if filename:
        render = False

    if ax is None:
        if not render or pause:
            plt.ioff()
        else:
            plt.ion()

        ax_supplied = False
        fig, ax = plt.subplots(**kwargs)

        try:
            if axis_index is None:
                axis_index = (0,) * ax.ndim

            ax = ax[axis_index]
        except (AttributeError, TypeError):
            pass  # only 1 axis, not an array
        except IndexError as exc:
            msg = f"axis_index={axis_index} is not compatible "
            msg += f"with arguments to subplots: {kwargs}"
            raise ValueError(msg) from exc
    else:
        fig = ax.get_figure()
        ax_supplied = True

        if not render:
            plt.ioff()

    yield fig, ax

    if ax is not None:
        ax.set_title(title)
        plt.tight_layout()
        ax.set_aspect("equal", anchor="C")

    # if no axes was supplied, finish the plot and return the figure and axes
    if ax_supplied:
        # if an axis was supplied, don't continue displaying or configuring the plot
        return

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


def draw_principal_axis(
    ax: matplotlib.axes.Axes,
    phi: float,
    cx: float,
    cy: float,
) -> None:
    """Draws the principal axis on a plot.

    Args:
        ax: Axes object on which to plot
        phi: Principal axis angle in radians
        cx: x-location of the centroid
        cy: y-location of the centroid
    """
    # get current axis limits
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    lims = (xmin, xmax, ymin, ymax)

    # form rotation matrix
    r = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])

    # get basis vectors in the directions of the principal axes
    x11_basis = r.dot(np.array([1, 0]))
    y22_basis = r.dot(np.array([0, 1]))

    def add_point(
        vec: list[list[float]],
        basis: np.ndarray,
        centroid: tuple[float, float],
        num: float,
        denom: float,
    ) -> None:
        """Adds a point to the list ``vec`` if there is an intersection.

        Args:
            vec: List of points to add to
            basis: Basis vector of principal axis
            centroid: Geometry centroid
            num: Numberator
            denom: Denominator
        """
        if denom != 0:
            point = basis * num / denom + centroid
            vec.append([point[0], point[1]])

    def get_principal_points(
        basis: np.ndarray,
        lims: tuple[float, float, float, float],
        centroid: tuple[float, float],
    ) -> np.ndarray:
        """Returns intersection points of prinicipal axis with bounding box.

        Determines the intersections of the principal axis with the four lines
        defining a bounding box around the limits of the cross-section. The middle two
        intersection points are returned for plotting.

        Args:
            basis: Basis (unit) vector in the direction of the principal axis
            lims: Tuple containing the axis limits ``(xmin, xmax, ymin, ymax)``
            centroid: Centroid ``(cx, cy)`` of the cross-section, through which the
                principal axis passes

        Return:
            List of intersection points
        """
        pts = []  # initialise list containing the intersection points

        # add intersection points to the list
        add_point(
            vec=pts,
            basis=basis,
            centroid=centroid,
            num=lims[0] - centroid[0],
            denom=basis[0],
        )
        add_point(
            vec=pts,
            basis=basis,
            centroid=centroid,
            num=lims[1] - centroid[0],
            denom=basis[0],
        )
        add_point(
            vec=pts,
            basis=basis,
            centroid=centroid,
            num=lims[2] - centroid[1],
            denom=basis[1],
        )
        add_point(
            vec=pts,
            basis=basis,
            centroid=centroid,
            num=lims[3] - centroid[1],
            denom=basis[1],
        )

        # sort point vector
        pts = np.array(pts)
        pts = pts[pts[:, 0].argsort()]  # stackoverflow sort numpy array by col

        # if there are four points, take the middle two points
        if len(pts) == 4:
            return pts[1:3, :]

        return pts

    # get intersection points for the 11 and 22 axes
    x11 = get_principal_points(
        basis=x11_basis,
        lims=lims,
        centroid=(cx, cy),
    )
    y22 = get_principal_points(
        basis=y22_basis,
        lims=lims,
        centroid=(cx, cy),
    )

    # plot the principal axis
    ax.plot(x11[:, 0], x11[:, 1], "k--", alpha=0.5, label="11-axis")
    ax.plot(y22[:, 0], y22[:, 1], "k-.", alpha=0.5, label="22-axis")


def print_results(
    cross_section: Section,
    fmt: str,
) -> None:
    """Prints the results that have been calculated to the terminal.

    Args:
        cross_section: Section object
        fmt: Number formatting string
    """
    if list(set(cross_section.materials)) != [DEFAULT_MATERIAL]:
        prefix = "E."
    else:
        prefix = ""

    table = Table(title="Section Properties")
    table.add_column("Property", justify="left", style="cyan", no_wrap=True)
    table.add_column("Value", justify="right", style="green")

    area = cross_section.get_area()
    if area is not None:
        table.add_row("A", "{:>{fmt}}".format(area, fmt=fmt))

    perimeter = cross_section.get_perimeter()
    if perimeter is not None:
        table.add_row("Perim.", "{:>{fmt}}".format(perimeter, fmt=fmt))

    if list(set(cross_section.materials)) != [DEFAULT_MATERIAL]:
        mass = cross_section.get_mass()
        if mass is not None:
            table.add_row("Mass", "{:>{fmt}}".format(mass, fmt=fmt))

    if list(set(cross_section.materials)) != [DEFAULT_MATERIAL]:
        ea = cross_section.get_ea()
        if ea is not None:
            table.add_row("E.A", "{:>{fmt}}".format(ea, fmt=fmt))

    (qx, qy) = cross_section.get_q()
    if qx is not None:
        table.add_row(prefix + "Qx", "{:>{fmt}}".format(qx, fmt=fmt))
        table.add_row(prefix + "Qy", "{:>{fmt}}".format(qy, fmt=fmt))

    (cx, cy) = cross_section.get_c()
    if cx is not None:
        table.add_row("cx", "{:>{fmt}}".format(cx, fmt=fmt))
        table.add_row("cy", "{:>{fmt}}".format(cy, fmt=fmt))

    (ixx_g, iyy_g, ixy_g) = cross_section.get_ig()
    if ixx_g is not None:
        table.add_row(prefix + "Ixx_g", "{:>{fmt}}".format(ixx_g, fmt=fmt))
        table.add_row(prefix + "Iyy_g", "{:>{fmt}}".format(iyy_g, fmt=fmt))
        table.add_row(prefix + "Ixy_g", "{:>{fmt}}".format(ixy_g, fmt=fmt))

    (ixx_c, iyy_c, ixy_c) = cross_section.get_ic()
    if ixx_c is not None:
        table.add_row(prefix + "Ixx_c", "{:>{fmt}}".format(ixx_c, fmt=fmt))
        table.add_row(prefix + "Iyy_c", "{:>{fmt}}".format(iyy_c, fmt=fmt))
        table.add_row(prefix + "Ixy_c", "{:>{fmt}}".format(ixy_c, fmt=fmt))

    (zxx_plus, zxx_minus, zyy_plus, zyy_minus) = cross_section.get_z()
    if zxx_plus is not None:
        table.add_row(prefix + "Zxx+", "{:>{fmt}}".format(zxx_plus, fmt=fmt))
        table.add_row(prefix + "Zxx-", "{:>{fmt}}".format(zxx_minus, fmt=fmt))
        table.add_row(prefix + "Zyy+", "{:>{fmt}}".format(zyy_plus, fmt=fmt))
        table.add_row(prefix + "Zyy-", "{:>{fmt}}".format(zyy_minus, fmt=fmt))

    (rx, ry) = cross_section.get_rc()
    if rx is not None:
        table.add_row("rx", "{:>{fmt}}".format(rx, fmt=fmt))
        table.add_row("ry", "{:>{fmt}}".format(ry, fmt=fmt))

    phi = cross_section.get_phi()
    (i11_c, i22_c) = cross_section.get_ip()
    if phi is not None:
        table.add_row("phi", "{:>{fmt}}".format(phi, fmt=fmt))
        table.add_row(prefix + "I11_c", "{:>{fmt}}".format(i11_c, fmt=fmt))
        table.add_row(prefix + "I22_c", "{:>{fmt}}".format(i22_c, fmt=fmt))

    (z11_plus, z11_minus, z22_plus, z22_minus) = cross_section.get_zp()
    if z11_plus is not None:
        table.add_row(prefix + "Z11+", "{:>{fmt}}".format(z11_plus, fmt=fmt))
        table.add_row(prefix + "Z11-", "{:>{fmt}}".format(z11_minus, fmt=fmt))
        table.add_row(prefix + "Z22+", "{:>{fmt}}".format(z22_plus, fmt=fmt))
        table.add_row(prefix + "Z22-", "{:>{fmt}}".format(z22_minus, fmt=fmt))

    (r11, r22) = cross_section.get_rp()
    if r11 is not None:
        table.add_row("r11", "{:>{fmt}}".format(r11, fmt=fmt))
        table.add_row("r22", "{:>{fmt}}".format(r22, fmt=fmt))

    if list(set(cross_section.materials)) != [DEFAULT_MATERIAL]:
        e_eff = cross_section.get_e_eff()
        g_eff = cross_section.get_g_eff()
        if e_eff is not None:
            table.add_row("E_eff", "{:>{fmt}}".format(e_eff, fmt=fmt))
            table.add_row("G_eff", "{:>{fmt}}".format(g_eff, fmt=fmt))

        nu_eff = cross_section.get_nu_eff()
        if nu_eff is not None:
            table.add_row("nu_eff", "{:>{fmt}}".format(nu_eff, fmt=fmt))

    j = cross_section.get_j()
    if j is not None:
        table.add_row("J", "{:>{fmt}}".format(j, fmt=fmt))

    gamma = cross_section.get_gamma()
    if gamma is not None:
        table.add_row("Iw", "{:>{fmt}}".format(gamma, fmt=fmt))

    (x_se, y_se) = cross_section.get_sc()
    if x_se is not None:
        table.add_row("x_se", "{:>{fmt}}".format(x_se, fmt=fmt))
        table.add_row("y_se", "{:>{fmt}}".format(y_se, fmt=fmt))

    (x_st, y_st) = cross_section.get_sc_t()
    if x_se is not None:
        table.add_row("x_st", "{:>{fmt}}".format(x_st, fmt=fmt))
        table.add_row("y_st", "{:>{fmt}}".format(y_st, fmt=fmt))

    (x1_se, y2_se) = cross_section.get_sc_p()
    if x1_se is not None:
        table.add_row("x1_se", "{:>{fmt}}".format(x1_se, fmt=fmt))
        table.add_row("y2_se", "{:>{fmt}}".format(y2_se, fmt=fmt))

    (a_sx, a_sy) = cross_section.get_As()
    if a_sx is not None:
        table.add_row(prefix + "A_sx", "{:>{fmt}}".format(a_sx, fmt=fmt))
        table.add_row(prefix + "A_sy", "{:>{fmt}}".format(a_sy, fmt=fmt))

    (a_s11, a_s22) = cross_section.get_As_p()
    if a_s11 is not None:
        table.add_row(prefix + "A_s11", "{:>{fmt}}".format(a_s11, fmt=fmt))
        table.add_row(prefix + "A_s22", "{:>{fmt}}".format(a_s22, fmt=fmt))

    (beta_x_plus, beta_x_minus, beta_y_plus, beta_y_minus) = cross_section.get_beta()
    if beta_x_plus is not None:
        table.add_row("betax+", "{:>{fmt}}".format(beta_x_plus, fmt=fmt))
        table.add_row("betax-", "{:>{fmt}}".format(beta_x_minus, fmt=fmt))
        table.add_row("betay+", "{:>{fmt}}".format(beta_y_plus, fmt=fmt))
        table.add_row("betay-", "{:>{fmt}}".format(beta_y_minus, fmt=fmt))

    (
        beta_11_plus,
        beta_11_minus,
        beta_22_plus,
        beta_22_minus,
    ) = cross_section.get_beta_p()
    if beta_x_plus is not None:
        table.add_row("beta11+", "{:>{fmt}}".format(beta_11_plus, fmt=fmt))
        table.add_row("beta11-", "{:>{fmt}}".format(beta_11_minus, fmt=fmt))
        table.add_row("beta22+", "{:>{fmt}}".format(beta_22_plus, fmt=fmt))
        table.add_row("beta22-", "{:>{fmt}}".format(beta_22_minus, fmt=fmt))

    (x_pc, y_pc) = cross_section.get_pc()
    if x_pc is not None:
        table.add_row("x_pc", "{:>{fmt}}".format(x_pc, fmt=fmt))
        table.add_row("y_pc", "{:>{fmt}}".format(y_pc, fmt=fmt))

    (sxx, syy) = cross_section.get_s()
    if sxx is not None:
        if list(set(cross_section.materials)) != [DEFAULT_MATERIAL]:
            table.add_row("M_p,xx", "{:>{fmt}}".format(sxx, fmt=fmt))
            table.add_row("M_p,yy", "{:>{fmt}}".format(syy, fmt=fmt))
        else:
            table.add_row("Sxx", "{:>{fmt}}".format(sxx, fmt=fmt))
            table.add_row("Syy", "{:>{fmt}}".format(syy, fmt=fmt))

    (sf_xx_plus, sf_xx_minus, sf_yy_plus, sf_yy_minus) = cross_section.get_sf()
    if sf_xx_plus is not None:
        table.add_row("SF_xx+", "{:>{fmt}}".format(sf_xx_plus, fmt=fmt))
        table.add_row("SF_xx-", "{:>{fmt}}".format(sf_xx_minus, fmt=fmt))
        table.add_row("SF_yy+", "{:>{fmt}}".format(sf_yy_plus, fmt=fmt))
        table.add_row("SF_yy-", "{:>{fmt}}".format(sf_yy_minus, fmt=fmt))

    (x11_pc, y22_pc) = cross_section.get_pc_p()
    if x_pc is not None:
        table.add_row("x11_pc", "{:>{fmt}}".format(x11_pc, fmt=fmt))
        table.add_row("y22_pc", "{:>{fmt}}".format(y22_pc, fmt=fmt))

    (s11, s22) = cross_section.get_sp()
    if s11 is not None:
        if list(set(cross_section.materials)) != [DEFAULT_MATERIAL]:
            table.add_row("M_p,11", "{:>{fmt}}".format(s11, fmt=fmt))
            table.add_row("M_p,22", "{:>{fmt}}".format(s22, fmt=fmt))
        else:
            table.add_row("S11", "{:>{fmt}}".format(s11, fmt=fmt))
            table.add_row("S22", "{:>{fmt}}".format(s22, fmt=fmt))

    (sf_11_plus, sf_11_minus, sf_22_plus, sf_22_minus) = cross_section.get_sf_p()
    if sf_11_plus is not None:
        table.add_row("SF_11+", "{:>{fmt}}".format(sf_11_plus, fmt=fmt))
        table.add_row("SF_11-", "{:>{fmt}}".format(sf_11_minus, fmt=fmt))
        table.add_row("SF_22+", "{:>{fmt}}".format(sf_22_plus, fmt=fmt))
        table.add_row("SF_22-", "{:>{fmt}}".format(sf_22_minus, fmt=fmt))

    console = Console()
    console.print(table)
    print("")
