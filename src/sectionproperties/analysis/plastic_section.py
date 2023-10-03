"""PlasticSection class for calculating plastic properties."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import brentq

import sectionproperties.analysis.fea as fea
import sectionproperties.pre.pre as pre


if TYPE_CHECKING:
    from rich.progress import Progress
    from scipy.optimize import RootResults

    from sectionproperties.analysis.section import Section
    from sectionproperties.pre.geometry import CompoundGeometry, Geometry


class PlasticSection:
    """Class for the plastic analysis of cross-sections.

    Stores the finite element geometry and provides methods to compute the plastic
    section properties.
    """

    def __init__(
        self,
        geometry: Geometry | CompoundGeometry,
    ) -> None:
        """Inits the PlasticSection class.

        Args:
            geometry: Section geometry object
        """
        self.geometry = geometry.align_center()
        self.geometry.compile_geometry()

        # initialize variables to be defined later within calculate_plastic_force
        self._c_top = [0.0, 0.0]
        self._c_bot = [0.0, 0.0]
        self._f_top = 0.0

    def calculate_plastic_properties(
        self,
        section: Section,
        verbose: bool,
        progress: Progress | None = None,
    ) -> None:
        """Calculates the plastic properties.

        Calculates the location of the plastic centroid with respect to the centroidal
        and principal bending axes, the plastic section moduli and shape factors.

        Stores the results in the
        :class:`~sectionproperties.post.post.SectionProperties` object belonging to the
        supplied :class:`~sectionproperties.analysis.section.Section` object.

        Args:
            section: Cross-section object containing the same geometry as this
                :class:`~sectionproperties.analysis.section.PlasticSection` object
            verbose: If set to True, prints convergence information to the terminal
            progress: Rich progress object

        Raises:
            RuntimeError: A geometric analysis has not yet been performed
        """
        if progress:
            task = progress.add_task(
                description="[red]Calculating plastic properties",
                total=4,
            )
        else:
            task = None

        # 1) Calculate plastic properties for centroidal axis
        # calculate distances to the extreme fibres
        fibres = self.calculate_extreme_fibres(angle=0)

        # 1a) Calculate x-axis plastic centroid
        y_pc, r = self.pc_algorithm(
            u=(1, 0),
            dlim=fibres[2:],  # fibres[2:] = ymin, ymax
            axis=1,
            verbose=verbose,
        )

        self.check_convergence(root_result=r, axis="x-axis")
        section.section_props.y_pc = y_pc
        section.section_props.sxx = self._f_top * abs(self._c_top[1] - self._c_bot[1])

        if verbose:
            self.print_verbose(d=y_pc, root_result=r, axis="x-axis")

        if progress and task is not None:
            progress.update(task_id=task, advance=1)

        # 1b) Calculate y-axis plastic centroid
        x_pc, r = self.pc_algorithm(
            u=(0, 1),
            dlim=fibres[0:2],  # fibres[0:2] = xmin, xmax
            axis=2,
            verbose=verbose,
        )

        self.check_convergence(root_result=r, axis="y-axis")
        section.section_props.x_pc = x_pc
        section.section_props.syy = self._f_top * abs(self._c_top[0] - self._c_bot[0])

        if verbose:
            self.print_verbose(d=x_pc, root_result=r, axis="y-axis")

        if progress and task is not None:
            progress.update(task_id=task, advance=1)

        # 2) Calculate plastic properties for principal axis
        # convert principal axis angle to radians
        if section.section_props.phi is not None:
            angle = section.section_props.phi * np.pi / 180
        else:
            raise RuntimeError("Run a geometric analysis prior to a plastic analysis.")

        # unit vectors in the axis directions
        ux = (np.cos(angle), np.sin(angle))
        uy = (-np.sin(angle), np.cos(angle))

        # calculate distances to the extreme fibres in the principal axis
        fibres = self.calculate_extreme_fibres(angle=section.section_props.phi)

        # 2a) Calculate 11-axis plastic centroid
        y22_pc, r = self.pc_algorithm(
            u=ux,
            dlim=fibres[2:],
            axis=1,
            verbose=verbose,
        )

        # calculate the centroids in the principal coordinate system
        c_top_p = fea.principal_coordinate(
            phi=section.section_props.phi, x=self._c_top[0], y=self._c_top[1]
        )
        c_bot_p = fea.principal_coordinate(
            phi=section.section_props.phi, x=self._c_bot[0], y=self._c_bot[1]
        )

        self.check_convergence(root_result=r, axis="11-axis")
        section.section_props.y22_pc = y22_pc
        section.section_props.s11 = self._f_top * abs(c_top_p[1] - c_bot_p[1])

        if verbose:
            self.print_verbose(d=y22_pc, root_result=r, axis="11-axis")

        if progress and task is not None:
            progress.update(task_id=task, advance=1)

        # 2b) Calculate 22-axis plastic centroid
        x11_pc, r = self.pc_algorithm(
            u=uy,
            dlim=fibres[0:2],
            axis=2,
            verbose=verbose,
        )

        # calculate the centroids in the principal coordinate system
        c_top_p = fea.principal_coordinate(
            phi=section.section_props.phi, x=self._c_top[0], y=self._c_top[1]
        )
        c_bot_p = fea.principal_coordinate(
            phi=section.section_props.phi, x=self._c_bot[0], y=self._c_bot[1]
        )

        self.check_convergence(root_result=r, axis="22-axis")
        section.section_props.x11_pc = x11_pc
        section.section_props.s22 = self._f_top * abs(c_top_p[0] - c_bot_p[0])

        if verbose:
            self.print_verbose(d=x11_pc, root_result=r, axis="22-axis")

        if progress and task is not None:
            progress.update(task_id=task, advance=1)

        # if there are no materials specified, calculate shape factors
        if list(set(section.materials)) == [pre.DEFAULT_MATERIAL]:
            if (
                section.section_props.zxx_plus
                and section.section_props.zxx_minus
                and section.section_props.zyy_plus
                and section.section_props.zyy_minus
            ):
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
            if (
                section.section_props.z11_plus
                and section.section_props.z11_minus
                and section.section_props.z22_plus
                and section.section_props.z22_minus
            ):
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

        if progress and task is not None:
            progress.update(
                task_id=task,
                description="[bold green]:white_check_mark: Plastic analysis complete",
            )

    def check_convergence(
        self,
        root_result: RootResults,
        axis: str,
    ) -> None:
        """Checks that the function solver converged and if not, raises a helpful error.

        Args:
            root_result: Result object from the root finder
            axis: Axis being considered by the function solver

        Raises:
            RuntimeError: If the function solver did not converge
        """
        if not root_result.converged:
            msg = f"Plastic centroid calculation about the {axis}"
            msg += " failed. Contact robbie.vanleeuwen@gmail.com with your"
            msg += f" analysis parameters. Termination flag: {root_result.flag}"
            raise RuntimeError(msg)

    def print_verbose(
        self,
        d: float,
        root_result: RootResults,
        axis: str,
    ) -> None:
        """Prints information related to the plastic solver convergence to the terminal.

        Args:
            d: Location of the plastic centroid axis
            root_result: Result object from the root finder
            axis: Axis being considered by the function solver
        """
        msg = f"---{axis} plastic centroid calculation converged at "
        msg += f"{d:.5e} in {root_result.iterations} iterations."
        print(msg)

    def calculate_extreme_fibres(
        self,
        angle: float,
    ) -> tuple[float, float, float, float]:
        """Calculates the section extreme fibres along and perpendicular to an axis.

        Args:
            angle: Angle (in degrees) along which to calculate the extreme fibre
                locations

        Returns:
            The location of the extreme fibres parallel (``u``) and perpendicular
            (``v``) to the axis, (``u_min``, ``u_max``, ``v_min``, ``v_max``)
        """
        # initialise mins and maxs
        pt = self.geometry.points[0]
        u, v = fea.principal_coordinate(phi=angle, x=pt[0], y=pt[1])
        u_min = u
        u_max = u
        v_min = v
        v_max = v

        # loop through all nodes in the mesh
        for pt in self.geometry.points[1:]:
            # determine the coordinate of the point wrt the axis
            u, v = fea.principal_coordinate(phi=angle, x=pt[0], y=pt[1])

            # update the mins and maxes where necessary
            u_min = min(u_min, u)
            u_max = max(u_max, u)
            v_min = min(v_min, v)
            v_max = max(v_max, v)

        return u_min, u_max, v_min, v_max

    def pc_algorithm(
        self,
        u: tuple[float, float],
        dlim: tuple[float, float],
        axis: int,
        verbose: bool,
    ) -> tuple[float, RootResults]:
        """Plastic centroid algorithm.

        An algorithm used for solving for the location of the plastic centroid. The
        algorithm searches for the location of the axis within the section depth
        (defined by unit vector ``u``) that satisfies force equilibrium.

        Args:
            u: Unit vector pointing in the direction of the axis
            dlim: Distances from the centroid to the extreme fibres perpendicular to the
                axis (``dmax``, ``dmin``)
            axis: The current axis direction, 1 (i.e. ``x`` or ``11``) or 2 (i.e. ``y``
                or ``22``)
            verbose: If set to True, prints convergence information to the terminal

        Returns:
            The perpendicular distance from the origin to the plastic centroid axis
            ``d``, and the ``scipy`` results object ``r`` - (``d``, ``r``)
        """
        # calculate vector perpendicular to u
        if axis == 1:
            u_p = np.array([-u[1], u[0]])
        else:
            u_p = np.array([u[1], -u[0]])

        return brentq(
            f=self.evaluate_force_eq,
            a=dlim[0],
            b=dlim[1],
            args=(u, u_p, verbose),
            full_output=True,
            disp=False,
            xtol=1e-6,
            rtol=1e-6,  # type:ignore
        )

    def evaluate_force_eq(
        self,
        d: float,
        u: tuple[float, float],
        u_p: tuple[float, float],
        verbose: bool,
    ) -> float:
        """Evaluates force equilibrium.

        Given a perpendicular distance ``d`` from the origin to an axis (defined by unit
        vector ``u``), calculates the force equilibrium between forces above and below
        that axis. The resultant force as a ratio of the total force, is returned.

        Args:
            d: Perpendicular distance from the origin to the current axis
            u: Unit vector defining the direction of the axis
            u_p: Unit vector perpendicular to the direction of the axis
            verbose: If set to True, prints convergence information to the terminal

        Returns:
            The force equilibrium norm
        """
        # p = finding a point on the axis by scaling the perpendicular
        p = d * u_p[0], d * u_p[1]

        # calculate force equilibrium
        f_top, f_bot = self.calculate_plastic_force(u=u, p=p)

        # calculate the force norm
        f_norm = (f_top - f_bot) / (f_top + f_bot)

        # print verbose results
        if verbose:
            print(f"d = {d}; f_norm = {f_norm}")

        # return the force norm (target for the root finding algorithm)
        return f_norm

    def calculate_plastic_force(
        self,
        u: tuple[float, float],
        p: tuple[float, float],
    ) -> tuple[float, float]:
        """Calculates the plastic force above and below an axis line.

        Sums the forces above and below the axis line defined by unit vector ``u`` and
        point ``p``.

        Args:
            u: Unit vector defining the direction of the axis line
            p: Point on the axis line

        Returns:
            Force in the above and below the axis line (``f_top``, ``f_bot``)
        """
        # initialise variables
        f_top, f_bot = 0, 0
        a_top, a_bot = 0, 0
        qx_top, qx_bot = 0, 0
        qy_top, qy_bot = 0, 0

        # split geometry above and below the line
        top_geoms, bot_geoms = self.geometry.split_section(point_i=p, vector=u)

        if top_geoms:
            # loop through all top geometries
            for top_geom in top_geoms:
                # get properties
                f_y = top_geom.material.yield_strength
                area_top = top_geom.calculate_area()
                cx, cy = top_geom.calculate_centroid()

                # sum properties
                a_top += area_top
                qx_top += cy * area_top
                qy_top += cx * area_top
                f_top += f_y * area_top

        if bot_geoms:
            # loop through all bottom geometries
            for bot_geom in bot_geoms:
                # get properties
                f_y = bot_geom.material.yield_strength
                area_bot = bot_geom.calculate_area()
                cx, cy = bot_geom.calculate_centroid()

                # sum properties
                a_bot += area_bot
                qx_bot += cy * area_bot
                qy_bot += cx * area_bot
                f_bot += f_y * area_bot

        # calculate centroid of force action
        try:
            self._c_top = [qy_top / a_top, qx_top / a_top]
            self._f_top = f_top
        except ZeroDivisionError:
            self._c_top = [0, 0]
            self._f_top = 0

        try:
            self._c_bot = [qy_bot / a_bot, qx_bot / a_bot]
        except ZeroDivisionError:
            self._c_bot = [0, 0]

        return f_top, f_bot
