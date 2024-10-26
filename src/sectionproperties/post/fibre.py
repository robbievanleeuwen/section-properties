"""Provides functionalities to export a section to a fibre section.

It can be used in `suanPan <https://github.com/TLCFEM/suanPan>`_ to perform further FEA.
"""

from __future__ import annotations

from abc import abstractmethod
from typing import TYPE_CHECKING

import sectionproperties.analysis.solver as solver
from sectionproperties.analysis.fea import Tri6, shape_function
from sectionproperties.analysis.section import Section
from sectionproperties.pre.geometry import CompoundGeometry, Geometry

if TYPE_CHECKING:
    import numpy as np
    import numpy.typing as npt


class Cell:
    """Holds the information of a fibre cell.

    Attributes:
        tag: The tag of the cell
        area: The area of the cell
        material: The material name of the cell
        y: The y-coordinate of the cell
        z: The z-coordinate of the cell
        omega: The warping function of the cell
        py: The derivative of the warping function with respect to y
        pz: The derivative of the warping function with respect to z
    """

    tag: int
    area: float
    material: str
    y: float
    z: float
    omega: float
    py: float
    pz: float

    def __init__(
        self,
        ele: Tri6,
        omega: npt.NDArray[np.float64] | None = None,
    ) -> None:
        """Converts a Tri6 element to a fibre cell.

        If `omega` is None, no warping is considered.

        Args:
            ele: The Tri6 element
            omega: The warping function. Defaults to ``None``.
        """
        n, dn, self.area, _, _ = shape_function(
            ele.coords, (0.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)
        )
        self.material = "-".join(ele.material.name.split())
        self.tag = ele.el_id + 1  # el_id is zero-based
        self.y, self.z = ele.coords[:, :3].mean(axis=1)
        if omega is not None:
            ele_omega = omega[ele.node_ids]
            self.omega = ele_omega.dot(n)
            self.py, self.pz = dn @ ele_omega
        else:
            self.omega = 0.0
            self.py = 0.0
            self.pz = 0.0

    @abstractmethod
    def export(self) -> str:
        """Export the cell to a string.

        Needs to be implemented by subclasses.

        Returns:
            str: The exported cell
        """
        raise NotImplementedError


class Cell2D(Cell):
    """A 2D cell section."""

    def export(self) -> str:
        """Export the cell to a 2D cell section.

        Returns:
            str: The exported cell
        """
        return f"section Cell2D {self.tag} {self.area:.8e} {self.material} {self.y:.8e}"


class Cell3D(Cell):
    """A 3D cell section."""

    def export(self) -> str:
        """Export the cell to a 3D cell section.

        Returns:
            str: The exported cell
        """
        return (
            f"section Cell3D {self.tag} {self.area:.8e} {self.material}"
            f" {self.y:.8e} {self.z:.8e}"
        )


class Cell3DOS(Cell):
    """A 3DOS cell section."""

    def export(self) -> str:
        """Export the cell to a 3DOS cell section.

        Returns:
            str: The exported cell
        """
        return (
            f"section Cell3DOS {self.tag} {self.area:.8e}"
            f" {self.omega:.8e} {self.py:.8e} {self.pz:.8e}"
            f" {self.material}"
            f" {self.y:.8e} {self.z:.8e}"
        )


def to_fibre_section(
    obj: Geometry | CompoundGeometry | Section,
    *,
    main_section_tag: int = 1,
    analysis_type: str = "3DOS",
    material_mapping: dict[str, int] | None = None,
    max_width: int = 160,
    save_to: str | None = None,
) -> str:
    """Export a section to the corresponding commands to create a fibre section.

    For a given geometry, this function computes necessary sectional properties and
    exports the corresponding commands to create a fibre/composite section that can be
    used in `suanPan <https://github.com/TLCFEM/suanPan>`_.

    Args:
        obj: The geometry/section to be exported
        main_section_tag: The tag of the main section. Defaults to ``1``.
        analysis_type: The type of analysis would be performed. Defaults to ``"3DOS"``.
        material_mapping: A dictionary mapping material names to material tags. Defaults
            to ``None``.
        max_width: The maximum width of a line in the output. Defaults to ``160``.
        save_to: The path to save the output to. Defaults to ``None``.

    Raises:
        TypeError: If `obj` is not a Geometry or Section
        ValueError: If `analysis_type` is not 2D, 3D or 3DOS

    Returns:
        str: The exported commands
    """
    if isinstance(obj, Geometry | CompoundGeometry):
        geometry = obj
    elif isinstance(obj, Section):  # pyright: ignore [reportUnnecessaryIsInstance]
        geometry = obj.geometry
    else:
        msg = f"Expected a Geometry or Section, got {type(obj).__name__}"
        raise TypeError(msg)

    analysis_type = analysis_type.upper()

    cell_class: type
    if analysis_type == "2D":
        cell_class = Cell2D
        fibre_class = "Fibre2D"
    elif analysis_type == "3D":
        cell_class = Cell3D
        fibre_class = "Fibre3D"
    elif analysis_type == "3DOS":
        cell_class = Cell3DOS
        fibre_class = "Fibre3DOS"
    else:
        msg = "Invalid analysis type, expected 2D, 3D or 3DOS"
        raise ValueError(msg)

    section = Section(geometry)

    if cell_class is Cell3DOS:
        stiffness, force = section.assemble_torsion()
        omega = solver.solve_direct_lagrange(stiffness, force)
    else:
        omega = None

    cells = [cell_class(ele, omega) for ele in section.elements]

    commands = """# This is generated by sectionproperties library.
# Please note the following:
#   1. The warping function and its derivatives (if present) are computed in the
#      local coordinate system (about origin).
#   2. Beware of the potential different orientations of beam section.
#   3. It may be necessary to manually adjust the material tags.
#   4. If uncertain, please validate the behaviour first.
"""

    part = f"section {fibre_class} {main_section_tag} "
    for cell in cells:
        cell.tag += main_section_tag
        part += f"{cell.tag} "
        if len(part) > max_width:
            commands += part + "\\\n"
            part = ""

    commands += part + "\n\n"
    commands += "\n".join(cell.export() for cell in cells)
    commands += "\n\n"

    if material_mapping is not None:
        for name, tag in material_mapping.items():
            commands = commands.replace("-".join(name.split()), str(tag))

    if isinstance(save_to, str):
        with open(save_to, "w") as f:
            f.write(commands)

    return commands
