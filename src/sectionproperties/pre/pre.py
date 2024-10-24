"""Classes and methods for generic pre-procesing in sectionproperties."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import cytriangle as triangle


@dataclass(eq=True, frozen=True)
class Material:
    """Class for structural materials.

    Provides a way of storing material properties related to a specific material. The
    color can be a multitude of different formats, refer to
    https://matplotlib.org/stable/api/colors_api.html and
    https://matplotlib.org/stable/gallery/color/named_colors.html for more information.

    Attributes:
        name: Material name
        elastic_modulus: Material modulus of elasticity
        poissons_ratio: Material Poisson's ratio
        yield_strength: Material yield strength
        density: Material density (mass per unit volume)
        color: Material color for rendering

    Example:
        The following example creates materials for concrete, steel and timber::

            from sectionproperties.pre import Material

            concrete = Material(
                name="Concrete",
                elastic_modulus=30.1e3,
                poissons_ratio=0.2,
                density=2.4e-6,
                yield_strength=32,
                color="lightgrey",
            )
            steel = Material(
                name="Steel",
                elastic_modulus=200e3,
                poissons_ratio=0.3,
                density=7.85e-6,
                yield_strength=500,
                color="grey",
            )
            timber = Material(
                name="Timber",
                elastic_modulus=8e3,
                poissons_ratio=0.35,
                density=6.5e-7,
                yield_strength=20,
                color="burlywood",
            )
    """

    name: str
    elastic_modulus: float
    poissons_ratio: float
    yield_strength: float
    density: float
    color: str

    @property
    def shear_modulus(self) -> float:
        """Returns the shear modulus of the material.

        Material shear modulus, derived from the elastic modulus and Poisson's ratio
        assuming an isotropic material.

        Returns:
            Shear modulus of the material
        """
        return self.elastic_modulus / (2 * (1 + self.poissons_ratio))


DEFAULT_MATERIAL = Material("default", 1, 0, 1, 1, "w")


def create_mesh(
    points: list[tuple[float, float]],
    facets: list[tuple[int, int]],
    holes: list[tuple[float, float]],
    control_points: list[tuple[float, float]],
    mesh_sizes: list[float] | float,
    min_angle: float,
    coarse: bool,
) -> dict[str, list[list[float]] | list[list[int]]]:
    """Generates a triangular mesh.

    Creates a quadratic triangular mesh using the ``CyTriangle`` module, which utilises
    the code ``Triangle``, by Jonathan Shewchuk.

    Args:
        points: List of points (``x``, ``y``) defining the vertices of the cross-section
        facets: List of point index pairs (``p1``, ``p2``) defining the edges of the
            cross-section
        holes: List of points (``x``, ``y``) defining the locations of holes within the
            cross-section. If there are no holes, provide an empty list [].
        control_points: A list of points (``x``, ``y``) that define different regions of
            the cross-section. A control point is an arbitrary point within a region
            enclosed by facets.
        mesh_sizes: List of maximum element areas for each region defined by a control
            point
        min_angle: The meshing algorithm adds vertices to the mesh to ensure that no
            angle smaller than the minimum angle (in degrees, rounded to 1 decimal
            place). Note that small angles between input segments cannot be eliminated.
            If the minimum angle is 20.7 deg or smaller, the triangulation algorithm is
            theoretically guaranteed to terminate (given sufficient precision). The
            algorithm often doesn't terminate for angles greater than 33 deg. Some
            meshes may require angles well below 20 deg to avoid problems associated
            with insufficient floating-point precision.
        coarse: If set to True, will create a coarse mesh (no area or quality
            constraints)

    Returns:
        Dictionary containing mesh data
    """
    if not isinstance(mesh_sizes, list):
        mesh_sizes = [mesh_sizes]

    tri: dict[str, Any] = {}  # create tri dictionary
    tri["vertices"] = points  # set point
    tri["segments"] = facets  # set facets

    if holes:
        tri["holes"] = holes  # set holes

    # prepare regions
    regions: list[list[float | int]] = []

    for i, cp in enumerate(control_points):
        rg = [cp[0], cp[1], i, mesh_sizes[i]]
        regions.append(rg)

    tri["regions"] = regions  # set regions

    # generate mesh
    if coarse:
        mesh = triangle.triangulate(tri, "pAo2")
    else:
        mesh = triangle.triangulate(tri, f"pq{min_angle:.1f}Aao2")

    return mesh
