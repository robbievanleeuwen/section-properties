from typing import Union, List
from dataclasses import dataclass
import numpy as np
import triangle


class GeometryError(Exception):
    """Exception raised when invalid geometry is found."""

    pass


@dataclass(eq=True, frozen=True)
class Material:
    """Class for structural materials.

    Provides a way of storing material properties related to a specific material. The color can be
    a multitude of different formats, refer to https://matplotlib.org/api/colors_api.html and
    https://matplotlib.org/examples/color/named_colors.html for more information.

    :param string name: Material name
    :param float elastic_modulus: Material modulus of elasticity
    :param float poissons_ratio: Material Poisson's ratio
    :param float yield_strength: Material yield strength
    :param float density: Material density (mass per unit volume)
    :param color: Material color for rendering
    :type color: :class:`matplotlib.colors`

    :cvar string name: Material name
    :cvar float elastic_modulus: Material modulus of elasticity
    :cvar float poissons_ratio: Material Poisson's ratio
    :cvar float shear_modulus: Material shear modulus, derived from the elastic modulus and
        Poisson's ratio assuming an isotropic material
    :cvar float density: Material density (mass per unit volume)
    :cvar float yield_strength: Material yield strength
    :cvar color: Material color for rendering
    :vartype color: :class:`matplotlib.colors`

    The following example creates materials for concrete, steel and timber::

        from sectionproperties.pre.pre import Material

        concrete = Material(
            name='Concrete', elastic_modulus=30.1e3, poissons_ratio=0.2, density=2.4e-6,
                yield_strength=32, color='lightgrey'
        )
        steel = Material(
            name='Steel', elastic_modulus=200e3, poissons_ratio=0.3, density=7.85e-6,
                yield_strength=500, color='grey'
        )
        timber = Material(
            name='Timber', elastic_modulus=8e3, poissons_ratio=0.35, density=6.5e-7,
                yield_strength=20, color='burlywood'
        )
    """

    name: str
    elastic_modulus: float
    poissons_ratio: float
    yield_strength: float
    density: float
    color: str

    @property
    def shear_modulus(self):
        return self.elastic_modulus / (2 * (1 + self.poissons_ratio))


DEFAULT_MATERIAL = Material("default", 1, 0, 1, 1, "w")


def create_mesh(
    points: List[List[float]],
    facets: List[List[float]],
    holes: List[List[float]],
    control_points: List[List[float]],
    mesh_sizes: Union[List[float], float],
):
    """Creates a quadratic triangular mesh using the triangle module, which utilises the code
    'Triangle', by Jonathan Shewchuk.

    :param points: List of points *(x, y)* defining the vertices of the cross-section
    :type points: list[list[float, float]]
    :param facets: List of point index pairs *(p1, p2)* defining the edges of the cross-section
    :type points: list[list[int, int]]
    :param holes: List of points *(x, y)* defining the locations of holes within the cross-section.
        If there are no holes, provide an empty list [].
    :type holes: list[list[float, float]]
    :param control_points: A list of points *(x, y)* that define different regions of the
        cross-section. A control point is an arbitrary point within a region enclosed by facets.
    :type control_points: list[list[float, float]]
    :param mesh_sizes: List of maximum element areas for each region defined by a control point
    :type mesh_sizes: list[float]

    :return: Dictionary containing mesh data
    :rtype: dict()
    """
    if not isinstance(mesh_sizes, list):
        mesh_sizes = [mesh_sizes]

    tri = {}  # create tri dictionary
    tri["vertices"] = points  # set point
    tri["segments"] = facets  # set facets

    if holes:
        tri["holes"] = holes  # set holes

    # prepare regions
    regions = []

    for (i, cp) in enumerate(control_points):
        regions.append([cp[0], cp[1], i, mesh_sizes[i]])

    tri["regions"] = regions  # set regions

    # generate mesh
    mesh = triangle.triangulate(tri, "pq30Aao2")

    return mesh
