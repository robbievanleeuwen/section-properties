from typing import Union, List
from dataclasses import dataclass
import numpy as np
import meshpy.triangle as triangle
import meshpy

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
    :param color: Material color for rendering
    :type color: :class:`matplotlib.colors`

    :cvar string name: Material name
    :cvar float elastic_modulus: Material modulus of elasticity
    :cvar float poissons_ratio: Material Poisson's ratio
    :cvar float shear_modulus: Material shear modulus, derived from the elastic modulus and
        Poisson's ratio assuming an isotropic material
    :cvar float yield_strength: Material yield strength
    :cvar color: Material color for rendering
    :vartype color: :class:`matplotlib.colors`

    The following example creates materials for concrete, steel and timber::

        from sectionproperties.pre.pre import Material

        concrete = Material(
            name='Concrete', elastic_modulus=30.1e3, poissons_ratio=0.2, yield_strength=32,
                color='lightgrey'
        )
        steel = Material(
            name='Steel', elastic_modulus=200e3, poissons_ratio=0.3, yield_strength=500,
                color='grey'
        )
        timber = Material(
            name='Timber', elastic_modulus=8e3, poissons_ratio=0.35, yield_strength=20,
                color='burlywood'
        )
    """
    name: str
    elastic_modulus: float
    poissons_ratio: float
    yield_strength: float
    color: str

    @property
    def shear_modulus(self):
        return self.elastic_modulus / (2 * (1 + self.poissons_ratio))

DEFAULT_MATERIAL = Material('default', 1, 0, 1, 'w')


def create_mesh(
    points: List[List[float]], 
    facets: List[List[float]], 
    holes: List[List[float]], 
    control_points: List[List[float]], 
    mesh_sizes: Union[List[float], float], 
    atol=1.0e-8,
    ) -> meshpy.triangle.MeshInfo:
    """Creates a quadratic triangular mesh using the meshpy module, which utilises the code
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
    :param atol: minimum permissable point distance from any section facet
    :type atol: float

    :return: Object containing generated mesh data
    :rtype: :class:`meshpy.triangle.MeshInfo`
    """

    # check_geometry(points, facets, holes, control_points, atol=atol)
    if not isinstance(mesh_sizes, list): mesh_sizes = [mesh_sizes]
    mesh = triangle.MeshInfo()  # create mesh info object
    mesh.set_points(points)  # set points
    mesh.set_facets(facets)  # set facets
    mesh.set_holes(holes)  # set holes

    # set regions
    mesh.regions.resize(len(control_points))  # resize regions list
    region_id = 0  # initialise region ID variable

    for (i, cp) in enumerate(control_points):
        mesh.regions[i] = [cp[0], cp[1], region_id, mesh_sizes[i]]
        region_id += 1

    mesh = triangle.build(
        mesh, min_angle=30, mesh_order=2, quality_meshing=True,
        attributes=True, volume_constraints=True)

    return mesh
