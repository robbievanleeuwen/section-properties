import pathlib
from typing import List, Union
from shapely.geometry.polygon import Polygon
from rhino_shapely_interop.importers import RhImporter


def load_3dm(r3dm_filepath: Union[pathlib.Path, str], **kwargs) -> List[Polygon]:
    """Load a Rhino `.3dm` file and import the single surface planer breps.

    :param r3dm_filepath:
        File path to the rhino `.3dm` file.
    :type r3dm_filepath: pathlib.Path or string
    :param \**kwargs:
        See below.
    :raises RuntimeError:
        A RuntimeError is raised if no polygons are found in the file.
        This is dependent on the keyword arguments.
        Try adjusting the keyword arguments if this error is raised.
    :return:
        List of Polygons found in the file.
    :rtype: List[shapely.geometry.Polygon]

    :Keyword Arguments:
        * *refine_num* (``int, optional``) --
            Bézier curve interpolation number. In Rhino a surface's edges are nurb based curves.
            Shapely does not support nurbs, so the individual Bézier curves are interpolated using straight lines.
            This parameter sets the number of straight lines used in the interpolation.
            Default is 1.
        * *vec1* (``numpy.ndarray, optional``) --
            A 3d vector in the Shapely plane. Rhino is a 3D geometry environment.
            Shapely is a 2D geometric library.
            Thus a 2D plane needs to be defined in Rhino that represents the Shapely coordinate system.
            `vec1` represents the 1st vector of this plane. It will be used as Shapely's x direction.
            Default is [1,0,0].
        * *vec2* (``numpy.ndarray, optional``) --
            Continuing from `vec1`, `vec2` is another vector to define the Shapely plane.
            It must not be [0,0,0] and it's only requirement is that it is any vector in the Shapely plane (but not equal to `vec1`).
            Default is [0,1,0].
        * *plane_distance* (``float, optional``) --
            The distance to the Shapely plane.
            Default is 0.
        * *project* (``boolean, optional``) --
            Controls if the breps are projected onto the plane in the direction of the Shapley plane's normal.
            Default is True.
        * *parallel* (``boolean, optional``) --
            Controls if only the rhino surfaces that have the same normal as the Shapely plane are yielded.
            If true, all non parallel surfaces are filtered out.
            Default is False.
    """
    rhi = RhImporter.from_file(str(r3dm_filepath))
    list_polygons = list(rhi.get_planer_brep(**kwargs))
    if len(list_polygons) == 0:
        raise RuntimeError(
            f"No shapely.Polygon objects found. "
            f"Consider adjusting the keyword arguments. "
            f"File name: {r3dm_filepath}. "
        )
    return list_polygons


def load_brep_encoding(brep: str, **kwargs) -> Polygon:
    """Load an encoded single surface planer brep.

    :param brep:
        Rhino3dm.Brep encoded as a string.
    :type brep: str
    :param \**kwargs:
        See below.
    :raises RuntimeError:
        A RuntimeError is raised if no polygons are found in the encoding.
        This is dependent on the keyword arguments.
        Try adjusting the keyword arguments if this error is raised.
    :return:
        The Polygons found in the encoding string.
    :rtype: shapely.geometry.Polygon

    :Keyword Arguments:
        * *refine_num* (``int, optional``) --
            Bézier curve interpolation number. In Rhino a surface's edges are nurb based curves.
            Shapely does not support nurbs, so the individual Bézier curves are interpolated using straight lines.
            This parameter sets the number of straight lines used in the interpolation.
            Default is 1.
        * *vec1* (``numpy.ndarray, optional``) --
            A 3d vector in the Shapely plane. Rhino is a 3D geometry environment.
            Shapely is a 2D geometric library.
            Thus a 2D plane needs to be defined in Rhino that represents the Shapely coordinate system.
            `vec1` represents the 1st vector of this plane. It will be used as Shapely's x direction.
            Default is [1,0,0].
        * *vec2* (``numpy.ndarray, optional``) --
            Continuing from `vec1`, `vec2` is another vector to define the Shapely plane.
            It must not be [0,0,0] and it's only requirement is that it is any vector in the Shapely plane (but not equal to `vec1`).
            Default is [0,1,0].
        * *plane_distance* (``float, optional``) --
            The distance to the Shapely plane.
            Default is 0.
        * *project* (``boolean, optional``) --
            Controls if the breps are projected onto the plane in the direction of the Shapley plane's normal.
            Default is True.
        * *parallel* (``boolean, optional``) --
            Controls if only the rhino surfaces that have the same normal as the Shapely plane are yielded.
            If true, all non parallel surfaces are filtered out.
            Default is False.
    """
    rhi = RhImporter.from_serialzed_brep(brep)
    geom = list(rhi.get_planer_brep(**kwargs))
    if len(geom) == 0:
        raise RuntimeError(f"No shapely.Polygon objects found for encoded object")
    return geom
