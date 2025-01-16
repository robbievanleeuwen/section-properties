"""Methods to load rhino files, interfacing with rhino-shapely-interop."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from rhino_shapely_interop.importers import RhImporter

if TYPE_CHECKING:
    import pathlib

    from shapely import Polygon


def load_3dm(
    r3dm_filepath: pathlib.Path | str,
    **kwargs: Any,
) -> list[Polygon]:
    """Load a Rhino ``.3dm`` file and import the single surface planer breps.

    Args:
        r3dm_filepath: File path to the rhino ``.3dm`` file.
        kwargs: See below.

    Keyword Args:
        refine_num (Optional[int]):  Bézier curve interpolation number. In Rhino a
            surface's edges are nurb based curves. Shapely does not support nurbs, so
            the individual Bézier curves are interpolated using straight lines. This
            parameter sets the number of straight lines used in the interpolation.
            Default is 1.
        vec1 (Optional[numpy.ndarray]): A 3d vector in the Shapely plane. Rhino is a 3D
            geometry environment. Shapely is a 2D geometric library. Thus a 2D plane
            needs to be defined in Rhino that represents the Shapely coordinate system.
            ``vec1`` represents the 1st vector of this plane. It will be used as
            Shapely's x direction. Default is [1,0,0].
        vec2 (Optional[numpy.ndarray]): Continuing from ``vec1``, ``vec2`` is another
            vector to define the Shapely plane. It must not be [0,0,0] and it's only
            requirement is that it is any vector in the Shapely plane (but not equal to
            ``vec1``). Default is [0,1,0].
        plane_distance (Optional[float]): The distance to the Shapely plane. Default is
            0.
        project (Optional[bool]): Controls if the breps are projected onto the plane in
            the direction of the Shapley plane's normal. Default is True.
        parallel (Optional[bool]): Controls if only the rhino surfaces that have the
            same normal as the Shapely plane are yielded. If true, all non parallel
            surfaces are filtered out. Default is False.

    Raises:
        RuntimeError: A RuntimeError is raised if no polygons are found in the file.
            This is dependent on the keyword arguments. Try adjusting the keyword
            arguments if this error is raised.

    Returns:
        List of Polygons found in the file.
    """
    rhi = RhImporter.from_file(str(r3dm_filepath))
    list_polygons = list(rhi.get_planer_brep(**kwargs))

    if len(list_polygons) == 0:
        msg = "No shapely.Polygon objects found. Consider adjusting the keyword "
        msg += f"arguments. File name: {r3dm_filepath}."
        raise RuntimeError(msg)

    return list_polygons


def load_brep_encoding(
    brep: str,
    **kwargs: Any,
) -> list[Polygon]:
    """Load an encoded single surface planer brep.

    Args:
        brep: Rhino3dm.Brep encoded as a string.
        kwargs: See below.

    Keyword Args:
        refine_num (Optional[int]):  Bézier curve interpolation number. In Rhino a
            surface's edges are nurb based curves. Shapely does not support nurbs, so
            the individual Bézier curves are interpolated using straight lines. This
            parameter sets the number of straight lines used in the interpolation.
            Default is 1.
        vec1 (Optional[numpy.ndarray]): A 3d vector in the Shapely plane. Rhino is a 3D
            geometry environment. Shapely is a 2D geometric library. Thus a 2D plane
            needs to be defined in Rhino that represents the Shapely coordinate system.
            ``vec1`` represents the 1st vector of this plane. It will be used as
            Shapely's x direction. Default is [1,0,0].
        vec2 (Optional[numpy.ndarray]): Continuing from ``vec1``, ``vec2`` is another
            vector to define the Shapely plane. It must not be [0,0,0] and it's only
            requirement is that it is any vector in the Shapely plane (but not equal to
            ``vec1``). Default is [0,1,0].
        plane_distance (Optional[float]): The distance to the Shapely plane. Default is
            0.
        project (Optional[bool]): Controls if the breps are projected onto the plane in
            the direction of the Shapley plane's normal. Default is True.
        parallel (Optional[bool]): Controls if only the rhino surfaces that have the
            same normal as the Shapely plane are yielded. If true, all non parallel
            surfaces are filtered out. Default is False.

    Raises:
        RuntimeError: A RuntimeError is raised if no polygons are found in the encoding.
            This is dependent on the keyword arguments. Try adjusting the keyword
            arguments if this error is raised.

    Returns:
        The Polygons found in the encoding string.
    """
    rhi = RhImporter.from_serialzed_brep(brep)
    geom = list(rhi.get_planer_brep(**kwargs))

    if len(geom) == 0:
        msg = "No shapely.Polygon objects found for encoded object."
        raise RuntimeError(msg)

    return geom
