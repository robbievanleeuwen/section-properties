from shapely import Polygon
from shapely.geometry.collection import GeometryCollection

class DxfImporter:
    polygons: list[Polygon]

    def __init__(
        self,
        filename: str,
    ) -> None: ...
    def process(
        self,
        spline_delta: float = ...,
        degrees_per_segment: float = ...,
    ) -> str: ...
    def cleanup(
        self,
        simplify: bool = ...,
        zip_length: float = ...,
        retry_with_zip: bool = ...,
    ) -> str: ...
    def polygonize(
        self,
        simplify: bool = True,
        force_zip: bool = False,
        zip_length: float = 0.000001,
        retry_with_zip: bool = True,
    ) -> tuple[
        GeometryCollection, GeometryCollection, GeometryCollection, GeometryCollection
    ]: ...
