from shapely import Polygon

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
