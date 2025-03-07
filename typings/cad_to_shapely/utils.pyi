from shapely.geometry import Polygon

def find_holes(polygons: list[Polygon]) -> Polygon: ...
def filter_polygons(
    polygons: list[Polygon],
    filter_flag: int = ...,
) -> list[Polygon]: ...
