import pathlib
from shapely.geometry import Polygon, MultiPolygon
import cad_to_shapely as c2s
from sectionproperties.pre.sections import Geometry, CompoundGeometry


def load_dxf(dxf_filepath: pathlib.Path):
    """
            Import any-old-shape in dxf format for analysis.
            Code by aegis1980 and connorferster
        """
    if not dxf_filepath.exists():
        raise ValueError(f"The filepath does not exist: {dxf_filepath}")

    my_dxf = c2s.dxf.DxfImporter(dxf_filepath)
    my_dxf.process()
    my_dxf.cleanup()

    polygons = my_dxf.polygons
    new_polygons = c2s.utils.find_holes(polygons)
    if isinstance(new_polygons, MultiPolygon):
        return CompoundGeometry(new_polygons)
    elif isinstance(new_polygons, Polygon):
        return Geometry(new_polygons)
    else:
        print(f"No shapely.Polygon objects found in file: {dxf_filepath}")
