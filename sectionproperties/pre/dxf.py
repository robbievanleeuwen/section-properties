import pathlib
from shapely.geometry import Polygon, MultiPolygon
from sectionproperties.pre.sections import Geometry, CompoundGeometry

def load_dxf(dxf_filepath: pathlib.Path):
        """
            Import any-old-shape in dxf format for analysis.
            Code by aegis1980 and connorferster
        """
        c2s = None
        try:
            import cad_to_shapely as c2s # type: ignore
        except ImportError as e:
            print (e.message)
            print("To use 'from_dxf(...)' you need to 'pip install cad_to_shapely'")
            return

        if not dxf_filepath.exists():
            raise ValueError(f"The filepath does not exist: {dxf_filepath}")

        #TODO avoid step of making a temp file locally
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

        # # hole_cps = [] # control points for holes
        # # for hole in new_polygons.interiors:
        # #     p = c2s.utils.point_in_polygon(shapely.geometry.Polygon(hole))
        # #     hole_cps.append([p.x,p.y])

        # # create the custom geometry object
        # points = [list(elem) for elem in list(zip(new_polygons.exterior.xy[0],new_polygons.exterior.xy[1]))]
        # cp = c2s.utils.point_in_polygon(new_polygons)

        # for hole_lr in new_polygons.interiors:
        #     pts = [list(elem) for elem in list(zip(hole_lr.xy[0],hole_lr.xy[1]))]
        #     points.extend(pts)

        # facets = c2s.utils.facets(new_polygons,inc_holes=True)

        # return CustomSection(points = points, facets=facets, holes =hole_cps, control_points=[[cp.x,cp.y]] )