from __future__ import annotations
from typing import List, Optional, Union, Tuple

import pathlib
import logging 
from icecream import ic
from IPython.display import display_svg
import more_itertools
import numpy as np
from shapely.geometry import Polygon, MultiPolygon, LinearRing, LineString, Point, GeometryCollection
from shapely.ops import split
import shapely
import matplotlib.pyplot as plt
import sectionproperties.pre.pre as pre
import sectionproperties.pre.bisect_section as bisect
import sectionproperties.post.post as post


log = logging.getLogger('shapely')
log_path = pathlib.Path("C:\\Users\\cferster\\Desktop\\sectionproperties logs\\shapley.log")
logging.basicConfig(filename=log_path, filemode='w', format="%(message)s", level=logging.DEBUG)

class Geometry:
    """Class for defining the geometry of a contiguous section of a single material.

    Provides an interface for the user to specify the geometry defining a section. A method
    is provided for generating a triangular mesh, transforming the section (e.g. translation,
    rotation, perimeter offset, mirroring), aligning the geometry to another geometry, and
    designating stress recovery points.

    :cvar geom: a Polygon object that defines the geometry
    :vartype geom: shapely.geometry.Polygon

    :cvar material: Optional, a Material to associate with this geometry
    :vartype material: Optional[sectionproperties.pre.Material]
    """
    def __init__(self, 
        geom: shapely.geometry.Polygon, 
        material: pre.Material = pre.DEFAULT_MATERIAL,
        ):
        """Inits the Geometry class.
        """
        if isinstance(geom, MultiPolygon):
            raise ValueError(f"Use CompoundGeometry(...) for a MultiPolygon object.")
        if not isinstance(geom, Polygon):
            raise ValueError(f"Argument is not a valid shapely.geometry.Polygon object: {geom}")
        self.geom = geom
        self.material = material
        self.control_points = []
        self.shift = []
        self.points = []
        self.facets = []
        self.holes = []
        self.perimeter = []
        self._recovery_points = []
        self.compile_geometry()

    def _repr_svg_(self):
        print("sectionproperties.pre.sections.Geometry")
        print(f"object at: {hex(id(self))}")
        return self.geom._repr_svg_()

    @staticmethod
    def from_points(
        points: List[List[float]], 
        facets: Optional[List[List[int]]],
        holes: Optional[List[List[float]]],
        control_points: Optional[List[List[float]]],
        ):
        """
        An interface for the creation of Geometry objects through the definition of points, 
        facets, and holes. 

        :cvar points: List of points *(x, y)* defining the vertices of the section geometry.
        If facets are not provided, it is a assumed the that the list of points are ordered
        around the perimeter, either clockwise or anti-clockwise
        :vartype points: list[list[float]]
        :cvar facets: Optional. A list of *(start, end)* indexes of vertices defining the edges
        of the section geoemtry. Can be used to define both external and internal perimeters of holes.
        Facets are assumed to be described in the order of exterior perimeter, interior perimeter 1,
        interior perimeter 2, etc.
        :vartype facets: list[list[int]]
        :cvar holes: Optional. A list of points *(x, y)* that define interior regions as
        being holes or voids. The point can be located anywhere within the hole region.
        Only one point is required per hole region.
        """
        if facets is None: return Geometry(Polygon(points))
        if holes is not None and facets is None:
            raise ValueError(
                "If holes coordinates are provided then facets must also be provided "
                "to distinguish between exterior and interior edges."
                )
        if len(control_points) > 1:
            raise ValueError(
                "A Geometry object can only have one contiguous region (with holes)."
                "Did you mean to use CompoundGeometry.from_points()?"
        )
        if holes is None: holes = []
        prev_facet = []

        # Initialize the total number of accumulators needed
        # Always an exterior, plus, a separate accumulator for each interior region
        exterior = []
        interiors = [[] for _ in holes] # initialize an empty facet list for every hole
        interior_counter = 0 # To keep track of interior regions
        active_list = exterior # The active_list is the list being accumulated on

        for facet in facets: # Loop through facets for graph connectivity
            i_idx, _ = facet
            if not prev_facet: # Add the first facet vertex to exterior and move on
                active_list.append(points[i_idx])
                prev_facet = facet
                continue

            prev_j_idx = prev_facet[1] # Look at the last j_idx to test for a break in the chain of edges
            if i_idx != prev_j_idx and holes: #If there is a break in the chain of edges...
                if active_list == exterior: # ...and we are still accumulating on the exterior...
                    active_list = interiors[interior_counter] # ... then move to the interior accumulator
                else: # ...or if we are already in the interior accumulator...
                    interior_counter += 1 # ...then start the next interior accumulator for a new hole.
                    active_list = interiors[interior_counter]
                active_list.append(points[i_idx])                
            else:
                active_list.append(points[i_idx]) # Only need i_idx b/c shapely auto-closes polygons
            prev_facet = facet
        
        exterior_geometry = Polygon(exterior)
        interior_polys = [Polygon(interior) for interior in interiors]
        interior_geometry = MultiPolygon(interior_polys)
        geometry = Geometry(exterior_geometry - interior_geometry)
        return geometry

    def create_facets_and_control_points(self):
        self.perimeter = None
        self.perimeter = list(range(len(self.geom.exterior.coords)))
        self.holes = []
        self.points = []
        self.facets = []
        self.points, self.facets = create_points_and_facets(self.geom)
        self.control_points = tuple(self.geom.representative_point().coords)

        for hole in self.geom.interiors:
            hole_polygon = Polygon(hole)
            self.holes += tuple(hole_polygon.representative_point().coords)
        return

    def compile_geometry(self): # Alias
        self.create_facets_and_control_points()


    def create_mesh(self, mesh_size: float):
        """Creates a quadratic triangular mesh from the Geometry object.

        :param mesh_size: A float describing the maximum mesh element area to be used
        within the Geometry-object finite-element mesh.
        :type mesh_size: float

        :return: Geometry-object with mesh data stored in .mesh attribute. Returned
        Geometry-object is self, not a new instance.
        :rtype: :class:`sectionproperties.pre.sections.Geometry`

        The following example creates a circular cross-section with a diameter of 50 with 64
        points, and generates a mesh with a maximum triangular area of 2.5::

            import sectionproperties.pre.sections as sections

            geometry = sections.CircularSection(d=50, n=64)
            geometry = geometry.create_mesh(mesh_size=2.5)

        ..  figure:: ../images/sections/circle_mesh.png
            :align: center
            :scale: 75 %

            Mesh generated from the above geometry.
        """
        self.mesh = pre.create_mesh(
            self.points, self.facets, self.holes, self.control_points, mesh_size
            )
        return self

    def align_left(self, align_to: Geometry, inner: bool = False):
        """
        Returns a new Geometry object, translated in x, so that the right-most point 
        of the new object will be aligned to left-most point of the other Geometry object.

        :param align_to: Another Geometry to align to.
        :type align_to: sectionproperties.pre.sections.Geometry

        :param inner: Default False. If True, align the left-most point of this
        object to the left-most point of 'align_to'. 
        :type align_to: bool

        :return: Geometry object translated to new alignment
        :rtype: :class:`sections.pre.sections.Geometry`
        """
        self_extents = self.calculate_extents()
        align_to_extents = align_to.calculate_extents()
        self_align_x = self_extents[1] # max x
        if inner: 
            self_align_x = self_extents[0] # min x
        align_to_min_x = align_to_extents[0]
        x_offset = align_to_min_x - self_align_x
        new_geom = self.shift_section(x_offset=x_offset)
        return new_geom

    def align_top(self, align_to: Geometry, inner: bool = False):
        """
        Returns a new Geometry object, tranlsated in y, so that the bottom-most point
        of the new object will be aligned to top-most point of the other Geometry object.

        :param align_to: Another Geometry to align to.
        :type align_to: sectionproperties.pre.sections.Geometry

        :param inner: Default False. If True, align the top-most point of this
        object to the top-most point of 'align_to'. 
        :type align_to: bool

        :return: Geometry object translated to new alignment
        :rtype: :class:`sections.pre.sections.Geometry`
        """
        self_extents = self.calculate_extents()
        align_to_extents = align_to.calculate_extents()
        self_align_y = self_extents[2] # min y
        if inner: 
            self_align_y = self_extents[3] # max y
        align_to_max_y = align_to_extents[3]
        y_offset = align_to_max_y - self_align_y
        new_geom = self.shift_section(y_offset=y_offset)
        return new_geom

    def align_right(self, align_to: Geometry, inner: bool = False):
        """
        Returns a new Geometry object, tranlsated in x, so that the left-most point
        of the new object will be aligned to right-most point of the other Geometry object.

        :param align_to: Another Geometry to align to.
        :type align_to: sectionproperties.pre.sections.Geometry

        :param inner: Default False. If True, align the top-most point of this
        object to the top-most point of 'align_to'. 
        :type align_to: bool

        :return: Geometry object translated to new alignment
        :rtype: :class:`sections.pre.sections.Geometry`
        """
        self_extents = self.calculate_extents()
        align_to_extents = align_to.calculate_extents()
        self_align_x = self_extents[0] # min x
        if inner: 
            self_align_x = self_extents[1] # max x
        align_to_max_x = align_to_extents[1]
        x_offset = align_to_max_x - self_align_x
        new_geom = self.shift_section(x_offset=x_offset)
        return new_geom

    def align_bottom(self, align_to: Geometry, inner: bool = False):
        """
        Returns a new Geometry object, tranlsated in y, so that the top-most point 
        of the new object will be aligned to bottom-most of the other Geometry object.

        :param align_to: Another Geometry to align to.
        :type align_to: sectionproperties.pre.sections.Geometry

        :param inner: Default False. If True, align the bottom-most point of this
        object to the bottom-most point of 'align_to'. 
        :type align_to: bool

        :return: Geometry object translated to new alignment
        :rtype: :class:`sections.pre.sections.Geometry`
        """
        self_extents = self.calculate_extents()
        align_to_extents = align_to.calculate_extents()
        self_align_y = self_extents[3] # max y
        if inner:
            self_align_y = self_extents[2] # min y
        align_to_min_y = align_to_extents[2]
        y_offset = align_to_min_y - self_align_y
        new_geom = self.shift_section(y_offset=y_offset)
        return new_geom

    def align_center(self, align_to: Optional[Geometry] = None):
        """
        Returns a new Geometry object, tranlsated in both x and y, so that the 
        center-point of the new object's centroid will be aligned with 
        centroid of the object in 'align_to'. If 'align_to' is None then the new
        object will be aligned with it's centroid at the origin.

        :param align_to: Another Geometry to align to or None (default is None)
        :type align_to: Optional[sectionproperties.pre.sections.Geometry]

        :return: Geometry object translated to new alignment
        :rtype: :class:`sections.pre.sections.Geometry`
        """
        cx, cy = list(self.geom.centroid.coords)[0]
        if align_to is None:
            shift_x, shift_y = -cx, -cy
        else:
            align_cx, align_cy = list(align_to.geom.centroid.coords)[0]
            shift_x = align_cx - cx
            shift_y = align_cy - cy
        new_geom = self.shift_section(x_offset=shift_x, y_offset=shift_y)
        return new_geom


    def shift_section(self, x_offset=0., y_offset=0.,):
        """
        Returns a new Geometry object translated by 'x_offset' and 'y_offset'.

        :param x_offset: Distance in x-direction by which to shift the geometry.
        :type x_offset: float
        :param y_offset: Distance in y-direction by which to shift the geometry.
        :type y_offset: float

        :return: New Geometry-object shifted by 'x_offset' and 'y_offset'
        :rtype: :class:`sections.pre.sections.Geometry`
        """
        new_geom = Geometry(shapely.affinity.translate(self.geom, x_offset, y_offset), self.material)
        return new_geom


    def rotate_section(self, angle: float, rot_point: Union[List[float], str] = "center", use_radians: bool=False):
        """Rotates the geometry and specified angle about a point. If the rotation point is not
        provided, rotates the section about the center of the geometry's bounding box.

        :param float angle: Angle (degrees by default) by which to rotate the section. A positive angle leads
            to a counter-clockwise rotation.
        :param rot_point: Optional. Point *(x, y)* about which to rotate the section. If not provided, will rotate
        about the center of the geometry's bounding box. Default = 'center'.
        :type rot_point: list[float, float]
        :param use_radians: Boolean to indicate whether 'angle' is in degrees or radians. If True, 'angle' is interpreted as radians.
        
        :return: New Geometry-object rotated by 'angle' about 'rot_point'
        :rtype: :class:`sections.pre.sections.Geometry`
        
        The following example rotates a 200UB25 section clockwise by 30 degrees::

            import sectionproperties.pre.sections as sections

            geometry = sections.i_section(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8)
            new_geometry = geometry.rotate_section(angle=-30)
        """
        new_geom = Geometry(shapely.affinity.rotate(self.geom, angle, rot_point, use_radians), self.material)
        return new_geom


    def mirror_section(self, axis: str ='x', mirror_point: Union[List[float], str] = 'center'):
        """Mirrors the geometry about a point on either the x or y-axis. 

        :param string axis: Axis about which to mirror the geometry, *'x'* or *'y'*
        :param mirror_point: Point about which to mirror the geometry *(x, y)*. 
        If no point is provided, mirrors the geometry about the centroid of the shape's bounding box.
        Default = 'center'.
        :type mirror_point: Union[list[float, float], str]

        :return: New Geometry-object mirrored on 'axis' about 'mirror_point'
        :rtype: :class:`sections.pre.sections.Geometry`

        The following example mirrors a 200PFC section about the y-axis and the point (0, 0)::

            import sectionproperties.pre.sections as sections

            geometry = sections.pfc_section(d=200, b=75, t_f=12, t_w=6, r=12, n_r=8)
            new_geometry = geometry.mirror_section(axis='y', mirror_point=[0, 0])
        """
        x_mirror = 1
        y_mirror = 1
        if axis == "x": x_mirror = -x_mirror
        elif axis == "y": y_mirror = -y_mirror
        new_geom = Geometry(shapely.affinity.scale(self.geom, x_mirror, y_mirror, origin=mirror_point), self.material)
        return new_geom


    def split_section(self, 
        point_i: Tuple[float, float], 
        point_j: Optional[Tuple[float, float]] = None,
        vector: Union[Optional[Tuple[float, float]], np.ndarray] = None,
        ) -> Tuple[List[Geometry], List[Geometry]]:
        """Splits, or bisects, the geometry about an infinite line, as defined by two points
        on the line or by one point on the line and a vector. Either 'point_j' or 'vector'
        must be given. If 'point_j' is given, 'vector' is ignored.

        Returns a tuple of two lists each containing new Geometry instances representing the 
        "top" and "bottom" portions, respectively, of the bisected geometry.

        If the line is a vertical line then the "right" and "left" portions, respectively, are 
        returned.

        :param point_i: A tuple of *(x, y)* coordinates to define a first point on the line
        :type point_i: Tuple[float, float]

        :param point_j: Optional. A tuple of *(x, y)* coordinates to define a second point on the line
        :type point_j: Tuple[float, float]

        :param vector: Optional. A tuple or numpy ndarray of *(x, y)* components to define the line direction.
        :type vector: Union[Tuple[float, float], numpy.ndarray]

        :return: A tuple of lists containing Geometry objects that are bisected about the 
        infinite line defined by the two given points. The first item in the tuple represents
        the geometries on the "top" of the line (or to the "right" of the line, if vertical) and
        the second item represents the geometries to the "bottom" of the line (or 
        to the "left" of the line, if vertical).

        :rtype: Tuple[List[Geometry], List[Geometry]]

        The following example splits a 200PFC section about the y-axis::

            import sectionproperties.pre.sections as sections
            from shapely.geometry import LineString

            geometry = sections.pfc_section(d=200, b=75, t_f=12, t_w=6, r=12, n_r=8)
            right_geom, left_geom = geometry.split_section((0, 0), (0, 1))
        """
        if point_j:
            vector = np.array([
                point_j[0] - point_i[0], 
                point_j[1] - point_i[1],
                ])
        elif vector is not None:
            vector = np.array(vector)
        elif not point_j and not vector:
            raise ValueError("Either a second point or a vector must be given to define the line.")
        bounds = self.calculate_extents()
        line_segment = bisect.create_line_segment(point_i, vector, bounds)
        top_right_polys, bottom_left_polys = bisect.group_top_and_bottom_polys(
            split(self.geom, line_segment), line_segment
        )

        # Create new Geometry instances from polys, preserve original material assignments
        top_right_geoms = [Geometry(poly, self.material) for poly in top_right_polys]
        bottom_left_geoms = [Geometry(poly, self.material) for poly in bottom_left_polys]

        return (top_right_geoms, bottom_left_geoms)


    def offset_section_perimeter(self, amount:float = 0, resolution: float = 12):
        """Dilates or erodes the section perimeter by a discrete amount. 

        :param amount: Distance to offset the section by. A -ve value "erodes" the section.
        A +ve value "dilates" the section.
        :type amount: float
        :param resolution: Number of segments used to approximate a quarter circle around a point
        :type resolution: float

        :return: Geometry object translated to new alignment
        :rtype: :class:`sections.pre.sections.Geometry`

        The following example erodes a 200PFC section by 3::

            import sectionproperties.pre.sections as sections

            geometry = sections.pfc_section(d=200, b=75, t_f=12, t_w=6, r=12, n_r=8)
            new_geometry = geometry.offset_section_perimeter(amount=-3)
        """
        new_geom = self.geom.buffer(
            distance=amount, 
            join_style=1, 
            resolution=resolution
            )
        if isinstance(new_geom, MultiPolygon):
            compound_geom = CompoundGeometry([Geometry(poly, self.material) for poly in new_geom])
            return compound_geom
        single_geom = Geometry(new_geom, self.material)
        return single_geom


    def plot_geometry(self, ax=None, pause=True, labels=False, perimeter=False):
        """Plots the geometry defined by the input section. If no axes object is supplied a new
        figure and axis is created.

        :param ax: Axes object on which the mesh is plotted
        :type ax: :class:`matplotlib.axes.Axes`
        :param bool pause: If set to true, the figure pauses the script until the window is closed.
            If set to false, the script continues immediately after the window is rendered.
        :param bool labels: If set to true, node and facet labels are displayed
        :param bool perimeter: If set to true, boldens the perimeter of the cross-section

        :return: Matplotlib figure and axes objects (fig, ax)
        :rtype: (:class:`matplotlib.figure.Figure`, :class:`matplotlib.axes`)

        The following example creates a CHS discretised with 64 points, with a diameter of 48 and
        thickness of 3.2, and plots the geometry::

            import sectionproperties.pre.sections as sections

            geometry = sections.chs(d=48, t=3.2, n=64)
            geometry.plot_geometry()

        ..  figure:: ../images/sections/chs_geometry.png
            :align: center
            :scale: 75 %

            Geometry generated by the above example.
        """
        # if no axes object is supplied, create and setup the plot
        fig = None
        if ax is None:
            ax_supplied = False
            (fig, ax) = plt.subplots()
            post.setup_plot(ax, pause)
        else:
            ax_supplied = True

        for (i, f) in enumerate(self.facets):
            if perimeter:
                if i in self.perimeter:
                    linewidth = 3
                else:
                    linewidth = 1.5
            else:
                linewidth = 1.5

            # plot the points and facets
            if i == 0:
                ax.plot([self.points[f[0]][0], self.points[f[1]][0]],
                        [self.points[f[0]][1], self.points[f[1]][1]],
                        'ko-', markersize=2, linewidth=linewidth, label='Points & Facets')
            else:
                ax.plot([self.points[f[0]][0], self.points[f[1]][0]],
                        [self.points[f[0]][1], self.points[f[1]][1]],
                        'ko-', markersize=2, linewidth=linewidth)

        for (i, h) in enumerate(self.holes):
            # plot the holes
            if i == 0:
                ax.plot(h[0], h[1], 'rx', markersize=5, label='Holes')
            else:
                ax.plot(h[0], h[1], 'rx', markersize=5)

        for (i, cp) in enumerate(self.control_points):
            # plot the control points
            if i == 0:
                ax.plot(cp[0], cp[1], 'bo', markersize=5,
                        label='Control Points')
            else:
                ax.plot(cp[0], cp[1], 'bo', markersize=5)


        # display the labels
        if labels:
            # plot control_point labels
            # With shapely, it will be useful to have numbered regions
            # to match with lists of mesh_sizes and Materials
            # With shapely, enumerated points and facets becomes less useful
            for (i, pt) in enumerate(self.control_points):
                ax.annotate(str(i), xy=pt, color='b')

            # for (i, pt) in enumerate(self.points):
            #     ax.annotate(str(i), xy=pt, color='r')

            # # plot facet labels
            # for (i, fct) in enumerate(self.facets):
            #     pt1 = self.points[fct[0]]
            #     pt2 = self.points[fct[1]]
            #     xy = [(pt1[0] + pt2[0]) / 2, (pt1[1] + pt2[1]) / 2]

            #     ax.annotate(str(i), xy=xy, color='b')

        # if no axes object is supplied, finish the plot
        if not ax_supplied:
            post.finish_plot(ax, pause, title='Cross-Section Geometry')
            return (fig, ax)
        return (fig, ax)

    def calculate_extents(self):
        """Calculates the minimum and maximum x and y-values amongst the list of points; 
        the points that describe the bounding box of the Geometry instance.

        :return: Minimum and maximum x and y-values *(x_min, x_max, y_min, y_max)*
        :rtype: tuple(float, float, float, float)
        """
        min_x, min_y, max_x, max_y = self.geom.bounds
        return (min_x, max_x, min_y, max_y)

    def calculate_perimeter(self):
        """Calculates the exterior perimeter of the geometry.

        :return: Geometry perimeter.
        :rtype: float
        """
        perimeter = self.geom.exterior.length
        return perimeter

    def calculate_area(self):
        """Calculates the area of the geometry.

        :return: Geometry area.
        :rtype: float
        """
        area = self.geom.area
        return area

    def calculate_centroid(self):
        """Calculates the centroid of the geometry as a tuple
        of *(x,y)* coordinates.

        :return: Geometry centroid.
        :rtype: Tuple[float, float]
        """
        cx, cy = self.geom.centroid.coords[0]
        return (cx, cy)

    @property
    def recovery_points(self):
        """
        Returns four stress recovery points for the section geometry. If the Geometry instance was
        created by a NASTRAN geometry function, e.g. sectionproperties.pre.nastran_sections.nastran_bar(),
        then the recovery points will be pre-set on the Geometry instance.
        """
        return self._recovery_points

    
    @recovery_points.setter
    def recovery_points(self, new_points: Union[List[list], List[tuple], List[Point]]) -> list:
        # The points are in the .geom polygon
        # intersection_exists = (MultiPoint(new_points) & self.geom).equals(MultiPoint(new_points))
        # only_four_points = len(new_points) == 4
        # if intersection_exists and only_four_points:
        self._recovery_points = new_points
        # elif intersection_exists and not only_four_points:
        #     raise ValueError("There must be exactly four recovery points")
        # elif not intersection_exists and only_four_points:
        #     raise ValueError("Not all of the points entered exist on the current geometry.")
        # else:
        #     raise ValueError("There must be exactly four recovery points and they must all exist on the current geometry.")

    @recovery_points.getter
    def recovery_points(self, new_points: Union[List[list], List[tuple], List[Point]]) -> list:
        return self._recovery_points

    def __or__(self, other):
        """
        Perform union on Geometry objects with the | operator
        """
        try:
            new_polygon = self.geom | other.geom
            if isinstance(new_polygon, MultiPolygon): 
                return CompoundGeometry([Geometry(polygon, self.material) for polygon in new_polygon.geoms])
            return Geometry(new_polygon, self.material)
        except:
            raise ValueError(
        f"Cannot perform 'union' on these two Geometry instances: {self} | {other}"
        )

    def __xor__(self, other):
        """
        Perform symmetric difference on Geometry objects with the ^ operator
        """
        try:
            new_polygon = self.geom ^ other.geom
            if isinstance(new_polygon, MultiPolygon): 
                return CompoundGeometry([Geometry(polygon, self.material) for polygon in new_polygon.geoms])
            return Geometry(new_polygon, self.material)
        except:
            raise ValueError(
        f"Cannot perform 'symmetric difference' on these two Geometry instances: {self} ^ {other}"
        )

    def __sub__(self, other):
        """
        Perform difference on Geometry objects with the - operator
        """
        try:
            new_polygon = self.geom - other.geom
            if isinstance(new_polygon, MultiPolygon): 
                return CompoundGeometry([Geometry(polygon, self.material) for polygon in new_polygon.geoms])
            return Geometry(new_polygon, self.material)
        except:
            raise ValueError(
        f"Cannot perform 'difference' on these two Geometry instances: {self} - {other}"
        )

    def __add__(self, other):
        """
        Combine Geometry objects into a CompoundGeometry using the + operator
        """
        try:
            return CompoundGeometry([self, other])
        except:
            raise ValueError(
            f"Cannot create new CompoundGeometry with these objects: {self} + {other}"
        )


    def __and__(self, other):
        """
        Perform intersection on Geometry objects with the & operator
        """
        try:
            new_polygon = self.geom & other.geom
            if isinstance(new_polygon, MultiPolygon): 
                return CompoundGeometry([Geometry(polygon, self.material) for polygon in new_polygon.geoms])
            return Geometry(new_polygon, self.material)
        except:
            raise ValueError(
        f"Cannot perform 'intersection' on these two Geometry instances: {self} & {other}"
        )

### 

# TODO: Create setter for adding Material to CompoundGeometry
class CompoundGeometry(Geometry):
    """Class for defining a geometry of multiple distinct regions, each potentially
    having different material properties.

    CompoundGeometry instances are composed of multiple Geometry objects. As with
    Geometry objects, CompoundGeometry objects have methods for generating a triangular
    mesh over all geometries, transforming the collection of geometries as though they
    were one (e.g. translation, rotation, and mirroring), and aligning the CompoundGeometry
    to another Geometry (or to another CompoundGeometry).

    CompoundGeometry objects can be created directly between two or more Geometry
    objects by using the + operator.

    :cvar geoms: either a list of Geometry objects or a shapely.geometry.MultiPolygon
    instance.
    :vartype geoms: List[sectionproperties.pre.sections.Geometry] or shapely.geometry.MultiPolygon
    """
    def __init__(self, geoms: Union[MultiPolygon, List[Geometry]]):
        if isinstance(geoms, MultiPolygon):
            self.geoms = [Geometry(poly) for poly in geoms.geoms]
            self.geom = geoms
        elif isinstance(geoms, list):
            processed_geoms = []
            for item in geoms:
                if isinstance(item, CompoundGeometry): 
                    # "Flatten" any CompoundGeometry objects in the 'geoms' list
                    processed_geoms += item.geoms 
                elif isinstance(item, Geometry):
                    processed_geoms.append(item)
            self.geoms = processed_geoms
            self.geom = MultiPolygon([geom.geom for geom in processed_geoms])

        self.control_points = []
        self.points = [] 
        self.facets = [] 
        self.holes = [] 
        self.perimeter = []
        self.compile_geometry()

        # self.mesh = None # Previously not a property

    def _repr_svg_(self):
        """
        Returns an svg representation of the CompoundGeometry.
        Wraps shapely.geometry.MultiPolygon._repr_svg_() by returning
        self.geom._repr_svg_()
        """
        print("sectionproperties.pre.sections.CompoundGeometry")
        print(f"object at: {hex(id(self))}")
        return self.geom._repr_svg_()

    @staticmethod
    def from_points(
        points: List[List[float]], 
        facets: Optional[List[List[int]]], 
        holes: Optional[List[List[float]]],
        control_points: Optional[List[List[float]]],
        ):
        """
        An interface for the creation of CompoundGeometry objects through the definition of points, 
        facets, and holes. 

        :cvar points: List of points *(x, y)* defining the vertices of the section geometry.
        If facets are not provided, it is a assumed the that the list of points are ordered
        around the perimeter, either clockwise or anti-clockwise
        :vartype points: list[list[float]]
        :cvar facets: Optional. A list of *(start, end)* indexes of vertices defining the edges
        of the section geoemtry. Can be used to define both external and internal perimeters of holes.
        Facets are assumed to be described in the order of exterior perimeter, interior perimeter 1,
        interior perimeter 2, etc.
        :vartype facets: list[list[int]]
        :cvar holes: Optional. A list of points *(x, y)* that define interior regions as
        being holes or voids. The point can be located anywhere within the hole region.
        Only one point is required per hole region.
        :vartype holes: list[list[float]]
        :cvar holes: Optional. A list of points *(x, y)* that define interior regions as
        being holes or voids. The point can be located anywhere within the hole region.
        Only one point is required per hole region.
        :vartype holes: list[list[float]]
        """
        # First, generate all invidual polygons from points and facets
        current_polygon_points = []
        all_polygons = []
        prev_facet = None
        for facet in facets:
            i_idx, j_idx = facet
            if not prev_facet: # Add the first facet vertex to exterior and move on
                current_polygon_points.append(points[i_idx])
                prev_facet = facet
                continue
            prev_j_idx = prev_facet[1]
            if i_idx != prev_j_idx: #If there is a break in the chain of edges...
                current_polygon_points.append(points[prev_j_idx])
                all_polygons.append(Polygon(current_polygon_points))
                current_polygon_points = [points[i_idx]]
            else:
                current_polygon_points.append(points[i_idx]) # Only need i_idx b/c shapely auto-closes polygons
            prev_facet = facet
        else:
            current_polygon_points.append(points[j_idx])
            all_polygons.append(Polygon(current_polygon_points))

        # Then classify all of the collected polygons as either "exterior" or "interior"
        exteriors = []
        interiors = []
        for polygon in all_polygons:
            hole_coord_in_polygon = any([polygon.contains(Point(hole_coord)) for hole_coord in holes])
            ctrl_coord_in_polygon = any([polygon.contains(Point(ctrl_coord)) for ctrl_coord in control_points])
            if hole_coord_in_polygon and not ctrl_coord_in_polygon:
                interiors.append(polygon)
            else:
                exteriors.append(polygon)

        # Create the holes by subtracting interior regions from exterior regions
        exterior_geometry = MultiPolygon(exteriors)
        if not interiors:
            return CompoundGeometry([Geometry(exterior) for exterior in exteriors])
        if len(interiors) == 1:
            interior_geometry = Polygon(interiors)
            return CompoundGeometry(exterior_geometry - interior_geometry)
        else:
            interior_geometry = MultiPolygon(interiors)
            return CompoundGeometry(exterior_geometry - interior_geometry)


    def shift_section(self, x_offset: float = 0, y_offset: float = 0):
        """
        Returns a new CompoundGeometry object translated by 'x_offset' and 'y_offset'.

        :param x_offset: Distance in x-direction by which to shift the geometry.
        :type x_offset: float
        :param y_offset: Distance in y-direction by which to shift the geometry.
        :type y_offset: float

        :return: CompoundGeometry object shifted by 'x_offset' and 'y_offset'
        :rtype: :class:`sections.pre.sections.CompoundGeometry`
        """
        geoms_acc = []
        for geom in self.geoms:
            geoms_acc.append(geom.shift_section(x_offset=x_offset, y_offset=y_offset))
        new_geom = CompoundGeometry(geoms_acc)
        return new_geom

    def rotate_section(self, angle, rot_point=None, use_radians=False):
        """Rotates the compound geometry and specified angle about a point. If the rotation point is not
        provided, rotates the section about the center of the compound geometry's bounding box.

        :param float angle: Angle (degrees by default) by which to rotate the section. A positive angle leads
            to a counter-clockwise rotation.
        :param rot_point: Optional. Point *(x, y)* about which to rotate the section. If not provided, will rotate
        about the center of the compound geometry's bounding box. Default = 'center'.
        :type rot_point: list[float, float]
        :param use_radians: Boolean to indicate whether 'angle' is in degrees or radians. If True, 'angle' is interpreted as radians.
        
        :return: CompoundGeometry object rotated by 'angle' about 'rot_point'
        :rtype: :class:`sections.pre.sections.CompoundGeometry`
        
        The following example rotates a 200UB25 section clockwise by 30 degrees::

            import sectionproperties.pre.sections as sections

            geometry_1 = sections.i_section(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8)
            geometry_2 = sections.rectangular_section(d=20, b=133)
            compound = geometry_2.align_center(geometry_1).align_top(geometry_1) + geometry_1
            new_compound = compound.rotate_section(angle=-30)
        """
        geoms_acc = []
        for geom in self.geoms:
            geoms_acc.append(geom.rotate_section(angle, rot_point, use_radians))
        new_geom = CompoundGeometry(geoms_acc)
        return new_geom

    def mirror_section(self, axis='x', mirror_point: Union[List[float], str] = 'center'):
        """Mirrors the geometry about a point on either the x or y-axis. 

        :param string axis: Axis about which to mirror the geometry, *'x'* or *'y'*
        :param mirror_point: Point about which to mirror the geometry *(x, y)*. 
        If no point is provided, mirrors the geometry about the centroid of the shape's bounding box.
        Default = 'center'.
        :type mirror_point: Union[list[float, float], str]

        :return: Geometry object mirrored on 'axis' about 'mirror_point'
        :rtype: :class:`sections.pre.sections.Geometry`

        The following example mirrors a 200PFC section about the y-axis and the point (0, 0)::

            import sectionproperties.pre.sections as sections

            geometry_1 = sections.i_section(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8)
            geometry_2 = sections.rectangular_section(d=20, b=133)
            compound = geometry_2.align_center(geometry_1).align_top(geometry_1) + geometry_1
            new_compound = compound.mirror_section(axis='y')
        """
        geoms_acc = []
        for geom in self.geoms:
            geoms_acc.append(geom.mirror_section(axis, mirror_point))
        new_geom = CompoundGeometry(geoms_acc)
        return new_geom


    def split_section(self, 
        point_i: Tuple[float, float], 
        point_j: Optional[Tuple[float, float]] = None,
        vector: Union[Optional[Tuple[float, float]], np.ndarray] = None,
        ) -> Tuple[List[Geometry], List[Geometry]]:
        """Splits, or bisects, the geometry about an infinite line, as defined by two points
        on the line or by one point on the line and a vector. Either 'point_j' or 'vector'
        must be given. If 'point_j' is given, 'vector' is ignored.

        Returns a tuple of two lists each containing new Geometry instances representing the 
        "top" and "bottom" portions, respectively, of the bisected geometry.

        If the line is a vertical line then the "right" and "left" portions, respectively, are 
        returned.

        :param point_i: A tuple of *(x, y)* coordinates to define a first point on the line
        :type point_i: Tuple[float, float]

        :param point_j: Optional. A tuple of *(x, y)* coordinates to define a second point on the line
        :type point_j: Tuple[float, float]

        :param vector: Optional. A tuple or numpy ndarray of *(x, y)* components to define the line direction.
        :type vector: Union[Tuple[float, float], numpy.ndarray]

        :return: A tuple of lists containing Geometry objects that are bisected about the 
        infinite line defined by the two given points. The first item in the tuple represents
        the geometries on the "top" of the line (or to the "right" of the line, if vertical) and
        the second item represents the geometries to the "bottom" of the line (or 
        to the "left" of the line, if vertical).

        :rtype: Tuple[List[Geometry], List[Geometry]]

        The following example splits a 200PFC section about the y-axis::

            import sectionproperties.pre.sections as sections
            from shapely.geometry import LineString

            geometry = sections.pfc_section(d=200, b=75, t_f=12, t_w=6, r=12, n_r=8)
            right_geom, left_geom = geometry.split_section((0, 0), (0, 1))
        """
        top_geoms_acc = []
        bottom_geoms_acc = []
        for geom in self.geoms:
            top_geoms, bottom_geoms = geom.split_section(point_i, point_j, vector)
            top_geoms_acc += top_geoms
            bottom_geoms_acc += bottom_geoms
        return (top_geoms_acc, bottom_geoms_acc)



    def offset_section_perimeter(self, amount:float = 0, resolution: float = 12):
        """Dilates or erodes perimeter of the individual geometries within the CompoundGeometry
        object by a discrete amount. Note, because the individual geometries have their own
        perimeters offset independently, sections don't "stick" as though they were a joined section.
        Any aligned geometries within the CompoundGeometry will need to be re-aligned after
        the perimeter offset.

        :param amount: Distance to offset the section by. A -ve value "erodes" the section. A +ve
        value "dilates" the section.
        :type amount: float
        :param resolution: Number of segments used to approximate a quarter circle around a point
        :type resolution: float

        :return: Geometry object translated to new alignment
        :rtype: :class:`sections.pre.sections.Geometry`

        The following example erodes a 200PFC section by 3::

            import sectionproperties.pre.sections as sections

            geometry_1 = sections.i_section(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8)
            geometry_2 = sections.rectangular_section(d=20, b=133)
            compound = geometry_2.align_center(geometry_1).align_top(geometry_1) + geometry_1
            new_geometry = compound.offset_section_perimeter(amount=-3)
        """
        geoms_acc = []
        for geom in self.geoms:
            geoms_acc.append(geom.offset_section_perimeter(amount, resolution))
        new_geom = CompoundGeometry(geoms_acc)
        return new_geom

    def compile_geometry(self):
        """
        Converts the shapely.geometry.Polygon objects stored in self.geoms into lists of 
        points, facets, control_points, and hole points.
        """
        self.points = []
        self.facets = []
        self.control_points = []
        self.holes = []

        # loop through all sections
        for geom in self.geoms:
            if not all([geom.points, geom.facets, geom.control_points]):
                geom.create_facets_and_control_points() # If not previously done

            # add points and count points
            # skip duplicate points
            for point in geom.points:
                if list(point) not in self.points:
                    self.points.append(list(point))

            # map facets from original Polygon points to collected MultiPolygon points
            # b/c points are not being duplicated, have to find the original point and
            # map the facet to the existing MultiPolygon point
            for facet in geom.facets:
                i_pnt, j_pnt = geom.points[facet[0]], geom.points[facet[1]]
                i_pnt_idx = self.points.index(i_pnt)
                j_pnt_idx = self.points.index(j_pnt)
                self.facets.append([i_pnt_idx, j_pnt_idx])

            # add holes
            for hole in geom.holes:
                self.holes.append(tuple(hole))

            # add control points
            for control_point in geom.control_points:
                self.control_points.append(tuple(control_point))

        # Check for holes created inadvertently from combined sections
        unionized_geometry = None
        for geom in self.geoms:
            if unionized_geometry is None:
                unionized_geometry = geom
                continue
            unionized_geometry = unionized_geometry | geom
        if len(unionized_geometry.holes) > len(self.holes):
            inadvertent_holes = [tuple(hole_coords) for hole_coords in unionized_geometry.holes if hole_coords not in self.holes]
            self.holes += inadvertent_holes
            #change x,y coords to tuples, only

    def calculate_perimeter(self):
        """
        Returns the length of the exterior convex hull of the CompoundGeometry.
        """
        return self.geom.convex_hull.exterior.length


### Helper functions for Geometry

def create_facets(points_list: list, connect_back: bool = False, offset: int = 0) -> list:
    """
    Returns a list of lists of integers representing the "facets" connecting
    the list of coordinates in 'loc'. It is assumed that 'loc' coordinates are 
    already in their order of connectivity.
    
    'loc': a list of coordinates
    'connect_back': if True, then the last facet pair will be [len(loc), offset]
    'offset': an integer representing the value that the facets should begin incrementing from.
    """
    idx_peeker = more_itertools.peekable([idx+offset for idx, _ in enumerate(points_list)])
    facets = [[item, idx_peeker.peek(offset)] for item in idx_peeker]
    if connect_back:
        return facets
    return facets[:-1]

def create_exterior_points(shape: Polygon) -> list:
    """
    Return a list of lists representing x,y pairs of the exterior
    perimeter of `polygon`.
    """
    acc = [list(coord) for coord in shape.exterior.coords]
    return acc

def create_interior_points(lr: LinearRing) -> list:
    """
    Return a list of lists representing x,y pairs of the exterior
    perimeter of `polygon`.
    """
    acc = [list(coord) for coord in lr.coords]
    return acc


def create_points_and_facets(shape: Polygon) -> tuple:
    """
    Return a list of lists representing x,y pairs of the exterior
    perimeter of `polygon`.
    """
    master_count = 0
    points = []
    facets = []
    
    # Shape perimeter
    for coords in list(shape.exterior.coords[:-1]): # The last point == first point (shapely)
        points.append(list(coords))
        master_count += 1
    facets += create_facets(points, connect_back=True)
    exterior_count = master_count
    
    # Holes
    for idx, hole in enumerate(shape.interiors):
        break_count = master_count
        int_points = []
        for coords in hole.coords[:-1]: # The last point == first point (shapely)
            int_points.append(list(coords))
            master_count += 1
        
        offset = break_count*(idx > 0) + exterior_count*(idx < 1) # (idx > 0) is like a 'step function'
        facets += create_facets(int_points, connect_back=True, offset=offset)
        points += int_points
        
    return points, facets


def draw_radius(pt: list, r: float, theta: float, n, ccw: bool = True): # Changed 'anti' to ccw to match shapely
    """Adds a quarter radius of points to the points list - centered at point *pt*, with radius
    *r*, starting at angle *theta*, with *n* points. If r = 0, adds pt only.

    :param pt: Centre of radius *(x,y)*
    :type pt: list[float, float]
    :param float r: Radius
    :param float theta: Initial angle
    :param int n: Number of points
    :param bool ccw: Counter-clockwise rotation?
    """
    points = []
    if r == 0:
        points.append(pt)
        return points

    if ccw:
        mult = 1
    else:
        mult = -1

    # calculate radius of points
    for i in range(n):
        # determine angle
        t = theta + mult * i * 1.0 / max(1, n - 1) * np.pi * 0.5

        x = pt[0] + r * np.cos(t)
        y = pt[1] + r * np.sin(t)
        points.append([x, y])
    return points

def rectangular_section(b, d):
    """Constructs a rectangular section with the bottom left corner at the origin *(0, 0)*, with
    depth *d* and width *b*.

    :param float d: Depth (y) of the rectangle
    :param float b: Width (x) of the rectangle

    The following example creates a rectangular cross-section with a depth of 100 and width of 50,
    and generates a mesh with a maximum triangular area of 5::

        import sectionproperties.pre.sections as sections

        geometry = sections.rectangular_section(d=100, b=50)
        mesh = geometry.create_mesh(mesh_sizes=[5])

    ..  figure:: ../images/sections/rectangle_geometry.png
        :align: center
        :scale: 75 %

        Rectangular section geometry.

    ..  figure:: ../images/sections/rectangle_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """
    # min_x = 0
    # min_y = 0
    # max_x = b
    # max_y = d

    # rectangle = box(min_x, min_y, max_x, max_y)
    points = [[0, 0], [b, 0], [b, d], [0, d]]
    rectangle = Polygon(points)
    return Geometry(rectangle)


def circular_section(d: float, n: int):
    """Constructs a solid circle centered at the origin *(0, 0)* with diameter *d* and using *n*
    points to construct the circle.

    :param float d: Diameter of the circle
    :param int n: Number of points discretising the circle

    The following example creates a circular geometry with a diameter of 50 with 64 points,
    and generates a mesh with a maximum triangular area of 2.5::

        import sectionproperties.pre.sections as sections

        geometry = sections.circular_section(d=50, n=64)
        mesh = geometry.create_mesh(mesh_sizes=[2.5])

    ..  figure:: ../images/sections/circle_geometry.png
        :align: center
        :scale: 75 %

        Circular section geometry.

    ..  figure:: ../images/sections/circle_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """
    x_off, y_off = (0, 0)
    points = []
    # loop through each point on the circle
    for i in range(n):
        # determine polar angle
        theta = i * 2 * np.pi * 1.0 / n

        # calculate location of the point
        x = 0.5 * d * np.cos(theta) + x_off
        y = 0.5 * d * np.sin(theta) + y_off

        # append the current point to the points list
        points.append([x, y])


    circle = Polygon(points)
    return Geometry(circle)


def circular_hollow_section(d: float, t: float, n: int):
    """Constructs a circular hollow section (CHS) centered at the origin *(0, 0)*, with diameter *d* and
    thickness *t*, using *n* points to construct the inner and outer circles.

    :param float d: Outer diameter of the CHS
    :param float t: Thickness of the CHS
    :param int n: Number of points discretising the inner and outer circles

    The following example creates a CHS discretised with 64 points, with a diameter of 48 and
    thickness of 3.2, and generates a mesh with a maximum triangular area of 1.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.circular_hollow_section(d=48, t=3.2, n=64)
        mesh = geometry.create_mesh(mesh_sizes=[1.0])

    ..  figure:: ../images/sections/chs_geometry.png
        :align: center
        :scale: 75 %

        CHS geometry.

    ..  figure:: ../images/sections/chs_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """
    points_inner = []
    points_outer = []
    # loop through each point of the CHS
    for i in range(n):
        # determine polar angle
        theta = i * 2 * np.pi * 1.0 / n

        # calculate location of outer and inner points
        x_outer = 0.5 * d * np.cos(theta)
        y_outer = 0.5 * d * np.sin(theta)
        x_inner = (0.5 * d - t) * np.cos(theta)
        y_inner = (0.5 * d - t) * np.sin(theta)

        # append the current points to the points list
        points_outer.append([x_outer, y_outer])
        points_inner.append([x_inner, y_inner])

    inner_circle = Polygon(points_inner)
    outer_circle = Polygon(points_outer)
    return Geometry(outer_circle - inner_circle)


def elliptical_section(d_y: float, d_x: float, n: int):
    """Constructs a solid ellipse centered at the origin *(0, 0)* with vertical diameter *d_y* and
    horizontal diameter *d_x*, using *n* points to construct the ellipse.

    :param float d_y: Diameter of the ellipse in the y-dimension
    :param float d_x: Diameter of the ellipse in the x-dimension
    :param int n: Number of points discretising the ellipse

    The following example creates an elliptical cross-section with a vertical diameter of 25 and
    horizontal diameter of 50, with 40 points, and generates a mesh with a maximum triangular area
    of 1.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.elliptical_section(d_y=25, d_x=50, n=40)
        mesh = geometry.create_mesh(mesh_sizes=[1.0])

    ..  figure:: ../images/sections/ellipse_geometry.png
        :align: center
        :scale: 75 %

        Elliptical section geometry.

    ..  figure:: ../images/sections/ellipse_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """
    points = []

    # loop through each point on the ellipse
    for i in range(n):
        # determine polar angle
        theta = i * 2 * np.pi * 1.0 / n

        # calculate location of the point
        x = 0.5 * d_x * np.cos(theta)
        y = 0.5 * d_y * np.sin(theta)

        # append the current point to the points list
        points.append([x, y])

    ellipse = Polygon(points)
    return Geometry(ellipse)    


def elliptical_hollow_section(d_y:float, d_x:float, t:float, n:int):
    """Constructs an elliptical hollow section (EHS) centered at the origin *(0, 0)*, with outer vertical
    diameter *d_y*, outer horizontal diameter *d_x*, and thickness *t*, using *n* points to
    construct the inner and outer ellipses.

    :param float d_y: Diameter of the ellipse in the y-dimension
    :param float d_x: Diameter of the ellipse in the x-dimension
    :param float t: Thickness of the EHS
    :param int n: Number of points discretising the inner and outer ellipses

    The following example creates a EHS discretised with 30 points, with a outer vertical diameter
    of 25, outer horizontal diameter of 50, and thickness of 2.0, and generates a mesh with a
    maximum triangular area of 0.5::

        import sectionproperties.pre.sections as sections

        geometry = sections.elliptical_hollow_section(d_y=25, d_x=50, t=2.0, n=64)
        mesh = geometry.create_mesh(mesh_sizes=[0.5])

    ..  figure:: ../images/sections/ehs_geometry.png
        :align: center
        :scale: 75 %

        EHS geometry.

    ..  figure:: ../images/sections/ehs_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """
    points_inner = []
    points_outer = []
    # loop through each point of the EHS
    for i in range(n):
        # determine polar angle
        theta = i * 2 * np.pi * 1.0 / n

        # calculate location of outer and inner points
        x_outer = 0.5 * d_x * np.cos(theta)
        y_outer = 0.5 * d_y * np.sin(theta)
        x_inner = ((0.5 * d_x) - t) * np.cos(theta)
        y_inner = ((0.5 * d_y) - t) * np.sin(theta)

        # append the current points to the points list
        points_outer.append([x_outer, y_outer])
        points_inner.append([x_inner, y_inner])

    outer = Polygon(points_outer)
    inner = Polygon(points_inner)
    return Geometry(outer - inner)



def rectangular_hollow_section(b: float, d: float, t: float, r_out: float, n_r: int):
    """Constructs a rectangular hollow section (RHS) centered at *(b/2, d/2)*, with depth *d*, width *b*,
    thickness *t* and outer radius *r_out*, using *n_r* points to construct the inner and outer
    radii. If the outer radius is less than the thickness of the RHS, the inner radius is set to
    zero.

    :param float d: Depth of the RHS
    :param float b: Width of the RHS
    :param float t: Thickness of the RHS
    :param float r_out: Outer radius of the RHS
    :param int n_r: Number of points discretising the inner and outer radii

    The following example creates an RHS with a depth of 100, a width of 50, a thickness of 6 and
    an outer radius of 9, using 8 points to discretise the inner and outer radii. A mesh is
    generated with a maximum triangular area of 2.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.rectangular_hollow_section(d=100, b=50, t=6, r_out=9, n_r=8)
        mesh = geometry.create_mesh(mesh_sizes=[2.0])

    ..  figure:: ../images/sections/rhs_geometry.png
        :align: center
        :scale: 75 %

        RHS geometry.

    ..  figure:: ../images/sections/rhs_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """
    points_inner = []
    points_outer = []
    # calculate internal radius
    r_in = max(r_out - t, 0)
    # construct the outer radius points
    points_outer += draw_radius([r_out, r_out], r_out, np.pi, n_r)
    points_outer += draw_radius([b - r_out, r_out], r_out, 1.5 * np.pi, n_r)
    points_outer += draw_radius([b - r_out, d - r_out], r_out, 0, n_r)
    points_outer += draw_radius([r_out, d - r_out], r_out, 0.5 * np.pi, n_r)

    points_inner += draw_radius([t + r_in, t + r_in], r_in, np.pi, n_r)
    points_inner += draw_radius([b - t - r_in, t + r_in], r_in, 1.5 * np.pi, n_r)
    points_inner += draw_radius([b - t - r_in, d - t - r_in], r_in, 0, n_r)
    points_inner += draw_radius([t + r_in, d - t - r_in], r_in, 0.5 * np.pi, n_r)

    outer = Polygon(points_outer)
    inner = Polygon(points_inner)
    return Geometry(outer-inner)


def i_section(d: float, b: float, t_f: float, t_w: float, r: float, n_r: int): # More specific description and less ambiguous? e.g. not an "S" section.
    """Constructs an I-section centered at *(b/2, d/2)*, with depth *d*, width *b*, flange
    thickness *t_f*, web thickness *t_w*, and root radius *r*, using *n_r* points to construct the
    root radius.

    :param float d: Depth of the I-section
    :param float b: Width of the I-section
    :param float t_f: Flange thickness of the I-section
    :param float t_w: Web thickness of the I-section
    :param float r: Root radius of the I-section
    :param int n_r: Number of points discretising the root radius

    The following example creates an I-section with a depth of 203, a width of 133, a flange
    thickness of 7.8, a web thickness of 5.8 and a root radius of 8.9, using 16 points to
    discretise the root radius. A mesh is generated with a maximum triangular area of 3.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.i_section(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=16)
        mesh = geometry.create_mesh(mesh_sizes=[3.0])

    ..  figure:: ../images/sections/isection_geometry.png
        :align: center
        :scale: 75 %

        I-section geometry.

    ..  figure:: ../images/sections/isection_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """
    points = []

    # add first three points
    points.append([0, 0])
    points.append([b, 0])
    points.append([b, t_f])

    # construct the bottom right radius
    pt = [b * 0.5 + t_w * 0.5 + r, t_f + r]
    points += draw_radius(pt, r, 1.5 * np.pi, n_r, False)

    # construct the top right radius
    pt = [b * 0.5 + t_w * 0.5 + r, d - t_f - r]
    points += draw_radius(pt, r, np.pi, n_r, False)

    # add the next four points
    points.append([b, d - t_f])
    points.append([b, d])
    points.append([0, d])
    points.append([0, d - t_f])

    # construct the top left radius
    pt = [b * 0.5 - t_w * 0.5 - r, d - t_f - r]
    points += draw_radius(pt, r, 0.5 * np.pi, n_r, False)

    # construct the bottom left radius
    pt = [b * 0.5 - t_w * 0.5 - r, t_f + r]
    points += draw_radius(pt, r, 0, n_r, False)

    # # add the last point
    points.append([0, t_f])
    i_section = Polygon(points)
    return Geometry(i_section)


def mono_i_section(d, b_t, b_b, t_fb, t_ft, t_w, r, n_r):
    """Constructs a monosymmetric I-section centered at *(max(b_t, b_b)/2, d/2)*, with depth *d*,
    top flange width *b_t*, bottom flange width *b_b*, top flange thickness *t_ft*, top flange
    thickness *t_fb*, web thickness *t_w*, and root radius *r*, using *n_r* points to construct the
    root radius.

    :param float d: Depth of the I-section
    :param float b_t: Top flange width
    :param float b_b: Bottom flange width
    :param float t_ft: Top flange thickness of the I-section
    :param float t_fb: Bottom flange thickness of the I-section
    :param float t_w: Web thickness of the I-section
    :param float r: Root radius of the I-section
    :param int n_r: Number of points discretising the root radius

    The following example creates a monosymmetric I-section with a depth of 200, a top flange width
    of 50, a top flange thickness of 12, a bottom flange width of 130, a bottom flange thickness of
    8, a web thickness of 6 and a root radius of 8, using 16 points to discretise the root radius.
    A mesh is generated with a maximum triangular area of 3.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.mono_i_section(
            d=200, b_t=50, b_b=130, t_ft=12, t_fb=8, t_w=6, r=8, n_r=16
        )
        mesh = geometry.create_mesh(mesh_sizes=[3.0])

    ..  figure:: ../images/sections/monoisection_geometry.png
        :align: center
        :scale: 75 %

        I-section geometry.

    ..  figure:: ../images/sections/monoisection_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """
    points = []
    # calculate central axis
    x_central = max(b_t, b_b) * 0.5

    # add first three points
    points.append([x_central - b_b * 0.5, 0])
    points.append([x_central + b_b * 0.5, 0])
    points.append([x_central + b_b * 0.5, t_fb])

    # construct the bottom right radius
    pt = [x_central + t_w * 0.5 + r, t_fb + r]
    points += draw_radius(pt, r, 1.5 * np.pi, n_r, False)

    # construct the top right radius
    pt = [x_central + t_w * 0.5 + r, d - t_ft - r]
    points += draw_radius(pt, r, np.pi, n_r, False)

    # add the next four points
    points.append([x_central + b_t * 0.5, d - t_ft])
    points.append([x_central + b_t * 0.5, d])
    points.append([x_central - b_t * 0.5, d])
    points.append([x_central - b_t * 0.5, d - t_ft])

    # construct the top left radius
    pt = [x_central - t_w * 0.5 - r, d - t_ft - r]
    points += draw_radius(pt, r, 0.5 * np.pi, n_r, False)

    # construct the bottom left radius
    pt = [x_central - t_w * 0.5 - r, t_fb + r]
    points += draw_radius(pt, r, 0, n_r, False)

    # add the last point
    points.append([x_central - b_b * 0.5, t_fb])

    polygon = Polygon(points)
    return Geometry(polygon)


def tapered_flange_i_section(d, b, t_f, t_w, r_r, r_f, alpha, n_r):
    """Constructs a Tapered Flange I-section centered at *(b/2, d/2)*, with depth *d*, width *b*,
    mid-flange thickness *t_f*, web thickness *t_w*, root radius *r_r*, flange radius *r_f* and
    flange angle *alpha*, using *n_r* points to construct the radii.

    :param float d: Depth of the Tapered Flange I-section
    :param float b: Width of the Tapered Flange I-section
    :param float t_f: Mid-flange thickness of the Tapered Flange I-section (measured at the point
        equidistant from the face of the web to the edge of the flange)
    :param float t_w: Web thickness of the Tapered Flange I-section
    :param float r_r: Root radius of the Tapered Flange I-section
    :param float r_f: Flange radius of the Tapered Flange I-section
    :param float alpha: Flange angle of the Tapered Flange I-section (degrees)
    :param int n_r: Number of points discretising the radii

    The following example creates a Tapered Flange I-section with a depth of 588, a width of 191, a
    mid-flange thickness of 27.2, a web thickness of 15.2, a root radius of 17.8, a flange radius
    of 8.9 and a flange angle of 8°, using 16 points to discretise the radii. A mesh is generated
    with a maximum triangular area of 20.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.tapered_flange_i_section(
            d=588, b=191, t_f=27.2, t_w=15.2, r_r=17.8, r_f=8.9, alpha=8, n_r=16
        )
        mesh = geometry.create_mesh(mesh_sizes=[20.0])

    ..  figure:: ../images/sections/taperedisection_geometry.png
        :align: center
        :scale: 75 %

        I-section geometry.

    ..  figure:: ../images/sections/taperedisection_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """
    points = []

    # calculate alpha in radians
    alpha_rad = np.pi * alpha / 180

    # calculate the height of the flange toe and dimensions of the straight
    x1 = b * 0.25 - t_w * 0.25 - r_f * (1 - np.sin(alpha_rad))
    y1 = x1 * np.tan(alpha_rad)
    x2 = b * 0.25 - t_w * 0.25 - r_r * (1 - np.sin(alpha_rad))
    y2 = x2 * np.tan(alpha_rad)
    y_t = t_f - y1 - r_f * np.cos(alpha_rad)

    # add first two points
    points.append([0, 0])
    points.append([b, 0])

    # construct the bottom right flange toe radius
    if r_f == 0:
        points.append([b, y_t])
    else:
        for i in range(n_r):
            # determine polar angle
            theta = i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)

            # calculate the locations of the radius points
            x = b - r_f + r_f * np.cos(theta)
            y = y_t + r_f * np.sin(theta)

            # append the current points to the points list
            points.append([x, y])

    # construct the bottom right root radius
    if r_r == 0:
        points.append([b * 0.5 + t_w * 0.5, t_f + y2])
    else:
        for i in range(n_r):
            # determine polar angle
            theta = (
                3.0 / 2 * np.pi - alpha_rad) - (i * 1.0 / max(1, n_r - 1) * (
                    np.pi * 0.5 - alpha_rad)
            )

            # calculate the locations of the radius points
            x = b * 0.5 + t_w * 0.5 + r_r + r_r * np.cos(theta)
            y = t_f + y2 + r_r * np.cos(alpha_rad) + r_r * np.sin(theta)

            # append the current points to the points list
            points.append([x, y])

    # construct the top right root radius
    if r_r == 0:
        points.append([b * 0.5 + t_w * 0.5, d - t_f - y2])
    else:
        for i in range(n_r):
            # determine polar angle
            theta = np.pi - i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)

            # calculate the locations of the radius points
            x = b * 0.5 + t_w * 0.5 + r_r + r_r * np.cos(theta)
            y = d - t_f - y2 - r_r * np.cos(alpha_rad) + r_r * np.sin(theta)

            # append the current points to the points list
            points.append([x, y])

    # construct the top right flange toe radius
    if r_f == 0:
        points.append([b, d - y_t])
    else:
        for i in range(n_r):
            # determine polar angle
            theta = (
                3.0 * np.pi / 2 + alpha_rad) + i * 1.0 / max(1, n_r - 1) * (
                np.pi * 0.5 - alpha_rad
            )

            # calculate the locations of the radius points
            x = b - r_f + r_f * np.cos(theta)
            y = d - y_t + r_f * np.sin(theta)

            # append the current points to the points list
            points.append([x, y])

    # add the next two points
    points.append([b, d])
    points.append([0, d])

    # construct the top left flange toe radius
    if r_f == 0:
        points.append([0, d - y_t])
    else:
        for i in range(n_r):
            # determine polar angle
            theta = np.pi + (i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad))

            # calculate the locations of the radius points
            x = r_f + r_f * np.cos(theta)
            y = d - y_t + r_f * np.sin(theta)

            # append the current points to the points list
            points.append([x, y])

    # construct the top left root radius
    if r_r == 0:
        points.append([b * 0.5 - t_w * 0.5, d - t_f - y2])
    else:
        for i in range(n_r):
            # determine polar angle
            theta = (
                np.pi * 0.5 - alpha_rad) - (i * 1.0 / max(1, n_r - 1) * (
                    np.pi * 0.5 - alpha_rad)
            )

            # calculate the locations of the radius points
            x = b * 0.5 - t_w * 0.5 - r_r + r_r * np.cos(theta)
            y = d - t_f - y2 - r_r * np.cos(alpha_rad) + r_r * np.sin(theta)

            # append the current points to the points list
            points.append([x, y])

    # construct the bottom left root radius
    if r_r == 0:
        points.append([b * 0.5 - t_w * 0.5, t_f + y2])
    else:
        for i in range(n_r):
            # determine polar angle
            theta = -i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)

            # calculate the locations of the radius points
            x = b * 0.5 - t_w * 0.5 - r_r + r_r * np.cos(theta)
            y = t_f + y2 + r_r * np.cos(alpha_rad) + r_r * np.sin(theta)

            # append the current points to the points list
            points.append([x, y])

    # construct the bottom left flange toe radius
    if r_f == 0:
        points.append([0, y_t])
    else:
        for i in range(n_r):
            # determine polar angle
            theta = (
                np.pi * 0.5 + alpha_rad) + (i * 1.0 / max(1, n_r - 1) * (
                    np.pi * 0.5 - alpha_rad)
            )

            # calculate the locations of the radius points
            x = r_f + r_f * np.cos(theta)
            y = y_t + r_f * np.sin(theta)

            # append the current points to the points list
            points.append([x, y])

    polygon = Polygon(points)
    return Geometry(polygon)


def channel_section(d, b, t_f, t_w, r, n_r):
    """Constructs a parallel-flange channel (PFC) section with the bottom left corner at the origin *(0, 0)*, with depth *d*,
    width *b*, flange thickness *t_f*, web  thickness *t_w* and root radius *r*, using *n_r* points
    to construct the root radius.

    :param float d: Depth of the PFC section
    :param float b: Width of the PFC section
    :param float t_f: Flange thickness of the PFC section
    :param float t_w: Web thickness of the PFC section
    :param float r: Root radius of the PFC section
    :param int n_r: Number of points discretising the root radius
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a PFC section with a depth of 250, a width of 90, a flange
    thickness of 15, a web thickness of 8 and a root radius of 12, using 8 points to discretise the
    root radius. A mesh is generated with a maximum triangular area of 5.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.channel_section(d=250, b=90, t_f=15, t_w=8, r=12, n_r=8)
        mesh = geometry.create_mesh(mesh_sizes=[5.0])

    ..  figure:: ../images/sections/pfc_geometry.png
        :align: center
        :scale: 75 %

        PFC geometry.

    ..  figure:: ../images/sections/pfc_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """
    points = []

    # add first three points
    points.append([0, 0])
    points.append([b, 0])
    points.append([b, t_f])

    # construct the bottom right radius
    pt = [t_w + r, t_f + r]
    points += draw_radius(pt, r, 1.5 * np.pi, n_r, False)

    # construct the top right radius
    pt = [t_w + r, d - t_f - r]
    points += draw_radius(pt, r, np.pi, n_r, False)

    # add last three points
    points.append([b, d - t_f])
    points.append([b, d])
    points.append([0, d])

    polygon = Polygon(points)
    return Geometry(polygon)


def tapered_flange_channel(d, b, t_f, t_w, r_r, r_f, alpha, n_r):
    """Constructs a Tapered Flange Channel section with the bottom left corner at the origin
    *(0, 0)*, with depth *d*, width *b*, mid-flange thickness *t_f*, web thickness *t_w*, root
    radius *r_r*, flange radius *r_f* and flange angle *alpha*, using *n_r* points to construct the
    radii.

    :param float d: Depth of the Tapered Flange Channel section
    :param float b: Width of the Tapered Flange Channel section
    :param float t_f: Mid-flange thickness of the Tapered Flange Channel section (measured at the
        point equidistant from the face of the web to the edge of the flange)
    :param float t_w: Web thickness of the Tapered Flange Channel section
    :param float r_r: Root radius of the Tapered Flange Channel section
    :param float r_f: Flange radius of the Tapered Flange Channel section
    :param float alpha: Flange angle of the Tapered Flange Channel section (degrees)
    :param int n_r: Number of points discretising the radii

    The following example creates a Tapered Flange Channel section with a depth of 10, a width of
    3.5, a mid-flange thickness of 0.575, a web thickness of 0.475, a root radius of 0.575, a
    flange radius of 0.4 and a flange angle of 8°, using 16 points to discretise the radii. A mesh
    is generated with a maximum triangular area of 0.02::

        import sectionproperties.pre.sections as sections

        geometry = sections.tapered_flange_channel(
            d=10, b=3.5, t_f=0.575, t_w=0.475, r_r=0.575, r_f=0.4, alpha=8, n_r=16
        )
        mesh = geometry.create_mesh(mesh_sizes=[0.02])

    ..  figure:: ../images/sections/taperedchannel_geometry.png
        :align: center
        :scale: 75 %

        Tapered flange channel geometry.

    ..  figure:: ../images/sections/taperedchannel_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """
    points = []

    # calculate alpha in radians
    alpha_rad = np.pi * alpha / 180

    # calculate the height of the flange toe and dimensions of the straight
    x1 = b * 0.5 - t_w * 0.5 - r_f * (1 - np.sin(alpha_rad))
    y1 = x1 * np.tan(alpha_rad)
    x2 = b * 0.5 - t_w * 0.5 - r_r * (1 - np.sin(alpha_rad))
    y2 = x2 * np.tan(alpha_rad)
    y_t = t_f - y1 - r_f * np.cos(alpha_rad)

    # add first two points
    points.append([0, 0])
    points.append([b, 0])

    # construct the bottom right flange toe radius
    if r_f == 0:
        points.append([b, y_t])
    else:
        for i in range(n_r):
            # determine polar angle
            theta = i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)

            # calculate the locations of the radius points
            x = b - r_f + r_f * np.cos(theta)
            y = y_t + r_f * np.sin(theta)

            # append the current points to the points list
            points.append([x, y])

    # construct the bottom right root radius
    if r_r == 0:
        points.append([t_w, t_f + y2])
    else:
        for i in range(n_r):
            # determine polar angle
            theta = (
                3.0 / 2 * np.pi - alpha_rad) - (i * 1.0 / max(1, n_r - 1) * (
                    np.pi * 0.5 - alpha_rad)
            )

            # calculate the locations of the radius points
            x = t_w + r_r + r_r * np.cos(theta)
            y = t_f + y2 + r_r * np.cos(alpha_rad) + r_r * np.sin(theta)

            # append the current points to the points list
            points.append([x, y])

    # construct the top right root radius
    if r_r == 0:
        points.append([t_w, d - t_f - y2])
    else:
        for i in range(n_r):
            # determine polar angle
            theta = np.pi - i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)

            # calculate the locations of the radius points
            x = t_w + r_r + r_r * np.cos(theta)
            y = d - t_f - y2 - r_r * np.cos(alpha_rad) + r_r * np.sin(
                theta)

            # append the current points to the points list
            points.append([x, y])

    # construct the top right flange toe radius
    if r_f == 0:
        points.append([b, d - y_t])
    else:
        for i in range(n_r):
            # determine polar angle
            theta = (
                3.0 * np.pi / 2 + alpha_rad) + (i * 1.0 / max(1, n_r - 1) * (
                    np.pi * 0.5 - alpha_rad)
            )

            # calculate the locations of the radius points
            x = b - r_f + r_f * np.cos(theta)
            y = d - y_t + r_f * np.sin(theta)

            # append the current points to the points list
            points.append([x, y])

    # add the final two points
    points.append([b, d])
    points.append([0, d])

    polygon = Polygon(points)
    return Geometry(polygon)


def tee_section(d, b, t_f, t_w, r, n_r):
    """Constructs a Tee section with the top left corner at *(0, d)*, with depth *d*, width *b*,
    flange thickness *t_f*, web thickness *t_w* and root radius *r*, using *n_r* points to
    construct the root radius.

    :param float d: Depth of the Tee section
    :param float b: Width of the Tee section
    :param float t_f: Flange thickness of the Tee section
    :param float t_w: Web thickness of the Tee section
    :param float r: Root radius of the Tee section
    :param int n_r: Number of points discretising the root radius

    The following example creates a Tee section with a depth of 200, a width of 100, a flange
    thickness of 12, a web thickness of 6 and a root radius of 8, using 8 points to discretise the
    root radius. A mesh is generated with a maximum triangular area of 3.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.tee_section(d=200, b=100, t_f=12, t_w=6, r=8, n_r=8)
        mesh = geometry.create_mesh(mesh_sizes=[3.0])

    ..  figure:: ../images/sections/tee_geometry.png
        :align: center
        :scale: 75 %

        Tee section geometry.

    ..  figure:: ../images/sections/tee_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """
    points = []
    # add first two points
    points.append([b * 0.5 - t_w * 0.5, 0])
    points.append([b * 0.5 + t_w * 0.5, 0])

    # construct the top right radius
    pt = [b * 0.5 + t_w * 0.5 + r, d - t_f - r]
    points += draw_radius(pt, r, np.pi, n_r, False)

    # add next four points
    points.append([b, d - t_f])
    points.append([b, d])
    points.append([0, d])
    points.append([0, d - t_f])

    # construct the top left radius
    pt = [b * 0.5 - t_w * 0.5 - r, d - t_f - r]
    points += draw_radius(pt, r, 0.5 * np.pi, n_r, False)

    polygon = Polygon(points)
    return Geometry(polygon)


def angle_section(d, b, t, r_r, r_t, n_r):
    """Constructs an angle section with the bottom left corner at the origin *(0, 0)*, with depth
    *d*, width *b*, thickness *t*, root radius *r_r* and toe radius *r_t*, using *n_r* points to
    construct the radii.

    :param float d: Depth of the angle section
    :param float b: Width of the angle section
    :param float t: Thickness of the angle section
    :param float r_r: Root radius of the angle section
    :param float r_t: Toe radius of the angle section
    :param int n_r: Number of points discretising the radii

    The following example creates an angle section with a depth of 150, a width of 100, a thickness
    of 8, a root radius of 12 and a toe radius of 5, using 16 points to discretise the radii. A
    mesh is generated with a maximum triangular area of 2.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.angle_section(d=150, b=100, t=8, r_r=12, r_t=5, n_r=16)
        mesh = geometry.create_mesh(mesh_sizes=[2.0])

    ..  figure:: ../images/sections/angle_geometry.png
        :align: center
        :scale: 75 %

        Angle section geometry.

    ..  figure:: ../images/sections/angle_mesh.png
        :align: center
        :scale: 75 %
    """

    points = []

    # add first two points
    points.append([0, 0])
    points.append([b, 0])

    # construct the bottom toe radius
    pt = [b - r_t, t - r_t]
    points += draw_radius(pt, r_t, 0, n_r)

    # construct the root radius
    pt = [t + r_r, t + r_r]
    points += draw_radius(pt, r_r, 1.5 * np.pi, n_r, False)

    # construct the top toe radius
    pt = [t - r_t, d - r_t]
    points += draw_radius(pt, r_t, 0, n_r)

    # add the next point
    points.append([0, d])

    polygon = Polygon(points)
    return Geometry(polygon)


def cee_section(d, b, l, t, r_out, n_r):
    """Constructs a Cee section (typical of cold-formed steel) with the bottom left corner at the 
    origin *(0, 0)*, with depth *d*, width *b*, lip *l*, thickness *t* and outer radius *r_out*,
    using *n_r* points to construct the radius. If the outer radius is less than the thickness 
    of the Cee Section, the inner radius is set to zero.

    :param float d: Depth of the Cee section
    :param float b: Width of the Cee section
    :param float l: Lip of the Cee section
    :param float t: Thickness of the Cee section
    :param float r_out: Outer radius of the Cee section
    :param int n_r: Number of points discretising the outer radius
    :raises Exception: Lip length must be greater than the outer radius

    The following example creates a Cee section with a depth of 125, a width of 50, a lip of 30, a
    thickness of 1.5 and an outer radius of 6, using 8 points to discretise the radius. A mesh is
    generated with a maximum triangular area of 0.25::

        import sectionproperties.pre.sections as sections

        geometry = sections.cee_section(d=125, b=50, l=30, t=1.5, r_out=6, n_r=8)
        mesh = geometry.create_mesh(mesh_sizes=[0.25])

    ..  figure:: ../images/sections/cee_geometry.png
        :align: center
        :scale: 75 %

        Cee section geometry.

    ..  figure:: ../images/sections/cee_mesh.png
        :align: center
        :scale: 75 %
    """
    # ensure the lip length is greater than the outer radius
    if l < r_out:
        raise Exception('Lip length must be greater than the outer radius')

    points = []

    # calculate internal radius
    r_in = max(r_out - t, 0)

    # construct the outer bottom left radius
    points += draw_radius([r_out, r_out], r_out, np.pi, n_r)

    # construct the outer bottom right radius
    points += draw_radius([b - r_out, r_out], r_out, 1.5 * np.pi, n_r)

    if r_out != l:
        # add next two points
        points.append([b, l])
        points.append([b - t, l])

    # construct the inner bottom right radius
    points += draw_radius([b - t - r_in, t + r_in], r_in, 0, n_r, False)

    # construct the inner bottom left radius
    points += draw_radius([t + r_in, t + r_in], r_in, 1.5 * np.pi, n_r, False)

    # construct the inner top left radius
    points += draw_radius([t + r_in, d - t - r_in], r_in, np.pi, n_r, False)

    # construct the inner top right radius
    points += draw_radius(
        [b - t - r_in, d - t - r_in], r_in, 0.5 * np.pi, n_r, False)

    if r_out != l:
        # add next two points
        points.append([b - t, d - l])
        points.append([b, d - l])

    # construct the outer top right radius
    points += draw_radius([b - r_out, d - r_out], r_out, 0, n_r)

    # construct the outer top left radius
    points += draw_radius([r_out, d - r_out], r_out, 0.5 * np.pi, n_r)

    polygon = Polygon(points)
    return Geometry(polygon)


def zed_section(d, b_l, b_r, l, t, r_out, n_r):
    """Constructs a zed section with the bottom left corner at the origin *(0, 0)*, with depth *d*,
    left flange width *b_l*, right flange width *b_r*, lip *l*, thickness *t* and outer radius
    *r_out*, using *n_r* points to construct the radius. If the outer radius is less than the
    thickness of the Zed Section, the inner radius is set to zero.

    :param float d: Depth of the zed section
    :param float b_l: Left flange width of the Zed section
    :param float b_r: Right flange width of the Zed section
    :param float l: Lip of the Zed section
    :param float t: Thickness of the Zed section
    :param float r_out: Outer radius of the Zed section
    :param int n_r: Number of points discretising the outer radius

    The following example creates a zed section with a depth of 100, a left flange width of 40, a
    right flange width of 50, a lip of 20, a thickness of 1.2 and an outer radius of 5, using 8
    points to discretise the radius. A mesh is generated with a maximum triangular area of 0.15::

        import sectionproperties.pre.sections as sections

        geometry = sections.zed_section(d=100, b_l=40, b_r=50, l=20, t=1.2, r_out=5, n_r=8)
        mesh = geometry.create_mesh(mesh_sizes=[0.15])

    ..  figure:: ../images/sections/zed_geometry.png
        :align: center
        :scale: 75 %

        zed section geometry.

    ..  figure:: ../images/sections/zed_mesh.png
        :align: center
        :scale: 75 %
    """
    # ensure the lip length is greater than the outer radius
    if l < r_out:
        raise Exception('Lip length must be greater than the outer radius')

    points = []

    # calculate internal radius
    r_in = max(r_out - t, 0)

    # construct the outer bottom left radius
    points += draw_radius([r_out, r_out], r_out, np.pi, n_r)

    # construct the outer bottom right radius
    points += draw_radius([b_r - r_out, r_out], r_out, 1.5 * np.pi, n_r)

    if r_out != l:
        # add next two points
        points.append([b_r, l])
        points.append([b_r - t, l])

    # construct the inner bottom right radius
    points += draw_radius([b_r - t - r_in, t + r_in], r_in, 0, n_r, False)

    # construct the inner bottom left radius
    points += draw_radius([t + r_in, t + r_in], r_in, 1.5 * np.pi, n_r, False)

    # construct the outer top right radius
    points += draw_radius([t - r_out, d - r_out], r_out, 0, n_r)

    # construct the outer top left radius
    points += draw_radius([t - b_l + r_out, d - r_out], r_out, 0.5 * np.pi, n_r)

    if r_out != l:
        # add the next two points
        points.append([t - b_l, d - l])
        points.append([t - b_l + t, d - l])

    # construct the inner top left radius
    points += draw_radius([2 * t - b_l + r_in, d - t - r_in], r_in, np.pi, n_r, False)

    # construct the inner top right radius
    points += draw_radius([-r_in, d - t - r_in], r_in, 0.5 * np.pi, n_r, False)
    polygon = Polygon(points)
    return Geometry(polygon)

zee_section = zed_section # An alias for our American friends (and friends who use "American English")

def cruciform_section(d, b, t, r, n_r):
    """Constructs a cruciform section centered at the origin *(0, 0)*, with depth *d*, width *b*,
    thickness *t* and root radius *r*, using *n_r* points to construct the root radius.

    :param float d: Depth of the cruciform section
    :param float b: Width of the cruciform section
    :param float t: Thickness of the cruciform section
    :param float r: Root radius of the cruciform section

    The following example creates a cruciform section with a depth of 250, a width of 175, a
    thickness of 12 and a root radius of 16, using 16 points to discretise the radius. A mesh is
    generated with a maximum triangular area of 5.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.cruciform_section(d=250, b=175, t=12, r=16, n_r=16)
        mesh = geometry.create_mesh(mesh_sizes=[5.0])

    ..  figure:: ../images/sections/cruciform_geometry.png
        :align: center
        :scale: 75 %

        Cruciform section geometry.

    ..  figure:: ../images/sections/cruciform_mesh.png
        :align: center
        :scale: 75 %
    """
    points = []

    # add first two points
    points.append([-t * 0.5, -d * 0.5])
    points.append([t * 0.5, -d * 0.5])

    # construct the bottom right radius
    pt = [0.5 * t + r, -0.5 * t - r]
    points += draw_radius(pt, r, np.pi, n_r, False)

    # add the next two points
    points.append([0.5 * b, -t * 0.5])
    points.append([0.5 * b, t * 0.5])

    # construct the top right radius
    pt = [0.5 * t + r, 0.5 * t + r]
    points += draw_radius(pt, r, 1.5 * np.pi, n_r, False)

    # add the next two points
    points.append([t * 0.5, 0.5 * d])
    points.append([-t * 0.5, 0.5 * d])

    # construct the top left radius
    pt = [-0.5 * t - r, 0.5 * t + r]
    points += draw_radius(pt, r, 0, n_r, False)

    # add the next two points
    points.append([-0.5 * b, t * 0.5])
    points.append([-0.5 * b, -t * 0.5])

    # construct the bottom left radius
    pt = [-0.5 * t - r, -0.5 * t - r]
    points += draw_radius(pt, r, 0.5 * np.pi, n_r, False)

    polygon = Polygon(points)
    return Geometry(polygon)

def polygon_section(d, t, n_sides, r_in=0, n_r=1, rot=0):
    """Constructs a regular hollow polygon section centered at *(0, 0)*, with a pitch circle
    diameter of bounding polygon *d*, thickness *t*, number of sides *n_sides* and an optional
    inner radius *r_in*, using *n_r* points to construct the inner and outer radii (if radii is
    specified).

    :param float d: Pitch circle diameter of the outer bounding polygon (i.e. diameter of circle
        that passes through all vertices of the outer polygon)
    :param float t: Thickness of the polygon section wall
    :param float r_in: Inner radius of the polygon corners. By default, if not specified, a polygon
        with no corner radii is generated.
    :param int n_r: Number of points discretising the inner and outer radii, ignored if no inner
        radii is specified
    :param rot: Initial counterclockwise rotation in degrees. By default bottom face is aligned
        with x axis.
    :raises Exception: Number of sides in polygon must be greater than or equal to 3

    The following example creates an Octagonal section (8 sides) with a diameter of 200, a
    thickness of 6 and an inner radius of 20, using 12 points to discretise the inner and outer
    radii. A mesh is generated with a maximum triangular area of 5::

        import sectionproperties.pre.sections as sections

        geometry = sections.polygon_section(d=200, t=6, n_sides=8, r_in=20, n_r=12)
        mesh = geometry.create_mesh(mesh_sizes=[5])

    ..  figure:: ../images/sections/polygon_geometry.png
        :align: center
        :scale: 75 %

        Octagonal section geometry.

    ..  figure:: ../images/sections/polygon_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """
    outer_points = []
    inner_points = []
    
    if n_sides < 3:
        msg = 'n_sides required to be greater than 3 for PolygonSection class'
        raise Exception(msg)

    # initial rotation
    rot = rot * np.pi / 180  # radians

    # determine triangular segment angle
    alpha = 2 * np.pi / n_sides  # radians

    # determine distance from origin to point perpendicular on face of side
    a_out = d / 2 * np.cos(alpha / 2)
    a_in = a_out - t

    # determine side length for outer & inner faces neglecting radii
    side_length_out = d * np.sin(alpha / 2)
    side_length_in = a_in / a_out * side_length_out

    # check limit on internal radii, if exceeded then radii merge to circle
    if r_in > a_in:
        r_in = a_in
        circle = True
    else:
        circle = False

    # calculate external radius, if r_in is zero, r_out also is zero
    if r_in == 0:
        r_out = 0
        n_r = 1
    else:
        r_out = r_in + t

    # equivalent side length of half the corner radii triangular segment
    c_out = r_out * (side_length_out / 2) / a_out
    c_in = r_in * (side_length_in / 2) / a_in

    # determine straight side length between corner radii (if present)
    side_length_straight_out = side_length_out - (2 * c_out)
    side_length_straight_in = side_length_in - (2 * c_in)

    # temp list for repeating geometry
    outer_base_points = []
    inner_base_points = []

    # start at bottom face, constructing one corner radii, then rotate by initial rotation +
    # alpha and repeat for n_side number of times to form full section perimeter

    # construct the first radius (bottom right)
    for i in range(n_r):
        # determine polar angle
        theta = 1 / 2 * np.pi + i * 1.0 / max(1, n_r - 1) * alpha

        # calculate location of inner and outer points
        x_outer = side_length_straight_out / 2 - r_out * np.cos(theta)
        y_outer = -a_out + r_out - r_out * np.sin(theta)
        x_inner = side_length_straight_in / 2 - r_in * np.cos(theta)
        y_inner = -a_in + r_in - r_in * np.sin(theta)

        # append the current temporary points to the temporary points list
        outer_base_points.append([x_outer, y_outer])
        inner_base_points.append([x_inner, y_inner])

    # if radii merged to circle with an outer diameter of a_out then skip last point as causes
    # overlapping end points which causes meshing issues if geometry is not cleaned by user
    # if circle:
    #     outer_base_points = base_points[0:-2]

    # iterate and add subsequent corner radii one point at a time for each side

    for i in range(n_sides):
        for point in outer_base_points:
            point_new = rotate(point, alpha * i + rot)
            outer_points.append(point_new)

        for point in inner_base_points:
            point_new = rotate(point, alpha * i + rot)
            inner_points.append(point_new)

    outer_polygon = Polygon(outer_points)
    inner_polygon = Polygon(inner_points)
    return Geometry(outer_polygon - inner_polygon)


def rotate(point, angle: float):
    """
    Rotate a point counterclockwise by a given angle around origin [0, 0]

    :param list point: Point coordinates to be rotated
    :param float angle: Angle to rotate point coordinates
    :return: Coordinates of rotated point
    :rtype: list[float, float]
    """

    pt_x, pt_y = point

    c = np.cos(angle)
    s = np.sin(angle)

    new_x = c * pt_x - s * pt_y
    new_y = s * pt_x + c * pt_y

    return [new_x, new_y]


def box_girder_section(d, b_t, b_b, t_ft, t_fb, t_w):
    """Constructs a box girder section centered at at *(max(b_t, b_b)/2, d/2)*, with depth *d*, top
    width *b_t*, bottom width *b_b*, top flange thickness *t_ft*, bottom flange thickness *t_fb*
    and web thickness *t_w*.

    :param float d: Depth of the Box Girder section
    :param float b_t: Top width of the Box Girder section
    :param float b_b: Bottom width of the Box Girder section
    :param float t_ft: Top flange thickness of the Box Girder section
    :param float t_fb: Bottom flange thickness of the Box Girder section
    :param float t_w: Web thickness of the Box Girder section

    The following example creates a Box Girder section with a depth of 1200, a top width of 1200, a
    bottom width of 400, a top flange thickness of 16, a bottom flange thickness of 12 and a web
    thickness of 8. A mesh is generated with a maximum triangular area of 5.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.box_girder_section(d=1200, b_t=1200, b_b=400, t_ft=100, t_fb=80, t_w=50)
        mesh = geometry.create_mesh(mesh_sizes=[200.0])

    ..  figure:: ../images/sections/box_girder_geometry.png
        :align: center
        :scale: 75 %

        Box Girder geometry.

    ..  figure:: ../images/sections/box_girder_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """
    outer_points = []
    inner_points = []

    # calculate central axis
    x_c = max(b_t, b_b) * 0.5

    # determine side wall angle
    if b_t < b_b:
        phi_b = np.arctan2(d, 0.5 * (b_b - b_t))
        phi_t = np.pi - phi_b
    else:
        phi_t = np.arctan2(d, 0.5 * (b_t - b_b))
        phi_b = np.pi - phi_t

    # determine inner wall x-offsets
    x_bot = t_fb / np.tan(np.pi - phi_b)
    x_top = t_ft / np.tan(np.pi - phi_t)
    web_x = abs(t_w / np.sin(np.pi - phi_b))

    # add outer points
    outer_points.append([x_c - 0.5 * b_b, 0])
    outer_points.append([x_c + 0.5 * b_b, 0])
    outer_points.append([x_c + 0.5 * b_t, d])
    outer_points.append([x_c - 0.5 * b_t, d])

    # add inner points
    inner_points.append([x_c - 0.5 * b_b - x_bot + web_x, t_fb])
    inner_points.append([x_c + 0.5 * b_b + x_bot - web_x, t_fb])
    inner_points.append([x_c + 0.5 * b_t + x_top - web_x, d - t_ft])
    inner_points.append([x_c - 0.5 * b_t - x_top + web_x, d - t_ft])

    outer_polygon = Polygon(outer_points)
    inner_polygon = Polygon(inner_points)

    return Geometry(outer_polygon - inner_polygon)

def dowel_array(
    b: int, 
    d: int, 
    bars_x: int, 
    bars_y: int, 
    cover: float, 
    diam: float, 
    n_r: int = 20, 
    perimeter_only: bool = False,
    ):
    """Constructs an array of circular sections within a space of width 'b' and
    depth, 'd' representing a section of cylindrical dowels.

    :param float b: Width of the total array space in the x direction
    :param float b: Depth of the total array space in the y direction
    :param float bars_x: Number of circles to place within array space distance 'b'
    :param float bars_y: Number of circles to place within array space distance 'd'
    :param float cover: A +ve distance that serves as a buffer between the edges of 'b' and 'd'
    and the tangent of the outer circles that would normally touch 'b' and 'd'.
    :param float diam: The diameter of each circle in the array
    :param float n_r: Number of points used to discretize each circular section
    :param bool perimeter_only: If True, only return circular sections around the perimeter
    of the array space. Default is False.

    The following example creates an array of circular sections within a space that has a depth
    of 1200, a top width of 800, has four circles in the x direction, six circles in the y
    direction, a cover distance of 50 (applies to all sides), each having a diameter of 25,
    and discretized with 20 points.

        import sectionproperties.pre.sections as sections

        geometry = sections.dowel_array(d=1200, b=800, bars_x=4, bars_y=6, cover=50, diam=25, n_r=20)
        mesh = geometry.create_mesh(mesh_sizes=[10.0])

    ..  figure:: ../images/sections/dowel_array_geometry.png
        :align: center
        :scale: 75 %

        Dowel array geometry.

    ..  figure:: ../images/sections/dowel_array_mesh.png
        :align: center
        :scale: 75 %

        Mesh generated from the above geometry.
    """
    total_x = b - 2*cover - diam
    total_y = d - 2*cover - diam
    spacing_x = total_x / (bars_x - 1)
    spacing_y = total_y / (bars_y - 1)
    bars_acc = []
    for x_pos in range(bars_x):
        for y_pos in range(bars_y):
            if perimeter_only: # Use circular section from Robbie with optional center location
                if (0 < x_pos < bars_x - 1) and (0 < y_pos < bars_y - 1):
                        continue
                else:
                    bars_acc.append(
                        circular_section(diam, n_r, [cover + diam/2 + spacing_x*x_pos, cover + diam/2 + spacing_y*y_pos])
                    )
            else:
                bars_acc.append(
                    circular_section(diam, n_r, [cover + diam/2 + spacing_x*x_pos, cover + diam/2 + spacing_y*y_pos])
                )
    return CompoundGeometry(bars_acc)
