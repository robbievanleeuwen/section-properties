from typing import List, Optional, Union
import more_itertools
import numpy as np
from shapely.geometry import Polygon, MultiPolygon, LinearRing, box, Point, MultiPoint
import shapely
import matplotlib.pyplot as plt
import sectionproperties.pre.pre as pre
import sectionproperties.post.post as post


class Geometry:
    """Class for defining the geometry of a contiguous section of a single material.

    Provides an interface for the user to specify the geometry defining a section. A method
    is provided for generating a triangular mesh, transforming the section (e.g. translation,
    rotation, perimeter offset, mirroring), aligning the geometry to another geometry, and
    designating stress recovery points.

    :cvar geom: a Polygon object that defines the geometry
    :vartype geom: shapely.geometry.Polygon
    :cvar points: List of points *(x, y)* defining the vertices of the section geometry. If geom
    is provided then points are ignored.
    :vartype points: list[list[float, float]]
    """
    def __init__(self, geom: shapely.geometry.Polygon = None):
        """Inits the Geometry class.
        Old args; control_points, shift
        """
        if geom is None and points is None:
            raise ValueError(f"Either geom")
        if isinstance(geom, MultiPolygon):
            raise ValueError(f"Use CompoundGeometry(...) for a MultiPolygon object.")
        if not isinstance(geom, Polygon):
            raise ValueError(f"Argument is not a valid shapely.geometry.Polygon object: {geom}")
        self.geom = geom
        self.control_points = [] # Given upon instantiation
        self.shift = [] # Given upon instantiation
        self.points = [] # Previously empty list
        self.facets = [] # Previously empty list
        self.holes = [] # Previously empty list
        self.perimeter = [] # Previously empty list
        self._recovery_points = []
        # self.mesh = None # Previously not a property


    def _repr_svg_(self):
        print("sectionproperties.pre.sections.Geometry")
        print(f"object at: {hex(id(self))}")
        return self.geom._repr_svg_()

    @staticmethod
    def from_points(
        points: List[List[float]], 
        facets: Optional[List[List[int]]] = None, 
        holes: Optional[List[List[float]]] = None,
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
        if holes is None and facets is not None:
            raise ValueError(
                "If holes coordinates are provided then facets must also be provided "
                "to distinguish between exterior and interior edges."
                )
        prev_facet = []
        exterior = []
        interiors = [[] for hole in holes] # initialize an empty facet list for every hole
        interior_counter = 0
        active_list = exterior # Like setting a pointer for the list we are accumulating on
        for facet in facets:
            i_idx, _ = facet
            if not prev_facet: # Add the first facet vertex to exterior and move on
                exterior.append(points[i_idx])
                prev_facet = facet
                continue
            if i_idx != prev_facet[1]: #If there is a break in the chain of edges...
                if active_list == exterior: # ...and we were still on the exterior...
                    active_list = interiors[interior_counter] # ... then move to interior
                else: # ...or if we are already in the interiors...
                    interior_counter += 1 # ...then start the next interior region.
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
        if not isinstance(self.geom, list):
            self.perimeter = list(range(len(self.geom.exterior.coords)))
        self.holes = []
        self.points = []
        self.facets = []
        self.points, self.facets = create_points_and_facets(self.geom)
        self.control_points = list(self.geom.representative_point().coords)

        for hole in self.geom.interiors:
            hole_polygon = Polygon(hole)
            self.holes += list(hole_polygon.representative_point().coords)
        return

    def compile_geometry(self): # Alias
        self.create_facets_and_control_points()


    def create_mesh(self, mesh_sizes: Union[float, list]):
        """Creates a quadratic triangular mesh from the Geometry object.

        :param mesh_sizes: A list of maximum element areas corresponding to each region within the
            cross-section geometry.
        :type mesh_size: list[float]

        :return: Object containing generated mesh data
        :rtype: :class:`mesh.triangle.MeshInfo`

        :raises AssertionError: If the number of mesh sizes does not match the number of regions

        The following example creates a circular cross-section with a diameter of 50 with 64
        points, and generates a mesh with a maximum triangular area of 2.5::

            import sectionproperties.pre.sections as sections

            geometry = sections.CircularSection(d=50, n=64)
            mesh = geometry.create_mesh(mesh_sizes=[2.5])

        ..  figure:: ../images/sections/circle_mesh.png
            :align: center
            :scale: 75 %

            Mesh generated from the above geometry.
        """
        self.compile_geometry()
        if isinstance(mesh_sizes, (float, int)): mesh_sizes = [mesh_sizes]*len(self.control_points)

        error_str = "Number of mesh_sizes ({0}), should match the number of regions ({1})".format(
            len(mesh_sizes), len(self.control_points)
        )
        assert(len(mesh_sizes) == len(self.control_points)), error_str

        self.mesh = pre.create_mesh(
            self.points, self.facets, self.holes, self.control_points, mesh_sizes)
        return self.mesh

    def align_left(self, align_to, inner: bool = False):
        """
        Aligns the "right-most" point of 'self' to the "left-most" point of 'align_to'
        """
        self_extents = self.calculate_extents()
        align_to_extents = align_to.calculate_extents()
        self_align_x = self_extents[1] # max x
        if inner: self_align_x = self_extents[0] # min x
        align_to_min_x = align_to_extents[0]
        x_offset = align_to_min_x - self_align_x
        return self.shift_section(x_offset=x_offset)


    def align_top(self, align_to, inner: bool = False):
        """
        Aligns the "bottom-most" point of 'self' to the "top-most" point of 'align_to'
        If 'inner' is True, aligns to the "inside" of the 'align_to' section
        """
        self_extents = self.calculate_extents()
        align_to_extents = align_to.calculate_extents()
        self_align_y = self_extents[2] # min y
        if inner: self_align_y = self_extents[3] # max y
        align_to_max_y = align_to_extents[3]
        y_offset = align_to_max_y - self_align_y
        return self.shift_section(y_offset=y_offset)


    def align_right(self, align_to, inner: bool = False):
        """
        Aligns the "left-most" point of 'self' to the "right-most" point of 'align_to'
        """
        self_extents = self.calculate_extents()
        align_to_extents = align_to.calculate_extents()
        self_align_x = self_extents[0] # min x
        if inner: self_align_x = self_extents[1] # max x
        align_to_max_x = align_to_extents[1]
        x_offset = align_to_max_x - self_align_x
        return self.shift_section(x_offset=x_offset)


    def align_bottom(self, align_to, inner: bool = False):
        """
        Aligns the "top-most" point of 'self' to the "bottom-most" point of 'align_to'
        """
        self_extents = self.calculate_extents()
        align_to_extents = align_to.calculate_extents()
        self_align_y = self_extents[3] # max y
        if inner: self_align_y = self_extents[2] # min y
        align_to_min_y = align_to_extents[2]
        y_offset = align_to_min_y - self_align_y
        return self.shift_section(y_offset=y_offset)

    def align_center(self, align_to):
        """
        Aligns the "bottom-most" point of 'self' to the "top-most" point of 'align_to'
        """
        self_extents = self.calculate_extents()
        align_to_extents = align_to.calculate_extents()
        self_center_x = (self_extents[1] - self_extents[0])/2 + self_extents[0]
        self_center_y = (self_extents[3] - self_extents[2])/2 + self_extents[2]
        align_to_center_x = (align_to_extents[1] - align_to_extents[0])/2 + align_to_extents[0]
        align_to_center_y = (align_to_extents[3] - align_to_extents[2])/2 + align_to_extents[2]
        x_offset = (align_to_center_x - self_center_x)
        y_offset = (align_to_center_y - self_center_y)
        return self.shift_section(x_offset=x_offset, y_offset=y_offset)


    def shift_section(self, x_offset=0., y_offset=0.,):
        """Shifts the cross-section parameters by the class variable vector *shift*."""

        new_geom = Geometry(shapely.affinity.translate(self.geom, x_offset, y_offset))
        # self.control_points = self.geom.representative_point() if self.control_points else None
        return new_geom


    def rotate_section(self, angle, rot_point=[], use_radians=False):
        """Rotates the geometry and specified angle about a point. If the rotation point is not
        provided, rotates the section about the first control point in the list of control points
        of the :class:`~sectionproperties.pre.sections.Geometry` object.

        :param float angle: Angle (degrees by default) by which to rotate the section. A positive angle leads
            to a counter-clockwise rotation.
        :param rot_point: Point *(x, y)* about which to rotate the section
        :type rot_point: list[float, float]
        :param use_radians: Boolean to indicate whether 'angle' is in degrees or radians. If True, 'angle' is interpreted as radians.

        The following example rotates a 200UB25 section clockwise by 30 degrees::

            import sectionproperties.pre.sections as sections

            geometry = sections.ISection(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8)
            geometry.rotate_section(angle=-30)
        """
        if rot_point == []: rot_point = "center"
        new_geom = Geometry(shapely.affinity.rotate(self.geom, angle, rot_point, use_radians))
        # self.control_points = self.geom.representative_point() if self.control_points else None
        return new_geom


    def mirror_section(self, axis='x', mirror_point: Union[List[float], str] = 'center'):
        """Mirrors the geometry about a point on either the x or y-axis. 
        
        Proposed change of behaviour to match shapely: 
        No longer: 
            If no point is provided,
            mirrors the geometry about the first control point in the list of control points of the
            :class:`~sectionproperties.pre.sections.Geometry` object.
        Instead:
            If no point is provided, mirrors the geometry about the centroid of the shape's bounding box.

        :param string axis: Axis about which to mirror the geometry, *'x'* or *'y'*
        :param mirror_point: Point about which to mirror the geometry *(x, y)*. 
        :type mirror_point: Union[list[float, float], str]

        The following example mirrors a 200PFC section about the y-axis and the point (0, 0)::

            import sectionproperties.pre.sections as sections

            geometry = sections.PfcSection(d=200, b=75, t_f=12, t_w=6, r=12, n_r=8)
            geometry.mirror_section(axis='y', mirror_point=[0, 0])
        """
        x_mirror = 1
        y_mirror = 1
        if axis == "x": x_mirror = -x_mirror
        elif axis == "y": y_mirror = -y_mirror
        new_geom = Geometry(shapely.affinity.scale(self.geom, x_mirror, y_mirror, mirror_point))
        return new_geom


    def offset_section_perimeter(self, amount:float = 0, resolution: float = 12):
        """Erodes the section (negative perimeter offset) by `amount` using. 

        :param amount: Distance to erode the section by. A -ve value "erodes" the section. A +ve
        value "dilates" the section.
        :type amount: float
        :param resolution: Number of segments used to approximate a quarter circle around a point
        :type resolution: float

        The following example erodes a 200PFC section by 3::

            import sectionproperties.pre.sections as sections

            geometry = sections.pfc_section(d=200, b=75, t_f=12, t_w=6, r=12, n_r=8)
            geometry.erode_section(amount=-3)
        """
        new_geom = self.geom.buffer(
            distance=amount, 
            join_style=1, 
            resolution=resolution
            )
        return Geometry(new_geom)


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

            geometry = sections.Chs(d=48, t=3.2, n=64)
            geometry.plot_geometry()

        ..  figure:: ../images/sections/chs_geometry.png
            :align: center
            :scale: 75 %

            Geometry generated by the above example.
        """
        # if no axes object is supplied, create and setup the plot
        self.compile_geometry()
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
                ax.plot(h[0], h[1], 'rx', markerSize=5, label='Holes')
            else:
                ax.plot(h[0], h[1], 'rx', markerSize=5)

        for (i, cp) in enumerate(self.control_points):
            # plot the control points
            if i == 0:
                ax.plot(cp[0], cp[1], 'bo', markerSize=5,
                        label='Control Points')
            else:
                ax.plot(cp[0], cp[1], 'bo', markerSize=5)


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
        """Calculates the minimum and maximum x and y-values amongst the list of points.

        :return: Minimum and maximum x and y-values *(x_min, x_max, y_min, y_max)*
        :rtype: tuple(float, float, float, float)
        """
        min_x, min_y, max_x, max_y = self.geom.bounds
        return (min_x, max_x, min_y, max_y)

    def calculate_perimeter(self):
        """Calculates the perimeter of the cross-section by summing the length of all facets in the
        ``perimeter`` class variable.

        :return: Cross-section perimeter, returns 0 if there is no perimeter defined
        :rtype: float
        """
        perimeter = self.geom.exterior.length
        return perimeter

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
        intersection_exists = MultiPoint(new_points) & self.geom == MultiPoint(new_points)
        only_four_points = len(new_points) == 4
        if intersection_exists and only_four_points:
            self._recovery_points = new_points
        elif intersection_exists and not only_four_points:
            raise ValueError("There must be exactly four recovery points")
        elif not intersection_exists and only_four_points:
            raise ValueError("Not all of the points entered exist on the current geometry.")
        else:
            raise ValueError("There must be exactly four recovery points and they must all exist on the current geometry.")

    @recovery_points.getter
    def recovery_points(self, new_points: Union[List[list], List[tuple], List[Point]]) -> list:
        return self._recovery_points

    def __or__(self, other):
        try:
            new_polygon = self.geom | other.geom
            if isinstance(new_polygon, MultiPolygon): 
                return CompoundGeometry([Geometry(polygon) for polygon in new_polygon.geoms])
            return Geometry(new_polygon)
        except:
            raise ValueError(
        f"Cannot perform 'union' on these two Geometry instances: {self} | {other}"
        )

    def __xor__(self, other):
        try:
            new_polygon = self.geom ^ other.geom
            if isinstance(new_polygon, MultiPolygon): 
                return CompoundGeometry([Geometry(polygon) for polygon in new_polygon.geoms])
            return Geometry(new_polygon)
        except:
            raise ValueError(
        f"Cannot perform 'symmetric difference' on these two Geometry instances: {self} ^ {other}"
        )

    def __sub__(self, other):
        try:
            new_polygon = self.geom - other.geom
            if isinstance(new_polygon, MultiPolygon): 
                return CompoundGeometry([Geometry(polygon) for polygon in new_polygon.geoms])
            return Geometry(new_polygon)
        except:
            raise ValueError(
        f"Cannot perform 'difference' on these two Geometry instances: {self} - {other}"
        )

    def __add__(self, other):
        try:
            return CompoundGeometry([self, other])
        except:
            raise ValueError(
            f"Cannot create new CompoundGeometry with these objects: {self} + {other}"
        )


    def __and__(self, other):
        try:
            new_polygon = self.geom & other.geom
            if isinstance(new_polygon, MultiPolygon): 
                return CompoundGeometry([Geometry(polygon) for polygon in new_polygon.geoms])
            return Geometry(new_polygon)
        except:
            raise ValueError(
        f"Cannot perform 'intersection' on these two Geometry instances: {self} & {other}"
        )

### 
class CompoundGeometry(Geometry):
    def __init__(self, geoms: Union[MultiPolygon, List[Geometry]]):
        if isinstance(geoms, MultiPolygon):
            self.geoms = [Geometry(geom) for geom in geoms.geoms]
            self.geom = geoms
        elif isinstance(geoms, list):
            processed_geoms = []
            for item in geoms:
                if isinstance(item, CompoundGeometry):
                    # Add the list of component Geometry objects to this instance
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
        # self.mesh = None # Previously not a property

    def _repr_svg_(self):
        print("sectionproperties.pre.sections.CompoundGeometry")
        print(f"object at: {hex(id(self))}")
        return self.geom._repr_svg_()

    def shift_section(self, x_offset: float = 0, y_offset: float = 0):
        geoms_acc = []
        for geom in self.geoms:
            geoms_acc.append(geom.shift_section(x_offset=x_offset, y_offset=y_offset))
        new_geom = CompoundGeometry(geoms_acc)
        return new_geom

    def rotate_section(self, angle, rot_point=None, use_radians=False):
        geoms_acc = []
        for geom in self.geoms:
            geoms_acc.append(geom.rotate_section(angle, rot_point, use_radians))
        new_geom = CompoundGeometry(geoms_acc)
        return new_geom

    def mirror_section(self, axis='x', mirror_point: Union[List[float], str] = 'center'):
        geoms_acc = []
        for geom in self.geoms:
            geoms_acc.append(geom.mirror_section(axis, mirror_point))
        new_geom = CompoundGeometry(geoms_acc)
        return new_geom

    def offset_section_perimeter(self, amount:float = 0, resolution: float = 12):
        geoms_acc = []
        for geom in self.geoms:
            geoms_acc.append(geom.offset_section_perimeter(amount, resolution))
        new_geom = CompoundGeometry(geoms_acc)
        return new_geom

    def compile_geometry(self):
        point_count = 0
        self.points = []
        self.facets = []
        self.control_points = []
        self.holes = []

        # loop through all sections
        for geom in self.geoms:
            if not all([geom.points, geom.facets, geom.control_points]):
                geom.create_facets_and_control_points() # If not previously done

            # add facets
            for facet in geom.facets:
                self.facets.append([facet[0] + point_count, facet[1] + point_count])

            # add points and count points
            for point in geom.points:
                self.points.append([point[0], point[1]])
                point_count += 1

            # add holes
            for hole in geom.holes:
                self.holes.append([hole[0], hole[1]])

            # add control points
            for control_point in geom.control_points:
                self.control_points.append([control_point[0], control_point[1]])

    def calculate_perimeter(self):
        return self.geom.convex_hull.exterior.length



    



### Helper functions for Geometry

def create_facets(loc: list, connect_back: bool = False, offset: int = 0) -> list:
    """
    Returns a list of lists of integers representing the "facets" connecting
    the list of coordinates in 'loc'. It is assumed that 'loc' coordinates are 
    already in their order of connectivity.
    
    'loc': a list of coordinates
    'connect_back': if True, then the last facet pair will be [len(loc), 0]
    'offset': an integer representing the value that the facets should begin incrementing from.
    """
    idx_peeker = more_itertools.peekable([idx+offset for idx, coords in enumerate(loc)])
    facets = [[item, idx_peeker.peek(0)] for item in idx_peeker]
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
    for coords in list(shape.exterior.coords):
        points.append(list(coords))
        master_count += 1
    facets += create_facets(points)
    exterior_count = master_count # Because increment after last iteration assumes another iteration
    
    # Holes
    for idx, hole in enumerate(shape.interiors):
        break_count = master_count
        int_points = []
        for coords in hole.coords:
            int_points.append(list(coords))
            master_count += 1
        
        offset = break_count*(idx > 0) + exterior_count*(idx < 1) # (idx > 0) is like a 'step function'
        facets += create_facets(int_points, offset = offset)
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
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a rectangular cross-section with a depth of 100 and width of 50,
    and generates a mesh with a maximum triangular area of 5::

        import sectionproperties.pre.sections as sections

        geometry = sections.RectangularSection(d=100, b=50)
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
    min_x = 0 - b/2
    min_y = 0 - d/2
    max_x = b/2
    max_y = d/2

    rectangle = box(min_x, min_y, max_x, max_y)
    return Geometry(rectangle)


def circular_section(d: float, n: int, center: List[float] = [0,0]):
    """Constructs a solid circle centered at the origin *(0, 0)* with diameter *d* and using *n*
    points to construct the circle.

    :param float d: Diameter of the circle
    :param int n: Number of points discretising the circle
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a circular cross-section with a diameter of 50 with 64 points,
    and generates a mesh with a maximum triangular area of 2.5::

        import sectionproperties.pre.sections as sections

        geometry = sections.CircularSection(d=50, n=64)
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
    x_off, y_off = center
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
    """Constructs a circular hollow section centered at the origin *(0, 0)*, with diameter *d* and
    thickness *t*, using *n* points to construct the inner and outer circles.

    :param float d: Outer diameter of the CHS
    :param float t: Thickness of the CHS
    :param int n: Number of points discretising the inner and outer circles
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a CHS discretised with 64 points, with a diameter of 48 and
    thickness of 3.2, and generates a mesh with a maximum triangular area of 1.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.Chs(d=48, t=3.2, n=64)
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
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates an elliptical cross-section with a vertical diameter of 25 and
    horizontal diameter of 50, with 40 points, and generates a mesh with a maximum triangular area
    of 1.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.EllipticalSection(d_y=25, d_x=50, n=40)
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
    """Constructs an elliptical hollow section centered at the origin *(0, 0)*, with outer vertical
    diameter *d_y*, outer horizontal diameter *d_x*, and thickness *t*, using *n* points to
    construct the inner and outer ellipses.

    :param float d_y: Diameter of the ellipse in the y-dimension
    :param float d_x: Diameter of the ellipse in the x-dimension
    :param float t: Thickness of the EHS
    :param int n: Number of points discretising the inner and outer ellipses
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a EHS discretised with 30 points, with a outer vertical diameter
    of 25, outer horizontal diameter of 50, and thickness of 2.0, and generates a mesh with a
    maximum triangular area of 0.5::

        import sectionproperties.pre.sections as sections

        geometry = sections.Ehs(d_y=25, d_x=50, t=2.0, n=64)
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
    """Constructs a rectangular hollow section centered at *(b/2, d/2)*, with depth *d*, width *b*,
    thickness *t* and outer radius *r_out*, using *n_r* points to construct the inner and outer
    radii. If the outer radius is less than the thickness of the RHS, the inner radius is set to
    zero.

    :param float d: Depth of the RHS
    :param float b: Width of the RHS
    :param float t: Thickness of the RHS
    :param float r_out: Outer radius of the RHS
    :param int n_r: Number of points discretising the inner and outer radii
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates an RHS with a depth of 100, a width of 50, a thickness of 6 and
    an outer radius of 9, using 8 points to discretise the inner and outer radii. A mesh is
    generated with a maximum triangular area of 2.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.Rhs(d=100, b=50, t=6, r_out=9, n_r=8)
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
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates an I-section with a depth of 203, a width of 133, a flange
    thickness of 7.8, a web thickness of 5.8 and a root radius of 8.9, using 16 points to
    discretise the root radius. A mesh is generated with a maximum triangular area of 3.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.ISection(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=16)
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
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a monosymmetric I-section with a depth of 200, a top flange width
    of 50, a top flange thickness of 12, a bottom flange width of 130, a bottom flange thickness of
    8, a web thickness of 6 and a root radius of 8, using 16 points to discretise the root radius.
    A mesh is generated with a maximum triangular area of 3.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.MonoISection(
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

        geometry = sections.TaperedFlangeISection(
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


def pfc_section(d, b, t_f, t_w, r, n_r):
    """Constructs a PFC section with the bottom left corner at the origin *(0, 0)*, with depth *d*,
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

        geometry = sections.PfcSection(d=250, b=90, t_f=15, t_w=8, r=12, n_r=8)
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
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a Tapered Flange Channel section with a depth of 10, a width of
    3.5, a mid-flange thickness of 0.575, a web thickness of 0.475, a root radius of 0.575, a
    flange radius of 0.4 and a flange angle of 8°, using 16 points to discretise the radii. A mesh
    is generated with a maximum triangular area of 0.02::

        import sectionproperties.pre.sections as sections

        geometry = sections.TaperedFlangeChannel(
            d=10, b=3.5, t_f=0.575, t_w=0.475, r_r=0.575, r_f=0.4, alpha=8, n_r=16
        )
        mesh = geometry.create_mesh(mesh_sizes=[0.02])

    ..  figure:: ../images/sections/taperedchannel_geometry.png
        :align: center
        :scale: 75 %

        I-section geometry.

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

        geometry = sections.TeeSection(d=200, b=100, t_f=12, t_w=6, r=8, n_r=8)
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
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates an angle section with a depth of 150, a width of 100, a thickness
    of 8, a root radius of 12 and a toe radius of 5, using 16 points to discretise the radii. A
    mesh is generated with a maximum triangular area of 2.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.AngleSection(d=150, b=100, t=8, r_r=12, r_t=5, n_r=16)
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
    """Constructs a Cee section with the bottom left corner at the origin *(0, 0)*, with depth *d*,
    width *b*, lip *l*, thickness *t* and outer radius *r_out*, using *n_r* points to construct the
    radius. If the outer radius is less than the thickness of the Cee Section, the inner radius is
    set to zero.

    :param float d: Depth of the Cee section
    :param float b: Width of the Cee section
    :param float l: Lip of the Cee section
    :param float t: Thickness of the Cee section
    :param float r_out: Outer radius of the Cee section
    :param int n_r: Number of points discretising the outer radius
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]
    :raises Exception: Lip length must be greater than the outer radius

    The following example creates a Cee section with a depth of 125, a width of 50, a lip of 30, a
    thickness of 1.5 and an outer radius of 6, using 8 points to discretise the radius. A mesh is
    generated with a maximum triangular area of 0.25::

        import sectionproperties.pre.sections as sections

        geometry = sections.CeeSection(d=125, b=50, l=30, t=1.5, r_out=6, n_r=8)
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
    """Constructs a Zed section with the bottom left corner at the origin *(0, 0)*, with depth *d*,
    left flange width *b_l*, right flange width *b_r*, lip *l*, thickness *t* and outer radius
    *r_out*, using *n_r* points to construct the radius. If the outer radius is less than the
    thickness of the Zed Section, the inner radius is set to zero.

    :param float d: Depth of the Zed section
    :param float b_l: Left flange width of the Zed section
    :param float b_r: Right flange width of the Zed section
    :param float l: Lip of the Zed section
    :param float t: Thickness of the Zed section
    :param float r_out: Outer radius of the Zed section
    :param int n_r: Number of points discretising the outer radius
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]
    :raises Exception: Lip length must be greater than the outer radius

    The following example creates a Zed section with a depth of 100, a left flange width of 40, a
    right flange width of 50, a lip of 20, a thickness of 1.2 and an outer radius of 5, using 8
    points to discretise the radius. A mesh is generated with a maximum triangular area of 0.15::

        import sectionproperties.pre.sections as sections

        geometry = sections.ZedSection(d=100, b_l=40, b_r=50, l=20, t=1.2, r_out=5, n_r=8)
        mesh = geometry.create_mesh(mesh_sizes=[0.15])

    ..  figure:: ../images/sections/zed_geometry.png
        :align: center
        :scale: 75 %

        Zed section geometry.

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
    :param int n_r: Number of points discretising the root radius
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a cruciform section with a depth of 250, a width of 175, a
    thickness of 12 and a root radius of 16, using 16 points to discretise the radius. A mesh is
    generated with a maximum triangular area of 5.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.CruciformSection(d=250, b=175, t=12, r=16, n_r=16)
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

        geometry = sections.PolygonSection(d=200, t=6, n_sides=8, r_in=20, n_r=12)
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
    if circle:
        base_points = base_points[0:-2]

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


def rotate(point, angle):
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
    """Constructs a Box Girder section centered at at *(max(b_t, b_b)/2, d/2)*, with depth *d*, top
    width *b_t*, bottom width *b_b*, top flange thickness *t_ft*, bottom flange thickness *t_fb*
    and web thickness *t_w*.

    :param float d: Depth of the Box Girder section
    :param float b_t: Top width of the Box Girder section
    :param float b_b: Bottom width of the Box Girder section
    :param float t_ft: Top lange thickness of the Box Girder section
    :param float t_fb: Bottom flange thickness of the Box Girder section
    :param float t_w: Web thickness of the Box Girder section
    :param shift: Vector that shifts the cross-section by *(x, y)*
    :type shift: list[float, float]

    The following example creates a Box Gider section with a depth of 1200, a top width of 1200, a
    bottom width of 400, a top flange thickness of 16, a bottom flange thickness of 12 and a web
    thickness of 8. A mesh is generated with a maximum triangular area of 5.0::

        import sectionproperties.pre.sections as sections

        geometry = sections.BoxGirderSection(d=1200, b_t=1200, b_b=400, t_ft=100, t_fb=80, t_w=50)
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

def dowel_array(b: int, d: int, bars_x: int, bars_y: int, cover: float, bar_diam: float, n_r: int = 20, perimeter_only: bool = False):
    total_x = b - 2*cover - bar_diam
    total_y = d - 2*cover - bar_diam
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
                        circular_section(bar_diam, n_r, [cover + bar_diam/2 + spacing_x*x_pos, cover + bar_diam/2 + spacing_y*y_pos])
                    )
            else:
                bars_acc.append(
                    circular_section(bar_diam, n_r, [cover + bar_diam/2 + spacing_x*x_pos, cover + bar_diam/2 + spacing_y*y_pos])
                )
    return CompoundGeometry(bars_acc)
