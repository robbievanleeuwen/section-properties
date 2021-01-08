import numpy as np
import meshpy.triangle as triangle


class GeometryError(Exception):
    """Exception raised when invalid geometry is found."""


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

    def __init__(self, name, elastic_modulus, poissons_ratio, yield_strength,
                 color='w'):
        """Inits the Material class"""

        self.name = name
        self.elastic_modulus = elastic_modulus
        self.poissons_ratio = poissons_ratio
        self.shear_modulus = elastic_modulus / (2 * (1 + poissons_ratio))
        self.yield_strength = yield_strength
        self.color = color


class GeometryCleanerMixin:
    """Class for cleaning :class:`~sectionproperties.pre.sections.Geometry` objects.

    Provides methods to clean various aspects of the geometry including:

    * Zipping nodes - Find nodes that are close together (relative and absolute tolerance) and
      deletes one of the nodes and rejoins the facets to the remaining node.
    * Removing zero length facets - Removes facets that start and end at the same point.
    * Remove duplicate facets - Removes facets that have the same starting and ending point as an
      existing facet.
    * Removing overlapping facets - Searches for facets that overlap each other, given a tolerance
      angle, and reconstructs a unique set of facets along the overlapping region.
    * Remove unused points - Removes points that are not connected to any facets.
    * Intersect facets - Searches for intersections between two facets and adds the intersection
      point to the points list and splits the intersected facets.

    Note that a geometry cleaning method is provided to all
    :class:`~sectionproperties.pre.sections.Geometry` objects.

    :cvar geometry: Geometry object to clean
    :vartype geometry: :class:`~sectionproperties.pre.sections.Geometry`
    :cvar bool verbose: If set to true, information related to the geometry cleaning process is
        printed to the terminal.

    The following example creates a back-to-back 200PFC geometry, rotates the geometry by 30
    degrees, and cleans the geometry before meshing::

        import sectionproperties.pre.sections as sections

        pfc_right = sections.PfcSection(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8)
        pfc_left = sections.PfcSection(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=8)
        pfc_left.mirror_section(axis='y', mirror_point=[0, 0])
        geometry = sections.MergedSection([pfc_left, pfc_right])
        geometry.rotate_section(angle=30)
        geometry.clean_geometry(verbose=True)
        mesh = geometry.create_mesh(mesh_sizes=[5, 5])

    ..  warning:: If the geometry were not cleaned in the previous example, the meshing algorithm
      would crash (most likely return a segment error). Cleaning the geometry is always recommended
      when creating a merged section, which may result in overlapping or intersecting facets, or
      duplicate nodes.
    """

    def zip_points(self, atol=1e-8, rtol=1e-5, verbose=False):
        """Zips points that are close to each other. Searches through the point list and merges two
        points if there are deemed to be sufficiently close. The average value of the coordinates
        is used for the new point. One of the points is deleted from the point list and the facet
        list is updated to remove references to the old points and renumber the remaining point
        indices in the facet list.

        :param float atol: Absolute tolerance for point zipping
        :param float rtol: Relative tolerance (to geometry extents) for point zipping
        """

        idx_to_remove = []

        # determine rtol
        (x_min, x_max, y_min, y_max) = self.calculate_extents()
        geom_range = max(x_max - x_min, y_max - y_min)
        rel_tol = rtol * geom_range

        # loop through the list of points
        for (i, pt1) in enumerate(self.points):
            # check all other points
            for (j, pt2) in enumerate(self.points[i + 1:]):
                # get point indices
                idx_1 = i
                idx_2 = i + j + 1

                # determine distance between two points
                dist = ((pt2[0] - pt1[0]) ** 2 + (pt2[1] - pt1[1]) ** 2) ** 0.5

                # if the points are close together and the point has not already been removed
                if (dist < atol or dist < rel_tol) and idx_2 not in idx_to_remove:
                    # update point1 (average of point1 + point2)
                    pt1[0] = 0.5 * (pt1[0] + pt2[0])
                    pt1[1] = 0.5 * (pt1[1] + pt2[1])

                    # join facets connected to pt2 to pt1 instead
                    self.replace_point_id(idx_2, idx_1)

                    # add pt2 to the list of points to remove
                    idx_to_remove.append(idx_2)

                    if verbose:
                        msg = "Zipped point {0} to point {1}".format(idx_2, idx_1)
                        print(msg)

        # sort list of indices to remove in reverse order so as not to compromise the indices
        idx_to_remove = sorted(idx_to_remove, reverse=True)

        for idx in idx_to_remove:
            self.remove_point_id(idx)

    def remove_zero_length_facets(self, verbose=False):
        """Searches through all facets and removes those that have the same starting and ending
        point."""

        idx_to_remove = []

        # loop through the list of facets
        for (idx, fct) in enumerate(self.facets):
            if fct[0] == fct[1]:
                idx_to_remove.append(idx)

        # sort list of indices to remove in reverse order so as not to compromise the indices
        idx_to_remove = sorted(idx_to_remove, reverse=True)

        for idx in idx_to_remove:
            self.facets.pop(idx)

            if verbose:
                print("Removed zero length facet {0}".format(idx))

    def remove_overlapping_facets(self, verbose=False):
        """Searches through all facet combinations and fixes facets that overlap within a
        tolerance."""

        cleaning = True

        while cleaning:
            # loop through the list of facets
            for (i, fct1) in enumerate(self.facets):
                broken = False

                # check all other facets
                for (j, fct2) in enumerate(self.facets[i + 1:]):
                    # get facet indices
                    idx_1 = i
                    idx_2 = i + j + 1

                    # get facets points
                    # facet 1: p -> p + r
                    p = np.array(self.points[fct1[0]])
                    r = self.points[fct1[1]] - p

                    # facet 2: q -> q + s
                    q = np.array(self.points[fct2[0]])
                    s = self.points[fct2[1]] - q

                    pts = self.is_overlap(p, q, r, s, fct1, fct2)

                    if pts is not None:
                        # delete both facets
                        idx_to_remove = sorted([idx_1, idx_2], reverse=True)
                        for idx in idx_to_remove:
                            self.facets.pop(idx)

                        # add new facets
                        for k in range(len(pts) - 1):
                            self.facets.append([pts[k], pts[k + 1]])

                        # remove duplicate facets
                        self.remove_duplicate_facets(verbose)

                        if verbose:
                            msg = "Removed overlapping facets {0}...".format(idx_to_remove)
                            msg += "Rebuilt with points: {0}".format(pts)
                            print(msg)

                        # break both loops and loop through all facets again
                        broken = True
                        break

                if broken:
                    break

            # if we've arrived at the end without detecting any overlaps
            if not broken:
                cleaning = False

    def remove_unused_points(self, verbose=False):
        """Searches through all facets and removes points that are not connected to any facets."""

        idx_to_remove = []
        facet_flattened = [i for fct in self.facets for i in fct]

        # loop through number of points
        for pt in range(len(self.points)):
            if pt not in facet_flattened:
                idx_to_remove.append(pt)

                if verbose:
                    print("Removed unused point {0}".format(pt))

        # sort list of indices to remove in reverse order so as not to compromise the indices
        idx_to_remove = sorted(idx_to_remove, reverse=True)

        for idx in idx_to_remove:
            self.remove_point_id(idx)

    def intersect_facets(self, verbose=False):
        """Searches through all facet combinations and finds facets that intersect each other. The
        intersection point is added and the facets rebuilt."""

        cleaning = True

        while cleaning:
            # loop through the list of facets
            for (i, fct1) in enumerate(self.facets):
                broken = False

                # check all other facets
                for (j, fct2) in enumerate(self.facets[i + 1:]):
                    # get facet indices
                    idx_1 = i
                    idx_2 = i + j + 1

                    # get facets points
                    # facet 1: p -> p + r
                    p = np.array(self.points[fct1[0]])
                    r = self.points[fct1[1]] - p

                    # facet 2: q -> q + s
                    q = np.array(self.points[fct2[0]])
                    s = self.points[fct2[1]] - q

                    pt = self.is_intersect(p, q, r, s)

                    if pt is not None:
                        # add point
                        self.points.append([pt[0], pt[1]])
                        pt_idx = len(self.points) - 1

                        # delete both facets
                        idx_to_remove = sorted([idx_1, idx_2], reverse=True)
                        for idx in idx_to_remove:
                            self.facets.pop(idx)

                        # rebuild facet 1
                        self.facets.append([fct1[0], pt_idx])
                        self.facets.append([pt_idx, fct1[1]])

                        # rebuild facet 2
                        self.facets.append([fct2[0], pt_idx])
                        self.facets.append([pt_idx, fct2[1]])

                        if verbose:
                            msg = "Intersected facets"
                            msg += " {0} and {1}".format(idx_1, idx_2)
                            msg += " at point: {0}".format(pt)
                            print(msg)

                        # break both loops and loop through all facets again
                        broken = True
                        break

                if broken:
                    break

            # if we've arrived at the end without detecting any overlaps
            if not broken:
                cleaning = False

    def replace_point_id(self, id_old, id_new):
        """Searches all facets and replaces references to point id_old with id_new.

        :param int id_old: Point index to be replaced
        :param int id_new: Point index to replace point id_old
        """

        # loop through all facets
        for (i, facet) in enumerate(self.facets):
            # loop through the point indices defining the facet
            for (j, point_id) in enumerate(facet):
                if point_id == id_old:
                    self.facets[i][j] = id_new

    def remove_point_id(self, point_id):
        """Removes point point_id from the points list and renumbers the references to points after
        point_id in the facet list.

        :param int point_id: Index of point to be removed
        """

        # remove index point_id from the points list
        self.points.pop(point_id)

        # renumber facet references to points after point_id
        for (i, facet) in enumerate(self.facets):
            # loop through the point indices defining the facet
            for (j, p_id) in enumerate(facet):
                # if the point index is greater the point to be deleted
                if p_id > point_id:
                    # decrement the point index
                    self.facets[i][j] -= 1

    @staticmethod
    def is_duplicate_facet(fct1, fct2):
        """Checks to see if to facets are duplicates.

        :param fct1: First facet to compare
        :type fct1: list[int, int]
        :param fct2: Second facet to compare
        :type fct2: list[int, int]
        :return: Whether or not the facets are identical
        :rtype: bool
        """

        # check for a facet duplicate
        return bool(fct1 == fct2 or fct1 == list(reversed(fct2)))

    @staticmethod
    def is_intersect(p, q, r, s):
        """Determines if the line segment p->p+r intersects q->q+s. Implements Gareth Rees's
        answer: https://stackoverflow.com/questions/563198.

        :param p: Starting point of the first line segment
        :type p: :class:`numpy.ndarray` [float, float]
        :param q: Starting point of the second line segment
        :type q: :class:`numpy.ndarray` [float, float]
        :param r: Vector of the first line segment
        :type r: :class:`numpy.ndarray` [float, float]
        :param s: Vector of the second line segment
        :type s: :class:`numpy.ndarray` [float, float]
        :returns: The intersection point of the line segments. If there is no intersection, returns
            None.
        :rtype: :class:`numpy.ndarray` [float, float]
        """

        if np.cross(r, s) != 0:
            # calculate t and u
            t = np.cross(q - p, s) / np.cross(r, s)
            u = np.cross(p - q, r) / np.cross(s, r)

            # modify from closed inequality (<=) to open (<) so end intersections are not picked up
            if (0 < t < 1) and (0 < u < 1):
                return p + t * r

        return None

    @staticmethod
    def is_overlap(p, q, r, s, fct1, fct2):
        """Determines if the line segment p->p+r overlaps q->q+s. Implements Gareth Rees's answer:
        https://stackoverflow.com/questions/563198.

        :param p: Starting point of the first line segment
        :type p: :class:`numpy.ndarray` [float, float]
        :param q: Starting point of the second line segment
        :type q: :class:`numpy.ndarray` [float, float]
        :param r: Vector of the first line segment
        :type r: :class:`numpy.ndarray` [float, float]
        :param s: Vector of the second line segment
        :type s: :class:`numpy.ndarray` [float, float]
        :param fct1: First facets
        :param fct2: Second facets
        :returns: A list containing the points required for facet rebuilding. If there is no
            rebuild to be done, returns None.
        :rtype: list[list[float, float]]
        """

        tol = 1e-3  # minimum angle tolerance (smaller is considered overlap)
        float_tol = 1e-12  # rounding error tolerance

        # relativise tolerance by length of smallest vector
        tol *= min(np.linalg.norm(r), np.linalg.norm(s))

        # are the line segments collinear?
        if abs(np.cross(r, s)) < tol:
            if abs(np.cross(q - p, r)) < tol:
                # CASE 1: two line segments are collinear
                # calculate end points of second segment in terms of the equation of the first line
                # segment (p + t * r)
                if np.dot(s, r) >= 0:
                    t0 = np.dot(q - p, r) / np.dot(r, r)
                    t1 = np.dot(q + s - p, r) / np.dot(r, r)
                else:
                    t0 = np.dot(q + s - p, r) / np.dot(r, r)
                    t1 = np.dot(q - p, r) / np.dot(r, r)

                # check interval [t0, t1] intersects (0, 1)
                if t0 < 1 - float_tol and float_tol < t1:
                    # recalculate t0 and t1 based on original assumptions
                    t0 = np.dot(q - p, r) / np.dot(r, r)
                    t1 = np.dot(q + s - p, r) / np.dot(r, r)

                    t = sorted(list(set([0.0, t0, 1.0, t1])))
                    idx_list = []

                    # loop through new points
                    for pt in t:
                        if pt == 0.0:
                            idx_list.append(fct1[0])
                        elif pt == 1.0:
                            idx_list.append(fct1[1])
                        elif pt == t0:
                            idx_list.append(fct2[0])
                        elif pt == t1:
                            idx_list.append(fct2[1])

                    return idx_list

                # else it is collinear and disjoint

        return None

    def remove_duplicate_facets(self, verbose=False):
        """Searches through all facets and removes facets that are duplicates, independent of the
        point order."""

        idx_to_remove = []

        # loop through the list of facets
        for (i, fct1) in enumerate(self.facets):
            # check all other facets
            for (j, fct2) in enumerate(self.facets[i + 1:]):
                # get facet indices
                idx_1 = i
                idx_2 = i + j + 1

                # check for a duplicate facet that has not already been deleted
                if (self.is_duplicate_facet(fct1, fct2)
                        and idx_2 not in idx_to_remove):
                    idx_to_remove.append(idx_2)

                    if verbose:
                        msg = "Removed duplicate facet: {0}".format(idx_2)
                        msg += " (identical to facet: {0})".format(idx_1)
                        print(msg)

        # sort list of indices to remove in reverse order so as not to compromise the indices
        idx_to_remove = sorted(idx_to_remove, reverse=True)

        for idx in idx_to_remove:
            self.facets.pop(idx)


def check_geometry(points, facets, holes, control_points, atol=1.0e-8):
    """Check that holes and control points do not lie on a facet line.

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
    :param atol: minimum permissable point distance from any section facet
    :type atol: float

    General method was adapted from:
    https://stackoverflow.com/a/54442561
    """
    # make sure points is a numpy array
    points = np.asarray(points)
    check_points = np.array(holes + control_points)

    # select the starting and ending points of each facet
    facet_points = points[facets, :]
    start = facet_points[:, 0, :]
    end = facet_points[:, 1, :]

    # create the vectors for each facet and normalize them to unit vectors
    facet_vectors = end - start
    lengths = np.sqrt(np.einsum('ij,ij->i', facet_vectors, facet_vectors))
    if np.any(lengths == 0.0):
        # if a length is zero, then a facet had the same starting and ending point
        facet_indices = np.nonzero(lengths == 0)[0]
        offending_facets = facet_points[facet_indices]
        facet_point_indices = np.array(facets)[facet_indices]
        message = ''
        for index, facet, point_indices in zip(
            facet_indices, offending_facets, facet_point_indices
        ):
            message += (
                f'\nFacet {index} with points {point_indices.tolist()} at ' f'{facet.tolist()}'
            )
        raise GeometryError(
            f'{offending_facets.shape[0]} facets have the same starting and '
            f'ending point:{message}'
        )
    unit_vectors = np.divide(facet_vectors, lengths.reshape(-1, 1))

    # vectors from each point to check to the start and end of each facet
    vec1 = np.array(check_points).reshape(1, -1, 2) - start.reshape(-1, 1, 2)
    vec2 = np.array(check_points).reshape(1, -1, 2) - end.reshape(-1, 1, 2)

    # signed parallel distance components
    signed_parallel_vec1 = np.einsum('ijk,ik->ij', -vec1, unit_vectors)
    signed_parallel_vec2 = np.einsum('ijk,ik->ij', vec2, unit_vectors)

    # clamped parallel distance
    clamped_distances = np.maximum.reduce(
        [signed_parallel_vec1, signed_parallel_vec2, np.zeros_like(signed_parallel_vec1)]
    )

    # calculate the cross product to get the perpendicular distance
    perpendicular_distances = np.cross(vec1, unit_vectors.reshape(-1, 1, 2))

    # get the minimum clamped distance
    minimum_distances = np.hypot(clamped_distances, perpendicular_distances)

    # make sure that the closest distance is not within the allowed tolerance
    if np.min(minimum_distances) < atol:
        # select the offending points and point type, then raise an exception with that information
        distance_check = (minimum_distances < atol).sum(axis=0)
        point_types = np.array(['Hole'] * len(holes) + ['Control Point'] * len(control_points))
        offending_types = point_types[distance_check > 0]
        offending_coords = check_points[distance_check > 0]
        message = ''
        for point_type, coordinates in zip(offending_types, offending_coords):
            message += f'\n{point_type} @ {coordinates}'
        raise GeometryError(
            f'{offending_coords.shape[0]} points are too close to a facet:' f'{message}'
        )


def create_mesh(points, facets, holes, control_points, mesh_sizes, atol=1.0e-8):
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

    check_geometry(points, facets, holes, control_points, atol=atol)

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
