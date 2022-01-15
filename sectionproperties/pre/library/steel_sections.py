import numpy as np
from shapely.geometry import Polygon
import sectionproperties.pre.geometry as geometry
import sectionproperties.pre.pre as pre
from sectionproperties.pre.library.utils import draw_radius, rotate


def circular_hollow_section(
    d: float, t: float, n: int, material: pre.Material = pre.DEFAULT_MATERIAL
) -> geometry.Geometry:
    """Constructs a circular hollow section (CHS) centered at the origin *(0, 0)*, with diameter *d* and
    thickness *t*, using *n* points to construct the inner and outer circles.

    :param float d: Outer diameter of the CHS
    :param float t: Thickness of the CHS
    :param int n: Number of points discretising the inner and outer circles
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a CHS discretised with 64 points, with a diameter of 48 and
    thickness of 3.2, and generates a mesh with a maximum triangular area of 1.0::

        from sectionproperties.pre.library.steel_sections import circular_hollow_section

        geometry = circular_hollow_section(d=48, t=3.2, n=64)
        geometry.create_mesh(mesh_sizes=[1.0])

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
    return geometry.Geometry(outer_circle - inner_circle, material)


def elliptical_hollow_section(
    d_y: float,
    d_x: float,
    t: float,
    n: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs an elliptical hollow section (EHS) centered at the origin *(0, 0)*, with outer vertical
    diameter *d_y*, outer horizontal diameter *d_x*, and thickness *t*, using *n* points to
    construct the inner and outer ellipses.

    :param float d_y: Diameter of the ellipse in the y-dimension
    :param float d_x: Diameter of the ellipse in the x-dimension
    :param float t: Thickness of the EHS
    :param int n: Number of points discretising the inner and outer ellipses
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a EHS discretised with 30 points, with a outer vertical diameter
    of 25, outer horizontal diameter of 50, and thickness of 2.0, and generates a mesh with a
    maximum triangular area of 0.5::

        from sectionproperties.pre.library.steel_sections import elliptical_hollow_section

        geometry = elliptical_hollow_section(d_y=25, d_x=50, t=2.0, n=64)
        geometry.create_mesh(mesh_sizes=[0.5])

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
    return geometry.Geometry(outer - inner, material)


def rectangular_hollow_section(
    b: float,
    d: float,
    t: float,
    r_out: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a rectangular hollow section (RHS) centered at *(b/2, d/2)*, with depth *d*, width *b*,
    thickness *t* and outer radius *r_out*, using *n_r* points to construct the inner and outer
    radii. If the outer radius is less than the thickness of the RHS, the inner radius is set to
    zero.

    :param float d: Depth of the RHS
    :param float b: Width of the RHS
    :param float t: Thickness of the RHS
    :param float r_out: Outer radius of the RHS
    :param int n_r: Number of points discretising the inner and outer radii
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates an RHS with a depth of 100, a width of 50, a thickness of 6 and
    an outer radius of 9, using 8 points to discretise the inner and outer radii. A mesh is
    generated with a maximum triangular area of 2.0::

        from sectionproperties.pre.library.steel_sections import rectangular_hollow_section

        geometry = rectangular_hollow_section(d=100, b=50, t=6, r_out=9, n_r=8)
        geometry.create_mesh(mesh_sizes=[2.0])

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
    return geometry.Geometry(outer - inner, material)


def polygon_hollow_section(
    d: float,
    t: float,
    n_sides: int,
    r_in: float = 0,
    n_r: int = 1,
    rot: float = 0,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
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
    :param float rot: Initial counterclockwise rotation in degrees. By default bottom face is
        aligned with x axis.
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry
    :raises Exception: Number of sides in polygon must be greater than or equal to 3

    The following example creates an Octagonal section (8 sides) with a diameter of 200, a
    thickness of 6 and an inner radius of 20, using 12 points to discretise the inner and outer
    radii. A mesh is generated with a maximum triangular area of 5::

        from sectionproperties.pre.library.steel_sections import polygon_hollow_section

        geometry = polygon_hollow_section(d=200, t=6, n_sides=8, r_in=20, n_r=12)
        geometry.create_mesh(mesh_sizes=[5])

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
        msg = "n_sides required to be greater than 3 for polygon_section()"
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

    for i in range(n_sides):
        for point in outer_base_points:
            point_new = rotate(point, alpha * i + rot)
            outer_points.append(point_new)

        for point in inner_base_points:
            point_new = rotate(point, alpha * i + rot)
            inner_points.append(point_new)

    outer_polygon = Polygon(outer_points)
    inner_polygon = Polygon(inner_points)
    return geometry.Geometry(outer_polygon - inner_polygon, material)


def i_section(
    d: float,
    b: float,
    t_f: float,
    t_w: float,
    r: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:  # More specific description and less ambiguous? e.g. not an "S" section.
    """Constructs an I Section centered at *(b/2, d/2)*, with depth *d*, width *b*, flange
    thickness *t_f*, web thickness *t_w*, and root radius *r*, using *n_r* points to construct the
    root radius.

    :param float d: Depth of the I Section
    :param float b: Width of the I Section
    :param float t_f: Flange thickness of the I Section
    :param float t_w: Web thickness of the I Section
    :param float r: Root radius of the I Section
    :param int n_r: Number of points discretising the root radius
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates an I Section with a depth of 203, a width of 133, a flange
    thickness of 7.8, a web thickness of 5.8 and a root radius of 8.9, using 16 points to
    discretise the root radius. A mesh is generated with a maximum triangular area of 3.0::

        from sectionproperties.pre.library.steel_sections import i_section

        geometry = i_section(d=203, b=133, t_f=7.8, t_w=5.8, r=8.9, n_r=16)
        geometry.create_mesh(mesh_sizes=[3.0])

    ..  figure:: ../images/sections/isection_geometry.png
        :align: center
        :scale: 75 %

        I Section geometry.

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
    return geometry.Geometry(i_section, material)


def mono_i_section(
    d: float,
    b_t: float,
    b_b: float,
    t_ft: float,
    t_fb: float,
    t_w: float,
    r: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a monosymmetric I Section centered at *(max(b_t, b_b)/2, d/2)*, with depth *d*,
    top flange width *b_t*, bottom flange width *b_b*, top flange thickness *t_ft*, top flange
    thickness *t_fb*, web thickness *t_w*, and root radius *r*, using *n_r* points to construct the
    root radius.

    :param float d: Depth of the I Section
    :param float b_t: Top flange width
    :param float b_b: Bottom flange width
    :param float t_ft: Top flange thickness of the I Section
    :param float t_fb: Bottom flange thickness of the I Section
    :param float t_w: Web thickness of the I Section
    :param float r: Root radius of the I Section
    :param int n_r: Number of points discretising the root radius
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a monosymmetric I Section with a depth of 200, a top flange width
    of 50, a top flange thickness of 12, a bottom flange width of 130, a bottom flange thickness of
    8, a web thickness of 6 and a root radius of 8, using 16 points to discretise the root radius.
    A mesh is generated with a maximum triangular area of 3.0::

        from sectionproperties.pre.library.steel_sections import mono_i_section

        geometry = mono_i_section(
            d=200, b_t=50, b_b=130, t_ft=12, t_fb=8, t_w=6, r=8, n_r=16
        )
        geometry.create_mesh(mesh_sizes=[3.0])

    ..  figure:: ../images/sections/monoisection_geometry.png
        :align: center
        :scale: 75 %

        I Section geometry.

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
    return geometry.Geometry(polygon, material)


def tapered_flange_i_section(
    d: float,
    b: float,
    t_f: float,
    t_w: float,
    r_r: float,
    r_f: float,
    alpha: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a Tapered Flange I Section centered at *(b/2, d/2)*, with depth *d*, width *b*,
    mid-flange thickness *t_f*, web thickness *t_w*, root radius *r_r*, flange radius *r_f* and
    flange angle *alpha*, using *n_r* points to construct the radii.

    :param float d: Depth of the Tapered Flange I Section
    :param float b: Width of the Tapered Flange I Section
    :param float t_f: Mid-flange thickness of the Tapered Flange I Section (measured at the point
        equidistant from the face of the web to the edge of the flange)
    :param float t_w: Web thickness of the Tapered Flange I Section
    :param float r_r: Root radius of the Tapered Flange I Section
    :param float r_f: Flange radius of the Tapered Flange I Section
    :param float alpha: Flange angle of the Tapered Flange I Section (degrees)
    :param int n_r: Number of points discretising the radii
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a Tapered Flange I Section with a depth of 588, a width of 191, a
    mid-flange thickness of 27.2, a web thickness of 15.2, a root radius of 17.8, a flange radius
    of 8.9 and a flange angle of 8°, using 16 points to discretise the radii. A mesh is generated
    with a maximum triangular area of 20.0::

        from sectionproperties.pre.library.steel_sections import tapered_flange_i_section

        geometry = tapered_flange_i_section(
            d=588, b=191, t_f=27.2, t_w=15.2, r_r=17.8, r_f=8.9, alpha=8, n_r=16
        )
        geometry.create_mesh(mesh_sizes=[20.0])

    ..  figure:: ../images/sections/taperedisection_geometry.png
        :align: center
        :scale: 75 %

        I Section geometry.

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
            theta = (3.0 / 2 * np.pi - alpha_rad) - (
                i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)
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
            theta = (3.0 * np.pi / 2 + alpha_rad) + i * 1.0 / max(1, n_r - 1) * (
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
            theta = (np.pi * 0.5 - alpha_rad) - (
                i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)
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
            theta = (np.pi * 0.5 + alpha_rad) + (
                i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)
            )

            # calculate the locations of the radius points
            x = r_f + r_f * np.cos(theta)
            y = y_t + r_f * np.sin(theta)

            # append the current points to the points list
            points.append([x, y])

    polygon = Polygon(points)
    return geometry.Geometry(polygon, material)


def channel_section(
    d: float,
    b: float,
    t_f: float,
    t_w: float,
    r: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
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
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a PFC section with a depth of 250, a width of 90, a flange
    thickness of 15, a web thickness of 8 and a root radius of 12, using 8 points to discretise the
    root radius. A mesh is generated with a maximum triangular area of 5.0::

        from sectionproperties.pre.library.steel_sections import channel_section

        geometry = channel_section(d=250, b=90, t_f=15, t_w=8, r=12, n_r=8)
        geometry.create_mesh(mesh_sizes=[5.0])

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
    return geometry.Geometry(polygon, material)


def tapered_flange_channel(
    d: float,
    b: float,
    t_f: float,
    t_w: float,
    r_r: float,
    r_f: float,
    alpha: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
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
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a Tapered Flange Channel section with a depth of 10, a width of
    3.5, a mid-flange thickness of 0.575, a web thickness of 0.475, a root radius of 0.575, a
    flange radius of 0.4 and a flange angle of 8°, using 16 points to discretise the radii. A mesh
    is generated with a maximum triangular area of 0.02::

        from sectionproperties.pre.library.steel_sections import tapered_flange_channel

        geometry = tapered_flange_channel(
            d=10, b=3.5, t_f=0.575, t_w=0.475, r_r=0.575, r_f=0.4, alpha=8, n_r=16
        )
        geometry.create_mesh(mesh_sizes=[0.02])

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
            theta = (3.0 / 2 * np.pi - alpha_rad) - (
                i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)
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
            y = d - t_f - y2 - r_r * np.cos(alpha_rad) + r_r * np.sin(theta)

            # append the current points to the points list
            points.append([x, y])

    # construct the top right flange toe radius
    if r_f == 0:
        points.append([b, d - y_t])
    else:
        for i in range(n_r):
            # determine polar angle
            theta = (3.0 * np.pi / 2 + alpha_rad) + (
                i * 1.0 / max(1, n_r - 1) * (np.pi * 0.5 - alpha_rad)
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
    return geometry.Geometry(polygon, material)


def tee_section(
    d: float,
    b: float,
    t_f: float,
    t_w: float,
    r: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a Tee section with the top left corner at *(0, d)*, with depth *d*, width *b*,
    flange thickness *t_f*, web thickness *t_w* and root radius *r*, using *n_r* points to
    construct the root radius.

    :param float d: Depth of the Tee section
    :param float b: Width of the Tee section
    :param float t_f: Flange thickness of the Tee section
    :param float t_w: Web thickness of the Tee section
    :param float r: Root radius of the Tee section
    :param int n_r: Number of points discretising the root radius
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a Tee section with a depth of 200, a width of 100, a flange
    thickness of 12, a web thickness of 6 and a root radius of 8, using 8 points to discretise the
    root radius. A mesh is generated with a maximum triangular area of 3.0::

        from sectionproperties.pre.library.steel_sections import tee_section

        geometry = tee_section(d=200, b=100, t_f=12, t_w=6, r=8, n_r=8)
        geometry.create_mesh(mesh_sizes=[3.0])

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
    return geometry.Geometry(polygon, material)


def angle_section(
    d: float,
    b: float,
    t: float,
    r_r: float,
    r_t: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs an angle section with the bottom left corner at the origin *(0, 0)*, with depth
    *d*, width *b*, thickness *t*, root radius *r_r* and toe radius *r_t*, using *n_r* points to
    construct the radii.

    :param float d: Depth of the angle section
    :param float b: Width of the angle section
    :param float t: Thickness of the angle section
    :param float r_r: Root radius of the angle section
    :param float r_t: Toe radius of the angle section
    :param int n_r: Number of points discretising the radii
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates an angle section with a depth of 150, a width of 100, a thickness
    of 8, a root radius of 12 and a toe radius of 5, using 16 points to discretise the radii. A
    mesh is generated with a maximum triangular area of 2.0::

        from sectionproperties.pre.library.steel_sections import angle_section

        geometry = angle_section(d=150, b=100, t=8, r_r=12, r_t=5, n_r=16)
        geometry.create_mesh(mesh_sizes=[2.0])

    ..  figure:: ../images/sections/angle_geometry.png
        :align: center
        :scale: 75 %

        Angle section geometry.

    ..  figure:: ../images/sections/angle_mesh.png
        :align: center
        :scale: 75 %
    """
    if r_t > t:
        raise ValueError(
            "The radius of the toe (r_t) cannot be larger than the toe thickness (t)."
        )

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
    return geometry.Geometry(polygon, material)


def cee_section(
    d: float,
    b: float,
    l: float,
    t: float,
    r_out: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
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
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a Cee section with a depth of 125, a width of 50, a lip of 30, a
    thickness of 1.5 and an outer radius of 6, using 8 points to discretise the radius. A mesh is
    generated with a maximum triangular area of 0.25::

        from sectionproperties.pre.library.steel_sections import cee_section

        geometry = cee_section(d=125, b=50, l=30, t=1.5, r_out=6, n_r=8)
        geometry.create_mesh(mesh_sizes=[0.25])

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
        raise Exception("Lip length must be greater than the outer radius")

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
    points += draw_radius([b - t - r_in, d - t - r_in], r_in, 0.5 * np.pi, n_r, False)

    if r_out != l:
        # add next two points
        points.append([b - t, d - l])
        points.append([b, d - l])

    # construct the outer top right radius
    points += draw_radius([b - r_out, d - r_out], r_out, 0, n_r)

    # construct the outer top left radius
    points += draw_radius([r_out, d - r_out], r_out, 0.5 * np.pi, n_r)

    polygon = Polygon(points)
    return geometry.Geometry(polygon, material)


def zed_section(
    d: float,
    b_l: float,
    b_r: float,
    l: float,
    t: float,
    r_out: float,
    n_r: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
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
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a zed section with a depth of 100, a left flange width of 40, a
    right flange width of 50, a lip of 20, a thickness of 1.2 and an outer radius of 5, using 8
    points to discretise the radius. A mesh is generated with a maximum triangular area of 0.15::

        from sectionproperties.pre.library.steel_sections import zed_section

        geometry = zed_section(d=100, b_l=40, b_r=50, l=20, t=1.2, r_out=5, n_r=8)
        geometry.create_mesh(mesh_sizes=[0.15])

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
        raise Exception("Lip length must be greater than the outer radius")

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
    return geometry.Geometry(polygon, material)


def box_girder_section(
    d: float,
    b_t: float,
    b_b: float,
    t_ft: float,
    t_fb: float,
    t_w: float,
    material: pre.Material = pre.DEFAULT_MATERIAL,
):
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

        from sectionproperties.pre.library.steel_sections import box_girder_section

        geometry = box_girder_section(d=1200, b_t=1200, b_b=400, t_ft=100, t_fb=80, t_w=50)
        geometry.create_mesh(mesh_sizes=[200.0])

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

    return geometry.Geometry(outer_polygon - inner_polygon, material)
