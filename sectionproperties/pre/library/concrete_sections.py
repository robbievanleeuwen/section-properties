from typing import Union, Optional
import numpy as np
import sectionproperties.pre.pre as pre
import sectionproperties.pre.geometry as geometry
import sectionproperties.pre.library.primitive_sections as primitive_sections


def concrete_rectangular_section(
    b: float,
    d: float,
    dia_top: float,
    n_top: int,
    dia_bot: float,
    n_bot: int,
    n_circle: int,
    cover: float,
    dia_side: Optional[float] = None,
    n_side: int = 0,
    area_top: Optional[float] = None,
    area_bot: Optional[float] = None,
    area_side: Optional[float] = None,
    conc_mat: pre.Material = pre.DEFAULT_MATERIAL,
    steel_mat: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.CompoundGeometry:
    """Constructs a concrete rectangular section of width *b* and depth *d*, with
    *n_top* top steel bars of diameter *dia_top*, *n_bot* bottom steel bars of diameter
    *dia_bot*, *n_side* left & right side steel bars of diameter *dia_side* discretised
    with *n_circle* points with equal side and top/bottom *cover* to the steel.

    :param float b: Concrete section width
    :param float d: Concrete section depth
    :param float dia_top: Diameter of the top steel reinforcing bars
    :param int n_top: Number of top steel reinforcing bars
    :param float dia_bot: Diameter of the bottom steel reinforcing bars
    :param int n_bot: Number of bottom steel reinforcing bars
    :param int n_circle: Number of points discretising the steel reinforcing bars
    :param float cover: Side and bottom cover to the steel reinforcing bars
    :param float dia_side: If provided, diameter of the side steel reinforcing bars
    :param int n_side: If provided, number of side bars either side of the section
    :param float area_top: If provided, constructs top reinforcing bars based on their
        area rather than diameter (prevents the underestimation of steel area due to
        circle discretisation)
    :param float area_bot: If provided, constructs bottom reinforcing bars based on
        their area rather than diameter (prevents the underestimation of steel area due
        to circle discretisation)
    :param float area_side: If provided, constructs side reinforcing bars based on
        their area rather than diameter (prevents the underestimation of steel area due
        to circle discretisation)
    :param conc_mat: Material to associate with the concrete
    :param steel_mat: Material to associate with the steel

    :raises ValueError: If the number of bars is not greater than or equal to 2 in an
        active layer

    The following example creates a 600D x 300W concrete beam with 3N20 bottom steel
    reinforcing bars and 30 mm cover::

        from sectionproperties.pre.library.concrete_sections import concrete_rectangular_section
        from sectionproperties.pre.pre import Material

        concrete = Material(
            name='Concrete', elastic_modulus=30.1e3, poissons_ratio=0.2, yield_strength=32,
            density=2.4e-6, color='lightgrey'
        )
        steel = Material(
            name='Steel', elastic_modulus=200e3, poissons_ratio=0.3, yield_strength=500,
            density=7.85e-6, color='grey'
        )

        geometry = concrete_rectangular_section(
            b=300, d=600, dia_top=20, n_top=0, dia_bot=20, n_bot=3, n_circle=24, cover=30,
            conc_mat=concrete, steel_mat=steel
        )
        geometry.create_mesh(mesh_sizes=[500])

    ..  figure:: ../images/sections/concrete_rectangular_section_geometry.png
        :align: center
        :scale: 50 %

        Concrete rectangular section geometry.

    ..  figure:: ../images/sections/concrete_rectangular_section_mesh.png
        :align: center
        :scale: 50 %

        Mesh generated from the above geometry.
    """

    if n_top == 1 or n_bot == 1:
        raise ValueError("If adding a reinforcing layer, provide 2 or more bars.")

    # create rectangular concrete geometry
    geom = primitive_sections.rectangular_section(b=b, d=d, material=conc_mat)

    # calculate reinforcing bar dimensions for top and bottom layers
    x_i_top = cover + dia_top / 2
    x_i_bot = cover + dia_bot / 2
    spacing_top = (b - 2 * cover - dia_top) / (n_top - 1)
    spacing_bot = (b - 2 * cover - dia_bot) / (n_bot - 1)

    # calculate reinforcing bar dimensions for side layers if specified
    if n_side != 0:
        x_i_side_left = cover + dia_side / 2
        x_i_side_right = b - x_i_side_left
        spacing_side = (d - 2 * cover - dia_top / 2 - dia_bot / 2) / (n_side + 1)

    if area_top is None:
        area_top = np.pi * dia_top**2 / 4
    if area_bot is None:
        area_bot = np.pi * dia_bot**2 / 4
    if area_side is None and dia_side is not None:
        area_side = np.pi * dia_side**2 / 4

    # add top bars
    for i in range(n_top):
        bar = primitive_sections.circular_section_by_area(
            area=area_top, n=n_circle, material=steel_mat
        )
        bar = bar.shift_section(
            x_offset=x_i_top + spacing_top * i, y_offset=d - cover - dia_top / 2
        )
        geom = (geom - bar) + bar

    # add bot bars
    for i in range(n_bot):
        bar = primitive_sections.circular_section_by_area(
            area=area_bot, n=n_circle, material=steel_mat
        )
        bar = bar.shift_section(
            x_offset=x_i_bot + spacing_bot * i, y_offset=cover + dia_bot / 2
        )
        geom = (geom - bar) + bar

    # add side bars if specified
    if n_side != 0:
        for i in range(n_side):
            bar_left = primitive_sections.circular_section_by_area(
                area=area_side, n=n_circle, material=steel_mat
            )
            bar_right = bar_left

            bar_left = bar_left.shift_section(
                x_offset=x_i_side_left,
                y_offset=cover + dia_bot / 2 + spacing_side * (i + 1),
            )
            bar_right = bar_right.shift_section(
                x_offset=x_i_side_right,
                y_offset=cover + dia_bot / 2 + spacing_side * (i + 1),
            )

            geom = (geom - bar_left - bar_right) + bar_left + bar_right

    return geom


def concrete_column_section(
    b: float,
    d: float,
    cover: float,
    n_bars_b: int,
    n_bars_d: int,
    dia_bar: float,
    bar_area: Optional[float] = None,
    conc_mat: pre.Material = pre.DEFAULT_MATERIAL,
    steel_mat: pre.Material = pre.DEFAULT_MATERIAL,
    filled: bool = False,
    n_circle: int = 4,
) -> geometry.CompoundGeometry:
    """Constructs a concrete rectangular section of width *b* and depth *d*, with
    steel bar reinforcing organized as an *n_bars_b* by *n_bars_d* array, discretised
    with *n_circle* points with equal sides and top/bottom *cover* to the steel which
    is taken as the clear cover (edge of bar to edge of concrete).

    :param float b: Concrete section width, parallel to the x-axis
    :param float d: Concrete section depth, parallel to the y-axis
    :param float cover: Clear cover, calculated as distance from edge of reinforcing bar to edge of section.
    :param int n_bars_b: Number of bars placed across the width of the section, minimum 2.
    :param int n_bars_d: Number of bars placed across the depth of the section, minimum 2.
    :param float dia_bar: Diameter of reinforcing bars. Used for calculating bar placement and,
        optionally, for calculating the bar area for section capacity calculations.
    :param float bar_area: Area of reinforcing bars. Used for section capacity calculations.
        If not provided, then dia_bar will be used to calculate the bar area.
    :param sectionproperties.pre.pre.Material conc_mat: Material to
        associate with the concrete
    :param sectionproperties.pre.pre.Material steel_mat: Material to
        associate with the reinforcing steel
    :param bool filled: When True, will populate the concrete section with an equally
        spaced 2D array of reinforcing bars numbering 'n_bars_b' by 'n_bars_d'.
        When False, only the bars around the perimeter of the array will be present.
    :param int n_circle: The number of points used to discretize the circle of the reinforcing
        bars. The bars themselves will have an exact area of 'bar_area' regardless of the
        number of points used in the circle. Useful for making the reinforcing bars look
        more circular when plotting the concrete section.

    :raises ValueError: If the number of bars in either 'n_bars_b' or 'n_bars_d' is not greater
        than or equal to 2.

    The following example creates a 600D x 300W concrete column with 25 mm diameter
    reinforcing bars each with 500 mm**2 area and 35 mm cover in a 3x6 array without
    the interior bars being filled::

        from sectionproperties.pre.library.concrete_sections import concrete_column_section
        from sectionproperties.pre.pre import Material

        concrete = Material(
            name='Concrete', elastic_modulus=30.1e3, poissons_ratio=0.2, yield_strength=32,
            density=2.4e-6, color='lightgrey'
        )
        steel = Material(
            name='Steel', elastic_modulus=200e3, poissons_ratio=0.3, yield_strength=500,
            density=7.85e-6, color='grey'
        )

        geometry = concrete_column_section(
            b=300, d=600, dia_bar=25, bar_area=500, cover=35, n_bars_b=3, n_bars_d=6,
            conc_mat=concrete, steel_mat=steel, filled=False, n_circle=4
        )
        geometry.create_mesh(mesh_sizes=[500])
    """

    concrete_geometry = primitive_sections.rectangular_section(b, d, material=conc_mat)
    bar_extents = concrete_geometry.offset_perimeter(
        -cover - dia_bar / 2
    ).calculate_extents()
    bar_x_min, bar_x_max, bar_y_min, bar_y_max = bar_extents

    b_edge_bars_x = np.linspace(bar_x_min, bar_x_max, n_bars_b)
    d_edge_bars_y = np.linspace(bar_y_min, bar_y_max, n_bars_d)

    if not filled:
        b_edge_bars_y1 = [bar_y_min] * n_bars_b
        b_edge_bars_y2 = [bar_y_max] * n_bars_b

        d_edge_bars_x1 = [bar_x_min] * n_bars_d
        d_edge_bars_x2 = [bar_x_max] * n_bars_d

        b_edge_bars_top = list(zip(b_edge_bars_x, b_edge_bars_y2))
        b_edge_bars_bottom = list(zip(b_edge_bars_x, b_edge_bars_y1))
        d_edge_bars_right = list(zip(d_edge_bars_x2, d_edge_bars_y))
        d_edge_bars_left = list(zip(d_edge_bars_x1, d_edge_bars_y))

        all_bar_coords = list(
            set(
                b_edge_bars_top
                + b_edge_bars_bottom
                + d_edge_bars_right
                + d_edge_bars_left
            )
        )
    if filled:
        xy = np.meshgrid(b_edge_bars_x, d_edge_bars_y)
        all_bar_coords = np.append(xy[0].reshape(-1, 1), xy[1].reshape(-1, 1), axis=1)

    if bar_area is None:
        bar_area = np.pi * dia_bar**2 / 4
    for bar_coord in all_bar_coords:
        concrete_geometry = add_bar(
            concrete_geometry,
            area=bar_area,
            material=steel_mat,
            x=bar_coord[0],
            y=bar_coord[1],
            n=n_circle,
        )
    return concrete_geometry


def concrete_tee_section(
    b: float,
    d: float,
    b_f: float,
    d_f: float,
    dia_top: float,
    n_top: int,
    dia_bot: float,
    n_bot: int,
    n_circle: int,
    cover: float,
    area_top: Optional[float] = None,
    area_bot: Optional[float] = None,
    conc_mat: pre.Material = pre.DEFAULT_MATERIAL,
    steel_mat: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.CompoundGeometry:
    """Constructs a concrete tee section of width *b*, depth *d*, flange width *b_f*
    and flange depth *d_f*, with *n_top* top steel bars of diameter *dia_top*, *n_bot*
    bottom steel bars of diameter *dia_bot*, discretised with *n_circle* points with
    equal side and top/bottom *cover* to the steel.

    :param float b: Concrete section width
    :param float d: Concrete section depth
    :param float b_f: Concrete section flange width
    :param float d_f: Concrete section flange depth
    :param float dia_top: Diameter of the top steel reinforcing bars
    :param int n_top: Number of top steel reinforcing bars
    :param float dia_bot: Diameter of the bottom steel reinforcing bars
    :param int n_bot: Number of bottom steel reinforcing bars
    :param int n_circle: Number of points discretising the steel reinforcing bars
    :param float cover: Side and bottom cover to the steel reinforcing bars
    :param float area_top: If provided, constructs top reinforcing bars based on their
        area rather than diameter (prevents the underestimation of steel area due to
        circle discretisation)
    :param float area_bot: If provided, constructs bottom reinforcing bars based on
        their area rather than diameter (prevents the underestimation of steel area due
        to circle discretisation)
    :param conc_mat: Material to associatewith the concrete
    :param steel_mat: Material toassociate with the steel

    :raises ValueErorr: If the number of bars is not greater than or equal to 2 in an
        active layer

    The following example creates a 900D x 450W concrete beam with a 1200W x 250D
    flange, with 5N24 steel reinforcing bars and 30 mm cover::

        from sectionproperties.pre.library.concrete_sections import concrete_tee_section
        from sectionproperties.pre.pre import Material

        concrete = Material(
            name='Concrete', elastic_modulus=30.1e3, poissons_ratio=0.2, yield_strength=32,
            density=2.4e-6, color='lightgrey'
        )
        steel = Material(
            name='Steel', elastic_modulus=200e3, poissons_ratio=0.3, yield_strength=500,
            density=7.85e-6, color='grey'
        )

        geometry = concrete_tee_section(
            b=450, d=900, b_f=1200, d_f=250, dia_top=24, n_top=0, dia_bot=24, n_bot=5,
            n_circle=24, cover=30, conc_mat=concrete, steel_mat=steel
        )
        geometry.create_mesh(mesh_sizes=[500])

    ..  figure:: ../images/sections/concrete_tee_section_geometry.png
        :align: center
        :scale: 50 %

        Concrete tee section geometry.

    ..  figure:: ../images/sections/concrete_tee_section_mesh.png
        :align: center
        :scale: 50 %

        Mesh generated from the above geometry.
    """

    if n_top == 1 or n_bot == 1:
        raise ValueError("If adding a reinforcing layer, provide 2 or more bars.")

    # generate tee-section
    geom = primitive_sections.rectangular_section(b=b, d=d - d_f, material=conc_mat)
    flange = primitive_sections.rectangular_section(b=b_f, d=d_f, material=conc_mat)
    geom += flange.align_center(align_to=geom).align_to(other=geom, on="top")

    # calculate reinforcing bar dimensions
    x_i_top = cover + dia_top / 2 + 0.5 * (b - b_f)
    x_i_bot = cover + dia_bot / 2
    spacing_top = (b_f - 2 * cover - dia_top) / (n_top - 1)
    spacing_bot = (b - 2 * cover - dia_bot) / (n_bot - 1)

    if area_top is None:
        area_top = np.pi * dia_top**2 / 4
    if area_bot is None:
        area_bot = np.pi * dia_bot**2 / 4

    # add top bars
    for i in range(n_top):
        bar = primitive_sections.circular_section_by_area(
            area=area_top, n=n_circle, material=steel_mat
        )
        bar = bar.shift_section(
            x_offset=x_i_top + spacing_top * i, y_offset=d - cover - dia_top / 2
        )
        geom = (geom - bar) + bar

    # add bot bars
    for i in range(n_bot):
        bar = primitive_sections.circular_section_by_area(
            area=area_bot, n=n_circle, material=steel_mat
        )
        bar = bar.shift_section(
            x_offset=x_i_bot + spacing_bot * i, y_offset=cover + dia_bot / 2
        )
        geom = (geom - bar) + bar

    return geom


def concrete_circular_section(
    d: float,
    n: int,
    dia: float,
    n_bar: int,
    n_circle: int,
    cover: float,
    area_conc: Optional[float] = None,
    area_bar: Optional[float] = None,
    conc_mat: pre.Material = pre.DEFAULT_MATERIAL,
    steel_mat: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.CompoundGeometry:
    """Constructs a concrete circular section of diameter *d* discretised with *n*
    points, with *n_bar* steel bars of diameter *dia*, discretised with *n_circle*
    points with equal side and bottom *cover* to the steel.

    :param float d: Concrete diameter
    :param float n: Number of points discretising the concrete section
    :param float dia: Diameter of the steel reinforcing bars
    :param int n_bar: Number of steel reinforcing bars
    :param int n_circle: Number of points discretising the steel reinforcing bars
    :param float cover: Side and bottom cover to the steel reinforcing bars
    :param float area_conc: If provided, constructs the concrete based on its area
        rather than diameter (prevents the underestimation of concrete area due to
        circle discretisation)
    :param float area_bar: If provided, constructs reinforcing bars based on their area
        rather than diameter (prevents the underestimation of steel area due to
    :param conc_mat: Material to associate with the concrete
    :param steel_mat: Material to associate with the steel

    :raises ValueErorr: If the number of bars is not greater than or equal to 2

    The following example creates a 450DIA concrete column with with 6N20 steel
    reinforcing bars and 45 mm cover::

        from sectionproperties.pre.library.concrete_sections import concrete_circular_section
        from sectionproperties.pre.pre import Material

        concrete = Material(
            name='Concrete', elastic_modulus=30.1e3, poissons_ratio=0.2, yield_strength=32,
            density=2.4e-6, color='lightgrey'
        )
        steel = Material(
            name='Steel', elastic_modulus=200e3, poissons_ratio=0.3, yield_strength=500,
            density=7.85e-6, color='grey'
        )

        geometry = concrete_circular_section(
            d=450, n=64, dia=20, n_bar=6, n_circle=24, cover=45, conc_mat=concrete, steel_mat=steel
        )
        geometry.create_mesh(mesh_sizes=[500])

    ..  figure:: ../images/sections/concrete_circular_section_geometry.png
        :align: center
        :scale: 50 %

        Concrete circular section geometry.

    ..  figure:: ../images/sections/concrete_circular_section_mesh.png
        :align: center
        :scale: 50 %

        Mesh generated from the above geometry.
    """

    if n_bar < 2:
        raise ValueError("Please provide 2 or more steel reinforcing bars.")

    # create circular geometry
    if area_conc:
        geom = primitive_sections.circular_section_by_area(
            area=area_conc, n=n, material=conc_mat
        )
    else:
        geom = primitive_sections.circular_section(d=d, n=n, material=conc_mat)

    # calculate bar geometry
    r = d / 2 - cover - dia / 2
    d_theta = 2 * np.pi / n_bar

    if area_bar is None:
        area_bar = np.pi * dia**2 / 4
    for i in range(n_bar):
        bar = primitive_sections.circular_section_by_area(
            area=area_bar, n=n_circle, material=steel_mat
        )
        bar = bar.shift_section(
            x_offset=r * np.cos(i * d_theta), y_offset=r * np.sin(i * d_theta)
        )
        geom = (geom - bar) + bar

    return geom


def add_bar(
    geometry: Union[geometry.Geometry, geometry.CompoundGeometry],
    area: float,
    material: pre.DEFAULT_MATERIAL,
    x: float,
    y: float,
    n: int = 4,
) -> geometry.CompoundGeometry:
    """Adds a reinforcing bar to a *sectionproperties* geometry.

    Bars are discretised by four points by default.

    :param geometry: Reinforced concrete geometry to which the new bar will be added
    :param area: Bar cross-sectional area
    :param material: Material object for the bar
    :param x: x-position of the bar
    :param y: y-position of the bar
    :param n: Number of points to discretise the bar circle

    :return: Reinforced concrete geometry with added bar
    """

    bar = primitive_sections.circular_section_by_area(
        area=area, n=n, material=material  # type: ignore
    ).shift_section(x_offset=x, y_offset=y)

    return (geometry - bar) + bar
