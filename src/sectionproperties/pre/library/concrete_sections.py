"""Concrete sections library."""

from __future__ import annotations

import numpy as np

import sectionproperties.pre.geometry as geometry
import sectionproperties.pre.library.primitive_sections as primitive_sections
import sectionproperties.pre.pre as pre


def concrete_rectangular_section(
    d: float,
    b: float,
    dia_top: float,
    area_top: float,
    n_top: int,
    c_top: float,
    dia_bot: float,
    area_bot: float,
    n_bot: int,
    c_bot: float,
    dia_side: float | None = None,
    area_side: float | None = None,
    n_side: int = 0,
    c_side: float = 0.0,
    n_circle: int = 4,
    conc_mat: pre.Material = pre.DEFAULT_MATERIAL,
    steel_mat: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.CompoundGeometry:
    """Constructs a reinforced concrete rectangular section.

    Constructs a reinforced concrete rectangular section of depth ``d`` and width ``b``.

    .. note::
        As the reinforcing bars are described by discretised circles, the area of each
        bar is required to ensure that the correct reinforcing area is provided.

    Args:
        d: Concrete section depth
        b: Concrete section width
        dia_top: Diameter of the top reinforcing bars, used for calculating bar
            placement
        area_top: Area of the top reinforcing bars
        n_top: Number of top, equally spaced reinforcing bars
        c_top: Clear cover to the top reinforcing bars
        dia_bot: Diameter of the bottom reinforcing bars, used for calculating bar
            placement
        area_bot: Area of the bottom reinforcing bars
        n_bot: Number of bottom, equally spaced reinforcing bars
        c_bot: Clear cover to the bottom reinforcing bars
        dia_side: Diameter of the side reinforcing bars, used for calculating bar
            placement. Defaults to ``None``.
        area_side: Area of the side reinforcing bars. Defaults to ``None``.
        n_side: Number of side, equally spaced reinforcing bars. Defaults to ``0``.
        c_side: Clear cover to the side reinforcing bars. Defaults to ``0.0``.
        n_circle: Number of points used to discretise the circular reinforcing bars.
            Defaults to ``4``.
        conc_mat: Material object to assign to the concrete area. Defaults to
            ``pre.DEFAULT_MATERIAL``.
        steel_mat: Material object to assign to the steel area. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        ValueError: Geometry generation failed

    Returns:
        Reinforced concrete rectangular section geometry

    Example:
        The following example creates a 600 mm deep x 300 mm wide concrete beam. The
        beam is reinforced with 3 x 16 mm top bars (20 mm cover), 3 x 20 mm bottom bars
        (30 mm cover), and 2 x 12 mm side bars (30 mm cover). A coarse finite element
        mesh is generated to show the different material regions:

        .. plot::
            :include-source: True
            :caption: Reinforced concrete rectangular section geometry

            from sectionproperties.pre import Material
            from sectionproperties.pre.library import concrete_rectangular_section
            from sectionproperties.analysis import Section

            concrete = Material(
                name="Concrete",
                elastic_modulus=30.1e3,
                poissons_ratio=0.2,
                yield_strength=32,
                density=2.4e-6,
                color="lightgrey",
            )
            steel = Material(
                name="Steel",
                elastic_modulus=200e3,
                poissons_ratio=0.3,
                yield_strength=500,
                density=7.85e-6,
                color="grey",
            )

            geom = concrete_rectangular_section(
                d=600,
                b=300,
                dia_top=16,
                area_top=200,
                n_top=3,
                c_top=20,
                dia_bot=20,
                area_bot=310,
                n_bot=3,
                c_bot=30,
                dia_side=12,
                area_side=110,
                n_side=2,
                c_side=30,
                n_circle=16,
                conc_mat=concrete,
                steel_mat=steel,
            )

            geom.create_mesh(mesh_sizes=[0])  # a size of zero creates a coarse mesh
            Section(geometry=geom).plot_mesh()
    """
    # create rectangular concrete geometry
    geom = primitive_sections.rectangular_section(b=b, d=d, material=conc_mat)

    # calculate reinforcing bar dimensions for top and bottom layers
    if n_top == 1:
        x_i_top = b / 2
        spacing_top = 0.0
    else:
        if c_side:
            x_i_top = c_side + dia_top / 2
            spacing_top = (b - 2 * c_side - dia_top) / (n_top - 1)
        else:
            x_i_top = c_top + dia_top / 2
            spacing_top = (b - 2 * c_top - dia_top) / (n_top - 1)

    if n_bot == 1:
        x_i_bot = b / 2
        spacing_bot = 0.0
    else:
        if c_side:
            x_i_bot = c_side + dia_bot / 2
            spacing_bot = (b - 2 * c_side - dia_bot) / (n_bot - 1)
        else:
            x_i_bot = c_bot + dia_bot / 2
            spacing_bot = (b - 2 * c_bot - dia_bot) / (n_bot - 1)

    # calculate reinforcing bar dimensions for side layers if specified
    if dia_side and n_side != 0:
        x_i_side_left = c_side + dia_side / 2
        x_i_side_right = b - x_i_side_left

        spacing_side = (d - c_top - c_bot - dia_top / 2 - dia_bot / 2) / (n_side + 1)
    else:
        x_i_side_left = 0.0
        x_i_side_right = 0.0
        spacing_side = 0.0

    # add top bars
    for idx in range(n_top):
        bar = primitive_sections.circular_section_by_area(
            area=area_top,
            n=n_circle,
            material=steel_mat,
        ).shift_section(
            x_offset=x_i_top + spacing_top * idx,
            y_offset=d - c_top - dia_top / 2,
        )
        geom = (geom - bar) + bar

    # add bot bars
    for i in range(n_bot):
        bar = primitive_sections.circular_section_by_area(
            area=area_bot,
            n=n_circle,
            material=steel_mat,
        ).shift_section(
            x_offset=x_i_bot + spacing_bot * i,
            y_offset=c_bot + dia_bot / 2,
        )
        geom = (geom - bar) + bar

    # add side bars
    if area_side:
        for i in range(n_side):
            bar_left = primitive_sections.circular_section_by_area(
                area=area_side,
                n=n_circle,
                material=steel_mat,
            ).shift_section(
                x_offset=x_i_side_left,
                y_offset=c_bot + dia_bot / 2 + spacing_side * (i + 1),
            )
            bar_right = bar_left.shift_section(x_offset=x_i_side_right - x_i_side_left)

            geom = (geom - bar_left - bar_right) + bar_left + bar_right

    if isinstance(geom, geometry.CompoundGeometry):
        return geom
    else:
        msg = "Concrete section generation failed."
        raise ValueError(msg)


def concrete_column_section(
    d: float,
    b: float,
    dia_bar: float,
    area_bar: float,
    n_x: int,
    n_y: int,
    cover: float,
    n_circle: int = 4,
    filled: bool = False,
    conc_mat: pre.Material = pre.DEFAULT_MATERIAL,
    steel_mat: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.CompoundGeometry:
    """Constructs a reinforced concrete column section.

    Constructs a reinforced concrete column section of depth ``d`` and width ``b``.

    .. note::
        As the reinforcing bars are described by discretised circles, the area of each
        bar is required to ensure that the correct reinforcing area is provided.

    Args:
        d: Concrete section depth, parallel to the y-axis
        b: Concrete section width, parallel to the x-axis
        dia_bar: Diameter of the reinforcing bars, used for calculating bar placement
        area_bar: Area of the reinforcing bars
        n_x: Number of bars placed across the width of the section, minimum 2
        n_y: Number of bars placed across the depth of the section, minimum 2
        cover: Clear cover to the reinforcing bars
        n_circle: Number of points used to discretise the circular reinforcing bars.
            Defaults to ``4``.
        filled:  When True, will populate the concrete section with an equally spaced
            2D array of reinforcing bars numbering 'n_x' by 'n_y'. When False, only the
            bars around the perimeter of the array will be present. Defaults to
            ``False``.
        conc_mat: Material object to assign to the concrete area. Defaults to
            ``pre.DEFAULT_MATERIAL``.
        steel_mat: Material object to assign to the steel area. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Reinforced concrete column section geometry

    Raises:
        ValueError: If the number of bars in either 'n_x' or 'n_y' is not greater than
            or equal to 2

    Example:
        The following example creates a 600 mm deep x 300 mm wide concrete column. The
        beam is reinforced with 25 mm diameter reinforcing bars with 35 mm cover. There
        are 3 bars in the x-direction and 6 bars in the y-direction. A coarse finite
        element mesh is generated to show the different material regions:

        .. plot::
            :include-source: True
            :caption: Reinforced concrete column section geometry

            from sectionproperties.pre import Material
            from sectionproperties.pre.library import concrete_column_section
            from sectionproperties.analysis import Section

            concrete = Material(
                name="Concrete",
                elastic_modulus=30.1e3,
                poissons_ratio=0.2,
                yield_strength=32,
                density=2.4e-6,
                color="lightgrey",
            )
            steel = Material(
                name="Steel",
                elastic_modulus=200e3,
                poissons_ratio=0.3,
                yield_strength=500,
                density=7.85e-6,
                color="grey",
            )

            geom = concrete_column_section(
                d=600,
                b=300,
                dia_bar=25,
                area_bar=500,
                n_x=3,
                n_y=6,
                cover=35,
                n_circle=4,
                filled=False,
                conc_mat=concrete,
                steel_mat=steel,
            )

            geom.create_mesh(mesh_sizes=[0])  # a size of zero creates a coarse mesh
            Section(geometry=geom).plot_mesh()
    """
    concrete_geometry = primitive_sections.rectangular_section(
        d=d, b=b, material=conc_mat
    )
    bar_extents = concrete_geometry.offset_perimeter(
        -cover - dia_bar / 2
    ).calculate_extents()
    bar_x_min, bar_x_max, bar_y_min, bar_y_max = bar_extents

    b_edge_bars_x: list[float] = np.linspace(bar_x_min, bar_x_max, n_x).tolist()
    d_edge_bars_y: list[float] = np.linspace(bar_y_min, bar_y_max, n_y).tolist()

    if not filled:
        b_edge_bars_y1 = [bar_y_min] * n_x
        b_edge_bars_y2 = [bar_y_max] * n_x

        d_edge_bars_x1 = [bar_x_min] * n_y
        d_edge_bars_x2 = [bar_x_max] * n_y

        b_edge_bars_top = list(zip(b_edge_bars_x, b_edge_bars_y2, strict=False))
        b_edge_bars_bottom = list(zip(b_edge_bars_x, b_edge_bars_y1, strict=False))
        d_edge_bars_right = list(zip(d_edge_bars_x2, d_edge_bars_y, strict=False))
        d_edge_bars_left = list(zip(d_edge_bars_x1, d_edge_bars_y, strict=False))

        all_bar_coords = list(
            set(
                b_edge_bars_top
                + b_edge_bars_bottom
                + d_edge_bars_right
                + d_edge_bars_left
            )
        )
    else:
        xy = np.meshgrid(b_edge_bars_x, d_edge_bars_y)
        all_bar_coords = np.append(xy[0].reshape(-1, 1), xy[1].reshape(-1, 1), axis=1)

    for bar_coord in all_bar_coords:
        concrete_geometry = add_bar(
            geometry=concrete_geometry,
            area=area_bar,
            material=steel_mat,
            x=bar_coord[0],
            y=bar_coord[1],
            n=n_circle,
        )

    if isinstance(concrete_geometry, geometry.CompoundGeometry):
        return concrete_geometry
    else:
        msg = "Concrete section generation failed."
        raise ValueError(msg)


def concrete_tee_section(
    d: float,
    b: float,
    d_f: float,
    b_f: float,
    dia_top: float,
    area_top: float,
    n_top: int,
    c_top: float,
    dia_bot: float,
    area_bot: float,
    n_bot: int,
    c_bot: float,
    dia_side: float | None = None,
    area_side: float | None = None,
    n_side: int = 0,
    c_side: float = 0.0,
    n_circle: int = 4,
    conc_mat: pre.Material = pre.DEFAULT_MATERIAL,
    steel_mat: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.CompoundGeometry:
    """Constructs a reinforced concrete tee section.

    Constructs a concrete tee section of depth ``d``, width ``b``, flange depth ``d_f``
    and flange width ``b_f``.

    Args:
        d: Concrete section depth
        b: Concrete section width
        d_f: Concrete section flange depth
        b_f: Concrete section flange width
        dia_top: Diameter of the top reinforcing bars, used for calculating bar
            placement
        area_top: Area of the top reinforcing bars
        n_top: Number of top, equally spaced reinforcing bars
        c_top: Clear cover to the top reinforcing bars
        dia_bot: Diameter of the bottom reinforcing bars, used for calculating bar
            placement
        area_bot: Area of the bottom reinforcing bars
        n_bot: Number of bottom, equally spaced reinforcing bars
        c_bot: Clear cover to the bottom reinforcing bars
        dia_side: Diameter of the side reinforcing bars, used for calculating bar
            placement. Defaults to ``None``.
        area_side: Area of the side reinforcing bars. Defaults to ``None``.
        n_side: Number of side, equally spaced reinforcing bars. Defaults to ``0``.
        c_side: Clear cover to the side reinforcing bars. Defaults to ``0.0``.
        n_circle: Number of points used to discretise the circular reinforcing bars.
            Defaults to ``4``.
        conc_mat: Material object to assign to the concrete area. Defaults to
            ``pre.DEFAULT_MATERIAL``.
        steel_mat: Material object to assign to the steel area. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        ValueError: Geometry generation failed

    Returns:
        Reinforced concrete tee section geometry

    Example:
        The following example creates a 900 mm deep x 450 mm wide concrete tee beam with
        a 150 mm deep x 1800 mm wide flange. The tee beam is reinforced with 12 x 24 mm
        top bars, 5 x 28 mm bottom bars and 4 x 16 mm side bars, with 42 mm cover all
        around (30 mm cover + 12 mm tie). A coarse finite element mesh is generated to
        show the different material regions:

        .. plot::
            :include-source: True
            :caption: Reinforced concrete tee section geometry

            from sectionproperties.pre import Material
            from sectionproperties.pre.library import concrete_tee_section
            from sectionproperties.analysis import Section

            concrete = Material(
                name="Concrete",
                elastic_modulus=30.1e3,
                poissons_ratio=0.2,
                yield_strength=32,
                density=2.4e-6,
                color="lightgrey",
            )
            steel = Material(
                name="Steel",
                elastic_modulus=200e3,
                poissons_ratio=0.3,
                yield_strength=500,
                density=7.85e-6,
                color="grey",
            )

            geom = concrete_tee_section(
                d=900,
                b=450,
                d_f=150,
                b_f=1800,
                dia_top=24,
                area_top=450,
                n_top=12,
                c_top=42,
                dia_bot=28,
                area_bot=620,
                n_bot=5,
                c_bot=42,
                dia_side=16,
                area_side=200,
                n_side=4,
                c_side=42,
                n_circle=16,
                conc_mat=concrete,
                steel_mat=steel,
            )

            geom.create_mesh(mesh_sizes=[0])  # a size of zero creates a coarse mesh
            Section(geometry=geom).plot_mesh()
    """
    # generate concrete geometry
    beam = primitive_sections.rectangular_section(b=b, d=d - d_f, material=conc_mat)
    flange = primitive_sections.rectangular_section(b=b_f, d=d_f, material=conc_mat)
    geom = beam + flange.align_to(other=beam, on="top").shift_section(
        x_offset=-(b_f / 2 - b / 2)
    )
    geom = geom.shift_section(x_offset=-b / 2)

    # calculate reinforcing bar dimensions for top and bottom layers
    if n_top == 1:
        x_i_top = 0.0
        spacing_top = 0.0
    else:
        if c_side:
            x_i_top = -b_f / 2 + c_side + dia_top / 2
            spacing_top = (b_f - 2 * c_side - dia_top) / (n_top - 1)
        else:
            x_i_top = -b_f / 2 + c_top + dia_top / 2
            spacing_top = (b_f - 2 * c_top - dia_top) / (n_top - 1)

    if n_bot == 1:
        x_i_bot = 0.0
        spacing_bot = 0.0
    else:
        if c_side:
            x_i_bot = -b / 2 + c_side + dia_bot / 2
            spacing_bot = (b - 2 * c_side - dia_bot) / (n_bot - 1)
        else:
            x_i_bot = -b / 2 + c_bot + dia_bot / 2
            spacing_bot = (b - 2 * c_bot - dia_bot) / (n_bot - 1)

    # calculate reinforcing bar dimensions for side layers if specified
    if dia_side and n_side != 0:
        x_i_side_left = -b / 2 + c_side + dia_side / 2
        x_i_side_right = -x_i_side_left

        spacing_side = (d - c_top - c_bot - dia_top / 2 - dia_bot / 2) / (n_side + 1)
    else:
        x_i_side_left = 0.0
        x_i_side_right = 0.0
        spacing_side = 0.0

    # add top bars
    for idx in range(n_top):
        bar = primitive_sections.circular_section_by_area(
            area=area_top,
            n=n_circle,
            material=steel_mat,
        ).shift_section(
            x_offset=x_i_top + spacing_top * idx,
            y_offset=d - c_top - dia_top / 2,
        )
        geom = (geom - bar) + bar

    # add bot bars
    for i in range(n_bot):
        bar = primitive_sections.circular_section_by_area(
            area=area_bot,
            n=n_circle,
            material=steel_mat,
        ).shift_section(
            x_offset=x_i_bot + spacing_bot * i,
            y_offset=c_bot + dia_bot / 2,
        )
        geom = (geom - bar) + bar

    # add side bars
    if area_side:
        for i in range(n_side):
            bar_left = primitive_sections.circular_section_by_area(
                area=area_side,
                n=n_circle,
                material=steel_mat,
            ).shift_section(
                x_offset=x_i_side_left,
                y_offset=c_bot + dia_bot / 2 + spacing_side * (i + 1),
            )
            bar_right = bar_left.shift_section(x_offset=x_i_side_right - x_i_side_left)

            geom = (geom - bar_left - bar_right) + bar_left + bar_right

    if isinstance(geom, geometry.CompoundGeometry):  # pyright: ignore [reportUnnecessaryIsInstance]
        return geom
    else:
        msg = "Concrete section generation failed."
        raise ValueError(msg)


def concrete_circular_section(
    d: float,
    area_conc: float,
    n_conc: int,
    dia_bar: float,
    area_bar: float,
    n_bar: int,
    cover: float,
    n_circle: int = 4,
    conc_mat: pre.Material = pre.DEFAULT_MATERIAL,
    steel_mat: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.CompoundGeometry:
    """Constructs a reinforced concrete circular section.

    Constructs a reinforced concrete circular section of diameter ``d``.

    .. note::
        As the concrete area and reinforcing bars are described by discretised circles,
        the area of the circle and each bar is required to ensure that the correct
        concrete and reinforcing area is provided.

    Args:
        d: Concrete diameter
        area_conc: Area of the circular concrete section
        n_conc: Number of points discretising the circular concrete section
        dia_bar: Diameter of the steel reinforcing bars
        area_bar: Area of the reinforcing bars
        n_bar: Number of steel reinforcing bars
        cover: Clear cover to the reinforcing bars
        n_circle: Number of points discretising the steel reinforcing bars. Defaults to
            ``4``.
        conc_mat: Material object to assign to the concrete area. Defaults to
            ``pre.DEFAULT_MATERIAL``.
        steel_mat: Material object to assign to the steel area. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Returns:
        Reinforced concrete circular section geometry

    Raises:
        ValueError: If the number of bars is not greater than or equal to 2

    Example:
        The following example creates a 450 mm diameter concrete column with 6 x 20 mm
        steel reinforcing bars with 45 mm cover. A coarse finite element mesh is
        generated to show the different material regions:

        .. plot::
            :include-source: True
            :caption: Reinforced concrete circular section geometry

            import numpy as np
            from sectionproperties.pre import Material
            from sectionproperties.pre.library import concrete_circular_section
            from sectionproperties.analysis import Section

            concrete = Material(
                name="Concrete",
                elastic_modulus=30.1e3,
                poissons_ratio=0.2,
                yield_strength=32,
                density=2.4e-6,
                color="lightgrey",
            )
            steel = Material(
                name="Steel",
                elastic_modulus=200e3,
                poissons_ratio=0.3,
                yield_strength=500,
                density=7.85e-6,
                color="grey",
            )

            geom = concrete_circular_section(
                d=450,
                area_conc=np.pi * 450 * 450 / 4,
                n_conc=64,
                dia_bar=20,
                area_bar=310,
                n_bar=6,
                cover=45,
                n_circle=24,
                conc_mat=concrete,
                steel_mat=steel,
            )

            geom.create_mesh(mesh_sizes=[0])  # a size of zero creates a coarse mesh
            Section(geometry=geom).plot_mesh()
    """
    if n_bar < 2:
        msg = "Please provide 2 or more steel reinforcing bars."
        raise ValueError(msg)

    # create circular geometry
    geom = primitive_sections.circular_section_by_area(
        area=area_conc, n=n_conc, material=conc_mat
    )

    # calculate bar geometry
    r = d / 2 - cover - dia_bar / 2
    d_theta = 2 * np.pi / n_bar

    for i in range(n_bar):
        bar = primitive_sections.circular_section_by_area(
            area=area_bar,
            n=n_circle,
            material=steel_mat,
        ).shift_section(
            x_offset=r * np.cos(i * d_theta),
            y_offset=r * np.sin(i * d_theta),
        )
        geom = (geom - bar) + bar

    if isinstance(geom, geometry.CompoundGeometry):
        return geom
    else:
        msg = "Concrete section generation failed."
        raise ValueError(msg)


def add_bar(
    geometry: geometry.Geometry | geometry.CompoundGeometry,
    area: float,
    material: pre.Material,
    x: float,
    y: float,
    n: int,
) -> geometry.CompoundGeometry:
    """Adds a reinforcing bar to a ``sectionproperties`` geometry.

    First removes the geometry through a subtraction operation, then adds the geometry
    on top of the newly created hole. This method avoids the doubling up of geometry.

    Args:
        geometry: Reinforced concrete geometry to which the new bar will be added
        area: Bar cross-sectional area
        material: Material object for the bar
        x: x-position of the bar
        y: y-position of the bar
        n: Number of points to discretise the bar circle

    Returns:
        Geometry object with added bar
    """
    bar = primitive_sections.circular_section_by_area(
        area=area,
        n=n,
        material=material,
    ).shift_section(x_offset=x, y_offset=y)

    return (geometry - bar) + bar
