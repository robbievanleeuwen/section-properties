"""Concrete sections library."""

from __future__ import annotations

from math import ceil

import numpy as np

import sectionproperties.pre.library.primitive_sections as primitive_sections
import sectionproperties.pre.pre as pre
from sectionproperties.pre.geometry import CompoundGeometry, Geometry


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
) -> CompoundGeometry:
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

    if isinstance(geom, CompoundGeometry):
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
) -> CompoundGeometry:
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

    if isinstance(concrete_geometry, CompoundGeometry):
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
) -> CompoundGeometry:
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

    if isinstance(geom, CompoundGeometry):  # pyright: ignore [reportUnnecessaryIsInstance]
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
) -> CompoundGeometry:
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

    if isinstance(geom, CompoundGeometry):
        return geom
    else:
        msg = "Concrete section generation failed."
        raise ValueError(msg)


def rectangular_wall(
    d: float,
    t: float,
    dia_bar: float,
    area_bar: float,
    spacing: float,
    cover: float,
    double: bool = True,
    n_circle: int = 4,
    conc_mat: pre.Material = pre.DEFAULT_MATERIAL,
    steel_mat: pre.Material = pre.DEFAULT_MATERIAL,
) -> CompoundGeometry:
    """Constructs a reinforced concrete rectangular wall section.

    Constructs a reinforced concrete rectangular wall section of depth ``d`` and
    thickness ``t``. The wall reinforcement has a maximum spacing of ``spacing`` and is
    doubly reinforced if ``double`` is set to ``True``, or singly reinforced if it is
    set to ``False``.

    .. note::
        As the reinforcing bars are described by discretised circles, the area of each
        bar is required to ensure that the correct reinforcing area is provided.

    Args:
        d: Concrete wall depth
        t: Concrete wall thickness
        dia_bar: Diameter of the reinforcing bars, used for calculating bar placement
        area_bar: Area of the reinforcing bars
        spacing: Maximum spacing of the reinforcement bars, the calculated spacing is
            equal to ``ceil(extent / spacing) + 1``
        cover: Clear cover to the reinforcing bars
        double: If set to ``True``, provides two layers of reinforcement to the wall. If
            set to ``False``, provides a single central layer of reinforcement to the
            wall. Defaults to ``True``.
        n_circle: Number of points used to discretise the circular reinforcing bars.
            Defaults to ``4``.
        conc_mat: Material object to assign to the concrete area. Defaults to
            ``pre.DEFAULT_MATERIAL``.
        steel_mat: Material object to assign to the steel area. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        ValueError: Geometry generation failed

    Returns:
        Reinforced concrete rectangular wall section geometry

    Example:
        The following example creates a 1000 mm long x 180 mm thick rectangular concrete
        wall, reinforced with a single layer of N12-200 with 30 mm cover. A coarse
        finite element mesh is generated to show the different material regions:

        .. plot::
            :include-source: True
            :caption: Reinforced concrete rectangular wall section geometry

            from sectionproperties.pre import Material
            from sectionproperties.pre.library import rectangular_wall
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

            geom = rectangular_wall(
                d=1000,
                t=180,
                dia_bar=12,
                area_bar=110,
                spacing=200,
                cover=30,
                double=False,
                n_circle=12,
                conc_mat=concrete,
                steel_mat=steel,
            )

            geom.create_mesh(mesh_sizes=[0])  # a size of zero creates a coarse mesh
            Section(geometry=geom).plot_mesh()
    """
    # create rectangular concrete geometry
    geom = primitive_sections.rectangular_section(d=d, b=t, material=conc_mat)

    # positions of bars to add
    xs: list[float] = []
    ys: list[float] = []

    # calculate number of bars along length of wall
    y_length = d - 2.0 * cover - dia_bar
    n_y = ceil(y_length / spacing) + 1

    # calculate y-position of bars
    start = cover + 0.5 * dia_bar
    stop = d - cover - 0.5 * dia_bar
    y_i = np.linspace(start=start, stop=stop, num=n_y)

    # calculate x-position of bars
    x_i = [cover + 0.5 * dia_bar, t - cover - 0.5 * dia_bar] if double else [0.5 * t]

    # add bars
    for y in y_i:
        for x in x_i:
            xs.append(x)
            ys.append(y)

    return add_bars(
        geometry=geom, area=area_bar, material=steel_mat, x=xs, y=ys, n=n_circle
    )


def cee_wall(
    d: float,
    b: float,
    t_f: float,
    t_w: float,
    dia_bar: float,
    area_bar: float,
    spacing: float,
    cover: float,
    double: bool = True,
    n_circle: int = 4,
    conc_mat: pre.Material = pre.DEFAULT_MATERIAL,
    steel_mat: pre.Material = pre.DEFAULT_MATERIAL,
) -> CompoundGeometry:
    """Constructs a reinforced concrete cee-shaped wall section.

    Constructs a reinforced concrete cee-shaped wall section of depth ``d``, width
    ``b``, flange thickness ``t_f`` and web thickness ``t_w``. The wall reinforcement
    has a maximum spacing of ``spacing`` and is doubly reinforced if ``double`` is set
    to ``True``, or singly reinforced if it is set to ``False``.

    .. note::
        As the reinforcing bars are described by discretised circles, the area of each
        bar is required to ensure that the correct reinforcing area is provided.

    Args:
        d: Concrete wall depth
        b: Concrete wall width
        t_f: Concrete wall flange thickness
        t_w: Concrete wall web thickness
        dia_bar: Diameter of the reinforcing bars, used for calculating bar placement
        area_bar: Area of the reinforcing bars
        spacing: Maximum spacing of the reinforcement bars, the calculated spacing is
            equal to ``ceil(extent / spacing) + 1``
        cover: Clear cover to the reinforcing bars
        double: If set to ``True``, provides two layers of reinforcement to the wall. If
            set to ``False``, provides a single central layer of reinforcement to the
            wall. Defaults to ``True``.
        n_circle: Number of points used to discretise the circular reinforcing bars.
            Defaults to ``4``.
        conc_mat: Material object to assign to the concrete area. Defaults to
            ``pre.DEFAULT_MATERIAL``.
        steel_mat: Material object to assign to the steel area. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        ValueError: Geometry generation failed

    Returns:
        Reinforced concrete cee-shaped wall section geometry

    Example:
        The following example creates a 2000 mm deep x 1500 mm wide cee-shaped concrete
        wall, with a 200 mm thick flange and 150 mm thick web. The wall is reinforced
        with a double layer of N16-150 with 30 mm cover. A coarse finite element mesh is
        generated to show the different material regions:

        .. plot::
            :include-source: True
            :caption: Reinforced concrete cee-shaped wall section geometry

            from sectionproperties.pre import Material
            from sectionproperties.pre.library import cee_wall
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

            geom = cee_wall(
                d=2000,
                b=1500,
                t_f=200,
                t_w=150,
                dia_bar=16,
                area_bar=200,
                spacing=150,
                cover=30,
                double=True,
                n_circle=12,
                conc_mat=concrete,
                steel_mat=steel,
            )

            geom.create_mesh(mesh_sizes=[0])  # a size of zero creates a coarse mesh
            Section(geometry=geom).plot_mesh()
    """
    # create cee concrete geometry
    geom_outer = primitive_sections.rectangular_section(d=d, b=b, material=conc_mat)
    geom_inner = (
        primitive_sections.rectangular_section(
            d=d - 2 * t_f, b=b - t_w, material=conc_mat
        )
        .align_center(align_to=geom_outer)
        .align_to(other=geom_outer, on="right", inner=True)
    )
    geom = geom_outer - geom_inner

    # positions of bars to add
    xs: list[float] = []
    ys: list[float] = []

    # calculate reinforcement positions
    # singly reinforced
    if not double:
        # calculate number of bars along length of wall
        x_length = b - cover - 0.5 * dia_bar - 0.5 * t_w
        n_x = ceil(x_length / spacing) + 1
        y_length = d - t_f
        n_y = ceil(y_length / spacing) + 1

        # calculate position of bars
        x_i = np.linspace(start=0.5 * t_w, stop=b - cover - 0.5 * dia_bar, num=n_x)
        y_i = np.linspace(start=0.5 * t_f, stop=d - 0.5 * t_f, num=n_y)

        # add bottom bars
        for x in x_i:
            xs.append(x)
            ys.append(0.5 * t_f)

        # add top bars
        for x in x_i:
            xs.append(x)
            ys.append(d - 0.5 * t_f)

        # add web bars
        for y in y_i[1:-1]:
            xs.append(0.5 * t_w)
            ys.append(y)
    # doubly reinforced
    else:
        # bottom/top outer bars
        x_length = b - 2.0 * cover - dia_bar
        n_x = ceil(x_length / spacing) + 1
        x_i = np.linspace(
            start=cover + 0.5 * dia_bar, stop=b - cover - 0.5 * dia_bar, num=n_x
        )

        # add bottom outer bars
        for x in x_i:
            xs.append(x)
            ys.append(cover + 0.5 * dia_bar)

        # add top outer bars
        for x in x_i:
            xs.append(x)
            ys.append(d - cover - 0.5 * dia_bar)

        # web outer bars
        y_length = d - 2 * cover - dia_bar
        n_y = ceil(y_length / spacing) + 1
        y_i = np.linspace(
            start=cover + 0.5 * dia_bar, stop=d - cover - 0.5 * dia_bar, num=n_y
        )

        # add web outer bars
        for y in y_i[1:-1]:
            xs.append(cover + 0.5 * dia_bar)
            ys.append(y)

        # bottom/top inner bars
        x_length = b - t_w
        n_x = ceil(x_length / spacing) + 1
        x_i = np.linspace(
            start=t_w - cover - 0.5 * dia_bar, stop=b - cover - 0.5 * dia_bar, num=n_x
        )

        # add bottom inner bars
        for x in x_i:
            xs.append(x)
            ys.append(t_f - cover - 0.5 * dia_bar)

        # add top inner bars
        for x in x_i:
            xs.append(x)
            ys.append(d - t_f + cover + 0.5 * dia_bar)

        # web inner bars
        y_length = d - 2 * t_f + 2 * cover + dia_bar
        n_y = ceil(y_length / spacing) + 1
        y_i = np.linspace(
            start=t_f - cover - 0.5 * dia_bar,
            stop=d - t_f + cover + 0.5 * dia_bar,
            num=n_y,
        )

        # add web inner bars
        for y in y_i[1:-1]:
            xs.append(t_w - cover - 0.5 * dia_bar)
            ys.append(y)

    return add_bars(
        geometry=geom, area=area_bar, material=steel_mat, x=xs, y=ys, n=n_circle
    )


def tee_wall(
    d: float,
    b: float,
    t_f: float,
    t_w: float,
    dia_bar: float,
    area_bar: float,
    spacing: float,
    cover: float,
    double: bool = True,
    n_circle: int = 4,
    conc_mat: pre.Material = pre.DEFAULT_MATERIAL,
    steel_mat: pre.Material = pre.DEFAULT_MATERIAL,
) -> CompoundGeometry:
    """Constructs a reinforced concrete tee-shaped wall section.

    Constructs a reinforced concrete tee-shaped wall section of depth ``d``, width
    ``b``, flange thickness ``t_f`` and web thickness ``t_w``. The wall reinforcement
    has a maximum spacing of ``spacing`` and is doubly reinforced if ``double`` is set
    to ``True``, or singly reinforced if it is set to ``False``.

    .. note::
        As the reinforcing bars are described by discretised circles, the area of each
        bar is required to ensure that the correct reinforcing area is provided.

    Args:
        d: Concrete wall depth
        b: Concrete wall width
        t_f: Concrete wall flange thickness
        t_w: Concrete wall web thickness
        dia_bar: Diameter of the reinforcing bars, used for calculating bar placement
        area_bar: Area of the reinforcing bars
        spacing: Maximum spacing of the reinforcement bars, the calculated spacing is
            equal to ``ceil(extent / spacing) + 1``
        cover: Clear cover to the reinforcing bars
        double: If set to ``True``, provides two layers of reinforcement to the wall. If
            set to ``False``, provides a single central layer of reinforcement to the
            wall. Defaults to ``True``.
        n_circle: Number of points used to discretise the circular reinforcing bars.
            Defaults to ``4``.
        conc_mat: Material object to assign to the concrete area. Defaults to
            ``pre.DEFAULT_MATERIAL``.
        steel_mat: Material object to assign to the steel area. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        ValueError: Geometry generation failed

    Returns:
        Reinforced concrete tee-shaped wall section geometry

    Example:
        The following example creates a 1200 mm deep x 1200 mm wide tee-shaped concrete
        wall, with a 150 mm thick flange and 120 mm thick web. The wall is reinforced
        with a single layer of N12-200 with 30 mm cover. A coarse finite element mesh is
        generated to show the different material regions:

        .. plot::
            :include-source: True
            :caption: Reinforced concrete tee-shaped wall section geometry

            from sectionproperties.pre import Material
            from sectionproperties.pre.library import tee_wall
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

            geom = tee_wall(
                d=1200,
                b=1200,
                t_f=150,
                t_w=120,
                dia_bar=12,
                area_bar=110,
                spacing=200,
                cover=30,
                double=False,
                n_circle=12,
                conc_mat=concrete,
                steel_mat=steel,
            )

            geom.create_mesh(mesh_sizes=[0])  # a size of zero creates a coarse mesh
            Section(geometry=geom).plot_mesh()
    """
    # create cee concrete geometry
    geom_outer = primitive_sections.rectangular_section(d=d, b=b, material=conc_mat)
    geom_inner_left = primitive_sections.rectangular_section(
        d=d - t_f, b=0.5 * (b - t_w)
    )
    geom_inner_right = geom_inner_left.mirror_section(
        axis="y", mirror_point=(0.5 * b, 0)
    )
    geom = geom_outer - geom_inner_left - geom_inner_right

    # positions of bars to add
    xs: list[float] = []
    ys: list[float] = []

    # calculate reinforcement positions
    # singly reinforced
    if not double:
        # calculate number of bars along length of wall
        x_length = b - 2 * cover - dia_bar
        n_x = ceil(x_length / spacing) + 1
        y_length = d - 0.5 * t_f - cover - 0.5 * dia_bar
        n_y = ceil(y_length / spacing) + 1

        # calculate position of bars
        x_i = np.linspace(
            start=cover + 0.5 * dia_bar, stop=b - cover - 0.5 * dia_bar, num=n_x
        )
        y_i = np.linspace(start=cover + 0.5 * dia_bar, stop=d - 0.5 * t_f, num=n_y)

        # add top bars
        for x in x_i:
            xs.append(x)
            ys.append(d - 0.5 * t_f)

        # add web bars
        for y in y_i[:-1]:
            xs.append(0.5 * b)
            ys.append(y)
    # doubly reinforced
    else:
        # top outer bars
        x_length = b - 2 * cover - dia_bar
        n_x = ceil(x_length / spacing) + 1
        x_i = np.linspace(
            start=cover + 0.5 * dia_bar, stop=b - cover - 0.5 * dia_bar, num=n_x
        )

        # add top outer bars
        for x in x_i:
            xs.append(x)
            ys.append(d - cover - 0.5 * dia_bar)

        # top inner bars
        x_length = 0.5 * (b - t_w)
        n_x = ceil(x_length / spacing) + 1
        x_i_l = np.linspace(
            start=cover + 0.5 * dia_bar,
            stop=0.5 * (b - t_w) + cover + 0.5 * dia_bar,
            num=n_x,
        )
        x_i_r = np.linspace(
            start=0.5 * (b + t_w) - cover - 0.5 * dia_bar,
            stop=b - cover - 0.5 * dia_bar,
            num=n_x,
        )
        x_i = x_i_l.tolist() + x_i_r.tolist()

        # add top inner bars
        for x in x_i:
            xs.append(x)
            ys.append(d - t_f + cover + 0.5 * dia_bar)

        # web bars
        y_length = d - t_f
        n_y = ceil(y_length / spacing) + 1
        y_i = np.linspace(
            start=cover + 0.5 * dia_bar, stop=d - t_f + cover + 0.5 * dia_bar, num=n_y
        )
        x_i = [
            0.5 * (b - t_w) + cover + 0.5 * dia_bar,
            0.5 * (b + t_w) - cover - 0.5 * dia_bar,
        ]

        # add web bars
        for y in y_i[:-1]:
            for x in x_i:
                xs.append(x)
                ys.append(y)

    return add_bars(
        geometry=geom, area=area_bar, material=steel_mat, x=xs, y=ys, n=n_circle
    )


def single_lift_core(
    d: float,
    b: float,
    t1: float,
    t2: float,
    a: float,
    dia_bar: float,
    area_bar: float,
    spacing: float,
    cover: float,
    double: bool = True,
    n_circle: int = 4,
    conc_mat: pre.Material = pre.DEFAULT_MATERIAL,
    steel_mat: pre.Material = pre.DEFAULT_MATERIAL,
) -> CompoundGeometry:
    """Constructs a reinforced concrete single lift core section.

    Constructs a reinforced concrete single lift core section of depth ``d``, width
    ``b``, top/bottom thickness ``t1``, left/right thickness ``t2`` and door opening
    width ``a``. The wall reinforcement has a maximum spacing of ``spacing`` and is
    doubly reinforced if ``double`` is set to ``True``, or singly reinforced if it is
    set to ``False``.

    .. note::
        As the reinforcing bars are described by discretised circles, the area of each
        bar is required to ensure that the correct reinforcing area is provided.

    Args:
        d: Concrete wall depth
        b: Concrete wall width
        t1: Top/bottom concrete wall thickness
        t2: Left/right concrete wall thickness
        a: Door opening width
        dia_bar: Diameter of the reinforcing bars, used for calculating bar placement
        area_bar: Area of the reinforcing bars
        spacing: Maximum spacing of the reinforcement bars, the calculated spacing is
            equal to ``ceil(extent / spacing) + 1``
        cover: Clear cover to the reinforcing bars
        double: If set to ``True``, provides two layers of reinforcement to the wall. If
            set to ``False``, provides a single central layer of reinforcement to the
            wall. Defaults to ``True``.
        n_circle: Number of points used to discretise the circular reinforcing bars.
            Defaults to ``4``.
        conc_mat: Material object to assign to the concrete area. Defaults to
            ``pre.DEFAULT_MATERIAL``.
        steel_mat: Material object to assign to the steel area. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        ValueError: Geometry generation failed

    Returns:
        Reinforced concrete single lift core section geometry

    Example:
        The following example creates a 2400 mm deep x 1800 mm wide concrete single lift
        core section, with top/bottom walls 150 mm thick, left/right walls 200 mm thick
        and a door opening of 1000 mm. The wall is reinforced with a double layer of
        N20-150 with 35 mm cover. A coarse finite element mesh is generated to show the
        different material regions:

        .. plot::
            :include-source: True
            :caption: Reinforced concrete single lift core section geometry

            from sectionproperties.pre import Material
            from sectionproperties.pre.library import single_lift_core
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

            geom = single_lift_core(
                d=2400,
                b=1800,
                t1=150,
                t2=200,
                a=1000,
                dia_bar=20,
                area_bar=310,
                spacing=150,
                cover=35,
                double=True,
                n_circle=12,
                conc_mat=concrete,
                steel_mat=steel,
            )

            geom.create_mesh(mesh_sizes=[0])  # a size of zero creates a coarse mesh
            Section(geometry=geom).plot_mesh()
    """
    # create lift core geometry
    geom_outer = primitive_sections.rectangular_section(d=d, b=b, material=conc_mat)
    geom_inner = primitive_sections.rectangular_section(
        d=d - 2 * t1, b=b - 2 * t2, material=conc_mat
    ).align_center(align_to=geom_outer)
    geom_door = (
        primitive_sections.rectangular_section(d=a, b=t2, material=conc_mat)
        .align_center(align_to=geom_outer)
        .align_to(other=geom_outer, on="right", inner=True)
    )
    geom = geom_outer - geom_inner - geom_door

    # positions of bars to add
    xs: list[float] = []
    ys: list[float] = []

    # calculate reinforcement positions
    # singly reinforced
    if not double:
        # calculate number of bars along length of wall
        x_length = b - t2
        n_x = ceil(x_length / spacing) + 1
        y_length_left = d - t1
        y_length_right = 0.5 * (d - a) - 0.5 * t1 - cover - 0.5 * dia_bar
        n_y_left = ceil(y_length_left / spacing) + 1
        n_y_right = ceil(y_length_right / spacing) + 1

        # calculate position of bars
        x_i = np.linspace(start=0.5 * t2, stop=b - 0.5 * t2, num=n_x)
        y_i_l = np.linspace(start=0.5 * t1, stop=d - 0.5 * t1, num=n_y_left)
        y_i_r_b = np.linspace(
            start=0.5 * t1, stop=0.5 * (d - a) - cover - 0.5 * dia_bar, num=n_y_right
        )
        y_i_r_t = np.linspace(
            start=0.5 * (d + a) + cover + 0.5 * dia_bar,
            stop=d - 0.5 * t1,
            num=n_y_right,
        )

        # add bot bars
        for x in x_i:
            xs.append(x)
            ys.append(0.5 * t1)

        # add top bars
        for x in x_i:
            xs.append(x)
            ys.append(d - 0.5 * t1)

        # add left bars
        for y in y_i_l[1:-1]:
            xs.append(0.5 * t2)
            ys.append(y)

        # add right bot bars
        for y in y_i_r_b[1:]:
            xs.append(b - 0.5 * t2)
            ys.append(y)

        # add right top bars
        for y in y_i_r_t[:-1]:
            xs.append(b - 0.5 * t2)
            ys.append(y)
    # doubly reinforced
    else:
        # top/bot outer bars
        x_length = b - 2 * cover - dia_bar
        n_x = ceil(x_length / spacing) + 1
        x_i = np.linspace(
            start=cover + 0.5 * dia_bar, stop=b - cover - 0.5 * dia_bar, num=n_x
        )
        y_i = [cover + 0.5 * dia_bar, d - cover - 0.5 * dia_bar]

        # add top/bot outer bars
        for x in x_i:
            for y in y_i:
                xs.append(x)
                ys.append(y)

        # top/bot inner bars
        x_length = b - 2 * t2 + 2 * cover + dia_bar
        n_x = ceil(x_length / spacing) + 1
        x_i = np.linspace(
            start=t2 - cover - 0.5 * dia_bar,
            stop=b - t2 + cover + 0.5 * dia_bar,
            num=n_x,
        )
        y_i = [t1 - cover - 0.5 * dia_bar, d - t1 + cover + 0.5 * dia_bar]

        # add top/bot inner bars
        for x in x_i:
            for y in y_i:
                xs.append(x)
                ys.append(y)

        # left outer bars
        y_length = d - 2 * cover - dia_bar
        n_y = ceil(y_length / spacing) + 1
        y_i = np.linspace(
            start=cover + 0.5 * dia_bar, stop=d - cover - 0.5 * dia_bar, num=n_y
        )

        # add left outer bars
        for y in y_i[1:-1]:
            xs.append(cover + 0.5 * dia_bar)
            ys.append(y)

        # left inner bars
        y_length = d - 2 * t1 + 2 * cover + dia_bar
        n_y = ceil(y_length / spacing) + 1
        y_i = np.linspace(
            start=t1 - cover - 0.5 * dia_bar,
            stop=d - t1 + cover + 0.5 * dia_bar,
            num=n_y,
        )

        # add left inner bars
        for y in y_i[1:-1]:
            xs.append(t2 - cover - 0.5 * dia_bar)
            ys.append(y)

        # right outer bars
        y_length = 0.5 * (d - a) - 2 * cover - dia_bar
        n_y = ceil(y_length / spacing) + 1
        y_i_b = np.linspace(
            start=cover + 0.5 * dia_bar,
            stop=0.5 * (d - a) - cover - 0.5 * dia_bar,
            num=n_y,
        )
        y_i_t = np.linspace(
            start=0.5 * (d + a) + cover + 0.5 * dia_bar,
            stop=d - cover - 0.5 * dia_bar,
            num=n_y,
        )
        y_i = y_i_b.tolist() + y_i_t.tolist()

        # add right outer bars
        for y in y_i[1:-1]:
            xs.append(b - cover - 0.5 * dia_bar)
            ys.append(y)

        # right inner bars
        y_length = 0.5 * (d - a) - t1
        n_y = ceil(y_length / spacing) + 1
        y_i_b = np.linspace(
            start=t1 - cover - 0.5 * dia_bar,
            stop=0.5 * (d - a) - cover - 0.5 * dia_bar,
            num=n_y,
        )
        y_i_t = np.linspace(
            start=0.5 * (d + a) + cover + 0.5 * dia_bar,
            stop=d - t1 + cover + 0.5 * dia_bar,
            num=n_y,
        )
        y_i = y_i_b.tolist() + y_i_t.tolist()

        # add right inner bars
        for y in y_i[1:-1]:
            xs.append(b - t2 + cover + 0.5 * dia_bar)
            ys.append(y)

    return add_bars(
        geometry=geom, area=area_bar, material=steel_mat, x=xs, y=ys, n=n_circle
    )


def double_lift_core_a(
    d: float,
    b: float,
    t1: float,
    t2: float,
    a1: float,
    a2: float,
    dia_bar: float,
    area_bar: float,
    spacing: float,
    cover: float,
    double: bool = True,
    n_circle: int = 4,
    conc_mat: pre.Material = pre.DEFAULT_MATERIAL,
    steel_mat: pre.Material = pre.DEFAULT_MATERIAL,
) -> CompoundGeometry:
    """Constructs a reinforced concrete double lift core (type A) section.

    Constructs a reinforced concrete double lift core (type A) section of depth ``d``,
    width ``b``, top/bottom thickness ``t1``, left/right thickness ``t2``, door opening
    width ``a1`` and door pier weidth ``a2``. The wall reinforcement has a maximum
    spacing of ``spacing`` and is doubly reinforced if ``double`` is set to ``True``, or
    singly reinforced if it is set to ``False``.

    Type A refers to the fact that there is no wall dividing the two lift cores.

    .. note::
        As the reinforcing bars are described by discretised circles, the area of each
        bar is required to ensure that the correct reinforcing area is provided.

    Args:
        d: Concrete wall depth
        b: Concrete wall width
        t1: Top/bottom concrete wall thickness
        t2: Left/right concrete wall thickness
        a1: Door opening width
        a2: Door pier width
        dia_bar: Diameter of the reinforcing bars, used for calculating bar placement
        area_bar: Area of the reinforcing bars
        spacing: Maximum spacing of the reinforcement bars, the calculated spacing is
            equal to ``ceil(extent / spacing) + 1``
        cover: Clear cover to the reinforcing bars
        double: If set to ``True``, provides two layers of reinforcement to the wall. If
            set to ``False``, provides a single central layer of reinforcement to the
            wall. Defaults to ``True``.
        n_circle: Number of points used to discretise the circular reinforcing bars.
            Defaults to ``4``.
        conc_mat: Material object to assign to the concrete area. Defaults to
            ``pre.DEFAULT_MATERIAL``.
        steel_mat: Material object to assign to the steel area. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        ValueError: Geometry generation failed
        ValueError: If the lift door opening geometry is invalid

    Returns:
        Reinforced concrete double lift core (type A) section geometry

    Example:
        The following example creates a 3200 mm deep x 1800 mm wide double lift core
        (type A) concrete section, with top/bottom walls 200 mm thick, left/right walls
        150 mm thick, a door opening of 900 mm and a dor pier width of 600 mm. The wall
        is reinforced with a double layer of N16-200 with 30 mm cover. A coarse finite
        element mesh is generated to show the different material regions:

        .. plot::
            :include-source: True
            :caption: Reinforced concrete double lift core (type A) section geometry

            from sectionproperties.pre import Material
            from sectionproperties.pre.library import double_lift_core_a
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

            geom = double_lift_core_a(
                d=3200,
                b=1800,
                t1=200,
                t2=150,
                a1=900,
                a2=600,
                dia_bar=16,
                area_bar=200,
                spacing=200,
                cover=30,
                double=True,
                n_circle=12,
                conc_mat=concrete,
                steel_mat=steel,
            )

            geom.create_mesh(mesh_sizes=[0])  # a size of zero creates a coarse mesh
            Section(geometry=geom).plot_mesh()
    """
    # check geometry
    if 2 * a1 + a2 > d - 2 * t1:
        msg = "(2 * a1 + a2) is larger than the lift opening."
        raise ValueError(msg)

    # create lift core geometry
    geom_outer = primitive_sections.rectangular_section(d=d, b=b, material=conc_mat)
    geom_inner = primitive_sections.rectangular_section(
        d=d - 2 * t1, b=b - 2 * t2, material=conc_mat
    ).align_center(align_to=geom_outer)
    geom_door = primitive_sections.rectangular_section(d=a1, b=t2, material=conc_mat)
    nib = 0.5 * (d - 2 * a1 - a2 - 2 * t1)
    geom_door_t = geom_door.shift_section(x_offset=b - t2, y_offset=d - t1 - nib - a1)
    geom_door_b = geom_door.shift_section(x_offset=b - t2, y_offset=t1 + nib)
    geom = geom_outer - geom_inner - geom_door_t - geom_door_b

    # positions of bars to add
    xs: list[float] = []
    ys: list[float] = []

    # calculate reinforcement positions
    # singly reinforced
    if not double:
        # calculate number of bars along length of wall
        x_length = b - t2
        n_x = ceil(x_length / spacing) + 1
        y_length = d - t1
        y_nib = nib + 0.5 * t1
        y_pier = a2 - 2 * cover - dia_bar
        n_y = ceil(y_length / spacing) + 1
        n_nib = ceil(y_nib / spacing) + 1
        n_pier = ceil(y_pier / spacing) + 1

        # calculate position of bars
        x_i = np.linspace(start=0.5 * t2, stop=b - 0.5 * t2, num=n_x)
        y_i = np.linspace(start=0.5 * t1, stop=d - 0.5 * t1, num=n_y)
        y_nib_b = np.linspace(
            start=0.5 * t1, stop=t1 + nib - cover - 0.5 * dia_bar, num=n_nib
        )
        y_nib_t = np.linspace(
            start=d - t1 - nib + cover + 0.5 * dia_bar, stop=d - 0.5 * t1, num=n_nib
        )
        y_pier = np.linspace(
            start=0.5 * (d - a2) + cover + 0.5 * dia_bar,
            stop=0.5 * (d + a2) - cover - 0.5 * dia_bar,
            num=n_pier,
        )

        # add bot bars
        for x in x_i:
            xs.append(x)
            ys.append(0.5 * t1)

        # add top bars
        for x in x_i:
            xs.append(x)
            ys.append(d - 0.5 * t1)

        # add left bars
        for y in y_i[1:-1]:
            xs.append(0.5 * t2)
            ys.append(y)

        # add nib bot bars
        for y in y_nib_b[1:]:
            xs.append(b - 0.5 * t2)
            ys.append(y)

        # add nib top bars
        for y in y_nib_t[:-1]:
            xs.append(b - 0.5 * t2)
            ys.append(y)

        # add pier bars
        for y in y_pier:
            xs.append(b - 0.5 * t2)
            ys.append(y)
    # doubly reinforced
    else:
        # top/bot outer bars
        x_length = b - 2 * cover - dia_bar
        n_x = ceil(x_length / spacing) + 1
        x_i = np.linspace(
            start=cover + 0.5 * dia_bar, stop=b - cover - 0.5 * dia_bar, num=n_x
        )
        y_i = [cover + 0.5 * dia_bar, d - cover - 0.5 * dia_bar]

        # add top/bot outer bars
        for x in x_i:
            for y in y_i:
                xs.append(x)
                ys.append(y)

        # top/bot inner bars
        x_length = b - 2 * t2 + 2 * cover + dia_bar
        n_x = ceil(x_length / spacing) + 1
        x_i = np.linspace(
            start=t2 - cover - 0.5 * dia_bar,
            stop=b - t2 + cover + 0.5 * dia_bar,
            num=n_x,
        )
        y_i = [t1 - cover - 0.5 * dia_bar, d - t1 + cover + 0.5 * dia_bar]

        # add top/bot inner bars
        for x in x_i:
            for y in y_i:
                xs.append(x)
                ys.append(y)

        # left outer bars
        y_length = d - 2 * cover - dia_bar
        n_y = ceil(y_length / spacing) + 1
        y_i = np.linspace(
            start=cover + 0.5 * dia_bar, stop=d - cover - 0.5 * dia_bar, num=n_y
        )

        # add left outer bars
        for y in y_i[1:-1]:
            xs.append(cover + 0.5 * dia_bar)
            ys.append(y)

        # left inner bars
        y_length = d - 2 * t1 + 2 * cover + dia_bar
        n_y = ceil(y_length / spacing) + 1
        y_i = np.linspace(
            start=t1 - cover - 0.5 * dia_bar,
            stop=d - t1 + cover + 0.5 * dia_bar,
            num=n_y,
        )

        # add left inner bars
        for y in y_i[1:-1]:
            xs.append(t2 - cover - 0.5 * dia_bar)
            ys.append(y)

        # nib outer bars
        y_nib = nib + t1 - 2 * cover - dia_bar
        n_nib = ceil(y_nib / spacing) + 1
        y_nib_b = np.linspace(
            start=cover + 0.5 * dia_bar,
            stop=t1 + nib - cover - 0.5 * dia_bar,
            num=n_nib,
        )
        y_nib_t = np.linspace(
            start=d - t1 - nib + cover + 0.5 * dia_bar,
            stop=d - cover - 0.5 * dia_bar,
            num=n_nib,
        )
        y_nib = y_nib_b.tolist() + y_nib_t.tolist()

        # add nib outer bars
        for y in y_nib[1:-1]:
            xs.append(b - cover - 0.5 * dia_bar)
            ys.append(y)

        # nib inner bars
        y_nib = nib
        n_nib = ceil(y_nib / spacing) + 1
        y_nib_b = np.linspace(
            start=t1 - cover - 0.5 * dia_bar,
            stop=t1 + nib - cover - 0.5 * dia_bar,
            num=n_nib,
        )
        y_nib_t = np.linspace(
            start=d - t1 - nib + cover + 0.5 * dia_bar,
            stop=d - t1 + cover + 0.5 * dia_bar,
            num=n_nib,
        )
        y_nib = y_nib_b.tolist() + y_nib_t.tolist()

        # add nib outer bars
        for y in y_nib[1:-1]:
            xs.append(b - t2 + cover + 0.5 * dia_bar)
            ys.append(y)

        # pier bars
        y_pier = a2 - 2 * cover - dia_bar
        n_pier = ceil(y_pier / spacing) + 1
        y_pier = np.linspace(
            start=0.5 * (d - a2) + cover + 0.5 * dia_bar,
            stop=0.5 * (d + a2) - cover - 0.5 * dia_bar,
            num=n_pier,
        )
        x_pier = [b - t2 + cover + 0.5 * dia_bar, b - cover - 0.5 * dia_bar]

        # add pier bars
        for y in y_pier:
            for x in x_pier:
                xs.append(x)
                ys.append(y)

    return add_bars(
        geometry=geom, area=area_bar, material=steel_mat, x=xs, y=ys, n=n_circle
    )


def double_lift_core_b(
    d: float,
    b: float,
    t1: float,
    t2: float,
    t3: float,
    a: float,
    dia_bar: float,
    area_bar: float,
    spacing: float,
    cover: float,
    double: bool = True,
    n_circle: int = 4,
    conc_mat: pre.Material = pre.DEFAULT_MATERIAL,
    steel_mat: pre.Material = pre.DEFAULT_MATERIAL,
) -> CompoundGeometry:
    """Constructs a reinforced concrete double lift core (type B) section.

    Constructs a reinforced concrete double lift core (type B) section of depth ``d``,
    width ``b``, top/bottom thickness ``t1``, left/right thickness ``t2``, central wall
    thickness ``t3`` and door opening width ``a``. The wall reinforcement has a maximum
    spacing of ``spacing`` and is doubly reinforced if ``double`` is set to ``True``, or
    singly reinforced if it is set to ``False``.

    Type B refers to the fact that there is a central wall dividing the two lift cores.

    .. note::
        As the reinforcing bars are described by discretised circles, the area of each
        bar is required to ensure that the correct reinforcing area is provided.

    Args:
        d: Concrete wall depth
        b: Concrete wall width
        t1: Top/bottom concrete wall thickness
        t2: Left/right concrete wall thickness
        t3: Central concrete wall thickness
        a: Door opening width
        dia_bar: Diameter of the reinforcing bars, used for calculating bar placement
        area_bar: Area of the reinforcing bars
        spacing: Maximum spacing of the reinforcement bars, the calculated spacing is
            equal to ``ceil(extent / spacing) + 1``
        cover: Clear cover to the reinforcing bars
        double: If set to ``True``, provides two layers of reinforcement to the wall. If
            set to ``False``, provides a single central layer of reinforcement to the
            wall. Defaults to ``True``.
        n_circle: Number of points used to discretise the circular reinforcing bars.
            Defaults to ``4``.
        conc_mat: Material object to assign to the concrete area. Defaults to
            ``pre.DEFAULT_MATERIAL``.
        steel_mat: Material object to assign to the steel area. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        ValueError: Geometry generation failed
        ValueError: If the lift door opening geometry is invalid

    Returns:
        Reinforced concrete double lift core (type B) section geometry

    Example:
        The following example creates a 6000 mm deep x 3000 mm wide double lift core
        (type B) concrete section, with top/bottom walls 250 mm thick, left/right walls
        220 mm thick, central wall 180 mm thick and a door opening of 1500 mm. The wall
        is reinforced with a single layer of N24-300 with 50 mm cover. A coarse finite
        element mesh is generated to show the different material regions:

        .. plot::
            :include-source: True
            :caption: Reinforced concrete double lift core (type B) section geometry

            from sectionproperties.pre import Material
            from sectionproperties.pre.library import double_lift_core_b
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

            geom = double_lift_core_b(
                d=6000,
                b=3000,
                t1=250,
                t2=220,
                t3=180,
                a=1500,
                dia_bar=24,
                area_bar=450,
                spacing=300,
                cover=50,
                double=True,
                n_circle=12,
                conc_mat=concrete,
                steel_mat=steel,
            )

            geom.create_mesh(mesh_sizes=[0])  # a size of zero creates a coarse mesh
            Section(geometry=geom).plot_mesh()
    """
    # check geometry
    if a > 0.5 * (d - 2 * t1 - t3):
        msg = "a is larger than a single lift opening."
        raise ValueError(msg)

    # calculate nib and a2
    nib = 0.25 * (d - 2 * t1 - 2 * a - t3)
    a2 = 2 * nib + t3

    # create lift core geometry
    geom_outer = primitive_sections.rectangular_section(d=d, b=b, material=conc_mat)
    geom_inner = primitive_sections.rectangular_section(
        d=0.5 * (d - 2 * t1 - t3), b=b - 2 * t2, material=conc_mat
    )
    geom_door = primitive_sections.rectangular_section(d=a, b=t2, material=conc_mat)
    geom = (
        geom_outer
        - geom_inner.shift_section(x_offset=t2, y_offset=t1)
        - geom_inner.shift_section(x_offset=t2, y_offset=0.5 * (d + t3))
        - geom_door.shift_section(x_offset=b - t2, y_offset=d - t1 - nib - a)
        - geom_door.shift_section(x_offset=b - t2, y_offset=t1 + nib)
    )

    # positions of bars to add
    xs: list[float] = []
    ys: list[float] = []

    # calculate reinforcement positions
    # singly reinforced
    if not double:
        # calculate number of bars along length of wall
        x_length = b - t2
        n_x = ceil(x_length / spacing) + 1
        y_length = d - t1
        y_nib = nib + 0.5 * t1
        y_pier = a2 - 2 * cover - dia_bar
        n_y = ceil(y_length / spacing) + 1
        n_nib = ceil(y_nib / spacing) + 1
        n_pier = ceil(y_pier / spacing) + 1

        # calculate position of bars
        x_i = np.linspace(start=0.5 * t2, stop=b - 0.5 * t2, num=n_x)
        y_i = np.linspace(start=0.5 * t1, stop=d - 0.5 * t1, num=n_y)
        y_nib_b = np.linspace(
            start=0.5 * t1, stop=t1 + nib - cover - 0.5 * dia_bar, num=n_nib
        )
        y_nib_t = np.linspace(
            start=d - t1 - nib + cover + 0.5 * dia_bar, stop=d - 0.5 * t1, num=n_nib
        )
        y_pier = np.linspace(
            start=0.5 * (d - a2) + cover + 0.5 * dia_bar,
            stop=0.5 * (d + a2) - cover - 0.5 * dia_bar,
            num=n_pier,
        )

        # add bot bars
        for x in x_i:
            xs.append(x)
            ys.append(0.5 * t1)

        # add mid bars
        for x in x_i[1:-1]:
            xs.append(x)
            ys.append(0.5 * d)

        # add top bars
        for x in x_i:
            xs.append(x)
            ys.append(d - 0.5 * t1)

        # add left bars
        for y in y_i[1:-1]:
            xs.append(0.5 * t2)
            ys.append(y)

        # add nib bot bars
        for y in y_nib_b[1:]:
            xs.append(b - 0.5 * t2)
            ys.append(y)

        # add nib top bars
        for y in y_nib_t[:-1]:
            xs.append(b - 0.5 * t2)
            ys.append(y)

        # add pier bars
        for y in y_pier:
            xs.append(b - 0.5 * t2)
            ys.append(y)
    # doubly reinforced
    else:
        # top/bot outer bars
        x_length = b - 2 * cover - dia_bar
        n_x = ceil(x_length / spacing) + 1
        x_i = np.linspace(
            start=cover + 0.5 * dia_bar, stop=b - cover - 0.5 * dia_bar, num=n_x
        )
        y_i = [cover + 0.5 * dia_bar, d - cover - 0.5 * dia_bar]

        # add top/bot outer bars
        for x in x_i:
            for y in y_i:
                xs.append(x)
                ys.append(y)

        # top/bot/central inner bars
        x_length = b - 2 * t2 + 2 * cover + dia_bar
        n_x = ceil(x_length / spacing) + 1
        x_i = np.linspace(
            start=t2 - cover - 0.5 * dia_bar,
            stop=b - t2 + cover + 0.5 * dia_bar,
            num=n_x,
        )
        y_i = [
            t1 - cover - 0.5 * dia_bar,
            0.5 * (d - t3) + cover + 0.5 * dia_bar,
            0.5 * (d + t3) - cover - 0.5 * dia_bar,
            d - t1 + cover + 0.5 * dia_bar,
        ]

        # add top/bot inner bars
        for x in x_i:
            for y in y_i:
                xs.append(x)
                ys.append(y)

        # left outer bars
        y_length = d - 2 * cover - dia_bar
        n_y = ceil(y_length / spacing) + 1
        y_i = np.linspace(
            start=cover + 0.5 * dia_bar, stop=d - cover - 0.5 * dia_bar, num=n_y
        )

        # add left outer bars
        for y in y_i[1:-1]:
            xs.append(cover + 0.5 * dia_bar)
            ys.append(y)

        # left inner bars
        y_length = 0.5 * (d - 2 * t1 - t3) + 2 * cover + dia_bar
        n_y = ceil(y_length / spacing) + 1
        y_i_b = np.linspace(
            start=t1 - cover - 0.5 * dia_bar,
            stop=0.5 * (d - t3) + cover + 0.5 * dia_bar,
            num=n_y,
        )
        y_i_t = np.linspace(
            start=0.5 * (d + t3) - cover - 0.5 * dia_bar,
            stop=d - t1 + cover + 0.5 * dia_bar,
            num=n_y,
        )
        y_i = y_i_b[1:-1].tolist() + y_i_t[1:-1].tolist()

        # add left inner bars
        for y in y_i:
            xs.append(t2 - cover - 0.5 * dia_bar)
            ys.append(y)

        # nib outer bars
        y_nib = nib + t1 - 2 * cover - dia_bar
        n_nib = ceil(y_nib / spacing) + 1
        y_nib_b = np.linspace(
            start=cover + 0.5 * dia_bar,
            stop=t1 + nib - cover - 0.5 * dia_bar,
            num=n_nib,
        )
        y_nib_t = np.linspace(
            start=d - t1 - nib + cover + 0.5 * dia_bar,
            stop=d - cover - 0.5 * dia_bar,
            num=n_nib,
        )
        y_nib = y_nib_b.tolist() + y_nib_t.tolist()

        # add nib outer bars
        for y in y_nib[1:-1]:
            xs.append(b - cover - 0.5 * dia_bar)
            ys.append(y)

        # nib inner bars
        y_nib = nib
        n_nib = ceil(y_nib / spacing) + 1
        y_nib_b = np.linspace(
            start=t1 - cover - 0.5 * dia_bar,
            stop=t1 + nib - cover - 0.5 * dia_bar,
            num=n_nib,
        )
        y_nib_t = np.linspace(
            start=d - t1 - nib + cover + 0.5 * dia_bar,
            stop=d - t1 + cover + 0.5 * dia_bar,
            num=n_nib,
        )
        y_nib = y_nib_b.tolist() + y_nib_t.tolist()

        # add nib outer bars
        for y in y_nib[1:-1]:
            xs.append(b - t2 + cover + 0.5 * dia_bar)
            ys.append(y)

        # outer pier bars
        y_pier = a2 - 2 * cover - dia_bar
        n_pier = ceil(y_pier / spacing) + 1
        y_pier = np.linspace(
            start=0.5 * (d - a2) + cover + 0.5 * dia_bar,
            stop=0.5 * (d + a2) - cover - 0.5 * dia_bar,
            num=n_pier,
        )

        # add outer pier bars
        for y in y_pier:
            xs.append(b - cover - 0.5 * dia_bar)
            ys.append(y)

        # inner pier bars
        y_pier = nib - 2 * cover - dia_bar
        n_pier = ceil(y_pier / spacing) + 1
        y_pier_b = np.linspace(
            start=0.5 * (d - a2) + cover + 0.5 * dia_bar,
            stop=0.5 * d - cover - 0.5 * dia_bar,
            num=n_pier,
        )
        y_pier_t = np.linspace(
            start=0.5 * d + cover + 0.5 * dia_bar,
            stop=0.5 * (d + a2) - cover - 0.5 * dia_bar,
            num=n_pier,
        )
        y_pier = y_pier_b[:-1].tolist() + y_pier_t[1:].tolist()

        # add outer pier bars
        for y in y_pier:
            xs.append(b - t2 + cover + 0.5 * dia_bar)
            ys.append(y)

    return add_bars(
        geometry=geom, area=area_bar, material=steel_mat, x=xs, y=ys, n=n_circle
    )


def stairwell(
    d: float,
    b: float,
    t1: float,
    t2: float,
    a: float,
    dia_bar: float,
    area_bar: float,
    spacing: float,
    cover: float,
    double: bool = True,
    n_circle: int = 4,
    conc_mat: pre.Material = pre.DEFAULT_MATERIAL,
    steel_mat: pre.Material = pre.DEFAULT_MATERIAL,
) -> CompoundGeometry:
    """Constructs a reinforced concrete stairwell section.

    Constructs a reinforced concrete stairwell section of depth ``d``, width ``b``,
    top/bottom thickness ``t1``, left/right thickness ``t2`` and door opening width
    ``a``. The wall reinforcement has a maximum spacing of ``spacing`` and is doubly
    reinforced if ``double`` is set to ``True``, or singly reinforced if it is set to
    ``False``.

    .. note::
        As the reinforcing bars are described by discretised circles, the area of each
        bar is required to ensure that the correct reinforcing area is provided.

    Args:
        d: Concrete wall depth
        b: Concrete wall width
        t1: Top/bottom concrete wall thickness
        t2: Left/right concrete wall thickness
        a: Door opening width
        dia_bar: Diameter of the reinforcing bars, used for calculating bar placement
        area_bar: Area of the reinforcing bars
        spacing: Maximum spacing of the reinforcement bars, the calculated spacing is
            equal to ``ceil(extent / spacing) + 1``
        cover: Clear cover to the reinforcing bars
        double: If set to ``True``, provides two layers of reinforcement to the wall. If
            set to ``False``, provides a single central layer of reinforcement to the
            wall. Defaults to ``True``.
        n_circle: Number of points used to discretise the circular reinforcing bars.
            Defaults to ``4``.
        conc_mat: Material object to assign to the concrete area. Defaults to
            ``pre.DEFAULT_MATERIAL``.
        steel_mat: Material object to assign to the steel area. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        ValueError: Geometry generation failed

    Returns:
        Reinforced concrete stairwell section geometry

    Example:
        The following example creates a 5000 mm deep x 2500 mm wide concrete stairwell
        section, with top/bottom walls 180 mm thick, left/right walls 200 mm thick and a
        door opening of 900 mm. The wall is reinforced with a double layer of N16-200
        with 30 mm cover. A coarse finite element mesh is generated to show the
        different material regions:

        .. plot::
            :include-source: True
            :caption: Reinforced concrete stairwell section geometry

            from sectionproperties.pre import Material
            from sectionproperties.pre.library import stairwell
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

            geom = stairwell(
                d=5000,
                b=2500,
                t1=180,
                t2=200,
                a=900,
                dia_bar=16,
                area_bar=200,
                spacing=200,
                cover=20,
                double=True,
                n_circle=12,
                conc_mat=concrete,
                steel_mat=steel,
            )

            geom.create_mesh(mesh_sizes=[0])  # a size of zero creates a coarse mesh
            Section(geometry=geom).plot_mesh()
    """
    # create stairwell geometry
    geom_outer = primitive_sections.rectangular_section(d=d, b=b, material=conc_mat)
    geom_inner = primitive_sections.rectangular_section(
        d=d - 2 * t1, b=b - 2 * t2, material=conc_mat
    ).align_center(align_to=geom_outer)
    geom_door = primitive_sections.rectangular_section(
        d=a, b=t2, material=conc_mat
    ).shift_section(x_offset=b - t2, y_offset=d - t1 - a)
    geom = geom_outer - geom_inner - geom_door

    # positions of bars to add
    xs: list[float] = []
    ys: list[float] = []

    # calculate reinforcement positions
    # singly reinforced
    if not double:
        # calculate number of bars along length of wall
        x_length = b - t2
        n_x = ceil(x_length / spacing) + 1
        y_length_left = d - t1
        y_length_right = d - 1.5 * t1 - a - cover - 0.5 * dia_bar
        n_y_left = ceil(y_length_left / spacing) + 1
        n_y_right = ceil(y_length_right / spacing) + 1

        # calculate position of bars
        x_i = np.linspace(start=0.5 * t2, stop=b - 0.5 * t2, num=n_x)
        y_i_l = np.linspace(start=0.5 * t1, stop=d - 0.5 * t1, num=n_y_left)
        y_i_r = np.linspace(
            start=0.5 * t1, stop=d - t1 - a - cover - 0.5 * dia_bar, num=n_y_right
        )

        # add bot bars
        for x in x_i:
            xs.append(x)
            ys.append(0.5 * t1)

        # add top bars
        for x in x_i:
            xs.append(x)
            ys.append(d - 0.5 * t1)

        # add left bars
        for y in y_i_l[1:-1]:
            xs.append(0.5 * t2)
            ys.append(y)

        # add right bars
        for y in y_i_r[1:]:
            xs.append(b - 0.5 * t2)
            ys.append(y)
    # doubly reinforced
    else:
        # top/bot outer bars
        x_length = b - 2 * cover - dia_bar
        n_x = ceil(x_length / spacing) + 1
        x_i = np.linspace(
            start=cover + 0.5 * dia_bar, stop=b - cover - 0.5 * dia_bar, num=n_x
        )
        y_i = [cover + 0.5 * dia_bar, d - cover - 0.5 * dia_bar]

        # add top/bot outer bars
        for x in x_i:
            for y in y_i:
                xs.append(x)
                ys.append(y)

        # top inner bars
        x_length = b - t2
        n_x = ceil(x_length / spacing) + 1
        x_i = np.linspace(
            start=t2 - cover - 0.5 * dia_bar,
            stop=b - cover - 0.5 * dia_bar,
            num=n_x,
        )

        # add top inner bars
        for x in x_i:
            xs.append(x)
            ys.append(d - t1 + cover + 0.5 * dia_bar)

        # bot inner bars
        x_length = b - 2 * t2 + 2 * cover + dia_bar
        n_x = ceil(x_length / spacing) + 1
        x_i = np.linspace(
            start=t2 - cover - 0.5 * dia_bar,
            stop=b - t2 + cover + 0.5 * dia_bar,
            num=n_x,
        )

        # add bot inner bars
        for x in x_i:
            xs.append(x)
            ys.append(t1 - cover - 0.5 * dia_bar)

        # left outer bars
        y_length = d - 2 * cover - dia_bar
        n_y = ceil(y_length / spacing) + 1
        y_i = np.linspace(
            start=cover + 0.5 * dia_bar, stop=d - cover - 0.5 * dia_bar, num=n_y
        )

        # add left outer bars
        for y in y_i[1:-1]:
            xs.append(cover + 0.5 * dia_bar)
            ys.append(y)

        # left inner bars
        y_length = d - 2 * t1 + 2 * cover + dia_bar
        n_y = ceil(y_length / spacing) + 1
        y_i = np.linspace(
            start=t1 - cover - 0.5 * dia_bar,
            stop=d - t1 + cover + 0.5 * dia_bar,
            num=n_y,
        )

        # add left inner bars
        for y in y_i[1:-1]:
            xs.append(t2 - cover - 0.5 * dia_bar)
            ys.append(y)

        # right inner bars
        y_length = d - 2 * t1 - a
        n_y = ceil(y_length / spacing) + 1
        y_i = np.linspace(
            start=t1 - cover - 0.5 * dia_bar,
            stop=d - t1 - a - cover - 0.5 * dia_bar,
            num=n_y,
        )

        # add right inner bars
        for y in y_i[1:]:
            xs.append(b - t2 + cover + 0.5 * dia_bar)
            ys.append(y)

        # right outer bars
        y_length = d - t1 - a - cover - 0.5 * dia_bar
        n_y = ceil(y_length / spacing) + 1
        y_i = np.linspace(
            start=cover + 0.5 * dia_bar,
            stop=d - t1 - a - cover - 0.5 * dia_bar,
            num=n_y,
        )

        # add right outer bars
        for y in y_i[1:]:
            xs.append(b - cover - 0.5 * dia_bar)
            ys.append(y)

    return add_bars(
        geometry=geom, area=area_bar, material=steel_mat, x=xs, y=ys, n=n_circle
    )


def add_bar(
    geometry: Geometry | CompoundGeometry,
    area: float,
    material: pre.Material,
    x: float,
    y: float,
    n: int,
) -> CompoundGeometry:
    """Adds a reinforcing bar to a ``sectionproperties`` geometry.

    First removes the geometry through a subtraction operation, then adds the geometry
    on top of the newly created hole. This method avoids the doubling up of

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


def add_bars(
    geometry: Geometry | CompoundGeometry,
    area: float,
    material: pre.Material,
    x: list[float],
    y: list[float],
    n: int,
) -> CompoundGeometry:
    """Adds a list of reinforcing bars to a ``sectionproperties`` geometry.

    First removes the geometry through a subtraction operation, then adds the geometry
    on top of the newly created hole. This method avoids the doubling up of

    Args:
        geometry: Reinforced concrete geometry to which the new bars will be added
        area: Bar cross-sectional area
        material: Material object for the bars
        x: x-positions of the bars
        y: y-positions of the bars
        n: Number of points to discretise the bar circle

    Returns:
        Geometry object with added bars
    """
    bars = CompoundGeometry(
        geoms=[
            primitive_sections.circular_section_by_area(
                area=area, n=n, material=material
            ).shift_section(x_offset=x_i, y_offset=y_i)
            for x_i, y_i in zip(x, y, strict=False)
        ]
    )

    return (geometry - bars) + bars
