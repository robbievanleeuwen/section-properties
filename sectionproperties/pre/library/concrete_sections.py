import numpy as np
from shapely.geometry import Polygon
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
    area_top: float = None,
    area_bot: float = None,
    conc_mat: pre.Material = pre.DEFAULT_MATERIAL,
    steel_mat: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.CompoundGeometry:
    """Constructs a concrete rectangular section of width *b* and depth *d*, with
    *n_top* top steel bars of diameter *dia_top*, *n_bot* bottom steel bars of diameter
    *dia_bot*, discretised with *n_circle* points with equal side and top/bottom
    *cover* to the steel.

    :param float b: Concrete section width
    :param float d: Concrete section depth
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
    :param Optional[sectionproperties.pre.pre.Material] conc_mat: Material to
        associate with the concrete
    :param Optional[sectionproperties.pre.pre.Material] steel_mat: Material to
        associate with the steel

    :raises ValueErorr: If the number of bars is not greater than or equal to 2 in an
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

    # calculate reinforcing bar dimensions
    x_i_top = cover + dia_top / 2
    x_i_bot = cover + dia_bot / 2
    spacing_top = (b - 2 * cover - dia_top) / (n_top - 1)
    spacing_bot = (b - 2 * cover - dia_bot) / (n_bot - 1)

    # add top bars
    for i in range(n_top):
        if area_top:
            bar = primitive_sections.circular_section_by_area(
                area=area_top, n=n_circle, material=steel_mat
            )
        else:
            bar = primitive_sections.circular_section(
                d=dia_top, n=n_circle, material=steel_mat
            )

        bar = bar.shift_section(
            x_offset=x_i_top + spacing_top * i, y_offset=d - cover - dia_top / 2
        )

        geom = (geom - bar) + bar

    # add bot bars
    for i in range(n_bot):
        if area_bot:
            bar = primitive_sections.circular_section_by_area(
                area=area_bot, n=n_circle, material=steel_mat
            )
        else:
            bar = primitive_sections.circular_section(
                d=dia_bot, n=n_circle, material=steel_mat
            )

        bar = bar.shift_section(
            x_offset=x_i_bot + spacing_bot * i, y_offset=cover + dia_bot / 2
        )

        geom = (geom - bar) + bar

    return geom


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
    area_top: float = None,
    area_bot: float = None,
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
    :param Optional[sectionproperties.pre.pre.Material] conc_mat: Material to associate
        with the concrete
    :param Optional[sectionproperties.pre.pre.Material] steel_mat: Material to
        associate with the steel

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

    # add top bars
    for i in range(n_top):
        if area_top:
            bar = primitive_sections.circular_section_by_area(
                area=area_top, n=n_circle, material=steel_mat
            )
        else:
            bar = primitive_sections.circular_section(
                d=dia_top, n=n_circle, material=steel_mat
            )

        bar = bar.shift_section(
            x_offset=x_i_top + spacing_top * i, y_offset=d - cover - dia_top / 2
        )

        geom = (geom - bar) + bar

    # add bot bars
    for i in range(n_bot):
        if area_bot:
            bar = primitive_sections.circular_section_by_area(
                area=area_bot, n=n_circle, material=steel_mat
            )
        else:
            bar = primitive_sections.circular_section(
                d=dia_bot, n=n_circle, material=steel_mat
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
    area_conc: float = None,
    area_bar: float = None,
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
    :param Optional[sectionproperties.pre.pre.Material] conc_mat: Material to associate
        with the concrete
    :param Optional[sectionproperties.pre.pre.Material] steel_mat: Material to
        associate with the steel

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

    for i in range(n_bar):
        if area_bar:
            bar = primitive_sections.circular_section_by_area(
                area=area_bar, n=n_circle, material=steel_mat
            )
        else:
            bar = primitive_sections.circular_section(
                d=dia, n=n_circle, material=steel_mat
            )

        bar = bar.shift_section(
            x_offset=r * np.cos(i * d_theta), y_offset=r * np.sin(i * d_theta)
        )

        geom = (geom - bar) + bar

    return geom
