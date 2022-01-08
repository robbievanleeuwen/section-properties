import numpy as np
from shapely.geometry import Polygon
import sectionproperties.pre.pre as pre
import sectionproperties.pre.geometry as geometry
import sectionproperties.pre.library.primitive_sections as primitive_sections


def concrete_rectangular_section(
    b: float,
    d: float,
    dia: float,
    n_bar: int,
    n_circle: int,
    cover: float,
    conc_mat: pre.Material = pre.DEFAULT_MATERIAL,
    steel_mat: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.CompoundGeometry:
    """Constructs a concrete rectangular section of width *b* and depth *d*, with *n_bar* steel bars
    of diameter *dia*, discretised with *n_circle* points with equal side and bottom *cover* to the
    steel.

    :param float b: Concrete section depth
    :param float d: Concrete section width
    :param float dia: Diameter of the steel reinforcing bars
    :param int n_bar: Number of steel reinforcing bars
    :param int n_circle: Number of points discretising the steel reinforcing bars
    :param float cover: Side and bottom cover to the steel reinforcing bars
    :param Optional[sectionproperties.pre.pre.Material] conc_mat: Material to associate with
        the concrete
    :param Optional[sectionproperties.pre.pre.Material] steel_mat: Material to associate with
        the steel

    :raises ValueErorr: If the number of bars is not greater than or equal to 2

    The following example creates a 600D x 300W concrete beam with 3N20 steel reinforcing bars and
    30 mm cover::

        from sectionproperties.pre.library.concrete_sections import concrete_rectangular_section
        from sectionproperties.pre.pre import Material

        concrete = Material(
            name='Concrete', elastic_modulus=30.1e3, poissons_ratio=0.2, yield_strength=32, color='lightgrey'
        )
        steel = Material(
            name='Steel', elastic_modulus=200e3, poissons_ratio=0.3, yield_strength=500, color='grey'
        )

        geometry = concrete_rectangular_section(
            b=300, d=600, dia=20, n_bar=3, n_circle=24, cover=30, conc_mat=concrete, steel_mat=steel
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

    if n_bar < 2:
        raise ValueError("Please provide 2 or more steel reinforcing bars.")

    geom = primitive_sections.rectangular_section(b=b, d=d, material=conc_mat)

    x_i = cover + dia / 2
    spacing = (b - 2 * cover - dia) / (n_bar - 1)

    for i in range(n_bar):
        bar = primitive_sections.circular_section(d=dia, n=n_circle, material=steel_mat)
        geom += bar.shift_section(x_offset=x_i + spacing * i, y_offset=cover + dia / 2)

    return geom


def concrete_tee_section(
    b: float,
    d: float,
    b_f: float,
    d_f: float,
    dia: float,
    n_bar: int,
    n_circle: int,
    cover: float,
    conc_mat: pre.Material = pre.DEFAULT_MATERIAL,
    steel_mat: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.CompoundGeometry:
    """Constructs a concrete tee section of width *b*, depth *d*, flange width *b_f* and flange
    depth *d_f* with *n_bar* steel bars of diameter *dia*, discretised with *n_circle* points with
    equal side and bottom *cover* to the steel.

    :param float b: Concrete section depth
    :param float d: Concrete section width
    :param float b_f: Concrete section flange depth
    :param float d_f: Concrete section flange width
    :param float dia: Diameter of the steel reinforcing bars
    :param int n_bar: Number of steel reinforcing bars
    :param int n_circle: Number of points discretising the steel reinforcing bars
    :param float cover: Side and bottom cover to the steel reinforcing bars
    :param Optional[sectionproperties.pre.pre.Material] conc_mat: Material to associate with
        the concrete
    :param Optional[sectionproperties.pre.pre.Material] steel_mat: Material to associate with
        the steel

    :raises ValueErorr: If the number of bars is not greater than or equal to 2

    The following example creates a 900D x 450W concrete beam with a 1200W x 250D flange, with 5N24
    steel reinforcing bars and 30 mm cover::

        from sectionproperties.pre.library.concrete_sections import concrete_tee_section
        from sectionproperties.pre.pre import Material

        concrete = Material(
            name='Concrete', elastic_modulus=30.1e3, poissons_ratio=0.2, yield_strength=32, color='lightgrey'
        )
        steel = Material(
            name='Steel', elastic_modulus=200e3, poissons_ratio=0.3, yield_strength=500, color='grey'
        )

        geometry = concrete_tee_section(
            b=450, d=900, b_f=1200, d_f=250, dia=24, n_bar=5, n_circle=24, cover=30,
            conc_mat=concrete, steel_mat=steel
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

    if n_bar < 2:
        raise ValueError("Please provide 2 or more steel reinforcing bars.")

    geom = primitive_sections.rectangular_section(b=b, d=d - d_f, material=conc_mat)
    flange = primitive_sections.rectangular_section(b=b_f, d=d_f, material=conc_mat)
    geom += flange.align_center(align_to=geom).align_to(other=geom, on="top")

    x_i = cover + dia / 2
    spacing = (b - 2 * cover - dia) / (n_bar - 1)

    for i in range(n_bar):
        bar = primitive_sections.circular_section(d=dia, n=n_circle, material=steel_mat)
        geom += bar.shift_section(x_offset=x_i + spacing * i, y_offset=cover + dia / 2)

    return geom


def concrete_circular_section(
    d: float,
    n: int,
    dia: float,
    n_bar: int,
    n_circle: int,
    cover: float,
    conc_mat: pre.Material = pre.DEFAULT_MATERIAL,
    steel_mat: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.CompoundGeometry:
    """Constructs a concrete circular section of diameter *d* discretised with *n* points, with
    *n_bar* steel bars of diameter *dia*, discretised with *n_circle* points with equal side and
    bottom *cover* to the steel.

    :param float d: Concrete diameter
    :param float n: Number of points discretising the concrete section
    :param float dia: Diameter of the steel reinforcing bars
    :param int n_bar: Number of steel reinforcing bars
    :param int n_circle: Number of points discretising the steel reinforcing bars
    :param float cover: Side and bottom cover to the steel reinforcing bars
    :param Optional[sectionproperties.pre.pre.Material] conc_mat: Material to associate with
        the concrete
    :param Optional[sectionproperties.pre.pre.Material] steel_mat: Material to associate with
        the steel

    :raises ValueErorr: If the number of bars is not greater than or equal to 2

    The following example creates a 450DIA concrete column with with 6N20 steel reinforcing bars
    and 45 mm cover::

        from sectionproperties.pre.library.concrete_sections import concrete_circular_section
        from sectionproperties.pre.pre import Material

        concrete = Material(
            name='Concrete', elastic_modulus=30.1e3, poissons_ratio=0.2, yield_strength=32, color='lightgrey'
        )
        steel = Material(
            name='Steel', elastic_modulus=200e3, poissons_ratio=0.3, yield_strength=500, color='grey'
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

    geom = primitive_sections.circular_section(d=d, n=n, material=conc_mat)

    r = d / 2 - cover - dia / 2
    d_theta = 2 * np.pi / n_bar

    for i in range(n_bar):
        bar = primitive_sections.circular_section(d=dia, n=n_circle, material=steel_mat)
        geom += bar.shift_section(
            x_offset=r * np.cos(i * d_theta), y_offset=r * np.sin(i * d_theta)
        )

    return geom
