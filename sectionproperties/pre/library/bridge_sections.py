import numpy as np
from shapely.geometry import Polygon
import sectionproperties.pre.geometry as geometry
import sectionproperties.pre.pre as pre


def super_t_girder_section(
    girder_type: int,
    girder_subtype: int = 2,
    w: float = 2100,
    t_w: float = None,
    t_f: float = 75,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a Super T Girder section to AS5100.5.

    :param int girder_type: Type of Super T (1 to 5)
    :param int girder_subtype: Era Super T (1: pre-2001, 2:contemporary)
    :param float w: Overall width of top flange
    :param float t_w: Web thickness of the Super-T section (defaults to those of AS5100.5 Tb D3(B))
    :param float t_f: Thickness of top flange (VIC (default) = 75 mm; NSW = 90 mm)
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    The following example creates a T5 Super-T section with a 180 mm overlay slab and assigns the
    different material properties::

        import sectionproperties.pre.library.bridge_sections as bridge_sections
        import sectionproperties.pre.library.standard_sections as standard_sections
        from sectionproperties.pre.pre import Material
        from sectionproperties.analysis.section import Section

        Dslab, w, t_f = 180, 2100, 75

        precast = Material(
            name="65 MPa",
            elastic_modulus=37.4e3,
            poissons_ratio=0.2,
            yield_strength=65,
            density=2.4e-6,
            color="grey",
        )
        insitu = Material(
            name="40 MPa",
            elastic_modulus=32.8e3,
            poissons_ratio=0.2,
            yield_strength=40,
            density=2.4e-6,
            color="lightgrey",
        )

        super_t = bridge_sections.super_t_girder_section(girder_type=5, w=w, material=precast)
        slab = standard_sections.rectangular_section(
            d=Dslab, b=w, material=insitu
        ).shift_section(-w / 2, t_f)

        geom = super_t + slab
        geom.plot_geometry()
        geom.create_mesh(mesh_sizes=[500])

        sec = Section(geom)
        sec.plot_mesh()

        sec.calculate_geometric_properties()
        sec.calculate_warping_properties()
        sec.display_results(fmt=".3f")


    Note that the properties are reported as ``modulus weighted`` properties (e.g. E.A) and can
    be normalized to the reference material by dividing by that elastic modulus::

        A_65 = section.get_ea() / precast.elastic_modulus

    The reported section centroids are already weighted.

    ..  figure:: ../images/sections/super_tee.png
        :align: center
        :scale: 40 %

        Super Tee Girder.
    """

    if girder_type < 1 or girder_type > 5:
        msg = "Super-T Girder Type must be between 1 and 5"
        raise Exception(msg)

    if girder_subtype not in [1, 2]:
        msg = "Only 2 subtypes: pre- and post-2001"
        raise Exception(msg)

    if girder_type == 5 and girder_subtype == 1:
        msg = "Only girders T1 to T4 before 2001"
        raise Exception(msg)

    # Some super-t constants, refer to AS5100.5, Figure D1(B)
    d_fillet, b_fillet = 75, 100
    d_recess, b_recess = 25, 25
    d_chamfer, b_chamfer = 13, 13
    web_slope = 10.556
    flg_slope = 5.347
    w_nom = 1027
    if girder_subtype == 1:
        flg_slope = 5.0
        w_nom = 920

    # Dims for the specific girder type
    d, t_b, t_w_nom = get_super_t_girder_dims(girder_type)

    # Overriding default web thickness?
    if t_w is None:
        t_w = t_w_nom

    # Origin is middle at level of bottom of top flange
    # Some geometrics of the slope
    web_hyp = np.sqrt(1 + web_slope ** 2)
    web_horiz = t_w * web_slope / web_hyp
    x_fillet = w_nom / 2 - d_fillet * 1 / web_hyp

    # initialise points variable
    points = []

    # Right recess
    pt_b = w_nom / 2 - web_horiz + (t_f - d_recess) / web_slope
    pt_a = pt_b + b_recess
    points.append([pt_b, t_f - d_recess])
    points.append([pt_a, t_f - d_recess])
    points.append([pt_a, t_f])

    # Right flange
    points.append([w / 2, t_f])
    points.append([w / 2, d_chamfer])
    points.append([w / 2 - b_chamfer, 0])
    points.append([w_nom / 2 + b_fillet, 0])
    points.append([x_fillet, -d_fillet])

    # Bottom outer
    btm_corner = x_fillet - d / web_slope
    points.append([btm_corner + d_chamfer / web_slope, -(d - d_chamfer)])
    points.append([btm_corner - b_chamfer, -d])
    points.append([-btm_corner + b_chamfer, -d])
    points.append([-btm_corner - d_chamfer / web_slope, -(d - d_chamfer)])

    # Left flange
    points.append([-x_fillet, -d_fillet])
    points.append([-w_nom / 2 - b_fillet, 0])
    points.append([-w / 2 + b_chamfer, 0])
    points.append([-w / 2, d_chamfer])
    points.append([-w / 2, t_f])

    # Left Recess
    points.append([-pt_a, t_f])
    points.append([-pt_a, t_f - d_recess])
    points.append([-pt_b, t_f - d_recess])

    # Bottom inner - find intersection point
    y10 = t_f - d_recess
    y20 = -(d - t_b)
    y_inner = (1 / (1 / web_slope - flg_slope)) * (
        y10 / web_slope - y20 * flg_slope - pt_b
    )
    x_inner = flg_slope * (y20 - y_inner)
    points.append([x_inner, y_inner])
    points.append([0, y20])
    points.append([-x_inner, y_inner])

    return geometry.Geometry(Polygon(points), material)


def i_girder_section(
    girder_type: int, material: pre.Material = pre.DEFAULT_MATERIAL
) -> geometry.Geometry:
    """Constructs a precast I girder section to AS5100.5.

    :param int girder_type: Type of I Girder (1 to 4)
    :param Optional[sectionproperties.pre.pre.Material]: Material to associate with this geometry

    As an example, replicate the table shown in AS5100.5 Fig. D1(A)::

        import pandas as pd
        import sectionproperties.pre.library.bridge_sections as bridge_sections
        from sectionproperties.analysis.section import Section

        df = pd.DataFrame(columns=["Ag", "Zt", "Zb", "I", "dy", "th"])

        for i in range(4):
            geom = bridge_sections.i_girder_section(girder_type=i + 1)
            dims = bridge_sections.get_i_girder_dims(girder_type=i + 1)
            d = sum(dims[-5:])
            geom.create_mesh(mesh_sizes=[200])
            geom.plot_geometry()
            sec = Section(geom)
            sec.plot_mesh()
            sec.calculate_geometric_properties()
            sec.calculate_warping_properties()

            A = sec.get_area()
            th = A / (sec.get_perimeter() / 2)

            df.loc[i] = [
                A,
                *(sec.get_z()[:2]),
                sec.get_ic()[0],
                d + sec.get_c()[1],
                th,
            ]

        print(df)

    Note that the section depth is obtained by summing the heights from the section dictionary in
    ``get_i_girder_dims()``.

    ..  figure:: ../images/sections/i_girder.png
        :align: center
        :scale: 40 %

        I Girder.
    """

    if girder_type < 1 or girder_type > 4:
        msg = "I Girder Type must be between 1 and 4"
        raise Exception(msg)

    b_tf, b_bf, b_w, h_tf, h_ts, h_w, h_bs, h_bf = get_i_girder_dims(girder_type)

    # Some section constants
    d = sum([h_tf, h_ts, h_w, h_bs, h_bf])
    inset_tf = (b_tf - b_w) / 2
    inset_bf = (b_bf - b_w) / 2

    # initialise points variable
    points = []

    # Origin at centre top; clockwise from top right corner
    points.append([b_tf / 2, 0])
    points.append([b_tf / 2, -h_tf])
    points.append([b_tf / 2 - inset_tf, -h_tf - h_ts])
    points.append([b_tf / 2 - inset_tf, -h_tf - h_ts - h_w])
    points.append([b_bf / 2, -d + h_bf])
    points.append([b_bf / 2, -d])
    points.append([-b_bf / 2, -d])
    points.append([-b_bf / 2, -d + h_bf])
    points.append([-b_tf / 2 + inset_tf, -h_tf - h_ts - h_w])
    points.append([-b_tf / 2 + inset_tf, -h_tf - h_ts])
    points.append([-b_tf / 2, -h_tf])
    points.append([-b_tf / 2, 0])

    return geometry.Geometry(Polygon(points), material)


def get_super_t_girder_dims(girder_type):
    """Returns a dictionary of Super-T dimensions, refer to AS5100.5, Appendix D

    :param int girder_type: Type of Super T (1 to 5)
    """

    girder_dims = {
        "T1": {"d": 675, "t_b": 240, "t_w": 100},
        "T2": {"d": 925, "t_b": 240, "t_w": 100},
        "T3": {"d": 1125, "t_b": 260, "t_w": 100},
        "T4": {"d": 1425, "t_b": 260, "t_w": 100},
        "T5": {"d": 1725, "t_b": 325, "t_w": 120},
    }

    key = f"T{girder_type}"

    # Rather delicious code that assigns the values to the keys as variables
    for key, val in girder_dims.items():
        exec(key + "=val")

    d = girder_dims[key]["d"]
    t_b = girder_dims[key]["t_b"]
    t_w = girder_dims[key]["t_w"]

    return d, t_b, t_w


def get_i_girder_dims(girder_type):
    """Returns a dictionary of I girder dimensions, refer to AS5100.5, Appendix D

    :param int girder_type: Type of I Girder (1 to 4)
    """

    girder_dims = {
        "T1": {
            "b_tf": 200,
            "b_bf": 300,
            "b_w": 120,
            "h_tf": 100,
            "h_ts": 40,
            "h_w": 420,
            "h_bs": 90,
            "h_bf": 100,
        },
        "T2": {
            "b_tf": 350,
            "b_bf": 450,
            "b_w": 150,
            "h_tf": 100,
            "h_ts": 100,
            "h_w": 450,
            "h_bs": 150,
            "h_bf": 100,
        },
        "T3": {
            "b_tf": 450,
            "b_bf": 500,
            "b_w": 150,
            "h_tf": 130,
            "h_ts": 150,
            "h_w": 545,
            "h_bs": 175,
            "h_bf": 150,
        },
        "T4": {
            "b_tf": 500,
            "b_bf": 650,
            "b_w": 150,
            "h_tf": 150,
            "h_ts": 175,
            "h_w": 650,
            "h_bs": 250,
            "h_bf": 175,
        },
    }

    key = f"T{girder_type}"

    b_tf = girder_dims[key]["b_tf"]
    b_bf = girder_dims[key]["b_bf"]
    b_w = girder_dims[key]["b_w"]
    h_tf = girder_dims[key]["h_tf"]
    h_ts = girder_dims[key]["h_ts"]
    h_w = girder_dims[key]["h_w"]
    h_bs = girder_dims[key]["h_bs"]
    h_bf = girder_dims[key]["h_bf"]

    return b_tf, b_bf, b_w, h_tf, h_ts, h_w, h_bs, h_bf
