"""Bridge sections library."""

from __future__ import annotations

import numpy as np
from shapely import Polygon

import sectionproperties.pre.geometry as geometry
import sectionproperties.pre.pre as pre


def super_t_girder_section(
    girder_type: int,
    girder_subtype: int = 2,
    w: float = 2100.0,
    t_w: float | None = None,
    t_f: float = 75.0,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a super T girder section to AS5100.5.

    Args:
        girder_type: Type of super T (1 to 5)
        girder_subtype: Era super T (1: pre-2001, 2: contemporary). Defaults to ``2``.
        w: Overall width of top flange. Defaults to ``2100.0``.
        t_w: Web thickness of the Super-T section, if ``None`` defaults to those of
            AS5100.5 Tb D3(B). Defaults to ``None``.
        t_f: Thickness of top flange (VIC = 75 mm; NSW = 90 mm). Defaults to ``75.0``.
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        ValueError: ``girder_type`` or ``girder_subtype`` are invalid

    Returns:
        Super T girder section geometry

    Example:
        The following example creates a T5 super T section with a 180 mm overlay slab:

        .. plot::
            :include-source: True
            :caption: Rectangular section geometry

            from sectionproperties.pre.library import super_t_girder_section
            from sectionproperties.pre.library import rectangular_section

            super_t = super_t_girder_section(girder_type=5, w=2100)
            slab = rectangular_section(d=180, b=2100).shift_section(
                x_offset=-1050, y_offset=75
            )
            (super_t + slab).plot_geometry()
    """
    if girder_type < 1 or girder_type > 5:
        msg = "Super-T Girder Type must be between 1 and 5"
        raise ValueError(msg)

    if girder_subtype not in (1, 2):
        msg = "Only 2 subtypes: pre- and post-2001"
        raise ValueError(msg)

    if girder_type == 5 and girder_subtype == 1:
        msg = "Only girders T1 to T4 before 2001"
        raise ValueError(msg)

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
    d, t_b, t_w_nom = get_super_t_girder_dims(girder_type=girder_type)

    # Overriding default web thickness?
    if t_w is None:
        t_w = t_w_nom

    # Origin is middle at level of bottom of top flange
    # Some geometrics of the slope
    web_hyp = np.sqrt(1 + web_slope**2)
    web_horiz = t_w * web_slope / web_hyp
    x_fillet = w_nom / 2 - d_fillet * 1 / web_hyp

    # initialise points variable
    points: list[tuple[float, float]] = []

    # Right recess
    pt_b = w_nom / 2 - web_horiz + (t_f - d_recess) / web_slope
    pt_a = pt_b + b_recess
    points.append((pt_b, t_f - d_recess))
    points.append((pt_a, t_f - d_recess))
    points.append((pt_a, t_f))

    # Right flange
    points.append((w / 2, t_f))
    points.append((w / 2, d_chamfer))
    points.append((w / 2 - b_chamfer, 0))
    points.append((w_nom / 2 + b_fillet, 0))
    points.append((x_fillet, -d_fillet))

    # Bottom outer
    btm_corner = x_fillet - d / web_slope
    points.append((btm_corner + d_chamfer / web_slope, -(d - d_chamfer)))
    points.append((btm_corner - b_chamfer, -d))
    points.append((-btm_corner + b_chamfer, -d))
    points.append((-btm_corner - d_chamfer / web_slope, -(d - d_chamfer)))

    # Left flange
    points.append((-x_fillet, -d_fillet))
    points.append((-w_nom / 2 - b_fillet, 0))
    points.append((-w / 2 + b_chamfer, 0))
    points.append((-w / 2, d_chamfer))
    points.append((-w / 2, t_f))

    # Left Recess
    points.append((-pt_a, t_f))
    points.append((-pt_a, t_f - d_recess))
    points.append((-pt_b, t_f - d_recess))

    # Bottom inner - find intersection point
    y10 = t_f - d_recess
    y20 = -(d - t_b)
    y_inner = (1 / (1 / web_slope - flg_slope)) * (
        y10 / web_slope - y20 * flg_slope - pt_b
    )
    x_inner = flg_slope * (y20 - y_inner)
    points.append((x_inner, y_inner))
    points.append((0, y20))
    points.append((-x_inner, y_inner))

    return geometry.Geometry(geom=Polygon(points), material=material)


def i_girder_section(
    girder_type: int,
    material: pre.Material = pre.DEFAULT_MATERIAL,
) -> geometry.Geometry:
    """Constructs a precast I girder section to AS5100.5.

    Args:
        girder_type: Type of I Girder (1 to 4)
        material: Material to associate with this geometry. Defaults to
            ``pre.DEFAULT_MATERIAL``.

    Raises:
        ValueError: ``girder_type`` is invalid

    Returns:
        I girder section geometry

    Example:
        The following example creates a type 1 precast I girder section with a 180 mm
        overlay slab:

        .. plot::
            :include-source: True
            :caption: Rectangular section geometry

            from sectionproperties.pre.library import i_girder_section
            from sectionproperties.pre.library import rectangular_section

            i_girder = i_girder_section(girder_type=1)
            slab = rectangular_section(d=180, b=1200).shift_section(x_offset=-1050)
            (i_girder + slab).plot_geometry()
    """
    if girder_type < 1 or girder_type > 4:
        msg = "I Girder Type must be between 1 and 4"
        raise ValueError(msg)

    b_tf, b_bf, b_w, h_tf, h_ts, h_w, h_bs, h_bf = get_i_girder_dims(girder_type)

    # Some section constants
    d = sum([h_tf, h_ts, h_w, h_bs, h_bf])
    inset_tf = (b_tf - b_w) / 2

    # initialise points variable
    points: list[tuple[float, float]] = []

    # Origin at centre top; clockwise from top right corner
    points.append((b_tf / 2, 0))
    points.append((b_tf / 2, -h_tf))
    points.append((b_tf / 2 - inset_tf, -h_tf - h_ts))
    points.append((b_tf / 2 - inset_tf, -h_tf - h_ts - h_w))
    points.append((b_bf / 2, -d + h_bf))
    points.append((b_bf / 2, -d))
    points.append((-b_bf / 2, -d))
    points.append((-b_bf / 2, -d + h_bf))
    points.append((-b_tf / 2 + inset_tf, -h_tf - h_ts - h_w))
    points.append((-b_tf / 2 + inset_tf, -h_tf - h_ts))
    points.append((-b_tf / 2, -h_tf))
    points.append((-b_tf / 2, 0))

    return geometry.Geometry(geom=Polygon(points), material=material)


def get_super_t_girder_dims(
    girder_type: int,
) -> tuple[int, int, int]:
    """Returns critical dimensions of super t girders, refer to AS5100.5, Appendix D.

    Args:
        girder_type: Type of Super T (1 to 5)

    Returns:
        Girder depth, base thickness and web thickness
    """
    girder_dims = {
        "T1": {"d": 675, "t_b": 240, "t_w": 100},
        "T2": {"d": 925, "t_b": 240, "t_w": 100},
        "T3": {"d": 1125, "t_b": 260, "t_w": 100},
        "T4": {"d": 1425, "t_b": 260, "t_w": 100},
        "T5": {"d": 1725, "t_b": 325, "t_w": 120},
    }

    key = f"T{girder_type}"

    d = girder_dims[key]["d"]
    t_b = girder_dims[key]["t_b"]
    t_w = girder_dims[key]["t_w"]

    return d, t_b, t_w


def get_i_girder_dims(
    girder_type: int,
) -> tuple[int, int, int, int, int, int, int, int]:
    """Returns critical I girder dimensions, refer to AS5100.5, Appendix D.

    Args:
        girder_type: Type of I Girder (1 to 4)

    Returns:
        Critical I girder dimensions
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
