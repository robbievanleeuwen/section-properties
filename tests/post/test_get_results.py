"""Validation of cross-section property get methods."""

from __future__ import annotations

import numpy as np
import pytest

from sectionproperties.analysis import Section
from sectionproperties.pre import Material
from sectionproperties.pre.library import rectangular_section


# sections setup
@pytest.fixture
def sec_no_mat_and_mat() -> tuple[Section, Section, Material]:
    """Creates two unit squares, one with materials, another without.

    Returns:
        Section objects and dummy material
    """
    geom_no_mat = rectangular_section(d=1, b=1)
    dummy_mat = Material(
        name="test",
        elastic_modulus=5,
        poissons_ratio=0,
        yield_strength=3,
        density=2,
        color="w",
    )
    geom_mat = rectangular_section(d=1, b=1, material=dummy_mat)
    geom_no_mat.create_mesh(mesh_sizes=0)
    geom_mat.create_mesh(mesh_sizes=0)

    return Section(geom_no_mat), Section(geom_mat), dummy_mat


def test_is_composite(sec_no_mat_and_mat):
    """Check whether section is composite or not."""
    rect_no_mat, rect_mat, _ = sec_no_mat_and_mat
    assert not rect_no_mat.is_composite()
    assert rect_mat.is_composite()

    # check built up geometric only section
    rect1 = rectangular_section(d=2, b=2)
    rect2 = rectangular_section(d=4, b=5).shift_section(x_offset=2)
    geom = rect1 + rect2
    geom.create_mesh(mesh_sizes=[0])
    sec = Section(geom)
    assert not sec.is_composite()

    # check section with hole geometric only section
    rect_out = rectangular_section(d=10, b=10)
    rect_in = rectangular_section(d=5, b=5).shift_section(x_offset=5, y_offset=5)
    geom = rect_out - rect_in
    geom.create_mesh(mesh_sizes=[0])
    sec = Section(geom)
    assert not sec.is_composite()


def test_get_e_ref(sec_no_mat_and_mat):
    """Check get_e_ref results."""
    rect_no_mat, _, dummy_mat = sec_no_mat_and_mat
    test_mat = Material(
        name="test",
        elastic_modulus=10,
        poissons_ratio=0,
        yield_strength=1,
        density=1,
        color="w",
    )

    assert rect_no_mat.get_e_ref(1) == 1
    assert rect_no_mat.get_e_ref(100) == 100
    assert rect_no_mat.get_e_ref(dummy_mat) == 5
    assert rect_no_mat.get_e_ref(test_mat) == 10


def test_no_analysis(sec_no_mat_and_mat):
    """Check errors when no analysis has been conducted or incorrect analysis."""
    rect_no_mat, rect_mat, _ = sec_no_mat_and_mat
    # check area
    with pytest.raises(RuntimeError):
        rect_no_mat.get_area()

    # check perimeter
    with pytest.raises(RuntimeError):
        rect_no_mat.get_perimeter()

    # check mass
    with pytest.raises(RuntimeError):
        rect_no_mat.get_mass()

    with pytest.raises(RuntimeError):
        rect_mat.get_mass()

    # check ea
    with pytest.raises(RuntimeError):
        rect_no_mat.get_mass()

    with pytest.raises(RuntimeError):
        rect_mat.get_mass()

    with pytest.raises(RuntimeError):
        rect_mat.get_mass()

    # check q
    with pytest.raises(RuntimeError):
        rect_no_mat.get_q()

    with pytest.raises(RuntimeError):
        rect_mat.get_q()

    # check eq
    with pytest.raises(RuntimeError):
        rect_no_mat.get_eq()

    with pytest.raises(RuntimeError):
        rect_mat.get_eq()

    # check ig
    with pytest.raises(RuntimeError):
        rect_no_mat.get_ig()

    with pytest.raises(RuntimeError):
        rect_mat.get_ig()

    # check eig
    with pytest.raises(RuntimeError):
        rect_no_mat.get_eig()

    with pytest.raises(RuntimeError):
        rect_mat.get_eig()

    # check c
    with pytest.raises(RuntimeError):
        rect_no_mat.get_c()

    # check ic
    with pytest.raises(RuntimeError):
        rect_no_mat.get_ic()

    with pytest.raises(RuntimeError):
        rect_mat.get_ic()

    # check eic
    with pytest.raises(RuntimeError):
        rect_no_mat.get_eic()

    with pytest.raises(RuntimeError):
        rect_mat.get_eic()

    # check z
    with pytest.raises(RuntimeError):
        rect_no_mat.get_z()

    with pytest.raises(RuntimeError):
        rect_mat.get_z()

    # check ez
    with pytest.raises(RuntimeError):
        rect_no_mat.get_ez()

    with pytest.raises(RuntimeError):
        rect_mat.get_ez()

    # check rc
    with pytest.raises(RuntimeError):
        rect_no_mat.get_rc()

    # check ip
    with pytest.raises(RuntimeError):
        rect_no_mat.get_ip()

    with pytest.raises(RuntimeError):
        rect_mat.get_ip()

    # check eip
    with pytest.raises(RuntimeError):
        rect_no_mat.get_eip()

    with pytest.raises(RuntimeError):
        rect_mat.get_eip()

    # check phi
    with pytest.raises(RuntimeError):
        rect_no_mat.get_phi()

    # check zp
    with pytest.raises(RuntimeError):
        rect_no_mat.get_zp()

    with pytest.raises(RuntimeError):
        rect_mat.get_zp()

    # check ezp
    with pytest.raises(RuntimeError):
        rect_no_mat.get_ezp()

    with pytest.raises(RuntimeError):
        rect_mat.get_ezp()

    # check rp
    with pytest.raises(RuntimeError):
        rect_no_mat.get_rc()

    # check nu_eff
    with pytest.raises(RuntimeError):
        rect_no_mat.get_nu_eff()

    with pytest.raises(RuntimeError):
        rect_mat.get_nu_eff()

    # check e_eff
    with pytest.raises(RuntimeError):
        rect_no_mat.get_e_eff()

    with pytest.raises(RuntimeError):
        rect_mat.get_e_eff()

    # check g_eff
    with pytest.raises(RuntimeError):
        rect_no_mat.get_g_eff()

    with pytest.raises(RuntimeError):
        rect_mat.get_g_eff()


def test_get_geometric_only(sec_no_mat_and_mat):
    """Check errors and results when a geometric analysis has been conducted."""
    rect_no_mat, rect_mat, dummy_mat = sec_no_mat_and_mat
    rect_no_mat.calculate_geometric_properties()
    rect_mat.calculate_geometric_properties()

    # check area
    assert rect_no_mat.get_area() == pytest.approx(1.0)
    assert rect_mat.get_area() == pytest.approx(1.0)

    # check perimeter
    assert rect_no_mat.get_perimeter() == pytest.approx(4.0)
    assert rect_mat.get_perimeter() == pytest.approx(4.0)

    # check mass
    with pytest.raises(RuntimeError):
        rect_no_mat.get_mass()

    assert rect_mat.get_mass() == pytest.approx(2.0)

    # check ea
    with pytest.raises(RuntimeError):
        rect_no_mat.get_ea()

    assert rect_mat.get_ea() == pytest.approx(5.0)
    assert rect_mat.get_ea(e_ref=2) == pytest.approx(2.5)
    assert rect_mat.get_ea(e_ref=dummy_mat) == pytest.approx(1.0)

    # check q
    qx, qy = rect_no_mat.get_q()
    assert qx == pytest.approx(0.5)
    assert qy == pytest.approx(0.5)

    with pytest.raises(RuntimeError):
        rect_mat.get_q()

    # check eq
    with pytest.raises(RuntimeError):
        rect_no_mat.get_eq()

    eqx, eqy = rect_mat.get_eq()
    assert eqx == pytest.approx(2.5)
    assert eqy == pytest.approx(2.5)
    eqx, eqy = rect_mat.get_eq(e_ref=2)
    assert eqx == pytest.approx(1.25)
    assert eqy == pytest.approx(1.25)
    eqx, eqy = rect_mat.get_eq(e_ref=dummy_mat)
    assert eqx == pytest.approx(0.5)
    assert eqy == pytest.approx(0.5)

    # check ig
    igxx, igyy, igxy = rect_no_mat.get_ig()
    assert igxx == pytest.approx(1 / 3.0)
    assert igyy == pytest.approx(1 / 3.0)
    assert igxy == pytest.approx(0.25)

    with pytest.raises(RuntimeError):
        rect_mat.get_ig()

    # check eig
    with pytest.raises(RuntimeError):
        rect_no_mat.get_eig()

    eigxx, eigyy, eigxy = rect_mat.get_eig()
    assert eigxx == pytest.approx(5.0 / 3)
    assert eigyy == pytest.approx(5.0 / 3)
    assert eigxy == pytest.approx(1.25)
    eigxx, eigyy, eigxy = rect_mat.get_eig(e_ref=2)
    assert eigxx == pytest.approx(5.0 / 6)
    assert eigyy == pytest.approx(5.0 / 6)
    assert eigxy == pytest.approx(0.625)
    eigxx, eigyy, eigxy = rect_mat.get_eig(e_ref=dummy_mat)
    assert eigxx == pytest.approx(1 / 3.0)
    assert eigyy == pytest.approx(1 / 3.0)
    assert eigxy == pytest.approx(0.25)

    # check c
    cx, cy = rect_no_mat.get_c()
    assert cx == pytest.approx(1 / 2.0)
    assert cy == pytest.approx(1 / 2.0)
    cx, cy = rect_mat.get_c()
    assert cx == pytest.approx(1 / 2.0)
    assert cy == pytest.approx(1 / 2.0)

    # check ic
    icxx, icyy, icxy = rect_no_mat.get_ic()
    assert icxx == pytest.approx(1 / 12.0)
    assert icyy == pytest.approx(1 / 12.0)
    assert icxy == pytest.approx(0.0)

    with pytest.raises(RuntimeError):
        rect_mat.get_ic()

    # check eic
    with pytest.raises(RuntimeError):
        rect_no_mat.get_eic()

    eicxx, eicyy, eicxy = rect_mat.get_eic()
    assert eicxx == pytest.approx(5.0 / 12)
    assert eicyy == pytest.approx(5.0 / 12)
    assert eicxy == pytest.approx(0.0)
    eicxx, eicyy, eicxy = rect_mat.get_eic(e_ref=2)
    assert eicxx == pytest.approx(5.0 / 24)
    assert eicyy == pytest.approx(5.0 / 24)
    assert eicxy == pytest.approx(0.0)
    eicxx, eicyy, eicxy = rect_mat.get_eic(e_ref=dummy_mat)
    assert eicxx == pytest.approx(1 / 12.0)
    assert eicyy == pytest.approx(1 / 12.0)
    assert eicxy == pytest.approx(0.0)

    # check z
    zxx_plus, zxx_minus, zyy_plus, zyy_minus = rect_no_mat.get_z()
    assert zxx_plus == pytest.approx(1 / 6.0)
    assert zxx_minus == pytest.approx(1 / 6.0)
    assert zyy_plus == pytest.approx(1 / 6.0)
    assert zyy_minus == pytest.approx(1 / 6.0)

    with pytest.raises(RuntimeError):
        rect_mat.get_z()

    # check ez
    with pytest.raises(RuntimeError):
        rect_no_mat.get_ez()

    ezxx_plus, ezxx_minus, ezyy_plus, ezyy_minus = rect_mat.get_ez()
    assert ezxx_plus == pytest.approx(5 / 6.0)
    assert ezxx_minus == pytest.approx(5 / 6.0)
    assert ezyy_plus == pytest.approx(5 / 6.0)
    assert ezyy_minus == pytest.approx(5 / 6.0)
    ezxx_plus, ezxx_minus, ezyy_plus, ezyy_minus = rect_mat.get_ez(e_ref=2)
    assert ezxx_plus == pytest.approx(5 / 12.0)
    assert ezxx_minus == pytest.approx(5 / 12.0)
    assert ezyy_plus == pytest.approx(5 / 12.0)
    assert ezyy_minus == pytest.approx(5 / 12.0)
    ezxx_plus, ezxx_minus, ezyy_plus, ezyy_minus = rect_mat.get_ez(e_ref=dummy_mat)
    assert ezxx_plus == pytest.approx(1 / 6.0)
    assert ezxx_minus == pytest.approx(1 / 6.0)
    assert ezyy_plus == pytest.approx(1 / 6.0)
    assert ezyy_minus == pytest.approx(1 / 6.0)

    # check rc
    rcx, rcy = rect_no_mat.get_rc()
    assert rcx == pytest.approx(np.sqrt(1 / 12.0))
    assert rcy == pytest.approx(np.sqrt(1 / 12.0))
    rcx, rcy = rect_mat.get_rc()
    assert rcx == pytest.approx(np.sqrt(1 / 12.0))
    assert rcy == pytest.approx(np.sqrt(1 / 12.0))

    # check ip
    i11, i22 = rect_no_mat.get_ip()
    assert i11 == pytest.approx(1 / 12.0)
    assert i22 == pytest.approx(1 / 12.0)

    with pytest.raises(RuntimeError):
        rect_mat.get_ip()

    # check eip
    with pytest.raises(RuntimeError):
        rect_no_mat.get_eip()

    ei11, ei22 = rect_mat.get_eip()
    assert ei11 == pytest.approx(5.0 / 12)
    assert ei22 == pytest.approx(5.0 / 12)
    ei11, ei22 = rect_mat.get_eip(e_ref=2)
    assert ei11 == pytest.approx(5.0 / 24)
    assert ei22 == pytest.approx(5.0 / 24)
    ei11, ei22 = rect_mat.get_eip(e_ref=dummy_mat)
    assert ei11 == pytest.approx(1 / 12.0)
    assert ei22 == pytest.approx(1 / 12.0)

    # check phi
    phi = rect_no_mat.get_phi()
    assert phi == pytest.approx(0.0)
    phi = rect_mat.get_phi()
    assert phi == pytest.approx(0.0)

    # check zp
    z11_plus, z11_minus, z22_plus, z22_minus = rect_no_mat.get_zp()
    assert z11_plus == pytest.approx(1 / 6.0)
    assert z11_minus == pytest.approx(1 / 6.0)
    assert z22_plus == pytest.approx(1 / 6.0)
    assert z22_minus == pytest.approx(1 / 6.0)

    with pytest.raises(RuntimeError):
        rect_mat.get_z()

    # check ezp
    with pytest.raises(RuntimeError):
        rect_no_mat.get_ezp()

    ez11_plus, ez11_minus, ez22_plus, ez22_minus = rect_mat.get_ezp()
    assert ez11_plus == pytest.approx(5 / 6.0)
    assert ez11_minus == pytest.approx(5 / 6.0)
    assert ez22_plus == pytest.approx(5 / 6.0)
    assert ez22_minus == pytest.approx(5 / 6.0)
    ez11_plus, ez11_minus, ez22_plus, ez22_minus = rect_mat.get_ezp(e_ref=2)
    assert ez11_plus == pytest.approx(5 / 12.0)
    assert ez11_minus == pytest.approx(5 / 12.0)
    assert ez22_plus == pytest.approx(5 / 12.0)
    assert ez22_minus == pytest.approx(5 / 12.0)
    ez11_plus, ez11_minus, ez22_plus, ez22_minus = rect_mat.get_ezp(e_ref=dummy_mat)
    assert ez11_plus == pytest.approx(1 / 6.0)
    assert ez11_minus == pytest.approx(1 / 6.0)
    assert ez22_plus == pytest.approx(1 / 6.0)
    assert ez22_minus == pytest.approx(1 / 6.0)

    # check rp
    r11, r22 = rect_no_mat.get_rp()
    assert r11 == pytest.approx(np.sqrt(1 / 12.0))
    assert r22 == pytest.approx(np.sqrt(1 / 12.0))
    r11, r22 = rect_mat.get_rp()
    assert r11 == pytest.approx(np.sqrt(1 / 12.0))
    assert r22 == pytest.approx(np.sqrt(1 / 12.0))

    # check nu_eff
    with pytest.raises(RuntimeError):
        rect_no_mat.get_nu_eff()

    assert rect_mat.get_nu_eff() == pytest.approx(0.0)

    # check e_eff
    with pytest.raises(RuntimeError):
        rect_no_mat.get_e_eff()

    assert rect_mat.get_e_eff() == pytest.approx(5.0)

    # check g_eff
    with pytest.raises(RuntimeError):
        rect_no_mat.get_g_eff()

    assert rect_mat.get_g_eff() == pytest.approx(2.5)

    # now check errors from performing no warping analysis
    # check j
    with pytest.raises(RuntimeError):
        rect_no_mat.get_j()

    with pytest.raises(RuntimeError):
        rect_mat.get_j()

    # check ej
    with pytest.raises(RuntimeError):
        rect_no_mat.get_ej()

    with pytest.raises(RuntimeError):
        rect_mat.get_ej()

    # check sc
    with pytest.raises(RuntimeError):
        rect_no_mat.get_sc()

    with pytest.raises(RuntimeError):
        rect_mat.get_sc()

    # check sc_p
    with pytest.raises(RuntimeError):
        rect_no_mat.get_sc_p()

    with pytest.raises(RuntimeError):
        rect_mat.get_sc_p()

    # check sc_t
    with pytest.raises(RuntimeError):
        rect_no_mat.get_sc_t()

    with pytest.raises(RuntimeError):
        rect_mat.get_sc_t()

    # check gamma
    with pytest.raises(RuntimeError):
        rect_no_mat.get_gamma()

    with pytest.raises(RuntimeError):
        rect_mat.get_gamma()

    # check egamma
    with pytest.raises(RuntimeError):
        rect_no_mat.get_egamma()

    with pytest.raises(RuntimeError):
        rect_mat.get_egamma()

    # check as
    with pytest.raises(RuntimeError):
        rect_no_mat.get_as()

    with pytest.raises(RuntimeError):
        rect_mat.get_as()

    # check eas
    with pytest.raises(RuntimeError):
        rect_no_mat.get_eas()

    with pytest.raises(RuntimeError):
        rect_mat.get_eas()

    # check as_p
    with pytest.raises(RuntimeError):
        rect_no_mat.get_as_p()

    with pytest.raises(RuntimeError):
        rect_mat.get_as_p()

    # check eas_p
    with pytest.raises(RuntimeError):
        rect_no_mat.get_eas_p()

    with pytest.raises(RuntimeError):
        rect_mat.get_eas_p()

    # check beta
    with pytest.raises(RuntimeError):
        rect_no_mat.get_beta()

    with pytest.raises(RuntimeError):
        rect_mat.get_beta()

    # check beta_p
    with pytest.raises(RuntimeError):
        rect_no_mat.get_beta_p()

    with pytest.raises(RuntimeError):
        rect_mat.get_beta_p()

    # check pc
    with pytest.raises(RuntimeError):
        rect_no_mat.get_pc()

    with pytest.raises(RuntimeError):
        rect_mat.get_pc()

    # check pc_p
    with pytest.raises(RuntimeError):
        rect_no_mat.get_pc_p()

    with pytest.raises(RuntimeError):
        rect_mat.get_pc_p()

    # check s
    with pytest.raises(RuntimeError):
        rect_no_mat.get_s()

    with pytest.raises(RuntimeError):
        rect_mat.get_s()

    # check mp
    with pytest.raises(RuntimeError):
        rect_no_mat.get_mp()

    with pytest.raises(RuntimeError):
        rect_mat.get_mp()

    # check sp
    with pytest.raises(RuntimeError):
        rect_no_mat.get_sp()

    with pytest.raises(RuntimeError):
        rect_mat.get_sp()

    # check mp_p
    with pytest.raises(RuntimeError):
        rect_no_mat.get_mp_p()

    with pytest.raises(RuntimeError):
        rect_mat.get_mp_p()

    # check sf
    with pytest.raises(RuntimeError):
        rect_no_mat.get_sf()

    with pytest.raises(RuntimeError):
        rect_mat.get_sf()

    # check sf_p
    with pytest.raises(RuntimeError):
        rect_no_mat.get_sf_p()

    with pytest.raises(RuntimeError):
        rect_mat.get_sf_p()


def test_get_warping(sec_no_mat_and_mat):
    """Check errors and results when a warping analysis has been conducted."""
    rect_no_mat, rect_mat, dummy_mat = sec_no_mat_and_mat
    rect_no_mat.calculate_geometric_properties()
    rect_no_mat.calculate_warping_properties()
    rect_mat.calculate_geometric_properties()
    rect_mat.calculate_warping_properties()

    # check j
    assert rect_no_mat.get_j() == pytest.approx(1 / 6.0)

    with pytest.raises(RuntimeError):
        rect_mat.get_j()

    # check ej
    with pytest.raises(RuntimeError):
        rect_no_mat.get_ej()

    assert rect_mat.get_ej() == pytest.approx(5 / 6.0)
    assert rect_mat.get_ej(e_ref=2) == pytest.approx(5 / 12.0)
    assert rect_mat.get_ej(e_ref=dummy_mat) == pytest.approx(1 / 6.0)

    # check sc
    x_se, y_se = rect_no_mat.get_sc()
    assert x_se == pytest.approx(0.5)
    assert y_se == pytest.approx(0.5)
    x_se, y_se = rect_mat.get_sc()
    assert x_se == pytest.approx(0.5)
    assert y_se == pytest.approx(0.5)

    # check sc_p
    x11_se, y22_se = rect_no_mat.get_sc_p()
    assert x11_se == pytest.approx(0.0)
    assert y22_se == pytest.approx(0.0)
    x11_se, y22_se = rect_mat.get_sc_p()
    assert x11_se == pytest.approx(0.0)
    assert y22_se == pytest.approx(0.0)

    # check sc_t
    x_st, y_st = rect_no_mat.get_sc_t()
    assert x_st == pytest.approx(0.5)
    assert y_st == pytest.approx(0.5)
    x_st, y_st = rect_mat.get_sc_t()
    assert x_st == pytest.approx(0.5)
    assert y_st == pytest.approx(0.5)

    # check gamma
    assert rect_no_mat.get_gamma() == pytest.approx(0.0)

    with pytest.raises(RuntimeError):
        rect_mat.get_gamma()

    # check egamma
    with pytest.raises(RuntimeError):
        rect_no_mat.get_egamma()

    assert rect_mat.get_egamma() == pytest.approx(0.0)
    assert rect_mat.get_egamma(e_ref=2) == pytest.approx(0.0)
    assert rect_mat.get_egamma(e_ref=dummy_mat) == pytest.approx(0.0)

    # check as
    a_sx, a_sy = rect_no_mat.get_as()
    assert a_sx == pytest.approx(0.9523810)
    assert a_sy == pytest.approx(0.9523810)

    with pytest.raises(RuntimeError):
        rect_mat.get_as()

    # check eas
    with pytest.raises(RuntimeError):
        rect_no_mat.get_eas()

    ea_sx, ea_sy = rect_mat.get_eas()
    assert ea_sx == pytest.approx(0.9523810 * 5)
    assert ea_sy == pytest.approx(0.9523810 * 5)
    ea_sx, ea_sy = rect_mat.get_eas(e_ref=2)
    assert ea_sx == pytest.approx(0.9523810 * 2.5)
    assert ea_sy == pytest.approx(0.9523810 * 2.5)
    ea_sx, ea_sy = rect_mat.get_eas(e_ref=dummy_mat)
    assert ea_sx == pytest.approx(0.9523810)
    assert ea_sy == pytest.approx(0.9523810)

    # check as_p
    a_s11, a_s22 = rect_no_mat.get_as_p()
    assert a_s11 == pytest.approx(0.9523810)
    assert a_s22 == pytest.approx(0.9523810)

    with pytest.raises(RuntimeError):
        rect_mat.get_as_p()

    # check eas_p
    with pytest.raises(RuntimeError):
        rect_no_mat.get_eas_p()

    ea_s11, ea_s22 = rect_mat.get_eas_p()
    assert ea_s11 == pytest.approx(0.9523810 * 5)
    assert ea_s22 == pytest.approx(0.9523810 * 5)
    ea_s11, ea_s22 = rect_mat.get_eas_p(e_ref=2)
    assert ea_s11 == pytest.approx(0.9523810 * 2.5)
    assert ea_s22 == pytest.approx(0.9523810 * 2.5)
    ea_s11, ea_s22 = rect_mat.get_eas_p(e_ref=dummy_mat)
    assert ea_s11 == pytest.approx(0.9523810)
    assert ea_s22 == pytest.approx(0.9523810)

    # check beta
    beta_x_plus, beta_x_minus, beta_y_plus, beta_y_minus = rect_no_mat.get_beta()
    assert beta_x_plus == pytest.approx(0.0)
    assert beta_x_minus == pytest.approx(0.0)
    assert beta_y_plus == pytest.approx(0.0)
    assert beta_y_minus == pytest.approx(0.0)
    beta_x_plus, beta_x_minus, beta_y_plus, beta_y_minus = rect_mat.get_beta()
    assert beta_x_plus == pytest.approx(0.0)
    assert beta_x_minus == pytest.approx(0.0)
    assert beta_y_plus == pytest.approx(0.0)
    assert beta_y_minus == pytest.approx(0.0)

    # check beta_p
    beta_11_plus, beta_11_minus, beta_22_plus, beta_22_minus = rect_no_mat.get_beta_p()
    assert beta_11_plus == pytest.approx(0.0)
    assert beta_11_minus == pytest.approx(0.0)
    assert beta_22_plus == pytest.approx(0.0)
    assert beta_22_minus == pytest.approx(0.0)
    beta_11_plus, beta_11_minus, beta_22_plus, beta_22_minus = rect_mat.get_beta_p()
    assert beta_11_plus == pytest.approx(0.0)
    assert beta_11_minus == pytest.approx(0.0)
    assert beta_22_plus == pytest.approx(0.0)
    assert beta_22_minus == pytest.approx(0.0)


def test_get_plastic(sec_no_mat_and_mat):
    """Check errors and results when a plastic analysis has been conducted."""
    rect_no_mat, rect_mat, _ = sec_no_mat_and_mat
    rect_no_mat.calculate_geometric_properties()
    rect_mat.calculate_geometric_properties()
    rect_no_mat.calculate_plastic_properties()
    rect_mat.calculate_plastic_properties()

    # check pc
    x_pc, y_pc = rect_no_mat.get_pc()
    assert x_pc == pytest.approx(0.5)
    assert y_pc == pytest.approx(0.5)
    x_pc, y_pc = rect_mat.get_pc()
    assert x_pc == pytest.approx(0.5)
    assert y_pc == pytest.approx(0.5)

    # check pc_p
    x_pc, y_pc = rect_no_mat.get_pc_p()
    assert x_pc == pytest.approx(0.5)
    assert y_pc == pytest.approx(0.5)
    x_pc, y_pc = rect_mat.get_pc_p()
    assert x_pc == pytest.approx(0.5)
    assert y_pc == pytest.approx(0.5)

    # check s
    sxx, syy = rect_no_mat.get_s()
    assert sxx == pytest.approx(0.25)
    assert syy == pytest.approx(0.25)

    with pytest.raises(RuntimeError):
        rect_mat.get_s()

    # check mp
    with pytest.raises(RuntimeError):
        rect_no_mat.get_mp()

    mp_xx, mp_yy = rect_mat.get_mp()
    assert mp_xx == pytest.approx(0.75)
    assert mp_yy == pytest.approx(0.75)

    # check sp
    sxx, syy = rect_no_mat.get_sp()
    assert sxx == pytest.approx(0.25)
    assert syy == pytest.approx(0.25)

    with pytest.raises(RuntimeError):
        rect_mat.get_sp()

    # check mp_p
    with pytest.raises(RuntimeError):
        rect_no_mat.get_mp_p()

    mp_xx, mp_yy = rect_mat.get_mp_p()
    assert mp_xx == pytest.approx(0.75)
    assert mp_yy == pytest.approx(0.75)

    # check sf
    sf_xx_plus, sf_xx_minus, sf_yy_plus, sf_yy_minus = rect_no_mat.get_sf()
    assert sf_xx_plus == pytest.approx(1.5)
    assert sf_xx_minus == pytest.approx(1.5)
    assert sf_yy_plus == pytest.approx(1.5)
    assert sf_yy_minus == pytest.approx(1.5)

    with pytest.raises(RuntimeError):
        rect_mat.get_sf()

    # check sf_p
    sf_11_plus, sf_11_minus, sf_22_plus, sf_22_minus = rect_no_mat.get_sf_p()
    assert sf_11_plus == pytest.approx(1.5)
    assert sf_11_minus == pytest.approx(1.5)
    assert sf_22_plus == pytest.approx(1.5)
    assert sf_22_minus == pytest.approx(1.5)

    with pytest.raises(RuntimeError):
        rect_mat.get_sf_p()


def test_get_effective_material():
    """Tests effective material properties."""
    mat1 = Material(
        name="test",
        elastic_modulus=1,
        poissons_ratio=0.5,  # g = 1 / 3
        yield_strength=1,
        density=1,
        color="w",
    )
    mat2 = Material(
        name="test",
        elastic_modulus=2,
        poissons_ratio=0.5,  # g = 2 / 3
        yield_strength=1,
        density=1,
        color="w",
    )
    mat3 = Material(
        name="test",
        elastic_modulus=1,
        poissons_ratio=0,  # g = 0.5
        yield_strength=1,
        density=1,
        color="w",
    )

    rect1 = rectangular_section(d=2, b=2, material=mat1)
    rect2 = rectangular_section(d=2, b=2, material=mat2)
    rect3 = rectangular_section(d=2, b=2, material=mat3).shift_section(x_offset=2)
    geom1 = rect1 + rect3
    geom2 = rect2 + rect3

    geom1.create_mesh(mesh_sizes=[0])
    geom2.create_mesh(mesh_sizes=[0])
    sec1 = Section(geom1)
    sec2 = Section(geom2)
    sec1.calculate_geometric_properties()
    sec2.calculate_geometric_properties()

    # test e_eff
    assert sec1.get_e_eff() == pytest.approx(1.0)
    assert sec2.get_e_eff() == pytest.approx((4 * 1 + 4 * 2) / 8)

    # test g_eff
    # g = e / (2 * (1 + nu))
    assert sec1.get_g_eff() == pytest.approx((0.5 + 1 / 3.0) / 2)
    assert sec2.get_g_eff() == pytest.approx((4 * 2 / 3.0 + 4 * 0.5) / 8)

    # test nu_eff
    # nu = ea / (2 * ga) - 1
    assert sec1.get_nu_eff() == pytest.approx(4 / (2 * 5 / 3.0) - 1)
    assert sec2.get_nu_eff() == pytest.approx(6 / (2 * 7 / 3.0) - 1)
