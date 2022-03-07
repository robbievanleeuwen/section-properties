import pytest_check as check
import numpy as np
import sectionproperties.pre.pre as pre
import sectionproperties.pre.library.concrete_sections as cs

r_tol = 1e-6
conc_mat = pre.Material(
    name="Concrete",
    elastic_modulus=32.8e3,
    poissons_ratio=0.2,
    density=2.4e-6,
    yield_strength=40,
    color="lightgrey",
)

steel_mat = pre.Material(
    name="Steel",
    elastic_modulus=200e3,
    poissons_ratio=0.3,
    density=7.85e-6,
    yield_strength=500,
    color="grey",
)


def test_concrete_rectangular_section():
    rect = cs.concrete_rectangular_section(
        b=300,
        d=600,
        dia_top=16,
        n_top=3,
        dia_bot=20,
        n_bot=3,
        n_circle=16,
        cover=30,
        area_top=200,
        area_bot=310,
        conc_mat=conc_mat,
        steel_mat=steel_mat,
    )

    # check geometry is created correctly
    conc_area = 0
    steel_area = 0

    for geom in rect.geoms:
        if geom.material == conc_mat:
            conc_area += geom.calculate_area()
        elif geom.material == steel_mat:
            steel_area += geom.calculate_area()
        else:
            raise ValueError(
                "Material {0} is not correctly assigned".format(geom.material)
            )

    net_area = 600 * 300
    actual_steel_area = 3 * (200 + 310)

    # check areas
    check.almost_equal(conc_area, net_area - actual_steel_area, rel=r_tol)
    check.almost_equal(steel_area, actual_steel_area, rel=r_tol)


def test_concrete_tee_section():
    rect = cs.concrete_tee_section(
        b=300,
        d=900,
        b_f=1200,
        d_f=200,
        dia_top=20,
        n_top=6,
        dia_bot=24,
        n_bot=3,
        n_circle=16,
        cover=30,
        area_top=310,
        area_bot=450,
        conc_mat=conc_mat,
        steel_mat=steel_mat,
    )

    # check geometry is created correctly
    conc_area = 0
    steel_area = 0

    for geom in rect.geoms:
        if geom.material == conc_mat:
            conc_area += geom.calculate_area()
        elif geom.material == steel_mat:
            steel_area += geom.calculate_area()
        else:
            raise ValueError(
                "Material {0} is not correctly assigned".format(geom.material)
            )

    net_area = 700 * 300 + 1200 * 200
    actual_steel_area = 6 * 310 + 3 * 450

    # check areas
    check.almost_equal(conc_area, net_area - actual_steel_area, rel=r_tol)
    check.almost_equal(steel_area, actual_steel_area, rel=r_tol)


def test_concrete_circular_section():
    rect = cs.concrete_circular_section(
        d=600,
        n=64,
        dia=20,
        n_bar=8,
        n_circle=16,
        cover=45,
        area_conc=np.pi * 600 * 600 / 4,
        area_bar=310,
        conc_mat=conc_mat,
        steel_mat=steel_mat,
    )

    # check geometry is created correctly
    conc_area = 0
    steel_area = 0

    for geom in rect.geoms:
        if geom.material == conc_mat:
            conc_area += geom.calculate_area()
        elif geom.material == steel_mat:
            steel_area += geom.calculate_area()
        else:
            raise ValueError(
                "Material {0} is not correctly assigned".format(geom.material)
            )

    net_area = np.pi * 600 * 600 / 4
    actual_steel_area = 8 * 310

    # check areas
    check.almost_equal(conc_area, net_area - actual_steel_area, rel=r_tol)
    check.almost_equal(steel_area, actual_steel_area, rel=r_tol)
