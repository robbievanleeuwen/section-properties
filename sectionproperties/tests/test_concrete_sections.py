import pytest_check as check
import numpy as np
import sectionproperties.pre.pre as pre
import sectionproperties.pre.library.concrete_sections as cs
import sectionproperties.pre.library.primitive_sections as ps

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


def test_concrete_column_section():

    concrete = pre.Material(
        name="Concrete",
        elastic_modulus=30.1e3,
        poissons_ratio=0.2,
        yield_strength=32,
        density=2.4e-6,
        color="lightgrey",
    )
    steel = pre.Material(
        name="Steel",
        elastic_modulus=200e3,
        poissons_ratio=0.3,
        yield_strength=500,
        density=7.85e-6,
        color="grey",
    )

    geometry = cs.concrete_column_section(
        b=300,
        d=600,
        dia_bar=40,
        bar_area=500,
        cover=40,
        n_bars_b=3,
        n_bars_d=6,
        conc_mat=concrete,
        steel_mat=steel,
        filled=False,
        n_circle=4,
    )  # NOTE: Bar diam and Bar area do not match. This is intentional.
    geometry.create_mesh(mesh_sizes=[500])

    # check geometry is created correctly
    conc_area = 0
    steel_area = 0

    for geom in geometry.geoms:
        if geom.material == concrete:
            conc_area += geom.calculate_area()
        elif geom.material == steel:
            steel_area += geom.calculate_area()
        else:
            raise ValueError(
                "Material {0} is not correctly assigned".format(geom.material)
            )

    net_area = 300 * 600
    actual_steel_area = 14 * 500.0

    check.almost_equal(conc_area, net_area - actual_steel_area, rel=r_tol)
    check.almost_equal(steel_area, actual_steel_area, rel=r_tol)

    bar_centroids = [tuple(geom.geom.centroid.coords[0]) for geom in geometry.geoms[1:]]

    from collections import Counter

    x_coords = Counter(round(coord[0], 0) for coord in bar_centroids)
    y_coords = Counter(round(coord[1], 0) for coord in bar_centroids)

    # Validate that we have 14 bars with the correct x-coordinates
    check.equal(x_coords.get(60), 6)
    check.equal(x_coords.get(150), 2)
    check.equal(x_coords.get(240), 6)

    # Validate that we have 14 bars with the correct y-coordinates
    check.equal(y_coords.get(60), 3)
    check.equal(y_coords.get(156), 2)
    check.equal(y_coords.get(252), 2)
    check.equal(y_coords.get(348), 2)
    check.equal(y_coords.get(444), 2)
    check.equal(y_coords.get(540), 3)


def test_add_bar():
    rect = ps.rectangular_section(b=400, d=600)
    steel = pre.Material(
        name="Steel",
        elastic_modulus=200e3,
        poissons_ratio=0.3,
        yield_strength=500,
        density=7.85e-6,
        color="grey",
    )
    rect = cs.add_bar(rect, area=500, x=100, y=100, material=steel)
    rect = cs.add_bar(rect, area=500, x=200, y=200, material=steel)
    rect_area = rect.geoms[0].geom.area
    steel_area = 2 * 500.0
    check.almost_equal(rect_area, 400 * 600 - steel_area)

    bar_1 = rect.geoms[1]
    bar_2 = rect.geoms[2]
    check.almost_equal(bar_1.calculate_centroid(), (100, 100))
    check.almost_equal(bar_2.calculate_centroid(), (200, 200))
