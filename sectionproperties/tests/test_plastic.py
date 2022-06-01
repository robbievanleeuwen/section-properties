import pytest

from sectionproperties.analysis.section import Section
import sectionproperties.pre.library.primitive_sections as sections
from sectionproperties.pre.pre import Material


def test_rectangle():
    fy = 500
    E = 200e3
    b = 50
    d = 100

    steel = Material(
        name="Steel",
        elastic_modulus=E,
        poissons_ratio=0.3,
        yield_strength=fy,
        density=8.05e-6,
        color="grey",
    )

    Sx = b * d * d / 4
    Mp = Sx * fy

    geom_mat = sections.rectangular_section(d=d, b=b, material=steel)
    geom_nomat = sections.rectangular_section(d=d, b=b)

    geom_mat.create_mesh(mesh_sizes=[2.5])
    geom_nomat.create_mesh(mesh_sizes=[2.5])

    sec_mat = Section(geom_mat)
    sec_nomat = Section(geom_nomat)

    sec_mat.calculate_geometric_properties()
    sec_mat.calculate_plastic_properties()

    sec_nomat.calculate_geometric_properties()
    sec_nomat.calculate_plastic_properties()

    assert sec_nomat.get_s()[0] == pytest.approx(Sx)
    assert sec_mat.get_s()[0] == pytest.approx(Mp)
    assert sec_mat.get_s()[0] / fy == sec_nomat.get_s()[0]
