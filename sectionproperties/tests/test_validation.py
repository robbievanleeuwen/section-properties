import pytest
from shapely.geometry import Polygon
import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import Section
from sectionproperties.tests.helper_functions import validate_properties


# Setup for angle section
angle = sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
angle.create_mesh(mesh_sizes=[2.5])
angle_section = Section(angle)
angle_section.calculate_geometric_properties()
angle_section.calculate_plastic_properties()
angle_section.calculate_warping_properties()

def test_angle_all_properties():
    assert angle_section.section_props.area == pytest.approx(2746.73, abs=2e-4)
    assert angle_section.section_props.perimeter == pytest.approx(471, abs=1e-3)
    assert angle_section.section_props.cx == pytest.approx(21.2255, abs=2e-4)
    assert angle_section.section_props.cy == pytest.approx(50.9893, abs=2e-4)
    assert angle_section.section_props.ixx_g == pytest.approx(1.3428e7, abs=2e-4)
    assert angle_section.section_props.iyy_g == pytest.approx(2.95629e6, abs=2e-4)
    assert angle_section.section_props.ixy_g == pytest.approx(1.08669e6, abs=2e-4)
    assert angle_section.section_props.ixx_c == pytest.approx(6.28678e6, abs=2e-4)
    assert angle_section.section_props.iyy_c == pytest.approx(1.71882e6, abs=3e-4)
    assert angle_section.section_props.ixy_c == pytest.approx(-1.88603e6, abs=3e-4)
    assert angle_section.section_props.zxx_plus == pytest.approx(63496.0, abs=2e-4)
    assert angle_section.section_props.zxx_minus == pytest.approx(123296, abs=2e-4)
    assert angle_section.section_props.zyy_plus == pytest.approx(24992.1, abs=3e-4)
    assert angle_section.section_props.zyy_minus == pytest.approx(80979.0, abs=2e-4)
    assert angle_section.section_props.rx == pytest.approx(47.8416, abs=2e-4)
    assert angle_section.section_props.ry == pytest.approx(25.0154, abs=2e-4)
    assert angle_section.section_props.i11_c == pytest.approx(6.96484e6, abs=2e-4)
    assert angle_section.section_props.i22_c == pytest.approx(1.04076e6, abs=2e-4)
    assert angle_section.section_props.phi == pytest.approx(19.7744 - 180, abs=2e-4)
    assert angle_section.section_props.z11_plus == pytest.approx(97751.9, abs=2e-4)
    assert angle_section.section_props.z11_minus == pytest.approx(69403.3, abs=2e-4)
    assert angle_section.section_props.z22_plus == pytest.approx(27959.0, abs=2e-4)
    assert angle_section.section_props.z22_minus == pytest.approx(20761.6, abs=3e-4)
    assert angle_section.section_props.r11 == pytest.approx(50.3556, abs=2e-4)
    assert angle_section.section_props.r22 == pytest.approx(19.4656, abs=2e-4)
    assert angle_section.section_props.sxx == pytest.approx(113541, abs=2e-4)
    assert angle_section.section_props.syy == pytest.approx(45724.6, abs=2e-4)
    assert angle_section.section_props.sf_xx_plus == pytest.approx(113541 / 63496.0, abs=2e-4)
    assert angle_section.section_props.sf_xx_minus == pytest.approx(113541 / 123296, abs=2e-4)
    assert angle_section.section_props.sf_yy_plus == pytest.approx(45724.6 / 24992.1, abs=3e-4)
    assert angle_section.section_props.sf_yy_minus == pytest.approx(45724.6 / 80979.0, abs=2e-4)
    assert angle_section.section_props.s11 == pytest.approx(121030, abs=2e-4)
    assert angle_section.section_props.s22 == pytest.approx(43760.6, abs=2e-4)
    assert angle_section.section_props.sf_11_plus == pytest.approx(121030 / 97751.9, abs=2e-4)
    assert angle_section.section_props.sf_11_minus == pytest.approx(121030 / 69403.3, abs=2e-4)
    assert angle_section.section_props.sf_22_plus == pytest.approx(43760.6 / 27959.0, abs=2e-4)
    assert angle_section.section_props.sf_22_minus == pytest.approx(43760.6 / 20761.6, abs=3e-4)
    assert angle_section.section_props.j == pytest.approx(135333, abs=1e-3)
    assert angle_section.section_props.gamma == pytest.approx(1.62288e8, abs=5e-4)
    assert angle_section.section_props.A_s11 == pytest.approx(885.444, abs=2e-4)
    assert angle_section.section_props.A_s22 == pytest.approx(1459.72, abs=4e-4)
    assert angle_section.section_props.x11_se == pytest.approx(28.719, abs=1e-3)
    assert angle_section.section_props.y22_se == pytest.approx(35.2348, abs=5e-4)


# Setup custom section
custom_geom_points = [
            [-10, 0], [110, 0], [100, 10], [55, 10], [55, 90], [100, 90], [110, 100], [110, 110],
            [-10, 110], [-10, 100], [0, 90], [45, 90], [45, 10], [-10, 10]
        ]
custom_geom = sections.Geometry(Polygon(custom_geom_points))
custom_geom.create_mesh(mesh_sizes=[5])
custom_section = Section(custom_geom)
custom_section.calculate_geometric_properties()
custom_section.calculate_plastic_properties()
custom_section.calculate_warping_properties()


def test_custom_section_all_properties():
    assert custom_section.section_props.area == pytest.approx(4250, abs=None)
    assert custom_section.section_props.cx == pytest.approx(49.3333, abs=None)
    assert custom_section.section_props.cy == pytest.approx(65.0196, abs=None)
    assert custom_section.section_props.ixx_g == pytest.approx(2.56725e7, abs=None)
    assert custom_section.section_props.iyy_g == pytest.approx(1.41858e7, abs=None)
    assert custom_section.section_props.ixy_g == pytest.approx(1.37979e7, abs=None)
    assert custom_section.section_props.ixx_c == pytest.approx(7.70542e6, abs=None)
    assert custom_section.section_props.iyy_c == pytest.approx(3.84228e6, abs=None)
    assert custom_section.section_props.ixy_c == pytest.approx(165472, abs=None)
    assert custom_section.section_props.zxx_plus == pytest.approx(171306, abs=None)
    assert custom_section.section_props.zxx_minus == pytest.approx(118509, abs=None)
    assert custom_section.section_props.zyy_plus == pytest.approx(63334.2, abs=None)
    assert custom_section.section_props.zyy_minus == pytest.approx(64757.5, abs=None)
    assert custom_section.section_props.rx == pytest.approx(42.5798, abs=None)
    assert custom_section.section_props.ry == pytest.approx(30.0677, abs=None)
    assert custom_section.section_props.phi == pytest.approx(177.552 - 180, abs=1e-4)
    assert custom_section.section_props.i11_c == pytest.approx(7.71249e6, abs=None)
    assert custom_section.section_props.i22_c == pytest.approx(3.8352e6, abs=None)
    assert custom_section.section_props.z11_plus == pytest.approx(162263, abs=None)
    assert custom_section.section_props.z11_minus == pytest.approx(114268, abs=None)
    assert custom_section.section_props.z22_plus == pytest.approx(60503.0, abs=None)
    assert custom_section.section_props.z22_minus == pytest.approx(62666.1, abs=None)
    assert custom_section.section_props.r11 == pytest.approx(42.5993, abs=None)
    assert custom_section.section_props.r22 == pytest.approx(30.04, abs=None)
    assert custom_section.section_props.sxx == pytest.approx(153196, abs=None)
    assert custom_section.section_props.syy == pytest.approx(101494, abs=None)
    assert custom_section.section_props.sf_xx_plus == pytest.approx(153196 / 171306, abs=None)
    assert custom_section.section_props.sf_xx_minus == pytest.approx(153196 / 118509, abs=None)
    assert custom_section.section_props.sf_yy_plus == pytest.approx(101494 / 63334.2, abs=None)
    assert custom_section.section_props.sf_yy_minus == pytest.approx(101494 / 64757.5, abs=None)
    assert custom_section.section_props.s11 == pytest.approx(153347, abs=None)
    assert custom_section.section_props.s22 == pytest.approx(101501, abs=None)
    assert custom_section.section_props.sf_11_plus == pytest.approx(153347 / 162263, abs=None)
    assert custom_section.section_props.sf_11_minus == pytest.approx(153347 / 114268, abs=None)
    assert custom_section.section_props.sf_22_plus == pytest.approx(101501 / 60503.0, abs=None)
    assert custom_section.section_props.sf_22_minus == pytest.approx(101501 / 62666.1, abs=None)
    assert custom_section.section_props.j == pytest.approx(347040, abs=5e-3)
    assert custom_section.section_props.gamma == pytest.approx(7.53539e9, abs=1e-3)
    assert custom_section.section_props.A_s11 == pytest.approx(2945.53, abs=5e-4)
    assert custom_section.section_props.A_s22 == pytest.approx(956.014, abs=5e-4)
    assert custom_section.section_props.x11_se == pytest.approx(1.9134, abs=5e-3)
    assert custom_section.section_props.y22_se == pytest.approx(3.02028, abs=5e-3)

# class TestValidation(unittest.TestCase):
#     def test_angle(self):
#         """Section properties are validated against results from the Strand7
#         beam section utility."""

#         geometry = sections.angle_section(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
#         geometry.compile_geometry()
#         geometry.create_mesh(mesh_sizes=[2.5])

#         section = Section(geometry)

#         val_list = []
        # val_list.append({"prop": "area", "val": 2746.73, "tol": 2e-4})
        # val_list.append({"prop": "perimeter", "val": 471, "tol": 1e-3})
        # val_list.append({"prop": "cx", "val": 21.2255, "tol": 2e-4})
        # val_list.append({"prop": "cy", "val": 50.9893, "tol": 2e-4})
        # val_list.append({"prop": "ixx_g", "val": 1.3428e7, "tol": 2e-4})
        # val_list.append({"prop": "iyy_g", "val": 2.95629e6, "tol": 2e-4})
        # val_list.append({"prop": "ixy_g", "val": 1.08669e6, "tol": 2e-4})
        # val_list.append({"prop": "ixx_c", "val": 6.28678e6, "tol": 2e-4})
        # val_list.append({"prop": "iyy_c", "val": 1.71882e6, "tol": 3e-4})
        # val_list.append({"prop": "ixy_c", "val": -1.88603e6, "tol": 3e-4})
        # val_list.append({"prop": "zxx_plus", "val": 63496.0, "tol": 2e-4})
        # val_list.append({"prop": "zxx_minus", "val": 123296, "tol": 2e-4})
        # val_list.append({"prop": "zyy_plus", "val": 24992.1, "tol": 3e-4})
        # val_list.append({"prop": "zyy_minus", "val": 80979.0, "tol": 2e-4})
        # val_list.append({"prop": "rx", "val": 47.8416, "tol": 2e-4})
        # val_list.append({"prop": "ry", "val": 25.0154, "tol": 2e-4})
        # val_list.append({"prop": "i11_c", "val": 6.96484e6, "tol": 2e-4})
        # val_list.append({"prop": "i22_c", "val": 1.04076e6, "tol": 2e-4})
        # val_list.append({"prop": "phi", "val": 19.7744 - 180, "tol": 2e-4})
        # val_list.append({"prop": "z11_plus", "val": 97751.9, "tol": 2e-4})
        # val_list.append({"prop": "z11_minus", "val": 69403.3, "tol": 2e-4})
        # val_list.append({"prop": "z22_plus", "val": 27959.0, "tol": 2e-4})
        # val_list.append({"prop": "z22_minus", "val": 20761.6, "tol": 3e-4})
        # val_list.append({"prop": "r11", "val": 50.3556, "tol": 2e-4})
        # val_list.append({"prop": "r22", "val": 19.4656, "tol": 2e-4})
        # val_list.append({"prop": "sxx", "val": 113541, "tol": 2e-4})
        # val_list.append({"prop": "syy", "val": 45724.6, "tol": 2e-4})
        # val_list.append({"prop": "sf_xx_plus", "val": 113541 / 63496.0, "tol": 2e-4})
        # val_list.append({"prop": "sf_xx_minus", "val": 113541 / 123296, "tol": 2e-4})
        # val_list.append({"prop": "sf_yy_plus", "val": 45724.6 / 24992.1, "tol": 3e-4})
        # val_list.append({"prop": "sf_yy_minus", "val": 45724.6 / 80979.0, "tol": 2e-4})
        # val_list.append({"prop": "s11", "val": 121030, "tol": 2e-4})
        # val_list.append({"prop": "s22", "val": 43760.6, "tol": 2e-4})
        # val_list.append({"prop": "sf_11_plus", "val": 121030 / 97751.9, "tol": 2e-4})
        # val_list.append({"prop": "sf_11_minus", "val": 121030 / 69403.3, "tol": 2e-4})
        # val_list.append({"prop": "sf_22_plus", "val": 43760.6 / 27959.0, "tol": 2e-4})
        # val_list.append({"prop": "sf_22_minus", "val": 43760.6 / 20761.6, "tol": 3e-4})
        # val_list.append({"prop": "j", "val": 135333, "tol": 1e-3})
        # val_list.append({"prop": "gamma", "val": 1.62288e8, "tol": 5e-4})
        # val_list.append({"prop": "A_s11", "val": 885.444, "tol": 2e-4})
        # val_list.append({"prop": "A_s22", "val": 1459.72, "tol": 4e-4})
        # val_list.append({"prop": "x11_se", "val": 28.719, "tol": 1e-3})
        # val_list.append({"prop": "y22_se", "val": 35.2348, "tol": 5e-4})

#         section.calculate_geometric_properties()
#         section.calculate_plastic_properties()
#         section.calculate_warping_properties()

#         validate_properties(self, val_list, section)

#     def test_custom(self):
#         """Section properties are validated against results from the Strand7
#         beam section utility."""

#         points = [
#             [-10, 0], [110, 0], [100, 10], [55, 10], [55, 90], [100, 90], [110, 100], [110, 110],
#             [-10, 110], [-10, 100], [0, 90], [45, 90], [45, 10], [-10, 10]
#         ]


#         geometry = sections.Geometry(Polygon(points))
#         mesh = geometry.create_mesh(mesh_sizes=[5])
#         section = Section(geometry, mesh)

#         val_list = []
#         val_list.append({"prop": "area", "val": 4250, "tol": None})
#         val_list.append({"prop": "cx", "val": 49.3333, "tol": None})
#         val_list.append({"prop": "cy", "val": 65.0196, "tol": None})
#         val_list.append({"prop": "ixx_g", "val": 2.56725e7, "tol": None})
#         val_list.append({"prop": "iyy_g", "val": 1.41858e7, "tol": None})
#         val_list.append({"prop": "ixy_g", "val": 1.37979e7, "tol": None})
#         val_list.append({"prop": "ixx_c", "val": 7.70542e6, "tol": None})
#         val_list.append({"prop": "iyy_c", "val": 3.84228e6, "tol": None})
#         val_list.append({"prop": "ixy_c", "val": 165472, "tol": None})
#         val_list.append({"prop": "zxx_plus", "val": 171306, "tol": None})
#         val_list.append({"prop": "zxx_minus", "val": 118509, "tol": None})
#         val_list.append({"prop": "zyy_plus", "val": 63334.2, "tol": None})
#         val_list.append({"prop": "zyy_minus", "val": 64757.5, "tol": None})
#         val_list.append({"prop": "rx", "val": 42.5798, "tol": None})
#         val_list.append({"prop": "ry", "val": 30.0677, "tol": None})
#         val_list.append({"prop": "phi", "val": 177.552 - 180, "tol": 1e-4})
#         val_list.append({"prop": "i11_c", "val": 7.71249e6, "tol": None})
#         val_list.append({"prop": "i22_c", "val": 3.8352e6, "tol": None})
#         val_list.append({"prop": "z11_plus", "val": 162263, "tol": None})
#         val_list.append({"prop": "z11_minus", "val": 114268, "tol": None})
#         val_list.append({"prop": "z22_plus", "val": 60503.0, "tol": None})
#         val_list.append({"prop": "z22_minus", "val": 62666.1, "tol": None})
#         val_list.append({"prop": "r11", "val": 42.5993, "tol": None})
#         val_list.append({"prop": "r22", "val": 30.04, "tol": None})
#         val_list.append({"prop": "sxx", "val": 153196, "tol": None})
#         val_list.append({"prop": "syy", "val": 101494, "tol": None})
#         val_list.append({"prop": "sf_xx_plus", "val": 153196 / 171306, "tol": None})
#         val_list.append({"prop": "sf_xx_minus", "val": 153196 / 118509, "tol": None})
#         val_list.append({"prop": "sf_yy_plus", "val": 101494 / 63334.2, "tol": None})
#         val_list.append({"prop": "sf_yy_minus", "val": 101494 / 64757.5, "tol": None})
#         val_list.append({"prop": "s11", "val": 153347, "tol": None})
#         val_list.append({"prop": "s22", "val": 101501, "tol": None})
#         val_list.append({"prop": "sf_11_plus", "val": 153347 / 162263, "tol": None})
#         val_list.append({"prop": "sf_11_minus", "val": 153347 / 114268, "tol": None})
#         val_list.append({"prop": "sf_22_plus", "val": 101501 / 60503.0, "tol": None})
#         val_list.append({"prop": "sf_22_minus", "val": 101501 / 62666.1, "tol": None})
#         val_list.append({"prop": "j", "val": 347040, "tol": 5e-3})
#         val_list.append({"prop": "gamma", "val": 7.53539e9, "tol": 1e-3})
#         val_list.append({"prop": "A_s11", "val": 2945.53, "tol": 5e-4})
#         val_list.append({"prop": "A_s22", "val": 956.014, "tol": 5e-4})
#         val_list.append({"prop": "x11_se", "val": 1.9134, "tol": 5e-3})
#         val_list.append({"prop": "y22_se", "val": 3.02028, "tol": 5e-3})

#         section.calculate_geometric_properties()
#         section.calculate_plastic_properties()
#         section.calculate_warping_properties()

#         validate_properties(self, val_list, section)


# if __name__ == "__main__":
#     unittest.main()
