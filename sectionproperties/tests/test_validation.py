import pytest_check as check
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

tol = 1e-3


def test_angle_all_properties():
    check.almost_equal(angle_section.section_props.area, 2747.059)
    check.almost_equal(angle_section.section_props.perimeter, 471.3501)
    check.almost_equal(angle_section.section_props.cx, 2.122282e1)
    check.almost_equal(angle_section.section_props.cy, 5.098127e1)
    check.almost_equal(angle_section.section_props.ixx_g, 1.342632e7)
    check.almost_equal(angle_section.section_props.iyy_g, 2.955753e6)
    check.almost_equal(angle_section.section_props.ixy_g, 1.086603e6)
    check.almost_equal(angle_section.section_props.ixx_c, 6.286470e6)
    check.almost_equal(angle_section.section_props.iyy_c, 1.718455e6)
    check.almost_equal(angle_section.section_props.ixy_c, -1.885622e6)
    check.almost_equal(angle_section.section_props.zxx_plus, 6.348769e4)
    check.almost_equal(angle_section.section_props.zxx_minus, 1.233094e5)
    check.almost_equal(angle_section.section_props.zyy_plus, 8.097207e4)
    check.almost_equal(angle_section.section_props.zyy_minus, 8.097207e4)
    check.almost_equal(angle_section.section_props.rx, 4.783761e1)
    check.almost_equal(angle_section.section_props.ry, 2.501124e1)
    check.almost_equal(angle_section.section_props.i11_c, 6.964263e6)
    check.almost_equal(angle_section.section_props.i22_c, 1.040662e6)
    check.almost_equal(angle_section.section_props.phi, -1.602289e2)
    check.almost_equal(angle_section.section_props.z11_plus, 9.775662e4)
    check.almost_equal(angle_section.section_props.z11_minus, 6.939239e4)
    check.almost_equal(angle_section.section_props.z22_plus, 2.796211e4)
    check.almost_equal(angle_section.section_props.z22_minus, 2.076613e4)
    check.almost_equal(angle_section.section_props.r11, 5.035048e1)
    check.almost_equal(angle_section.section_props.r22, 1.946350e1)
    check.almost_equal(angle_section.section_props.sxx, 1.135392e5)
    check.almost_equal(angle_section.section_props.syy, 4.572269e4)
    check.almost_equal(angle_section.section_props.sf_xx_plus, 1.788366)
    check.almost_equal(angle_section.section_props.sf_xx_minus, 9.207672e-1)
    check.almost_equal(angle_section.section_props.sf_yy_plus, 1.829944)
    check.almost_equal(angle_section.section_props.sf_yy_minus, 5.646723e-1)
    check.almost_equal(angle_section.section_props.s11, 1.210275e5)
    check.almost_equal(angle_section.section_props.s22, 4.376054e4)
    check.almost_equal(angle_section.section_props.sf_11_plus, 1.238049)
    check.almost_equal(angle_section.section_props.sf_11_minus, 1.744103)
    check.almost_equal(angle_section.section_props.sf_22_plus, 1.564994)
    check.almost_equal(angle_section.section_props.sf_22_minus, 2.107303)
    check.almost_equal(angle_section.section_props.j, 1.354663e5)
    check.almost_equal(angle_section.section_props.gamma, 162220735.49)
    check.almost_equal(angle_section.section_props.A_s11, 8.855951e2)
    check.almost_equal(angle_section.section_props.A_s22, 1.460240e3)
    check.almost_equal(angle_section.section_props.x11_se, 2.870404e1)
    check.almost_equal(angle_section.section_props.y22_se, 3.522141e1)


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
    check.almost_equal(custom_section.section_props.area, 4250)
    check.almost_equal(custom_section.section_props.cx, 4.933333e1)
    check.almost_equal(custom_section.section_props.cy, 6.501961e1)
    check.almost_equal(custom_section.section_props.ixx_g, 2.564250e7)
    check.almost_equal(custom_section.section_props.iyy_g, 1.418583e7)
    check.almost_equal(custom_section.section_props.ixy_g, 1.379792e7)
    check.almost_equal(custom_section.section_props.ixx_c, 7.705415e6)
    check.almost_equal(custom_section.section_props.iyy_c, 3.842278e6)
    check.almost_equal(custom_section.section_props.ixy_c, 1.654722e5)
    check.almost_equal(custom_section.section_props.zxx_plus, 1.713061e5)
    check.almost_equal(custom_section.section_props.zxx_minus, 1.185091e5)
    check.almost_equal(custom_section.section_props.zyy_plus, 6.333425e4)
    check.almost_equal(custom_section.section_props.zyy_minus, 6.475749e4)
    check.almost_equal(custom_section.section_props.rx, 4.25797e1)
    check.almost_equal(custom_section.section_props.ry, 3.006768e1)
    check.almost_equal(custom_section.section_props.phi, -2.448209)
    check.almost_equal(custom_section.section_props.i11_c, 7.712490e6)
    check.almost_equal(custom_section.section_props.i22_c, 3.835203e6)
    check.almost_equal(custom_section.section_props.z11_plus, 1.622630e5)
    check.almost_equal(custom_section.section_props.z11_minus, 1.142680e5)
    check.almost_equal(custom_section.section_props.z22_plus, 6.050295e4)
    check.almost_equal(custom_section.section_props.z22_minus, 6.266613e4)
    check.almost_equal(custom_section.section_props.r11, 4.259934e1)
    check.almost_equal(custom_section.section_props.r22, 3.003998e1)
    check.almost_equal(custom_section.section_props.sxx, 1.531971e5)
    check.almost_equal(custom_section.section_props.syy, 1.014943e5)
    check.almost_equal(custom_section.section_props.sf_xx_plus, 8.942884e-1)
    check.almost_equal(custom_section.section_props.sf_xx_minus, 8.942884e-1)
    check.almost_equal(custom_section.section_props.sf_yy_plus, 1.602519)
    check.almost_equal(custom_section.section_props.sf_yy_minus, 1.527298)
    check.almost_equal(custom_section.section_props.s11, 1.533463e5)
    check.almost_equal(custom_section.section_props.s22, 1.015010e5)
    check.almost_equal(custom_section.section_props.sf_11_plus, 9.450478e-1)
    check.almost_equal(custom_section.section_props.sf_11_minus, 1.341988)
    check.almost_equal(custom_section.section_props.sf_22_plus, 1.677621)
    check.almost_equal(custom_section.section_props.sf_22_minus, 1.619711)
    check.almost_equal(custom_section.section_props.j, 3.477399e5)
    check.almost_equal(custom_section.section_props.gamma, 7.532929e9)
    check.almost_equal(custom_section.section_props.A_s11, 2.945692e3)
    check.almost_equal(custom_section.section_props.A_s22, 9.564143e2)
    check.almost_equal(custom_section.section_props.x11_se, 1.916270)
    check.almost_equal(custom_section.section_props.y22_se, 3.017570)

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


# if __name__, "__main__":
#     unittest.main()
