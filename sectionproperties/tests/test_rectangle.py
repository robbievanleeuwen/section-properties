import unittest
import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection
from sectionproperties.tests.helper_functions import validate_properties


class TestRectangle(unittest.TestCase):
    """Section properties are mostly validated against hand calcs. Results from
    Strand7 and Roark's Formulas for Stress and Strain are used for some
    warping properties."""

    def setUp(self):
        self.geometry = sections.RectangularSection(d=100, b=50)
        self.mesh = self.geometry.create_mesh(mesh_sizes=[10])
        self.section = CrossSection(self.geometry, self.mesh)

    def test_geometric(self):
        self.section.calculate_geometric_properties()

        val_list = []
        val_list.append({"prop": "area", "val": 100 * 50, "tol": None})
        val_list.append({"prop": "perimeter", "val": 100 * 2 + 50 * 2, "tol": None})
        val_list.append({"prop": "ea", "val": 1 * 100 * 50, "tol": None})
        val_list.append({"prop": "qx", "val": 100 * 50 * 50, "tol": None})
        val_list.append({"prop": "qy", "val": 100 * 50 * 25, "tol": None})
        val_list.append({"prop": "ixx_g", "val": 50 * 100 ** 3 / 3, "tol": None})
        val_list.append({"prop": "iyy_g", "val": 100 * 50 ** 3 / 3, "tol": None})
        val_list.append({"prop": "ixy_g", "val": 100 * 50 * 50 * 25, "tol": None})
        val_list.append({"prop": "cx", "val": 50 / 2, "tol": None})
        val_list.append({"prop": "cy", "val": 100 / 2, "tol": None})
        val_list.append({"prop": "ixx_c", "val": 50 * 100 ** 3 / 12, "tol": None})
        val_list.append({"prop": "iyy_c", "val": 100 * 50 ** 3 / 12, "tol": None})
        val_list.append({"prop": "ixy_c", "val": 0, "tol": None})
        val_list.append({"prop": "zxx_plus", "val": 50 * 100 ** 2 / 6, "tol": None})
        val_list.append({"prop": "zxx_minus", "val": 50 * 100 ** 2 / 6, "tol": None})
        val_list.append({"prop": "zyy_plus", "val": 100 * 50 ** 2 / 6, "tol": None})
        val_list.append({"prop": "zyy_minus", "val": 100 * 50 ** 2 / 6, "tol": None})
        val_list.append({"prop": "rx", "val": (50 * 100 ** 3 / 12 / 100 / 50) ** 0.5, "tol": None})
        val_list.append({"prop": "ry", "val": (100 * 50 ** 3 / 12 / 100 / 50) ** 0.5, "tol": None})
        val_list.append({"prop": "i11_c", "val": 50 * 100 ** 3 / 12, "tol": None})
        val_list.append({"prop": "i22_c", "val": 100 * 50 ** 3 / 12, "tol": None})
        val_list.append({"prop": "phi", "val": 0, "tol": None})
        val_list.append({"prop": "z11_plus", "val": 50 * 100 ** 2 / 6, "tol": None})
        val_list.append({"prop": "z11_minus", "val": 50 * 100 ** 2 / 6, "tol": None})
        val_list.append({"prop": "z22_plus", "val": 100 * 50 ** 2 / 6, "tol": None})
        val_list.append({"prop": "z22_minus", "val": 100 * 50 ** 2 / 6, "tol": None})
        val_list.append(
            {"prop": "r11", "val": (50 * 100 ** 3 / 12 / 100 / 50) ** 0.5, "tol": None}
        )
        val_list.append(
            {"prop": "r22", "val": (100 * 50 ** 3 / 12 / 100 / 50) ** 0.5, "tol": None}
        )

        validate_properties(self, val_list, self.section)

    def test_plastic(self):
        self.section.calculate_geometric_properties()
        self.section.calculate_plastic_properties()

        val_list = []
        val_list.append({"prop": "x_pc", "val": 50 / 2, "tol": None})
        val_list.append({"prop": "y_pc", "val": 100 / 2, "tol": None})
        val_list.append({"prop": "x11_pc", "val": 50 / 2, "tol": None})
        val_list.append({"prop": "y22_pc", "val": 100 / 2, "tol": None})
        val_list.append({"prop": "sxx", "val": 50 * 100 ** 2 / 4, "tol": None})
        val_list.append({"prop": "syy", "val": 100 * 50 ** 2 / 4, "tol": None})
        val_list.append({"prop": "s11", "val": 50 * 100 ** 2 / 4, "tol": None})
        val_list.append({"prop": "s22", "val": 100 * 50 ** 2 / 4, "tol": None})
        val_list.append({"prop": "sf_xx_plus", "val": 1.5, "tol": None})
        val_list.append({"prop": "sf_xx_minus", "val": 1.5, "tol": None})
        val_list.append({"prop": "sf_yy_plus", "val": 1.5, "tol": None})
        val_list.append({"prop": "sf_yy_minus", "val": 1.5, "tol": None})
        val_list.append({"prop": "sf_11_plus", "val": 1.5, "tol": None})
        val_list.append({"prop": "sf_11_minus", "val": 1.5, "tol": None})
        val_list.append({"prop": "sf_22_plus", "val": 1.5, "tol": None})
        val_list.append({"prop": "sf_22_minus", "val": 1.5, "tol": None})

        validate_properties(self, val_list, self.section)

    def test_warping(self):
        self.section.calculate_geometric_properties()
        self.section.calculate_warping_properties()

        val_list = []
        val_list.append({"prop": "j", "val": 2861002, "tol": 0.04})  # roark's
        val_list.append({"prop": "j", "val": 2.85852e6, "tol": 2e-5})  # st7
        val_list.append({"prop": "gamma", "val": 3.17542e8, "tol": None})  # st7
        val_list.append({"prop": "x_se", "val": 50 / 2, "tol": None})
        val_list.append({"prop": "y_se", "val": 100 / 2, "tol": None})
        val_list.append({"prop": "x11_se", "val": 0, "tol": None})
        val_list.append({"prop": "y22_se", "val": 0, "tol": None})
        val_list.append({"prop": "x_st", "val": 50 / 2, "tol": None})
        val_list.append({"prop": "y_st", "val": 100 / 2, "tol": None})
        val_list.append({"prop": "A_sx", "val": 5 / 6 * 100 * 50, "tol": None})
        val_list.append({"prop": "A_sy", "val": 5 / 6 * 100 * 50, "tol": None})
        val_list.append({"prop": "A_s11", "val": 5 / 6 * 100 * 50, "tol": None})
        val_list.append({"prop": "A_s22", "val": 5 / 6 * 100 * 50, "tol": None})


if __name__ == "__main__":
    unittest.main()
