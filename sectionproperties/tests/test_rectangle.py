import unittest
import numpy as np
import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection


class TestRectangle(unittest.TestCase):
    def setUp(self):
        self.geometry = sections.RectangularSection(d=100, b=50)
        self.mesh = self.geometry.create_mesh(mesh_sizes=[10])
        self.section = CrossSection(self.geometry, self.mesh)

    def test_geometric(self):
        self.section.calculate_geometric_properties()

        self.assertAlmostEqual(self.section.get_area(), 100 * 50)
        self.assertAlmostEqual(self.section.get_ea(), 1 * 100 * 50)
        (qx, qy) = self.section.get_q()
        self.assertAlmostEqual(qx, 100 * 50 * 50)
        self.assertAlmostEqual(qy, 100 * 50 * 25)
        (ixx_g, iyy_g, ixy_g) = self.section.get_ig()
        self.assertAlmostEqual(ixx_g, 50 * 100 ** 3 / 3)
        self.assertAlmostEqual(iyy_g, 100 * 50 ** 3 / 3)
        self.assertAlmostEqual(ixy_g, 100 * 50 * 50 * 25)
        (cx, cy) = self.section.get_c()
        self.assertAlmostEqual(cx, 50 / 2)
        self.assertAlmostEqual(cy, 100 / 2)
        (ixx_c, iyy_c, ixy_c) = self.section.get_ic()
        self.assertAlmostEqual(ixx_c, 50 * 100 ** 3 / 12)
        self.assertAlmostEqual(iyy_c, 100 * 50 ** 3 / 12)
        self.assertAlmostEqual(ixy_c, 0)
        (zxx_plus, zxx_minus, zyy_plus, zyy_minus) = self.section.get_z()
        self.assertAlmostEqual(zxx_plus, 50 * 100 ** 2 / 6)
        self.assertAlmostEqual(zxx_minus, 50 * 100 ** 2 / 6)
        self.assertAlmostEqual(zyy_plus, 100 * 50 ** 2 / 6)
        self.assertAlmostEqual(zyy_minus, 100 * 50 ** 2 / 6)
        (rx, ry) = self.section.get_rc()
        self.assertAlmostEqual(rx, (50 * 100 ** 3 / 12 / 100 / 50) ** 0.5)
        self.assertAlmostEqual(ry, (100 * 50 ** 3 / 12 / 100 / 50) ** 0.5)
        (i11_c, i22_c) = self.section.get_ip()
        self.assertAlmostEqual(i11_c, 50 * 100 ** 3 / 12)
        self.assertAlmostEqual(i22_c, 100 * 50 ** 3 / 12)
        phi = self.section.get_phi()
        self.assertAlmostEqual(phi, 0)
        (z11_plus, z11_minus, z22_plus, z22_minus) = self.section.get_zp()
        self.assertAlmostEqual(z11_plus, 50 * 100 ** 2 / 6)
        self.assertAlmostEqual(z11_minus, 50 * 100 ** 2 / 6)
        self.assertAlmostEqual(z22_plus, 100 * 50 ** 2 / 6)
        self.assertAlmostEqual(z22_minus, 100 * 50 ** 2 / 6)
        (r11, r22) = self.section.get_rp()
        self.assertAlmostEqual(r11, (50 * 100 ** 3 / 12 / 100 / 50) ** 0.5)
        self.assertAlmostEqual(r22, (100 * 50 ** 3 / 12 / 100 / 50) ** 0.5)

    def test_plastic(self):
        self.section.calculate_geometric_properties()
        self.section.calculate_plastic_properties()

        (x_pc, y_pc) = self.section.get_pc()
        self.assertTrue(np.isclose(x_pc, 50 / 2))
        self.assertTrue(np.isclose(y_pc, 100 / 2))
        (x11_pc, y22_pc) = self.section.get_pc_p()
        self.assertTrue(np.isclose(x11_pc, 50 / 2))
        self.assertTrue(np.isclose(y22_pc, 100 / 2))
        (sxx, syy) = self.section.get_s()
        self.assertTrue(np.isclose(sxx, 50 * 100 ** 2 / 4))
        self.assertTrue(np.isclose(syy, 100 * 50 ** 2 / 4))
        (s11, s22) = self.section.get_sp()
        self.assertTrue(np.isclose(s11, 50 * 100 ** 2 / 4))
        self.assertTrue(np.isclose(s22, 100 * 50 ** 2 / 4))
        (sf_xx_plus, sf_xx_minus, sf_yy_plus,
         sf_yy_minus) = self.section.get_sf()
        self.assertTrue(np.isclose(sf_xx_plus, 1.5))
        self.assertTrue(np.isclose(sf_xx_minus, 1.5))
        self.assertTrue(np.isclose(sf_yy_plus, 1.5))
        self.assertTrue(np.isclose(sf_yy_minus, 1.5))
        (sf_11_plus, sf_11_minus, sf_22_plus,
         sf_22_minus) = self.section.get_sf_p()
        self.assertTrue(np.isclose(sf_11_plus, 1.5))
        self.assertTrue(np.isclose(sf_11_minus, 1.5))
        self.assertTrue(np.isclose(sf_22_plus, 1.5))
        self.assertTrue(np.isclose(sf_22_minus, 1.5))

    def test_warping(self):
        self.section.calculate_geometric_properties()
        self.section.calculate_warping_properties()

        j = self.section.get_j()
        self.assertTrue(np.isclose(j, 2861002, rtol=0.05, atol=0))  # roarks
        # TODO: add Strand7 check and more resources?
        (x_se, y_se) = self.section.get_sc()
        self.assertTrue(np.isclose(x_se, 50 / 2))
        self.assertTrue(np.isclose(y_se, 100 / 2))
        (x11_se, y22_se) = self.section.get_sc_p()
        self.assertTrue(np.isclose(x11_se, 50 / 2))
        self.assertTrue(np.isclose(y22_se, 100 / 2))
        (x_st, y_st) = self.section.get_sc_t()
        self.assertTrue(np.isclose(x_st, 50 / 2))
        self.assertTrue(np.isclose(y_st, 100 / 2))
        gamma = self.section.get_gamma()
        # TODO: check gamma
        (A_sx, A_sy) = self.section.get_As()
        self.assertTrue(np.isclose(A_sx, 5 / 6 * 100 * 50))
        self.assertTrue(np.isclose(A_sy, 5 / 6 * 100 * 50))
        (A_s11, A_s22) = self.section.get_As_p()
        self.assertTrue(np.isclose(A_s11, 5 / 6 * 100 * 50))
        self.assertTrue(np.isclose(A_s22, 5 / 6 * 100 * 50))


if __name__ == "__main__":
    unittest.main()
