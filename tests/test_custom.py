"""Testing for a custom section."""

import pytest

import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection


@pytest.mark.slow
class TestCustom:
    """Test a custom section section."""

    @pytest.fixture
    def helpers(self, helpers):
        self.validate_properties = helpers.validate_properties

    @pytest.fixture(scope='class')
    def section(self):
        points = [
            [-10, 0],
            [110, 0],
            [100, 10],
            [55, 10],
            [55, 90],
            [100, 90],
            [110, 100],
            [110, 110],
            [-10, 110],
            [-10, 100],
            [0, 90],
            [45, 90],
            [45, 10],
            [-10, 10],
        ]
        facets = [
            [0, 1],
            [1, 2],
            [2, 3],
            [3, 4],
            [4, 5],
            [5, 6],
            [6, 7],
            [7, 8],
            [8, 9],
            [9, 10],
            [10, 11],
            [11, 12],
            [12, 13],
            [13, 0],
        ]
        holes = []
        control_points = [[0, 5]]

        geometry = sections.CustomSection(points, facets, holes, control_points)
        mesh = geometry.create_mesh(mesh_sizes=[5])
        section = CrossSection(geometry, mesh)
        yield section

    def test_geometric(self, section, benchmark):
        benchmark(section.calculate_geometric_properties)

        val_list = []
        val_list.append({'prop': 'area', 'val': 4250})
        val_list.append({'prop': 'cx', 'val': 49.3333})
        val_list.append({'prop': 'cy', 'val': 65.0196})
        val_list.append({'prop': 'ixx_g', 'val': 2.56725e7})
        val_list.append({'prop': 'iyy_g', 'val': 1.41858e7})
        val_list.append({'prop': 'ixy_g', 'val': 1.37979e7})
        val_list.append({'prop': 'ixx_c', 'val': 7.70542e6})
        val_list.append({'prop': 'iyy_c', 'val': 3.84228e6})
        val_list.append({'prop': 'ixy_c', 'val': 165472})
        val_list.append({'prop': 'zxx_plus', 'val': 171306})
        val_list.append({'prop': 'zxx_minus', 'val': 118509})
        val_list.append({'prop': 'zyy_plus', 'val': 63334.2})
        val_list.append({'prop': 'zyy_minus', 'val': 64757.5})
        val_list.append({'prop': 'rx', 'val': 42.5798})
        val_list.append({'prop': 'ry', 'val': 30.0677})
        val_list.append({'prop': 'phi', 'val': 177.552 - 180, 'tol': 1e-4})
        val_list.append({'prop': 'i11_c', 'val': 7.71249e6})
        val_list.append({'prop': 'i22_c', 'val': 3.8352e6})
        val_list.append({'prop': 'z11_plus', 'val': 162263})
        val_list.append({'prop': 'z11_minus', 'val': 114268})
        val_list.append({'prop': 'z22_plus', 'val': 60503.0})
        val_list.append({'prop': 'z22_minus', 'val': 62666.1})
        val_list.append({'prop': 'r11', 'val': 42.5993})
        val_list.append({'prop': 'r22', 'val': 30.04})

        self.validate_properties(val_list, section)

    def test_plastic(self, section, benchmark):
        benchmark(section.calculate_plastic_properties)

        val_list = []
        val_list.append({'prop': 'sxx', 'val': 153196})
        val_list.append({'prop': 'syy', 'val': 101494})
        val_list.append({'prop': 's11', 'val': 153347})
        val_list.append({'prop': 's22', 'val': 101501})
        val_list.append({'prop': 'sf_xx_plus', 'val': 153196 / 171306})
        val_list.append({'prop': 'sf_xx_minus', 'val': 153196 / 118509})
        val_list.append({'prop': 'sf_yy_plus', 'val': 101494 / 63334.2})
        val_list.append({'prop': 'sf_yy_minus', 'val': 101494 / 64757.5})
        val_list.append({'prop': 'sf_11_plus', 'val': 153347 / 162263})
        val_list.append({'prop': 'sf_11_minus', 'val': 153347 / 114268})
        val_list.append({'prop': 'sf_22_plus', 'val': 101501 / 60503.0})
        val_list.append({'prop': 'sf_22_minus', 'val': 101501 / 62666.1})

        self.validate_properties(val_list, section)

    def test_warping(self, section, benchmark):
        benchmark(section.calculate_warping_properties)

        val_list = []
        val_list.append({'prop': 'j', 'val': 347040, 'tol': 5e-3})
        val_list.append({'prop': 'gamma', 'val': 7.53539e9, 'tol': 1e-3})
        val_list.append({'prop': 'A_s11', 'val': 2945.53, 'tol': 5e-4})
        val_list.append({'prop': 'A_s22', 'val': 956.014, 'tol': 5e-4})
        val_list.append({'prop': 'x11_se', 'val': 1.9134, 'tol': 5e-3})
        val_list.append({'prop': 'y22_se', 'val': 3.02028, 'tol': 5e-3})

        self.validate_properties(val_list, section)
