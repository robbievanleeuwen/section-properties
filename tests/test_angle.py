"""Testing for an angle section."""

import pytest

import sectionproperties.pre.sections as sections
from sectionproperties.analysis.cross_section import CrossSection


@pytest.mark.slow
class TestAngle:
    """Test an angle section."""

    @pytest.fixture
    def helpers(self, helpers):
        self.validate_properties = helpers.validate_properties

    @pytest.fixture(scope='class')
    def section(self):
        geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
        mesh = geometry.create_mesh(mesh_sizes=[2.5])
        section = CrossSection(geometry, mesh)
        yield section

    def test_geometric(self, section, benchmark):
        benchmark(section.calculate_geometric_properties)

        val_list = []
        val_list.append({'prop': 'area', 'val': 2746.73, 'tol': 2e-4})
        val_list.append({'prop': 'perimeter', 'val': 471, 'tol': 1e-3})
        val_list.append({'prop': 'cx', 'val': 21.2255, 'tol': 2e-4})
        val_list.append({'prop': 'cy', 'val': 50.9893, 'tol': 2e-4})
        val_list.append({'prop': 'ixx_g', 'val': 1.3428e7, 'tol': 2e-4})
        val_list.append({'prop': 'iyy_g', 'val': 2.95629e6, 'tol': 2e-4})
        val_list.append({'prop': 'ixy_g', 'val': 1.08669e6, 'tol': 2e-4})
        val_list.append({'prop': 'ixx_c', 'val': 6.28678e6, 'tol': 2e-4})
        val_list.append({'prop': 'iyy_c', 'val': 1.71882e6, 'tol': 3e-4})
        val_list.append({'prop': 'ixy_c', 'val': -1.88603e6, 'tol': 3e-4})
        val_list.append({'prop': 'zxx_plus', 'val': 63496.0, 'tol': 2e-4})
        val_list.append({'prop': 'zxx_minus', 'val': 123296, 'tol': 2e-4})
        val_list.append({'prop': 'zyy_plus', 'val': 24992.1, 'tol': 3e-4})
        val_list.append({'prop': 'zyy_minus', 'val': 80979.0, 'tol': 2e-4})
        val_list.append({'prop': 'rx', 'val': 47.8416, 'tol': 2e-4})
        val_list.append({'prop': 'ry', 'val': 25.0154, 'tol': 2e-4})
        val_list.append({'prop': 'i11_c', 'val': 6.96484e6, 'tol': 2e-4})
        val_list.append({'prop': 'i22_c', 'val': 1.04076e6, 'tol': 2e-4})
        val_list.append({'prop': 'phi', 'val': 19.7744 - 180, 'tol': 2e-4})
        val_list.append({'prop': 'z11_plus', 'val': 97751.9, 'tol': 2e-4})
        val_list.append({'prop': 'z11_minus', 'val': 69403.3, 'tol': 2e-4})
        val_list.append({'prop': 'z22_plus', 'val': 27959.0, 'tol': 2e-4})
        val_list.append({'prop': 'z22_minus', 'val': 20761.6, 'tol': 3e-4})
        val_list.append({'prop': 'r11', 'val': 50.3556, 'tol': 2e-4})
        val_list.append({'prop': 'r22', 'val': 19.4656, 'tol': 2e-4})

        self.validate_properties(val_list, section)

    def test_plastic(self, section, benchmark):
        benchmark(section.calculate_plastic_properties)

        val_list = []
        val_list.append({'prop': 'sxx', 'val': 113541, 'tol': 2e-4})
        val_list.append({'prop': 'syy', 'val': 45724.6, 'tol': 2e-4})
        val_list.append({'prop': 's11', 'val': 121030, 'tol': 2e-4})
        val_list.append({'prop': 's22', 'val': 43760.6, 'tol': 2e-4})
        val_list.append({'prop': 'sf_xx_plus', 'val': 113541 / 63496.0, 'tol': 2e-4})
        val_list.append({'prop': 'sf_xx_minus', 'val': 113541 / 123296, 'tol': 2e-4})
        val_list.append({'prop': 'sf_yy_plus', 'val': 45724.6 / 24992.1, 'tol': 3e-4})
        val_list.append({'prop': 'sf_yy_minus', 'val': 45724.6 / 80979.0, 'tol': 2e-4})
        val_list.append({'prop': 'sf_11_plus', 'val': 121030 / 97751.9, 'tol': 2e-4})
        val_list.append({'prop': 'sf_11_minus', 'val': 121030 / 69403.3, 'tol': 2e-4})
        val_list.append({'prop': 'sf_22_plus', 'val': 43760.6 / 27959.0, 'tol': 2e-4})
        val_list.append({'prop': 'sf_22_minus', 'val': 43760.6 / 20761.6, 'tol': 3e-4})

        self.validate_properties(val_list, section)

    def test_warping(self, section, benchmark):
        benchmark(section.calculate_warping_properties)

        val_list = []
        val_list.append({'prop': 'j', 'val': 135333, 'tol': 1e-3})
        val_list.append({'prop': 'gamma', 'val': 1.62288e8, 'tol': 5e-4})
        val_list.append({'prop': 'A_s11', 'val': 885.444, 'tol': 2e-4})
        val_list.append({'prop': 'A_s22', 'val': 1459.72, 'tol': 4e-4})
        val_list.append({'prop': 'x11_se', 'val': 28.719, 'tol': 1e-3})
        val_list.append({'prop': 'y22_se', 'val': 35.2348, 'tol': 5e-4})

        self.validate_properties(val_list, section)
