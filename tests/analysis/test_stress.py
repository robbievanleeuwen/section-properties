"""Tests for cross-section stress calculation."""

from __future__ import annotations

import math

import pytest

import sectionproperties.pre.library.primitive_sections as primitive_sections
from sectionproperties.analysis.section import Section

# define geometry
rect = primitive_sections.rectangular_section(b=50, d=100)
rect.create_mesh(mesh_sizes=0)  # coarse mesh
sec = Section(geometry=rect)


# test stress runtime errors
def test_stress_runtime_errors():
    """Tests for raising RuntimeErrors."""
    # check runtime error with no analysis performed
    with pytest.raises(RuntimeError):
        sec.calculate_stress()

    sec.calculate_geometric_properties()

    # check runtime errors with shear/torsion applied, no warping analysis
    with pytest.raises(RuntimeError):
        sec.calculate_stress(vx=1)

    with pytest.raises(RuntimeError):
        sec.get_stress_at_points(pts=[(10, 10)], vx=1)

    with pytest.raises(RuntimeError):
        sec.calculate_stress(vy=1)

    with pytest.raises(RuntimeError):
        sec.get_stress_at_points(pts=[(10, 10)], vy=1)

    with pytest.raises(RuntimeError):
        sec.calculate_stress(mzz=1)

    with pytest.raises(RuntimeError):
        sec.get_stress_at_points(pts=[(10, 10)], mzz=1)

    # check no runtime errors with no shear/torsion applied
    sec.calculate_stress(n=1)
    sec.calculate_stress(mxx=1)
    sec.calculate_stress(myy=1)
    sec.calculate_stress(m11=1)
    sec.calculate_stress(m22=1)

    sec.calculate_warping_properties()

    # check no runtime errors after warping analysis
    sec.calculate_stress(n=1)
    sec.calculate_stress(mxx=1)
    sec.calculate_stress(myy=1)
    sec.calculate_stress(m11=1)
    sec.calculate_stress(m22=1)
    sec.calculate_stress(vx=1)
    sec.calculate_stress(vy=1)
    sec.calculate_stress(mzz=1)
    sec.get_stress_at_points(pts=[(10, 10)], mzz=1)


def test_rectangle():
    """Tests bending stresses for a rectangle."""
    mxx = 7
    sy = 50.0 * 100.0**2 / 6.0
    sig_max = mxx / sy
    sig_0, sig_1, sig_2 = sec.get_stress_at_points(
        pts=[(25, 50), (25, 75), (25, 100)], mxx=mxx
    )
    assert sig_0 == pytest.approx((0, 0, 0))
    assert sig_1 == pytest.approx((sig_max / 2.0, 0, 0))
    assert sig_2 == pytest.approx((sig_max, 0, 0))


def test_rotated_rectangle():
    """Tests bending stresses for a rotated rectangle."""
    b = 50
    d = 100
    angle = math.atan(100 / 50)
    cx = b / 2 * math.cos(angle) - d / 2 * math.sin(angle)
    cy = b / 2 * math.sin(angle) + d / 2 * math.cos(angle)
    sy = b * d / 6.0 * cy
    mxx = 7
    sig_max = mxx / sy
    rot_rect = (
        primitive_sections.rectangular_section(b=b, d=d)
        .shift_section(x_offset=-b / 2, y_offset=-d / 2)
        .rotate_section(angle=angle, use_radians=True)
    )
    rot_rect.create_mesh(mesh_sizes=0)  # coarse mesh
    rot_sec = Section(geometry=rot_rect)
    rot_sec.calculate_geometric_properties()
    rot_sec.calculate_warping_properties()
    _, sig_1 = rot_sec.get_stress_at_points(pts=[(cx, 0), (cx, cy)], mxx=mxx)
    assert sig_1 == pytest.approx((sig_max, 0, 0))
