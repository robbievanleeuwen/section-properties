"""Tests for the steel sections library."""

from __future__ import annotations

import warnings

import pytest

import sectionproperties.pre.library.steel_sections as ss
from sectionproperties.analysis import Section


def test_angle_section_toe_thickness():
    """Tests the angle section when the toe thickness equals the thickness."""
    geom = ss.angle_section(d=20, b=20, t=3, r_r=3.5, r_t=3, n_r=16)
    geom.create_mesh(mesh_sizes=1.0)
    sec = Section(geom)

    # check for no additional warnings
    with pytest.warns() as record:
        warnings.warn("user", UserWarning, stacklevel=1)
        sec.calculate_geometric_properties()
        sec.calculate_warping_properties()

    assert len(record) == 1
