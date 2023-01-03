![Logo Dark](docs/source/images/logo_dark.png#gh-dark-mode-only)
![Logo Light](docs/source/images/logo.png#gh-light-mode-only)

[![Run Tests](https://github.com/robbievanleeuwen/section-properties/actions/workflows/tests.yml/badge.svg)](https://github.com/robbievanleeuwen/section-properties/actions/workflows/tests.yml) [![Lint with Black](https://github.com/robbievanleeuwen/section-properties/actions/workflows/black.yml/badge.svg)](https://github.com/robbievanleeuwen/section-properties/actions/workflows/black.yml) [![Documentation Status](https://readthedocs.org/projects/sectionproperties/badge/?version=latest)](https://sectionproperties.readthedocs.io/en/latest/?badge=latest)
 [![codecov](https://codecov.io/gh/robbievanleeuwen/section-properties/branch/master/graph/badge.svg?token=QCH9J4SG6P)](https://codecov.io/gh/robbievanleeuwen/section-properties) [![PyPI version](https://badge.fury.io/py/sectionproperties.svg)](https://badge.fury.io/py/sectionproperties) [![Python versions](https://img.shields.io/badge/python-3.8%20%7C%203.9%20%7C%203.10-blue?style=flat&logo=python)](https://badge.fury.io/py/sectionproperties) [![GitHub license](https://img.shields.io/github/license/robbievanleeuwen/section-properties)](https://github.com/robbievanleeuwen/section-properties/blob/master/LICENSE.md)

A python package for the analysis of arbitrary cross-sections using the finite element method. *sectionproperties* can be used to determine section properties to be used in structural design and visualise cross-sectional stresses resulting from combinations of applied forces and bending moments.

[Subscribe](http://eepurl.com/dMMUeg) to the mailing list!

## Installation:

For more detailed installation instructions, refer to the [documentation](https://sectionproperties.readthedocs.io/en/latest/rst/installation.html).

```
$ pip install sectionproperties
```

## Documentation:

*sectionproperties* is fully documented including a user walkthrough, examples, background theory and an API guide. The documentation can found at [https://sectionproperties.readthedocs.io/](https://sectionproperties.readthedocs.io/).

## Current Capabilities:

### Pre-Processor:
- [x] Python API
- [x] Geometry manipulation by Shapely
- [x] Common section geometry functions
- [x] Custom section geometry input
- [x] Rhino .3dm import
- [x] .dxf import
- [x] Perimeter offset tool
- [x] Quadratic triangular mesh generation
- [x] Composite material definition

### Cross-Section Analysis:
- [x] Global axis geometric section properties:
  - [x] Area
  - [x] Perimeter
  - [x] Mass
  - [x] First moments of area
  - [x] Second moments of area
  - [x] Elastic centroid
- [x] Centroidal axis geometric section properties:
  - [x] Second moments of area
  - [x] Elastic section moduli
  - [ ] Yield moment
  - [x] Radii of gyration
  - [x] Plastic centroid
  - [x] Plastic section moduli
  - [x] Shape factors
- [x] Principal axis geometric section properties:
  - [x] Second moments of area
  - [x] Elastic section moduli
  - [ ] Yield moment
  - [x] Radii of gyration
  - [x] Plastic centroid
  - [x] Plastic section moduli
  - [x] Shape factors
- [x] Warping section properties:
  - [x] Torsion constant
  - [x] Warping constant
  - [x] Monosymmetry constants
- [x] Shear section properties:
  - [x] Shear centre (elastic method)
  - [x] Shear centre (Trefftz's method)
  - [x] Shear areas (global axis)
  - [x] Shear areas (principal axis)
- [x] Cross-section stress analysis
- [x] Mohr's circles for stresses at a point

### Solver:
- [x] Direct solver
- [x] CGS iterative solver
- [x] Sparse matrices

### Post-Processor:
- [x] Plot geometry
- [x] Plot mesh
- [x] Plot centroids
- [x] Plot cross-section stresses
- [x] Retrieve cross-section stresses
- [ ] Generate cross-section report

## Disclaimer:

*sectionproperties* is an open source engineering tool that continues to benefit from the collaboration of many contributors. Although efforts have been made to ensure the that relevant engineering theories have been correctly implemented, it remains the user's responsibility to confirm and accept the output. Refer to the [license](LICENSE.md) for clarification of the conditions of use.
