<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/_static/logo-dark-mode.png">
  <source media="(prefers-color-scheme: light)" srcset="docs/_static/logo-light-mode.png">
  <img alt="sectionproperties logo" src="docs/_static/logo-light-mode.png">
</picture>

[![PyPI](https://img.shields.io/pypi/v/sectionproperties.svg)][pypi_]
[![Status](https://img.shields.io/pypi/status/sectionproperties.svg)][status]
[![Python Version](https://img.shields.io/pypi/pyversions/sectionproperties)][python version]
[![License](https://img.shields.io/pypi/l/sectionproperties)][license]
[![Read the documentation at https://sectionproperties.readthedocs.io/](https://img.shields.io/readthedocs/sectionproperties/latest.svg?label=Read%20the%20Docs)][read the docs]
[![Tests](https://github.com/robbievanleeuwen/section-properties/workflows/Tests/badge.svg)][tests]
[![Codecov](https://codecov.io/gh/robbievanleeuwen/section-properties/branch/master/graph/badge.svg)][codecov]
[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)][pre-commit]
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)][black]

[pypi_]: https://pypi.org/project/sectionproperties/
[status]: https://pypi.org/project/sectionproperties/
[python version]: https://pypi.org/project/sectionproperties
[read the docs]: https://sectionproperties.readthedocs.io/
[tests]: https://github.com/robbievanleeuwen/section-properties/actions?workflow=Tests
[codecov]: https://app.codecov.io/gh/robbievanleeuwen/section-properties
[pre-commit]: https://github.com/pre-commit/pre-commit
[black]: https://github.com/psf/black

A python package for the analysis of arbitrary cross-sections using the finite element method. _sectionproperties_ can be used to determine section properties to be used in structural design and visualise cross-sectional stresses resulting from combinations of applied forces and bending moments.

[Subscribe](http://eepurl.com/dMMUeg) to the mailing list!

## Installation

You can install _sectionproperties_ via [pip] from [PyPI]:

```console
> pip install sectionproperties
```

## Documentation

_sectionproperties_ is fully documented including a user walkthrough, examples, background theory and an API guide. The documentation can found at [https://sectionproperties.readthedocs.io/](https://sectionproperties.readthedocs.io/).

## Features

### Pre-Processor

- [x] Python API
- [x] Geometry manipulation by Shapely
- [x] Common section geometry functions
- [x] Custom section geometry input
- [x] Rhino .3dm import
- [x] .dxf import
- [x] Perimeter offset tool
- [x] Quadratic triangular mesh generation
- [x] Composite material definition

### Cross-Section Analysis

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

### Solver

- [x] Direct solver
- [x] CGS iterative solver
- [x] Sparse matrices

### Post-Processor

- [x] Plot geometry
- [x] Plot mesh
- [x] Plot centroids
- [x] Plot cross-section stresses
- [x] Retrieve cross-section stresses
- [ ] Generate cross-section report

## Contributing

Contributions are very welcome.
To learn more, see the [Contributor Guide].

## License

Distributed under the terms of the [MIT license][license],
_sectionproperties_ is free and open source software.

## Issues

If you encounter any problems,
please [file an issue] along with a detailed description.

## Disclaimer

_sectionproperties_ is an open source engineering tool that continues to benefit from the collaboration of many contributors. Although efforts have been made to ensure the that relevant engineering theories have been correctly implemented, it remains the user's responsibility to confirm and accept the output. Refer to the [license](LICENSE.md) for clarification of the conditions of use.

## Credits

This project was generated from [@cjolowicz]'s [Hypermodern Python Cookiecutter] template.

[@cjolowicz]: https://github.com/cjolowicz
[pypi]: https://pypi.org/
[hypermodern python cookiecutter]: https://github.com/cjolowicz/cookiecutter-hypermodern-python
[file an issue]: https://github.com/robbievanleeuwen/section-properties/issues
[pip]: https://pip.pypa.io/

<!-- github-only -->

[license]: https://github.com/robbievanleeuwen/section-properties/blob/master/LICENSE
[contributor guide]: https://github.com/robbievanleeuwen/section-properties/blob/master/CONTRIBUTING.md
