![alt text](logo.png "sectionproperties")

[![Build Status](https://travis-ci.com/robbievanleeuwen/section-properties.svg?branch=master)](https://travis-ci.com/robbievanleeuwen/section-properties) [![Documentation Status](https://readthedocs.org/projects/sectionproperties/badge/?version=latest)](https://sectionproperties.readthedocs.io/en/latest/?badge=latest)

A python package for the analysis of arbitrary cross-sections using the finite element method written by Robbie van Leeuwen. *sectionproperties* can be used to determine section properties to be used in structural design and visualise cross-sectional stresses resulting from combinations of applied forces and bending moments.

[Subscribe](http://eepurl.com/dMMUeg) to the mailing list!

## Installation:

For more detailed installation instructions, refer to the [documentation](https://sectionproperties.readthedocs.io/).

### UNIX (MacOS/Linux):

```
$ pip install sectionproperties
```

### Windows

Install *meshpy* by downloading the appropriate [installation wheel](https://www.lfd.uci.edu/~gohlke/pythonlibs/#meshpy).

Navigate to the location of the downloaded wheel and install using pip:

```
$ cd Downloads
$ pip install MeshPy‑2018.2.1‑cp36‑cp36m‑win_amd64.whl
```

Once *meshpy* has been installed, *sectionproperties* can be installed:

```
$ pip install sectionproperties
```

## Documentation:

*sectionproperties* has a fully documented python API which you can find at [https://sectionproperties.readthedocs.io/](https://sectionproperties.readthedocs.io/). To read more about the theory behind the program, its implementation and some more examples, check out my blog at [https://robbievanleeuwen.github.io/](https://robbievanleeuwen.github.io/).

## Current Capabilities:

### Pre-Processor:
- [x] Python API
- [x] Custom section geometry input
- [x] Common section geometry generators
- [x] Multiple geometry merging
- [x] Perimeter offset tool
- [x] Geometry cleaning
- [ ] JSON input file
- [ ] .dxf import
- [x] Quadratic triangular mesh generation
- [x] Composite material properties

### Cross-Section Analysis:
- [x] Global axis geometric section properties:
  - [x] Area
  - [x] Perimeter
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
- [x] Shear section properties:
  - [x] Shear centre (elastic method)
  - [x] Shear centre (Trefftz's method)
  - [x] Shear areas (global axis)
  - [x] Shear areas (principal axis)
- [x] Cross-section stresses

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
- [ ] Export to Paraview

### Additional Modules:
- [ ] Optimisation
- [ ] Reinforced Concrete
- [ ] Steel

## Disclaimer:

*sectionproperties* is an open source engineering tool that continues to benefit from the collaboration of many contributors. Although efforts have been made to ensure the that relevant engineering theories have been correctly implemented, it remains the user's responsibility to confirm and accept the output. Refer to the [license](LICENSE) for clarification of the conditions of use.
