# sectionproperties
a python package for the analysis of arbitrary cross-sections using the finite element method written by Robbie van Leeuwen. *sectionproperties* can be used to determine section properties to be used in structural design and visualise cross-sectional stresses resulting from combinations of applied forces and bending moments.

## Documentation:

*sectionproperties* has a fully documented python API which you can find at [https://sectionproperties.readthedocs.io/](https://sectionproperties.readthedocs.io/). To read more about the theory behind the program, its implementation and some more examples, check out my blog at [https://robbievanleeuwen.github.io/](https://robbievanleeuwen.github.io/).

## Current Capabilities:

*Note: this branch of sectionproperties is still in development and thus the calculations may not be correct until the branch is merged back to the master.*

### Pre-Processor:
- [x] Python API
- [x] Custom section geometry input
- [x] Common section geometry generators
- [x] Multiple geometry merging
- [x] Geometry cleaning
- [ ] JSON input file
- [ ] .dxf import
- [x] Quadratic triangular mesh generation
- [x] Composite material properties

### Cross-Section Analysis:
- [x] Global axis geometric section properties:
  - [x] Area
  - [x] First moments of area
  - [x] Second moments of area
  - [x] Elastic centroid
- [x] Centroidal axis geometric section properties:
  - [x] Second moments of area
  - [x] Elastic section moduli
  - [x] Radii of gyration
  - [ ] Plastic centroid
  - [ ] Plastic section moduli
  - [ ] Shape factors
- [x] Principal axis geometric section properties:
  - [x] Second moments of area
  - [x] Elastic section moduli
  - [x] Radii of gyration
  - [ ] Plastic centroid
  - [ ] Plastic section moduli
  - [ ] Shape factors
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
- [ ] Generate cross-section report
- [ ] Export to Paraview

### Additional Modules:
- [ ] Optimisation
- [ ] Reinforced Concrete
- [ ] Steel

## Change Log:
xxx
