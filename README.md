# sectionproperties
a python package for the analysis of arbitrary cross-sections using the finite element method.

## Current Capabilities:

### Pre-Processor:
- [x] Python API
- [x] Custom section geometry input
- [x] Common section geometry generators
- [x] Multiple geometry merging
- [ ] JSON input file
- [ ] .dxf import
- [x] Quadratic triangular mesh generation
- [ ] Composite material properties

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
- [ ] Cross-section stresses

### Solver:
- [x] Direct solver
- [x] CGS iterative solver
- [x] Sparse matrices

### Post-Processor:
- [x] Plot geometry
- [x] Plot mesh
- [x] Plot centroids
- [ ] Plot cross-section stresses
- [ ] Generate cross-section report
- [ ] Export to Paraview

### Additional Modules:
- [ ] Optimisation
- [ ] Reinforced Concrete
- [ ] Steel
