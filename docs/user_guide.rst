User Guide
==========

The ``sectionproperties`` user guide provides an overview of each step in the
``sectionproperties`` workflow, details the theoretical background for the analysis and
includes a section on validation.

.. toctree::
    :caption: Contents
    :maxdepth: 1

    user_guide/overview
    user_guide/materials
    user_guide/geometry
    user_guide/meshing
    user_guide/analysis
    user_guide/results
    user_guide/theory
    user_guide/validation

.. _label-features:

Features
--------

Pre-Processor
^^^^^^^^^^^^^

* ☑ Python API
* ☑ Geometry manipulation by Shapely
* ☑ Common section geometry functions
* ☑ Custom section geometry input
* ☑ Rhino .3dm import
* ☑ .dxf import
* ☑ Perimeter offset tool
* ☑ Quadratic triangular mesh generation
* ☑ Composite material definition

Cross-Section Analysis
^^^^^^^^^^^^^^^^^^^^^^

* ☑ Global axis geometric section properties

  * ☑ Area
  * ☑ Perimeter
  * ☑ Mass
  * ☑ First moments of area
  * ☑ Second moments of area
  * ☑ Elastic centroid

* ☑ Centroidal axis geometric section properties

  * ☑ Second moments of area
  * ☑ Elastic section moduli
  * ☐ Yield moment
  * ☑ Radii of gyration
  * ☑ Plastic centroid
  * ☑ Plastic section moduli
  * ☑ Shape factors

* ☑ Principal axis geometric section properties

  * ☑ Second moments of area
  * ☑ Elastic section moduli
  * ☐ Yield moment
  * ☑ Radii of gyration
  * ☑ Plastic centroid
  * ☑ Plastic section moduli
  * ☑ Shape factors

* ☑ Warping section properties

  * ☑ Torsion constant
  * ☑ Warping constant
  * ☑ Monosymmetry constants

* ☑ Shear section properties

  * ☑ Shear centre (elastic method)
  * ☑ Shear centre (Trefftz's method)
  * ☑ Shear areas (global axis)
  * ☑ Shear areas (principal axis)

* ☑ Cross-section stress analysis
* ☑ Mohr's circles for stresses at a point

Solver
^^^^^^

* ☑ Direct solver
* ☑ CGS iterative solver
* ☑ Sparse matrices

Post-Processor
^^^^^^^^^^^^^^

* ☑ Plot geometry
* ☑ Plot mesh
* ☑ Plot centroids
* ☑ Plot cross-section stresses
* ☑ Retrieve cross-section stresses
* ☐ Generate cross-section report
