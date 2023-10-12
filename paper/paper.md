---
title: "sectionproperties: A Python package for the analysis of arbitrary cross-sections using the finite element method"
tags:
  - Python
  - finite element method
  - cross-section
  - stress analysis
  - engineering
authors:
  - given-names: Robbie
    non-dropping-particle: van
    surname: Leeuwen
    orcid: 0009-0004-8056-3977
    affiliation: 1
    corresponding: true
  - name: Connor Ferster
    orcid: 0000-0000-0000-0000 # update
    affiliation: 2
affiliations:
  - name: Independent Researcher, Australia
    index: 1
  - name: Institution Name, Country
    index: 2
date: 12 October 2023
bibliography: paper.bib
---

# Summary

Properties of plane cross-sections are often required in engineering research, analysis,
and design. For example, cross-sectional properties are used to determine the
displacements, natural frequencies, and stresses within beams under complex loading.
`sectionproperties` is a Python package for the analysis of _arbitrary_ cross-sections
using the finite element method. `sectionproperties` can be used to determine
geometric and warping properties, as well as visualising cross-sectional stresses
resulting from combinations of applied loads. `sectionproperties` aims to provide a
pre-processor, analysis engine, and post-processor in a single open source and
accessible package, that can be used by researchers, practising engineers, and students.

# Statement of Need

Obtaining the geometric properties of simple shapes is a classical engineering problem
with well-defined analytical solutions. However, obtaining warping properties (e.g. for
torsion and shear analyses) involves solving partial differential equations [@pilkey].
While some analytical solutions exist for a small subset of geometries, the method for
obtaining these results is not able to be generalised to shapes commonly used in
engineering pratice. Further, the analysis of arbitrary composite geometries, in which a
cross-section could consist of any shape with any number of internal holes, and be made
up of any number of materials, complicates both geometric and warping computations.

To the best of our knowledge, there is no open source software available for the
computation of both geometric and warping propreties for composite, arbitary
cross-sections. While there are several commercial programs available, e.g.
[`RSECTION 1`](https://www.dlubal.com/en/products/cross-section-properties-software/rsection),
[`ShapeDesigner SaaS`](http://mechatools.com/en/shapedesigner.html), or
[`CADRE Profiler`](https://www.cadreanalytic.com/profiler.htm), none of these are
open source or provide an application programming interface (API) that would enable
these programs to be used for research. As a result, `sectionproperties` supports both
engineering practice and research, by implementing an open source solution to the
complex modelling problem of composite, arbitrary geometric and warping analyses.

IMPLEMENTATION

# Availability

The `sectionproperties` package is available on [GitHub](https://github.com/robbievanleeuwen/section-properties),
where the source code, issue tracker, CI workflow, and discussion board can be found.
The [documentation](https://sectionproperties.readthedocs.io) includes
[installation instructions](https://sectionproperties.readthedocs.io/en/stable/installation.html),
a detailed [user guide](https://sectionproperties.readthedocs.io/en/stable/user_guide.html),
a list of [examples](https://sectionproperties.readthedocs.io/en/stable/examples.html),
and an [API reference](https://sectionproperties.readthedocs.io/en/stable/api.html).

# Acknowledgements

We acknowledge the contributions from all the
[contributors](https://github.com/robbievanleeuwen/section-properties/graphs/contributors)
to `sectionproperties`.

# References
