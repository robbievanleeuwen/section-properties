---
title: "sectionproperties: A Python package for the analysis of arbitrary cross-sections using the finite element method"
tags:
  - python
  - computational mechanics
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
  - name: Connor Ferster
    orcid: 0009-0005-0861-2428
    affiliation: 2
affiliations:
  - name: Independent Researcher, Australia
    index: 1
  - name: Independent Researcher, Canada # update
    index: 2
date: 13 October 2023
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
pre-processor, analysis engine, and post-processor, in a single open source and
accessible package, that can be used by researchers, practising engineers, and students.

# Statement of Need

Obtaining the geometric properties of simple shapes is a classical engineering problem
with well-defined analytical solutions. However, obtaining warping properties, e.g. for
torsion and shear analyses, involves solving partial differential equations [@pilkey].
While some analytical solutions exist for a small subset of geometries, the method for
obtaining these results is not able to be generalised to shapes commonly used in
engineering pratice. Further, the analysis of arbitrary composite geometries, in which a
cross-section could consist of any shape with any number of internal holes, and be made
from any number of materials, complicates both geometric and warping computations.

To the best of our knowledge, there is no open source software available for the
computation of both geometric and warping propreties for composite, arbitary
cross-sections. While there are several commercial solutions available, e.g.
[`RSECTION 1`](https://www.dlubal.com/en/products/cross-section-properties-software/rsection),
[`ShapeDesigner SaaS`](http://mechatools.com/en/shapedesigner.html), or
[`CADRE Profiler`](https://www.cadreanalytic.com/profiler.htm), none of these are
open source or provide an application programming interface (API) that would enable
these programs to be used for research. As a result, `sectionproperties` supports both
engineering practice and research, by implementing an open source solution to the
complex modelling problem that is arbitrary composite geometric and warping analyses.

# Implementation

`sectionproperties` harnesses the power of Shapely [@shapely] to streamline geometry
generation, and triangle [@triangle] (a python port of Triangle [@shewchuck]) to produce
a triangular mesh of six-noded quadratic elements. The finite element method is used to
solve for the geometric and warping properties, the latter involving the solution of
partial differential equations and boundary value problems [@pilkey]. For example, the
Saint-Venant torsion constant ($J$) is obtained by solving for the warping function,
$\omega$ [@pilkey]:

$$
\nabla^2 \omega = 0
$$

subject to the boundary condition:

$$
\frac{\partial \omega}{\partial x} n_x + \frac{\partial \omega}{\partial y} n_y = y n_x - x n_y
$$

Using the finite element method, this problem is reduced to a set of linear equations of
the form:

$$
\textbf{K} \boldsymbol{\omega} = \textbf{F}
$$

where the stiffness matrix and load vector at the element level are defined as:

$$
\textbf{k}^e = \sum_{i=1}^6 w_i \textbf{B}_i^{\rm T} \textbf{B}_i J_i
$$

$$
\textbf{f}^e = \sum_{i=1}^6 w_i \textbf{B}_i^{\rm T}
\begin{bmatrix}
  \textbf{N}_i \textbf{y}_e \\
  -\textbf{N}_i \textbf{x}_e \\
\end{bmatrix} J_i
$$

In the above, $\textbf{N}$ and $\textbf{B}$ are the shape functions and their
derivatives, and $w_i$ and $J_i$ are the weights and Jacobians of the current
integration point. The boundary conditions neccesitate the inversion of a nearly
singular global stiffness matrix. As such, the Lagrangian multiplier method is used to
solve the set of linear equations of the form $\textbf{K} \textbf{u} = \textbf{F}$ by
introducing an extra constraint on the solution vector, whereby the mean value is equal
to zero [@larson].

$$
\begin{bmatrix}
  \textbf{K} & \textbf{C}^{\rm{T}} \\
  \textbf{C} & 0 \\
\end{bmatrix}
\begin{bmatrix}
  \textbf{u} \\
  \lambda \\
\end{bmatrix} =
\begin{bmatrix}
  \textbf{F} \\
  0 \\
\end{bmatrix}
$$

where $\textbf{C}$ is the assembly of $\sum_{i} w_i \textbf{N}_i J_i$, and $\lambda$
may be though of as a relatively small force acting to enforce the constraints. Once the
warping function has been evaluated, the Saint-Venant torsion constant can be calculated
as follows:

$$
J = I_{xx} + I_{yy} - \boldsymbol{\omega}^{\rm T} \textbf{K} \boldsymbol{\omega}
$$

The calculation of plastic properties is meshless, and is conducted using an iterative
method to enforce plastic equilibrium, yielding the plastic centroids. A full
description of the theoretical background underpinning `sectionproperties` can
be found in the
[documentation](https://sectionproperties.readthedocs.io/en/stable/user_guide/theory.html).

An example of some of the visualisation generated by `sectionproperties` can be seen in
\autoref{fig:example} below.

![Plot of the centroids and torsion stress distribution for a bulb-section modelled in `sectionproperties`.\label{fig:example}](figures/example.png)

# Software Development

The `sectionproperties` package is available on [GitHub](https://github.com/robbievanleeuwen/section-properties),
where the source code, issue tracker, CI workflow, and discussion board can be found.
Pre-commit hooks are used to ensure code quality and style is consistent across all
contributions. There is an extensive testing and validation suite used to ensure that
the output produced by `sectionproperties` is verified and repeatable, including a set
of benchmarking tests. `sectionproperties` has an actively maintained and complete
[documentation](https://sectionproperties.readthedocs.io), including
[installation instructions](https://sectionproperties.readthedocs.io/en/stable/installation.html),
a detailed [user guide](https://sectionproperties.readthedocs.io/en/stable/user_guide.html),
a list of [examples](https://sectionproperties.readthedocs.io/en/stable/examples.html),
and an [API reference](https://sectionproperties.readthedocs.io/en/stable/api.html).

# Conclusion

In this paper we have described `sectionproperties`, a Python package that calculates
the section properties of arbitrary sections. It is our hope that this project is used
by researchers and practising engineers to improve their experimental and analysis
workflows.

# Acknowledgements

We acknowledge the contributions from all the
[contributors](https://github.com/robbievanleeuwen/section-properties/graphs/contributors)
to `sectionproperties`.

# References
