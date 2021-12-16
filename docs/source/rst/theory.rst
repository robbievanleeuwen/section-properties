Theoretical Background
======================

Introduction
------------

intro

Mesh Generation
---------------

mesh


Finite Element Preliminaries
----------------------------

Element Type
^^^^^^^^^^^^

element

Isoparametric Representation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

iso

Shape Functions
"""""""""""""""

shape

Cartesian Partial Derivatives
"""""""""""""""""""""""""""""

cartesian

Numerical Integration
^^^^^^^^^^^^^^^^^^^^^

int

Extrapolation to Nodes
^^^^^^^^^^^^^^^^^^^^^^

extrap

Lagrangian Multiplier
^^^^^^^^^^^^^^^^^^^^^

lagrange

Calculation of Cross-Section Properties
---------------------------------------

Cross-Sectional Area
^^^^^^^^^^^^^^^^^^^^

area

First Moments of Area
^^^^^^^^^^^^^^^^^^^^^

area

Centroids
^^^^^^^^^

centroid

Second Moments of Area
^^^^^^^^^^^^^^^^^^^^^^

area

Radii of Gyration
^^^^^^^^^^^^^^^^^

radii


Elastic Section Moduli
^^^^^^^^^^^^^^^^^^^^^^

moduli

Plastic Section Moduli
^^^^^^^^^^^^^^^^^^^^^^

moduli

Principal Axis Properties
^^^^^^^^^^^^^^^^^^^^^^^^^

prinicpal

Torsion Constant
^^^^^^^^^^^^^^^^

torsion

Shear Properties
^^^^^^^^^^^^^^^^

shear

Shear Centre
""""""""""""

shear

Shear Deformation Coefficients
""""""""""""""""""""""""""""""

shear

Warping Constant
""""""""""""""""

warp

Cross-Section Stresses
----------------------

stresses

Axial Stresses
^^^^^^^^^^^^^^

axial

Bending Stresses
^^^^^^^^^^^^^^^^

Global Axis Bending
"""""""""""""""""""

bend

Principal Axis Bending
""""""""""""""""""""""

bend

Torsion Stresses
^^^^^^^^^^^^^^^^

twist

Shear Stresses
^^^^^^^^^^^^^^

shear

Principal Stresses
^^^^^^^^^^^^^^^^^^

principal

von Mises Stresses
^^^^^^^^^^^^^^^^^^

vm

Mohr's Circle
^^^^^^^^^^^^^

mohr

.. _label-theory-composite:

Composite Cross-Sections
------------------------

composite

..  Mention that Poisson's ratios should be relatively close as if the Poisson's
    ratio is largely variable, the basic contention that sig_x = sig_y = sig_xy = 0
    ceases to be applicable.

..  In this program an effective Poisson's ratio is calculated by determining a
    weighted E and G and then deriving a nu that is effective for the entire
    cross-section.
