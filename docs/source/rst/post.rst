.. _label-post:

Viewing the Results
===================

.. _label-print-properties:

Printing a List of the Section Properties
-----------------------------------------

A list of section properties that have been calculated by various analyses can
be printed to the terminal using the :func:`~sectionproperties.analysis.section.Section.display_results`
method that belongs to every
:class:`~sectionproperties.analysis.section.Section` object.

..  automethod:: sectionproperties.analysis.section.Section.display_results
    :noindex:

.. _label-get-methods:

Getting Specific Section Properties
-----------------------------------

Alternatively, there are a number of methods that can be called on the
:class:`~sectionproperties.analysis.section.Section` object to return
a specific section property:

Section Area
^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.section.Section.get_area
    :noindex:

Section Perimeter
^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.section.Section.get_perimeter
    :noindex:

Section Mass
^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.section.Section.get_mass
    :noindex:

Axial Rigidity
^^^^^^^^^^^^^^

If material properties have been specified, returns the axial rigidity of the
section.

..  automethod:: sectionproperties.analysis.section.Section.get_ea
    :noindex:

First Moments of Area
^^^^^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.section.Section.get_q
    :noindex:

Second Moments of Area
^^^^^^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.section.Section.get_ig
    :noindex:

..  automethod:: sectionproperties.analysis.section.Section.get_ic
    :noindex:

..  automethod:: sectionproperties.analysis.section.Section.get_ip
    :noindex:

Elastic Centroid
^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.section.Section.get_c
    :noindex:


Section Moduli
^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.section.Section.get_z
    :noindex:

..  automethod:: sectionproperties.analysis.section.Section.get_zp
    :noindex:

Radii of Gyration
^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.section.Section.get_rc
    :noindex:

..  automethod:: sectionproperties.analysis.section.Section.get_rp
    :noindex:


Principal Axis Angle
^^^^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.section.Section.get_phi
    :noindex:

Effective Material Properties
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.section.Section.get_e_eff
    :noindex:

..  automethod:: sectionproperties.analysis.section.Section.get_g_eff
    :noindex:

..  automethod:: sectionproperties.analysis.section.Section.get_nu_eff
    :noindex:


Torsion Constant
^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.section.Section.get_j
    :noindex:

Shear Centre
^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.section.Section.get_sc
    :noindex:

..  automethod:: sectionproperties.analysis.section.Section.get_sc_p
    :noindex:

Trefftz's Shear Centre
^^^^^^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.section.Section.get_sc_t
    :noindex:

Warping Constant
^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.section.Section.get_gamma
    :noindex:

Shear Area
^^^^^^^^^^

..  automethod:: sectionproperties.analysis.section.Section.get_As
    :noindex:

..  automethod:: sectionproperties.analysis.section.Section.get_As_p
    :noindex:

Monosymmetry Constants
^^^^^^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.section.Section.get_beta
    :noindex:

..  automethod:: sectionproperties.analysis.section.Section.get_beta_p
    :noindex:

Plastic Centroid
^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.section.Section.get_pc
    :noindex:

..  automethod:: sectionproperties.analysis.section.Section.get_pc_p
    :noindex:

Plastic Section Moduli
^^^^^^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.section.Section.get_s
    :noindex:

..  automethod:: sectionproperties.analysis.section.Section.get_sp
    :noindex:


Shape Factors
^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.section.Section.get_sf
    :noindex:

..  automethod:: sectionproperties.analysis.section.Section.get_sf_p
    :noindex:


.. _label-material-results:

How Material Properties Affect Results
--------------------------------------

If a :class:`~sectionproperties.pre.geometry.Geometry` containing a user defined
:class:`~sectioproperties.pre.pre.Material` is used to build a
:class:`~sectionproperties.analysis.section.Section`, *sectionproperties* will assume you
are performing a **composite analysis** and this will affect the way some of the results are
stored and presented.

In general, the calculation of gross composite section properties takes into account the elastic
modulus, Poisson's ratio and yield strength of each material in the section. Unlike many design
codes, *sectionproperties* is 'material property agnostic' and does not transform sections based on
a defined material property, e.g. in reinforced concrete analysis it is commonplace to transform
the reinforcing steel area based on the ratio between the elastic moduli,
:math:`n = E_{steel} / E_{conc}`. *sectionproperties* instead calculates the gross material
weighted properties, which is analogous to transforming with respect to a material property with
elastic modulus, :math:`E = 1`.

Using the example of a reinforced concrete section, *sectionproperties* will calculate the gross
section bending stiffness, :math:`(EI)_g`, rather than an effective concrete second moment of area,
:math:`I_{c,eff}`:

.. math::
  (EI)_g = E_s \times I_s + E_c \times I_c

If the user wanted to obtain the effective concrete second moment of area for a code calculation,
they could simply divide the gross bending stiffness by the elastic modulus for concrete:

.. math::
  I_{c,eff} = \frac{(EI)_g}{E_c}

With reference to the ``get`` methods described in :ref:`label-print-properties`, a
**composite analysis** will modify the following properties:

* First moments of area :func:`~sectionproperties.analysis.Section.get_q` - returns elastic
  modulus weighted first moments of area :math:`E.Q`
* Second moments of area :func:`~sectionproperties.analysis.Section.get_ig`,
  :func:`~sectionproperties.analysis.Section.get_ic`,
  :func:`~sectionproperties.analysis.Section.get_ip` - return elastic modulus weighted second
  moments of area :math:`E.I`
* Section moduli :func:`~sectionproperties.analysis.Section.get_z`,
  :func:`~sectionproperties.analysis.Section.get_zp` - return elastic modulus weighted section
  moduli :math:`E.Z`
* Torsion constant :func:`~sectionproperties.analysis.Section.get_j` - returns elastic
  modulus weighted torsion constant :math:`E.J`
* Warping constant :func:`~sectionproperties.analysis.Section.get_gamma` - returns elastic
  modulus weighted warping constant :math:`E.\Gamma`
* Shear areas :func:`~sectionproperties.analysis.Section.get_As`,
  :func:`~sectionproperties.analysis.Section.get_As_p` - return elastic modulus weighted shear
  areas :math:`E.A_s`
* Plastic section moduli :func:`~sectionproperties.analysis.Section.get_s`,
  :func:`~sectionproperties.analysis.Section.get_sp` - return yield strength weighted plastic
  section moduli, i.e. plastic moments :math:`M_p = f_y.S`

A **composite analysis** will also enable the user to retrieve effective gross section
area-weighted material properties:

* Effective elastic modulus :math:`E_{eff}` - :func:`~sectionproperties.analysis.Section.get_e_eff`
* Effective shear modulus :math:`G_{eff}` - :func:`~sectionproperties.analysis.Section.get_g_eff`
* Effective Poisson's ratio :math:`\nu_{eff}` -
  :func:`~sectionproperties.analysis.Section.get_nu_eff`

These values may be used to transform composite properties output by *sectionproperties* for
practical use, e.g. to calculate torsional rigidity:

.. math::
  (GJ)_g = \frac{G_{eff}}{E_{eff}} (EJ)_g

For further information, see the theoretical background to the calculation of
:ref:`label-theory-composite`.


Section Property Centroids Plots
--------------------------------

A plot of the centroids (elastic, plastic and shear centre) can be produced with
the finite element mesh in the background:

..  automethod:: sectionproperties.analysis.section.Section.plot_centroids
    :noindex:


Plotting Section Stresses
-------------------------

There are a number of methods that can be called from a :class:`~sectionproperties.analysis.section.StressResult`
object to plot the various cross-section stresses. These methods take the following form:

  :class:`~sectionproperties.analysis.section.StressResult`.plot_(*stress/vector*)_(*action*)_(*stresstype*)

where:

- *stress* denotes a contour plot and *vector* denotes a vector plot.
- *action* denotes the type of action causing the stress e.g. *mxx* for bending moment about the x-axis. Note that the action is omitted for stresses caused by the application of all actions.
- *stresstype* denotes the type of stress that is being plotted e.g. *zx* for the *x*-component of shear stress.

The examples shown in the methods below are performed on a 150x90x12 UA
(unequal angle) section. The :class:`~sectionproperties.analysis.section.Section`
object is created below::

  import sectionproperties.pre.library.steel_sections as steel_sections
  from sectionproperties.analysis.section import Section

  geometry = steel_sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
  mesh = geometry.create_mesh(mesh_sizes=[2.5])
  section = Section(geometry, mesh)

Primary Stress Plots
^^^^^^^^^^^^^^^^^^^^

Axial Stress (:math:`\sigma_{zz,N}`)
""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_n_zz
    :noindex:

Bending Stress (:math:`\sigma_{zz,Mxx}`)
""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_mxx_zz
    :noindex:

Bending Stress (:math:`\sigma_{zz,Myy}`)
""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_myy_zz
    :noindex:

Bending Stress (:math:`\sigma_{zz,M11}`)
""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_m11_zz
    :noindex:

Bending Stress (:math:`\sigma_{zz,M22}`)
""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_m22_zz
    :noindex:

Bending Stress (:math:`\sigma_{zz,\Sigma M}`)
"""""""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_m_zz
    :noindex:

Torsion Stress (:math:`\sigma_{zx,Mzz}`)
""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_mzz_zx
    :noindex:

Torsion Stress (:math:`\sigma_{zy,Mzz}`)
""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_mzz_zy
    :noindex:

Torsion Stress (:math:`\sigma_{zxy,Mzz}`)
"""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_mzz_zxy
    :noindex:

..  automethod:: sectionproperties.analysis.section.StressPost.plot_vector_mzz_zxy
    :noindex:

Shear Stress (:math:`\sigma_{zx,Vx}`)
"""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_vx_zx
    :noindex:

Shear Stress (:math:`\sigma_{zy,Vx}`)
"""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_vx_zy
    :noindex:

Shear Stress (:math:`\sigma_{zxy,Vx}`)
""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_vx_zxy
    :noindex:

..  automethod:: sectionproperties.analysis.section.StressPost.plot_vector_vx_zxy
    :noindex:

Shear Stress (:math:`\sigma_{zx,Vy}`)
"""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_vy_zx
    :noindex:

Shear Stress (:math:`\sigma_{zy,Vy}`)
"""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_vy_zy
    :noindex:

Shear Stress (:math:`\sigma_{zxy,Vy}`)
""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_vy_zxy
    :noindex:

..  automethod:: sectionproperties.analysis.section.StressPost.plot_vector_vy_zxy
    :noindex:

Shear Stress (:math:`\sigma_{zx,\Sigma V}`)
"""""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_v_zx
    :noindex:

Shear Stress (:math:`\sigma_{zy,\Sigma V}`)
"""""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_v_zy
    :noindex:

Shear Stress (:math:`\sigma_{zxy,\Sigma V}`)
""""""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_v_zxy
    :noindex:

..  automethod:: sectionproperties.analysis.section.StressPost.plot_vector_v_zxy
    :noindex:

Combined Stress Plots
^^^^^^^^^^^^^^^^^^^^^

Normal Stress (:math:`\sigma_{zz}`)
"""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_zz
    :noindex:

Shear Stress (:math:`\sigma_{zx}`)
""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_zx
    :noindex:

Shear Stress (:math:`\sigma_{zy}`)
""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_zy
    :noindex:

Shear Stress (:math:`\sigma_{zxy}`)
"""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_zxy
    :noindex:

..  automethod:: sectionproperties.analysis.section.StressPost.plot_vector_zxy
    :noindex:

Major Principal Stress (:math:`\sigma_{1}`)
"""""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_1
    :noindex:

Minor Principal Stress (:math:`\sigma_{3}`)
"""""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_3
    :noindex:

von Mises Stress (:math:`\sigma_{vM}`)
"""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_stress_vm
    :noindex:

Mohr's Circles for Stresses at a Point
""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.section.StressPost.plot_mohrs_circles
    :noindex:

Retrieving Section Stress
-------------------------

All cross-section stresses can be recovered using the :func:`~sectionproperties.analysis.section.StressPost.get_stress`
method that belongs to every
:class:`~sectionproperties.analysis.section.StressPost` object:

..  automethod:: sectionproperties.analysis.section.StressPost.get_stress
    :noindex:
