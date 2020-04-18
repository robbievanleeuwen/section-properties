.. _label-post:

Viewing the Results
===================

Printing a List of the Section Properties
-----------------------------------------

A list of section properties that have been calculated by various analyses can
be printed to the terminal using the :func:`~sectionproperties.analysis.cross_section.CrossSection.display_results`
method that belongs to every
:class:`~sectionproperties.analysis.cross_section.CrossSection` object.

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.display_results
    :noindex:

Getting Specific Section Properties
-----------------------------------

Alternatively, there are a number of methods that can be called on the
:class:`~sectionproperties.analysis.cross_section.CrossSection` object to return
a specific section property:

Cross-Section Area
^^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_area
    :noindex:

Cross-Section Perimeter
^^^^^^^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_perimeter
    :noindex:

Axial Rigidity
^^^^^^^^^^^^^^

If material properties have been specified, returns the axial rigidity of the
section.

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_ea
    :noindex:

First Moments of Area
^^^^^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_q
    :noindex:

Second Moments of Area
^^^^^^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_ig
    :noindex:

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_ic
    :noindex:

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_ip
    :noindex:

Elastic Centroid
^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_c
    :noindex:


Section Moduli
^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_z
    :noindex:

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_zp
    :noindex:

Radii of Gyration
^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_rc
    :noindex:

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_rp
    :noindex:


Principal Axis Angle
^^^^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_phi
    :noindex:


Torsion Constant
^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_j
    :noindex:

Shear Centre
^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_sc
    :noindex:

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_sc_p
    :noindex:

Trefftz's Shear Centre
^^^^^^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_sc_t
    :noindex:

Warping Constant
^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_gamma
    :noindex:

Shear Area
^^^^^^^^^^

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_As
    :noindex:

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_As_p
    :noindex:

Monosymmetry Constants
^^^^^^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_beta
    :noindex:

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_beta_p
    :noindex:

Plastic Centroid
^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_pc
    :noindex:

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_pc_p
    :noindex:

Plastic Section Moduli
^^^^^^^^^^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_s
    :noindex:

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_sp
    :noindex:


Shape Factors
^^^^^^^^^^^^^

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_sf
    :noindex:

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.get_sf_p
    :noindex:


Section Property Centroids Plots
--------------------------------

A plot of the centroids (elastic, plastic and shear centre) can be produced with
the finite element mesh in the background:

..  automethod:: sectionproperties.analysis.cross_section.CrossSection.plot_centroids
    :noindex:


Plotting Cross-Section Stresses
-------------------------------

There are a number of methods that can be called from a :class:`~sectionproperties.analysis.cross_section.StressResult`
object to plot the various cross-section stresses. These methods take the following form:

  :class:`~sectionproperties.analysis.cross_section.StressResult`.plot_(*stress/vector*)_(*action*)_(*stresstype*)

where:

- *stress* denotes a contour plot and *vector* denotes a vector plot.
- *action* denotes the type of action causing the stress e.g. *mxx* for bending moment about the x-axis. Note that the action is omitted for stresses caused by the application of all actions.
- *stresstype* denotes the type of stress that is being plotted e.g. *zx* for the *x*-component of shear stress.

The examples shown in the methods below are performed on a 150x90x12 UA
(unequal angle) section. The :class:`~sectionproperties.analysis.cross_section.CrossSection`
object is created below::

  import sectionproperties.pre.sections as sections
  from sectionproperties.analysis.cross_section import CrossSection

  geometry = sections.AngleSection(d=150, b=90, t=12, r_r=10, r_t=5, n_r=8)
  mesh = geometry.create_mesh(mesh_sizes=[2.5])
  section = CrossSection(geometry, mesh)

Primary Stress Plots
^^^^^^^^^^^^^^^^^^^^

Axial Stress (:math:`\sigma_{zz,N}`)
""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_n_zz
    :noindex:

Bending Stress (:math:`\sigma_{zz,Mxx}`)
""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_mxx_zz
    :noindex:

Bending Stress (:math:`\sigma_{zz,Myy}`)
""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_myy_zz
    :noindex:

Bending Stress (:math:`\sigma_{zz,M11}`)
""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_m11_zz
    :noindex:

Bending Stress (:math:`\sigma_{zz,M22}`)
""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_m22_zz
    :noindex:

Bending Stress (:math:`\sigma_{zz,\Sigma M}`)
"""""""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_m_zz
    :noindex:

Torsion Stress (:math:`\sigma_{zx,Mzz}`)
""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_mzz_zx
    :noindex:

Torsion Stress (:math:`\sigma_{zy,Mzz}`)
""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_mzz_zy
    :noindex:

Torsion Stress (:math:`\sigma_{zxy,Mzz}`)
"""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_mzz_zxy
    :noindex:

..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_vector_mzz_zxy
    :noindex:

Shear Stress (:math:`\sigma_{zx,Vx}`)
"""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_vx_zx
    :noindex:

Shear Stress (:math:`\sigma_{zy,Vx}`)
"""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_vx_zy
    :noindex:

Shear Stress (:math:`\sigma_{zxy,Vx}`)
""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_vx_zxy
    :noindex:

..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_vector_vx_zxy
    :noindex:

Shear Stress (:math:`\sigma_{zx,Vy}`)
"""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_vy_zx
    :noindex:

Shear Stress (:math:`\sigma_{zy,Vy}`)
"""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_vy_zy
    :noindex:

Shear Stress (:math:`\sigma_{zxy,Vy}`)
""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_vy_zxy
    :noindex:

..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_vector_vy_zxy
    :noindex:

Shear Stress (:math:`\sigma_{zx,\Sigma V}`)
"""""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_v_zx
    :noindex:

Shear Stress (:math:`\sigma_{zy,\Sigma V}`)
"""""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_v_zy
    :noindex:

Shear Stress (:math:`\sigma_{zxy,\Sigma V}`)
""""""""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_v_zxy
    :noindex:

..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_vector_v_zxy
    :noindex:

Combined Stress Plots
^^^^^^^^^^^^^^^^^^^^^

Normal Stress (:math:`\sigma_{zz}`)
""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_zz
    :noindex:

Shear Stress (:math:`\sigma_{zx}`)
"""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_zx
    :noindex:

Shear Stress (:math:`\sigma_{zy}`)
"""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_zy
    :noindex:

Shear Stress (:math:`\sigma_{zxy}`)
""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_zxy
    :noindex:

..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_vector_zxy
    :noindex:

von Mises Stress (:math:`\sigma_{vM}`)
"""""""""""""""""""""""""""""""""""""""
..  automethod:: sectionproperties.analysis.cross_section.StressPost.plot_stress_vm
    :noindex:


Retrieving Cross-Section Stress
-------------------------------

All cross-section stresses can be recovered using the :func:`~sectionproperties.analysis.cross_section.StressPost.get_stress`
method that belongs to every
:class:`~sectionproperties.analysis.cross_section.StressPost` object:

..  automethod:: sectionproperties.analysis.cross_section.StressPost.get_stress
    :noindex:
