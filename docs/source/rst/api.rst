Python API Documentation
========================

Pre-Processor Package
---------------------

.. _label-sections-module:

*sections* Module
^^^^^^^^^^^^^^^^^

Geometry Class
""""""""""""""
..  autoclass:: sectionproperties.pre.sections.Geometry
    :members:

CompoundGeometry Class
""""""""""""""""""""""
..  autoclass:: sectionproperties.pre.sections.CompoundGeometry
    :show-inheritance:
    
*sections* Functions
^^^^^^^^^^^^^^^^^^^^^

rectangular_section
"""""""""""""""""""
..  autofunction:: sectionproperties.pre.sections.rectangular_section
    :noindex:

circular_section
""""""""""""""""
..  autofunction:: sectionproperties.pre.sections.circular_section
    :noindex:

circular_hollow_section
"""""""""""""""""""""""
..  autofunction:: sectionproperties.pre.sections.circular_hollow_section
    :noindex:

elliptical_section
""""""""""""""""""
..  autofunction:: sectionproperties.pre.sections.elliptical_section
    :noindex:

elliptical_hollow_section
"""""""""""""""""""""""""
..  autofunction:: sectionproperties.pre.sections.elliptical_hollow_section
    :noindex:

rectangular_hollow_section
""""""""""""""""""""""""""
..  autofunction:: sectionproperties.pre.sections.rectangular_hollow_section
    :noindex:

i_section
"""""""""
  ..  autofunction:: sectionproperties.pre.sections.i_section
      :noindex:

mono_i_section
""""""""""""""
  ..  autofunction:: sectionproperties.pre.sections.mono_i_section
      :noindex:

tapered_flange_i_section
""""""""""""""""""""""""
  ..  autofunction:: sectionproperties.pre.sections.tapered_flange_i_section
      :noindex:

channel_section
"""""""""""""""
  ..  autofunction:: sectionproperties.pre.sections.channel_section
      :noindex:

tapered_flange_channel
""""""""""""""""""""""
  ..  autofunction:: sectionproperties.pre.sections.tapered_flange_channel
      :noindex:

tee_section
"""""""""""
  ..  autofunction:: sectionproperties.pre.sections.tee_section
      :noindex:

angle_section
""""""""""""""
  ..  autofunction:: sectionproperties.pre.sections.angle_section
      :noindex:

cee_section
""""""""""""
  ..  autofunction:: sectionproperties.pre.sections.cee_section
      :noindex:

zed_section
"""""""""""
  ..  autofunction:: sectionproperties.pre.sections.zed_section
      :noindex:

cruciform_section
""""""""""""""""""
  ..  autofunction:: sectionproperties.pre.sections.cruciform_section
      :noindex:

polygon_hollow_section
"""""""""""""""""""""""
  ..  autofunction:: sectionproperties.pre.sections.polygon_hollow_section
      :noindex:

box_girder_section
"""""""""""""""""""
  ..  autofunction:: sectionproperties.pre.sections.box_girder_section
      :noindex:

*pre* Module
^^^^^^^^^^^^

Material Class
""""""""""""""

  ..  autoclass:: sectionproperties.pre.pre.Material
      :show-inheritance:
      :members:

pre Functions
"""""""""""""

create_mesh
""""""""""""
..  autofunction:: sectionproperties.pre.pre.create_mesh


*nastran_sections* Module
^^^^^^^^^^^^^^^^^^^^^^^^^
This module contains sections as defined by Nastran and Nastran-based programs,
such as MYSTRAN and ASTROS.

nastran_bar
""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_bar
    :noindex:

nastran_box
""""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_box
    :noindex:

nastran_box1
"""""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_box1
    :noindex:

nastran_chan
"""""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_chan
    :noindex:

nastran_chan1
""""""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_chan1
    :noindex:

nastran_chan2
""""""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_chan2
    :noindex:

nastran_cross
"""""""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_cross
    :noindex:

nastran_dbox
"""""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_dbox
    :noindex:

nastran_fcross
"""""""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_fcross
    :noindex:

nastran_gbox
"""""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_gbox
    :noindex:

nastran_h
""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_h
    :noindex:

nastran_hat
""""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_hat
    :noindex:

nastran_hat1
"""""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_hat1
    :noindex:

nastran_hexa
"""""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_hexa
    :noindex:

nastran_i
"""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_i
    :noindex:

nastran_i1
"""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_i1
    :noindex:

nastran_l
""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_l
    :noindex:

nastran_rod
""""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_rod
    :noindex:

nastran_tee
""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_tee
    :noindex:

nastran_tee1
"""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_tee1
    :noindex:

nastran_tee2
"""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_tee2
    :noindex:

nastran_tube
"""""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_tube
    :noindex:

nastran_tube2
""""""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_tube2
    :noindex:

nastran_zed
""""""""""""""
..  autofunction:: sectionproperties.pre.nastran_sections.nastran_zed
    :noindex:

References
""""""""""
.. [1]  MSC Nastran Quick Reference Guide 2012,
        PBEAML - Simple Beam Cross-Section Property, pp. 2890-2894
        https://simcompanion.mscsoftware.com/infocenter/index?page=content&id=DOC10351
.. [2]  Siemens NX Nastran 12 Quick Reference Guide,
        PBEAML, pp. 16-59 - 16-62
        https://docs.plm.automation.siemens.com/data_services/resources/nxnastran/12/help/tdoc/en_US/pdf/QRG.pdf
.. [3]  AutoDesk Nastran Online Documentation, Nastran Reference Guide,
        Section 4 - Bulk Data, PBEAML
        http://help.autodesk.com/view/NSTRN/2018/ENU/?guid=GUID-B7044BA7-3C26-49DA-9EE7-DA7505FD4B2C
.. [4]  Users Reference Manual for the MYSTRAN General Purpose Finite Element Structural Analysis Computer Program,
        Jan. 2019, Section 6.4.1.53 - PBARL, pp. 154-156
        https://www.mystran.com/Executable/MYSTRAN-Users-Manual.pdf
.. [5]  Astros Enhancements - Volume III - Astros Theoretical Manual,
        Section 5.1.3.2, pp. 56
        https://apps.dtic.mil/dtic/tr/fulltext/u2/a308134.pdf

Analysis Package
----------------

*cross_section* Module
^^^^^^^^^^^^^^^^^^^^^^

Section Class
""""""""""""""""""

  ..  autoclass:: sectionproperties.analysis.cross_section.Section
      :show-inheritance:
      :members:

PlasticSection Class
""""""""""""""""""""

  ..  autoclass:: sectionproperties.analysis.cross_section.PlasticSection
      :show-inheritance:
      :members:

StressPost Class
""""""""""""""""

..  autoclass:: sectionproperties.analysis.cross_section.StressPost
    :show-inheritance:
    :members:

MaterialGroup Class
"""""""""""""""""""

  ..  autoclass:: sectionproperties.analysis.cross_section.MaterialGroup
      :show-inheritance:
      :members:

StressResult Class
""""""""""""""""""

..  autoclass:: sectionproperties.analysis.cross_section.StressResult
    :show-inheritance:
    :members:

SectionProperties Class
"""""""""""""""""""""""

..  autoclass:: sectionproperties.analysis.cross_section.SectionProperties
    :show-inheritance:
    :members:

*fea* Module
^^^^^^^^^^^^

Tri6 Class
""""""""""

..  autoclass:: sectionproperties.analysis.fea.Tri6
    :show-inheritance:
    :members:

fea Functions
"""""""""""""

..  autofunction:: sectionproperties.analysis.fea.gauss_points
..  autofunction:: sectionproperties.analysis.fea.shape_function
..  autofunction:: sectionproperties.analysis.fea.extrapolate_to_nodes
..  autofunction:: sectionproperties.analysis.fea.principal_coordinate
..  autofunction:: sectionproperties.analysis.fea.global_coordinate
..  autofunction:: sectionproperties.analysis.fea.point_above_line

*solver* Module
^^^^^^^^^^^^^^^

solver Functions
""""""""""""""""

..  autofunction:: sectionproperties.analysis.solver.solve_cgs
..  autofunction:: sectionproperties.analysis.solver.solve_cgs_lagrange
..  autofunction:: sectionproperties.analysis.solver.solve_direct
..  autofunction:: sectionproperties.analysis.solver.solve_direct_lagrange
..  autofunction:: sectionproperties.analysis.solver.function_timer

Post-Processor Package
----------------------

*post* Module
^^^^^^^^^^^^^

post Functions
""""""""""""""

..  autofunction:: sectionproperties.post.post.setup_plot
..  autofunction:: sectionproperties.post.post.finish_plot
..  autofunction:: sectionproperties.post.post.draw_principal_axis
..  autofunction:: sectionproperties.post.post.print_results
