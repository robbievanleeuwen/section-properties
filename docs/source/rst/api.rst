Python API Reference
====================

Pre-Processor Package
---------------------

.. _label-sections-module:

*geometry* Module
^^^^^^^^^^^^^^^^^

Geometry Class
""""""""""""""
..  autoclass:: sectionproperties.pre.geometry.Geometry
    :members:

CompoundGeometry Class
""""""""""""""""""""""
..  autoclass:: sectionproperties.pre.geometry.CompoundGeometry
    :show-inheritance:

load_dxf
""""""""
..  autofunction:: sectionproperties.pre.geometry.load_dxf

create_facets
"""""""""""""
..  autofunction:: sectionproperties.pre.geometry.create_facets

create_exterior_points
""""""""""""""""""""""
..  autofunction:: sectionproperties.pre.geometry.create_exterior_points

create_interior_points
""""""""""""""""""""""
..  autofunction:: sectionproperties.pre.geometry.create_interior_points

create_points_and_facets
""""""""""""""""""""""""
..  autofunction:: sectionproperties.pre.geometry.create_points_and_facets


*pre* Module
^^^^^^^^^^^^

Material Class
""""""""""""""

  ..  autoclass:: sectionproperties.pre.pre.Material
      :show-inheritance:
      :members:

create_mesh
"""""""""""
..  autofunction:: sectionproperties.pre.pre.create_mesh


*rhino* Module
^^^^^^^^^^^^^^

load_3dm
""""""""
..  autofunction:: sectionproperties.pre.rhino.load_3dm

load_brep_encoding
""""""""""""""""""
..  autofunction:: sectionproperties.pre.rhino.load_brep_encoding


*bisect_section* Module
^^^^^^^^^^^^^^^^^^^^^^^

create_line_segment
"""""""""""""""""""
..  autofunction:: sectionproperties.pre.bisect_section.create_line_segment

group_top_and_bottom_polys
""""""""""""""""""""""""""
..  autofunction:: sectionproperties.pre.bisect_section.group_top_and_bottom_polys

line_mx_plus_b
""""""""""""""
..  autofunction:: sectionproperties.pre.bisect_section.line_mx_plus_b

perp_mx_plus_b
""""""""""""""
..  autofunction:: sectionproperties.pre.bisect_section.perp_mx_plus_b

line_intersection
"""""""""""""""""
..  autofunction:: sectionproperties.pre.bisect_section.line_intersection

sum_poly_areas
""""""""""""""
..  autofunction:: sectionproperties.pre.bisect_section.sum_poly_areas


*primitive_sections* Module
^^^^^^^^^^^^^^^^^^^^^^^^^^^

rectangular_section
"""""""""""""""""""
..  autofunction:: sectionproperties.pre.library.primitive_sections.rectangular_section

circular_section
""""""""""""""""
..  autofunction:: sectionproperties.pre.library.primitive_sections.circular_section

circular_section_by_area
""""""""""""""""""""""""
..  autofunction:: sectionproperties.pre.library.primitive_sections.circular_section_by_area

elliptical_section
""""""""""""""""""
..  autofunction:: sectionproperties.pre.library.primitive_sections.elliptical_section

triangular_section
""""""""""""""""""
..  autofunction:: sectionproperties.pre.library.primitive_sections.triangular_section

triangular_radius_section
"""""""""""""""""""""""""
..  autofunction:: sectionproperties.pre.library.primitive_sections.triangular_radius_section

cruciform_section
""""""""""""""""""
  ..  autofunction:: sectionproperties.pre.library.primitive_sections.cruciform_section


*steel_sections* Module
^^^^^^^^^^^^^^^^^^^^^^^

circular_hollow_section
"""""""""""""""""""""""
..  autofunction:: sectionproperties.pre.library.steel_sections.circular_hollow_section

elliptical_hollow_section
"""""""""""""""""""""""""
..  autofunction:: sectionproperties.pre.library.steel_sections.elliptical_hollow_section

rectangular_hollow_section
""""""""""""""""""""""""""
..  autofunction:: sectionproperties.pre.library.steel_sections.rectangular_hollow_section

polygon_hollow_section
"""""""""""""""""""""""
  ..  autofunction:: sectionproperties.pre.library.steel_sections.polygon_hollow_section

i_section
"""""""""
  ..  autofunction:: sectionproperties.pre.library.steel_sections.i_section

mono_i_section
""""""""""""""
  ..  autofunction:: sectionproperties.pre.library.steel_sections.mono_i_section

tapered_flange_i_section
""""""""""""""""""""""""
  ..  autofunction:: sectionproperties.pre.library.steel_sections.tapered_flange_i_section

channel_section
"""""""""""""""
  ..  autofunction:: sectionproperties.pre.library.steel_sections.channel_section

tapered_flange_channel
""""""""""""""""""""""
  ..  autofunction:: sectionproperties.pre.library.steel_sections.tapered_flange_channel

tee_section
"""""""""""
  ..  autofunction:: sectionproperties.pre.library.steel_sections.tee_section

angle_section
""""""""""""""
  ..  autofunction:: sectionproperties.pre.library.steel_sections.angle_section

cee_section
""""""""""""
  ..  autofunction:: sectionproperties.pre.library.steel_sections.cee_section

zed_section
"""""""""""
  ..  autofunction:: sectionproperties.pre.library.steel_sections.zed_section

box_girder_section
"""""""""""""""""""
  ..  autofunction:: sectionproperties.pre.library.steel_sections.box_girder_section

bulb_section
""""""""""""
  ..  autofunction:: sectionproperties.pre.library.steel_sections.bulb_section


*concrete_sections* Module
^^^^^^^^^^^^^^^^^^^^^^^^^^

concrete_rectangular_section
""""""""""""""""""""""""""""
..  autofunction:: sectionproperties.pre.library.concrete_sections.concrete_rectangular_section

concrete_column_section
"""""""""""""""""""""""
..  autofunction:: sectionproperties.pre.library.concrete_sections.concrete_column_section

concrete_tee_section
""""""""""""""""""""
..  autofunction:: sectionproperties.pre.library.concrete_sections.concrete_tee_section

concrete_circular_section
"""""""""""""""""""""""""
..  autofunction:: sectionproperties.pre.library.concrete_sections.concrete_circular_section

add_bar
"""""""
..  autofunction:: sectionproperties.pre.library.concrete_sections.add_bar


*bridge_sections* Module
^^^^^^^^^^^^^^^^^^^^^^^^

super_t_girder_section
""""""""""""""""""""""
..  autofunction:: sectionproperties.pre.library.bridge_sections.super_t_girder_section

i_girder_section
""""""""""""""""
..  autofunction:: sectionproperties.pre.library.bridge_sections.i_girder_section

get_super_t_girder_dims
"""""""""""""""""""""""
..  autofunction:: sectionproperties.pre.library.bridge_sections.get_super_t_girder_dims

get_i_girder_dims
"""""""""""""""""""""""
..  autofunction:: sectionproperties.pre.library.bridge_sections.get_i_girder_dims


.. _label-nastran-sections:

*nastran_sections* Module
^^^^^^^^^^^^^^^^^^^^^^^^^
This module contains sections as defined by Nastran and Nastran-based programs,
such as MYSTRAN and ASTROS.

nastran_bar
"""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_bar

nastran_box
"""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_box

nastran_box1
""""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_box1

nastran_chan
""""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_chan

nastran_chan1
"""""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_chan1

nastran_chan2
"""""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_chan2

nastran_cross
"""""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_cross

nastran_dbox
""""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_dbox

nastran_fcross
""""""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_fcross

nastran_gbox
""""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_gbox

nastran_h
""""""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_h

nastran_hat
"""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_hat

nastran_hat1
""""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_hat1

nastran_hexa
""""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_hexa

nastran_i
"""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_i

nastran_i1
""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_i1

nastran_l
"""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_l

nastran_rod
"""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_rod

nastran_tee
"""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_tee

nastran_tee1
""""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_tee1

nastran_tee2
""""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_tee2

nastran_tube
""""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_tube

nastran_tube2
"""""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_tube2

nastran_zed
"""""""""""
..  autofunction:: sectionproperties.pre.library.nastran_sections.nastran_zed

Nastran References
""""""""""""""""""
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

*section* Module
^^^^^^^^^^^^^^^^^^^^^^

Section Class
"""""""""""""

  ..  autoclass:: sectionproperties.analysis.section.Section
      :show-inheritance:
      :members:

PlasticSection Class
""""""""""""""""""""

  ..  autoclass:: sectionproperties.analysis.section.PlasticSection
      :show-inheritance:
      :members:

StressPost Class
""""""""""""""""

..  autoclass:: sectionproperties.analysis.section.StressPost
    :show-inheritance:
    :members:

MaterialGroup Class
"""""""""""""""""""

  ..  autoclass:: sectionproperties.analysis.section.MaterialGroup
      :show-inheritance:
      :members:

StressResult Class
""""""""""""""""""

..  autoclass:: sectionproperties.analysis.section.StressResult
    :show-inheritance:
    :members:

SectionProperties Class
"""""""""""""""""""""""

..  autoclass:: sectionproperties.analysis.section.SectionProperties
    :show-inheritance:
    :members:

*fea* Module
^^^^^^^^^^^^

Tri6 Class
""""""""""

..  autoclass:: sectionproperties.analysis.fea.Tri6
    :show-inheritance:
    :members:

gauss_points
""""""""""""
..  autofunction:: sectionproperties.analysis.fea.gauss_points

shape_function
""""""""""""""
..  autofunction:: sectionproperties.analysis.fea.shape_function

extrapolate_to_nodes
""""""""""""""""""""
..  autofunction:: sectionproperties.analysis.fea.extrapolate_to_nodes

principal_coordinate
""""""""""""""""""""
..  autofunction:: sectionproperties.analysis.fea.principal_coordinate

global_coordinate
"""""""""""""""""
..  autofunction:: sectionproperties.analysis.fea.global_coordinate

point_above_line
""""""""""""""""
..  autofunction:: sectionproperties.analysis.fea.point_above_line


*solver* Module
^^^^^^^^^^^^^^^

solve_cgs
"""""""""
..  autofunction:: sectionproperties.analysis.solver.solve_cgs

solve_cgs_lagrange
""""""""""""""""""
..  autofunction:: sectionproperties.analysis.solver.solve_cgs_lagrange

solve_direct
""""""""""""
..  autofunction:: sectionproperties.analysis.solver.solve_direct

solve_direct_lagrange
"""""""""""""""""""""
..  autofunction:: sectionproperties.analysis.solver.solve_direct_lagrange


Post-Processor Package
----------------------

*post* Module
^^^^^^^^^^^^^

plotting_context
""""""""""""""""
..  autofunction:: sectionproperties.post.post.plotting_context

draw_principal_axis
"""""""""""""""""""
..  autofunction:: sectionproperties.post.post.draw_principal_axis

print_results
"""""""""""""
..  autofunction:: sectionproperties.post.post.print_results
