.. _label-section-library:

====================================================
Creating Section Geometries from the Section Library
====================================================

In order to make your life easier, there are a number of built-in functions that generate typical
structural cross-sections, resulting in :class:`~sectionproperties.pre.geometry.Geometry` objects.
These typical cross-sections reside in the ``sectionproperties.pre.library`` module.


.. _label-primitive-library:

Primitive Sections Library
==========================

Rectangular Section
-------------------
..  autofunction:: sectionproperties.pre.library.primitive_sections.rectangular_section
    :noindex:

Circular Section
----------------
..  autofunction:: sectionproperties.pre.library.primitive_sections.circular_section
    :noindex:

Circular Section By Area
------------------------
..  autofunction:: sectionproperties.pre.library.primitive_sections.circular_section_by_area
    :noindex:

Elliptical Section
------------------
..  autofunction:: sectionproperties.pre.library.primitive_sections.elliptical_section
    :noindex:

Triangular Section
------------------
..  autofunction:: sectionproperties.pre.library.primitive_sections.triangular_section
    :noindex:

Triangular Radius Section
-------------------------
..  autofunction:: sectionproperties.pre.library.primitive_sections.triangular_radius_section
    :noindex:

Cruciform Section
-----------------
  ..  autofunction:: sectionproperties.pre.library.primitive_sections.cruciform_section
      :noindex:


Steel Sections Library
======================

Circular Hollow Section (CHS)
-----------------------------
..  autofunction:: sectionproperties.pre.library.steel_sections.circular_hollow_section
    :noindex:

Elliptical Hollow Section (EHS)
-------------------------------
..  autofunction:: sectionproperties.pre.library.steel_sections.elliptical_hollow_section
    :noindex:

Rectangular Hollow Section (RHS)
--------------------------------
..  autofunction:: sectionproperties.pre.library.steel_sections.rectangular_hollow_section
    :noindex:

Polygon Hollow Section
----------------------
  ..  autofunction:: sectionproperties.pre.library.steel_sections.polygon_hollow_section
      :noindex:

I Section
---------
  ..  autofunction:: sectionproperties.pre.library.steel_sections.i_section
      :noindex:

Monosymmetric I Section
-----------------------
  ..  autofunction:: sectionproperties.pre.library.steel_sections.mono_i_section
      :noindex:

Tapered Flange I Section
------------------------
  ..  autofunction:: sectionproperties.pre.library.steel_sections.tapered_flange_i_section
      :noindex:

Parallel Flange Channel (PFC) Section
-------------------------------------
  ..  autofunction:: sectionproperties.pre.library.steel_sections.channel_section
      :noindex:

Tapered Flange Channel Section
------------------------------
  ..  autofunction:: sectionproperties.pre.library.steel_sections.tapered_flange_channel
      :noindex:

Tee Section
-----------
  ..  autofunction:: sectionproperties.pre.library.steel_sections.tee_section
      :noindex:

Angle Section
-------------
  ..  autofunction:: sectionproperties.pre.library.steel_sections.angle_section
      :noindex:

Cee Section
-----------
  ..  autofunction:: sectionproperties.pre.library.steel_sections.cee_section
      :noindex:

Zed Section
-----------
  ..  autofunction:: sectionproperties.pre.library.steel_sections.zed_section
      :noindex:

Box Girder Section
------------------
  ..  autofunction:: sectionproperties.pre.library.steel_sections.box_girder_section
      :noindex:

Bulb Section
------------
  ..  autofunction:: sectionproperties.pre.library.steel_sections.bulb_section
      :noindex:


.. _label-concrete-library:

Concrete Sections Library
=========================

Concrete Rectangular Section
----------------------------
..  autofunction:: sectionproperties.pre.library.concrete_sections.concrete_rectangular_section
    :noindex:

Concrete Column Section
-----------------------
..  autofunction:: sectionproperties.pre.library.concrete_sections.concrete_column_section
    :noindex:

Concrete Tee Section
--------------------
..  autofunction:: sectionproperties.pre.library.concrete_sections.concrete_tee_section
    :noindex:

Concrete Circular Section
-------------------------
..  autofunction:: sectionproperties.pre.library.concrete_sections.concrete_circular_section
    :noindex:

Add Bar
-------
..  autofunction:: sectionproperties.pre.library.concrete_sections.add_bar
    :noindex:


.. _label-bridge-library:

Bridge Sections Library
=======================

Super Tee Girder Section
------------------------
..  autofunction:: sectionproperties.pre.library.bridge_sections.super_t_girder_section
    :noindex:

I Girder Section
----------------
..  autofunction:: sectionproperties.pre.library.bridge_sections.i_girder_section
    :noindex:


Nastran Sections Library
========================

See :ref:`label-nastran-sections`.
