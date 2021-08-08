.. _label-geom_mesh:

====================================================
Creating Geometries, Meshes, and Material Properties
====================================================

Before performing a cross-section analysis, the geometry of the cross-section and a finite element
mesh must be created. Optionally, material properties can be applied to different regions of the
cross-section. If materials are not applied, then a default material is used to report the geometric
properties.


Section Geometry
================

**New in v2.0.0**
There are two types of geometry objects in sectionproperties:

* The :class:`~sectionproperties.pre.sections.Geometry` class, for section geometries with a single, contiguous region
* The :class:`~sectionproperties.pre.sections.CompoundGeometry` class, comprised of two or more `Geometry` objects

A :class:`~sectionproperties.pre.sections.Geometry` is instantiated with two arguments: 

#. A shapely.geometry.Polygon object
#. An optional :class:`~sectionproperties.pre.pre.Material` object

.. note::
   If a Material is not given, then the default material is assigned to the `Geometry.material` attribute

A :class:`~sectionproperties.pre.sections.CompoundGeometry` is instantiated with **one argument** which can be one of **two types**:
#. A list of :class:`~sectionproperties.pre.sections.Geometry`; or
#. A shapely.geometry.MultiPolygon

.. note::
   A CompoundGeometry does not have a `.material` attribute and a :class:`~sectionproperties.pre.pre.Material`
   cannot be assigned to a CompoundGeometry directly. Since a CompoundGeometry is simply a combination of Geometry objects,
   the Material(s) should be assigned to the individual Geometry objects that comprise the CompoundGeometry. CompoundGeoemtry
   objects created from a `shapely.geometry.MultiPolygon` will have its constituent Geometry objects assigned the default material.


Defining Material Properties
----------------------------

Materials are defined in *sectionproperties* by creating a :class:`~sectionproperties.pre.pre.Material` object:

..  autoclass:: sectionproperties.pre.pre.Material
    :noindex:
    :show-inheritance:

Each :class:`~sectionproperties.pre.sections.Geometry` contains its own material definition,
which is stored in the :attr:`~sectionproperties.pre.sections.Geometry.material` attribute. A geometry's material
may be altered at any time by simply assigning a new :class:`~sectionproperties.pre.pre.Material` to the ``.material`` attribute.

.. note::
   A :class:`~sectionproperties.pre.sections.CompoundGeometry` is composed of multiple :class:`~sectionproperties.pre.sections.Geometry` objects,
   each of which have their own 

Creating Section Geometries
===========================

In addition to creating geometries directly from `shapely.geometry.Polygon` and/or `shapely.geometry.MultiPolygon` objects
directly, there are other ways to create geometries for analysis:

#. From lists of points, facets, hole regions, and control regions
#. From DXF files
#. From `sectionproperties`'s "factory functions" of common shapes
#. Using transformation methods on existing geometries and/or by applying set operations (e.g. `|`, `+`, `-`, `&`, `^`)

For the first three approaches, an optional `materials` parameter can be passed containing a :class:`~sectionproperties.pre.pre.Material`
(or a list of `Material` objects) to associate with the newly created geometry(ies). The material attribute can be altered afterward in a `Geometry`
object at any time by simply assigning a different :class:`~sectionproperties.pre.pre.Material` to the `.material` attribute


Geometry from points, facets, holes, and control_points
-------------------------------------------------------

In sectionproperties v1.x.x, geometries were created by specifying lists of `points`, `facets`, `holes`, and `control_points` and
this functionality has been preserved in v2.0.0 by using `:attr:`~sectionproperties.pre.sections.Geometry.from_points()`

.. class:: 
   sectionproperties.pre.sections.Geometry 
   .. method:: from_points()

For simple geometries (i.e. single-region shapes without holes), if the points are an ordered sequence of coordinates, only the `points` 
argument is required (`facets`, `holes`, and `control_points` are optional). If the geometry has holes, then all arguments are required.

If the geometry has multiple regions, then `:attr:`~sectionproperties.pre.sections.CompoundGeometry.from_points()` class method must be used.

.. class:: 
   sectionproperties.pre.sections.CompoundGeometry 
   .. method:: from_points()


Geometry from DXF Files
-----------------------

Geometries can now be created from DXF files using the :attr:`~sectionproperties.pre.sections.Geometry.from_dxf()` method. The returned geometry will either be a `Geometry`
or `CompoundGeometry` object depending on the geometry in the file (depending on the number of contiguous regions).

.. class:: 
   sectionproperties.pre.sections.Geometry 
   .. method:: from_dxf()

.. class:: 
   sectionproperties.pre.sections.CompoundGeometry 
   .. method:: from_dxf()

Creating Common Structural Geometries from sectionproperties "factory functions"
--------------------------------------------------------------------------------

In order to make your life easier, there are a number of built-in functions that generate typical
structural cross-sections that create :class:`~sectionproperties.pre.sections.Geometry` objects.


Rectangular Section
^^^^^^^^^^^^^^^^^^^
..  autofunction:: sectionproperties.pre.sections.rectangular_section
    :noindex:

Circular Section
^^^^^^^^^^^^^^^^
..  autofunction:: sectionproperties.pre.sections.circular_section
    :noindex:

Circular Hollow Section (CHS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
..  autofunction:: sectionproperties.pre.sections.circular_hollow_section
    :noindex:

Elliptical Section
^^^^^^^^^^^^^^^^^^
..  autofunction:: sectionproperties.pre.sections.elliptical_section
    :noindex:

Elliptical Hollow Section (EHS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
..  autofunction:: sectionproperties.pre.sections.elliptical_hollow_section
    :noindex:

Rectangular Hollow Section (RHS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
..  autofunction:: sectionproperties.pre.sections.rectangular_hollow_section
    :noindex:

I-Section
^^^^^^^^^
  ..  autofunction:: sectionproperties.pre.sections.i_section
      :noindex:

Monosymmetric I-Section
^^^^^^^^^^^^^^^^^^^^^^^
  ..  autofunction:: sectionproperties.pre.sections.mono_i_section
      :noindex:

Tapered Flange I-Section
^^^^^^^^^^^^^^^^^^^^^^^^
  ..  autofunction:: sectionproperties.pre.sections.tapered_flange_i_section
      :noindex:

Parallel Flange Channel (PFC) Section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  ..  autofunction:: sectionproperties.pre.sections.channel_section
      :noindex:

Tapered Flange Channel Section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  ..  autofunction:: sectionproperties.pre.sections.tapered_flange_channel
      :noindex:

Tee Section
^^^^^^^^^^^
  ..  autofunction:: sectionproperties.pre.sections.tee_section
      :noindex:

Angle Section
^^^^^^^^^^^^^
  ..  autofunction:: sectionproperties.pre.sections.angle_section
      :noindex:

Cee Section
^^^^^^^^^^^
  ..  autofunction:: sectionproperties.pre.sections.cee_section
      :noindex:

Zed Section
^^^^^^^^^^^
  ..  autofunction:: sectionproperties.pre.sections.zed_section
      :noindex:

Cruciform Section
^^^^^^^^^^^^^^^^^
  ..  autofunction:: sectionproperties.pre.sections.cruciform_section
      :noindex:

Polygon Hollow Section
^^^^^^^^^^^^^^^^^^^^^^
  ..  autofunction:: sectionproperties.pre.sections.polygon_hollow_section
      :noindex:

Box Girder Section
^^^^^^^^^^^^^^^^^^
  ..  autofunction:: sectionproperties.pre.sections.box_girder_section
      :noindex:  

Combining Geometries using set operations `+`, `-`, `|`, `^`, `&`
-----------------------------------------------------------------

Both `Geometry` and `CompoundGeometry` objects can be manipulated using Python's set operators:
- `|`  Bitwise OR - Performs a union on the two geometries
- `-`  Bitwise DIFFERENCE - Performs a subtraction, subtracting the second geometry from the first
- `&`  Bitwise AND - Performs an intersection operation, returning the regions of geometry common to both
- `^`  Bitwise XOR - Performs a symmetric difference operation, returning the regions of geometry that are not overlapping
- `+`  Addition - Combines two geometries into a `CompoundGeometry`

Manipulating Geometries
=======================

Each geometry instance is able to be manipulated in 2D space for the purpose of creating novel, custom section geometries
that the user may require. 

.. note::   
   Operations on geometries are _non-destructive_. For each operation, a new geometry object is returned.

   This gives sectionproperties geoemtries a _fluent API_ meaning that transformation methods can be
   chained together. Please see :doc:`./advanced_geom` for examples.

Shifting, Aligning, Mirroring, Rotating, etc.
---------------------------------------------

  .. automethod:: sectionproperties.pre.sections.Geometry.align_center
     :noindex:

  .. automethod:: sectionproperties.pre.sections.Geometry.align_to
     :noindex:
   
  .. automethod:: sectionproperties.pre.sections.Geometry.mirror_section
     :noindex:

  .. automethod:: sectionproperties.pre.sections.Geometry.offset_perimeter
     :noindex:

  .. automethod:: sectionproperties.pre.sections.Geometry.rotate_section
     :noindex:
 
  .. automethod:: sectionproperties.pre.sections.Geometry.shift_points
     :noindex:

  .. automethod:: sectionproperties.pre.sections.Geometry.shift_section
     :noindex:

  .. automethod:: sectionproperties.pre.sections.Geometry.split_section
     :noindex:

Visualising the Geometry
------------------------

Visualization of geometry objects is best performed in the Jupyter computing environment,
however, most visualization can also be done in any environment which supports display of
matplotlib plots.

There are generally two ways to visualize geometry objects:

#. In the Jupyter computing environment, geometry objects utilize their underlying
   ``shapely.geometry.Polygon`` object's ``_repr_svg_`` method to show the geometry
   as it's own representation.
#. By using the :attr:`~sectionproperties.pre.sections.Geometry.plot_geometry()` method

..  automethod:: sectionproperties.pre.sections.Geometry.plot_geometry
    :noindex:

.. note::
   You can also use ``.plot_geometry()`` with ``CompoundGeometry`` objects

Generating a Mesh
-----------------

A finite element mesh is required to perform a cross-section analysis. After a geometry has been created,
a finite element mesh can then be created for the geometry by using the 
:attr:`~sectionproperties.pre.sections.Geometry.create_mesh()` method:

..  automethod:: sectionproperties.pre.sections.Geometry.create_mesh
    :noindex:

..  warning:: The length of ``mesh_sizes`` must match the number of regions
  in the geometry object.

Once the mesh has been created, it is stored within the geometry object and the geometry object
can then be passed to :class:`~sectionproperties.analysis.cross_section.Section` for analysis.

Please see :doc:`./analysis` for further information on performing analyses.
