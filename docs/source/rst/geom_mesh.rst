.. _label-geom_mesh:

Creating a Geometry, Mesh and Material Properties
=================================================

Before performing a cross-section analysis, the geometry of the cross-section and a finite element
mesh must be created. Optionally, material properties can be applied to different regions of the
cross-section. If materials are not applied, then a default material is used to report the geometric
properties.


Section Geometry
""""""""""""""""

There are two types of geometry objects in sectionproperties:

* The :class:`~sectionproperties.pre.sections.Geometry` class, for section geometries with a single, contiguous region
* The :class:`~sectionproperties.pre.sections.CompoundGeometry` class, for section geometries with multiple distinct regions

**New in v2.0.0**
sectionproperties geometry classes are now powered by shapely. A :class:``~sectionproperties.pre.sections.Geometry`` object
takes a shapely Polygon object. A :class:`~sectionproperties.pre.sections.CompoundGeometry` object takes a shapely MultiPolygon
object. 

In both cases, the original shapely geometry is stored in the geometry class in the `.geom` attribute.

Geometry for sectionproperties analysis can be created in several ways: 

1. From shapely Polygon and MultiPolygon objects
2. From points, facets, hole regions, and control_point regions
3. From DXF files
4. From sectionproperties "factory functions" of common shapes
5. By applying set operations (e.g. `|`, `+`, `-`, `&`, `^`) and transformation methods on existing geometries

For the first four methods, an optional `materials` parameter can be passed containing a :class:`~sectionproperties.pre.pre.Material` object
(or a list of `Material` objects) to associate with the newly created geometry(ies). The material attribute can be altered in a `Geometry`
object at any time by simply assigning a material to the `.material` attribute.

Geometry from shapely Polygon and MultiPolygon objects
------------------------------------------------------
Using shapely, any 2D shape can be created from an ordered sequence of coordinates. A Polygon is created 
from a coordinate sequence with one enclosed perimeter. A MultiPolygon can be created from multiple Polygons.

* :class:`~sectionproperties.pre.sections.Geometry` objects are instantiated with a shapely Polygon object and an optional :class:`~sectionproperties.pre.pre.Material`
* :class:`~sectionproperties.pre.sections.CompoundGeometry` objects are instantiated with a shapely MultiPolygon object -OR- a list of Polygon objects and an optional :class:`~sectionproperties.pre.pre.Material`

The geometry objects in sectionproperties make use of the shapely geometry's `_repr_svg_` method for rich display 
of geometries in Jupyter environments.


Geometry from points, facets, holes, and control_points
-------------------------------------------------------

In sectionproperties v1.x.x, geometries were created by specifying lists of `points`, `facets`, `holes`, and `control_points`.

To maintain backward compatability, geometries can be created this same way using the `Geometry.from_points()` class method
and the `CompoundGeometry.from_points()` class method. 

For simple geometries (i.e. single-region shapes without holes), if the points are an ordered sequence of coordinates, only the `points` 
argument is required (`facets`, `holes`, and `control_points` are optional). If the geometry has holes, then all arguments are required.

If the geometry has multiple regions, then `CompoundGeometry.from_points()` class method must be used.


Geometry from DXF Files
-----------------------

Geometries can now be created from DXF files using the `Geometry.from_dxf()` method. The returned geometry will either be a `Geometry`
or `CompoundGeometry` object depending on the geometry in the file (one or more contiguous regions).


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

Operations on geometries are non-destructive: for each operation, a new `Geometry` (or `CompoundGeometry`) instance is created.
No geometries are "altered in place".

Manipulating Geometries
"""""""""""""""""""""""

Each geometry instance is able to be manipulated in 2D space for the purpose of creating novel, custom section geometries
that the user may require. 

Shifting, Aligning, Mirroring, Rotating, etc.
---------------------------------------------

These manipulation methods are all non-destructive: they all return a new geometry instance
instead of altering the geometry in place.

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

Geometry objects can be visualised by using the
:func:`~sectionproperties.pre.sections.Geometry.plot_geometry` method:

..  automethod:: sectionproperties.pre.sections.Geometry.plot_geometry
    :noindex:

Generating a Mesh
-----------------

A finite element mesh is required to perform a cross-section analysis. A finite element mesh can
be created by using the :func:`~sectionproperties.pre.sections.Geometry.create_mesh` method:

..  automethod:: sectionproperties.pre.sections.Geometry.create_mesh
    :noindex:

..  warning:: The length of ``mesh_sizes`` must match the number of regions
  in the geometry object.

Defining Material Properties
----------------------------

Composite cross-sections can be analysed by specifying different material properties for each
section of the mesh. Materials are defined in *sectionproperties* by creating a
:class:`~sectionproperties.pre.pre.Material` object:

..  autoclass:: sectionproperties.pre.pre.Material
    :noindex:
    :show-inheritance:
