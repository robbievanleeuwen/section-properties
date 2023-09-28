.. _label-geom_mesh:

====================================================
Creating Geometries, Meshes, and Material Properties
====================================================

Before performing a cross-section analysis, the geometry of the cross-section and a finite element
mesh must be created. Optionally, material properties can be applied to different regions of the
cross-section. If materials are not applied, then a default material is used to report the
geometric properties.


Section Geometry
================

**New in** ``v2.0.0`` - There are two types of geometry objects in sectionproperties:

* The :class:`~sectionproperties.pre.geometry.Geometry` class, for section geometries with a
  single, contiguous region
* The :class:`~sectionproperties.pre.geometry.CompoundGeometry` class, comprised of two or more
  `Geometry` objects


Geometry Class
--------------

..  autoclass:: sectionproperties.pre.geometry.Geometry
   :noindex:

A :class:`~sectionproperties.pre.geometry.Geometry` is instantiated with two arguments:

#. A :class:`shapely.geometry.Polygon` object
#. An optional :class:`~sectionproperties.pre.pre.Material` object

.. note::
   If a :class:`~sectionproperties.pre.pre.Material` is not given, then the default material is
   assigned to the ``Geometry.material`` attribute. The default material has an elastic modulus of
   1, a Poisson's ratio of 0 and a yield strength of 1.


CompoundGeometry Class
----------------------

..  autoclass:: sectionproperties.pre.geometry.CompoundGeometry
   :noindex:

A :class:`~sectionproperties.pre.geometry.CompoundGeometry` is instantiated with **one argument**
which can be one of **two types**:

#. A list of :class:`~sectionproperties.pre.geometry.Geometry` objects; or
#. A :class:`shapely.geometry.MultiPolygon`

.. note::
   A :class:`~sectionproperties.pre.geometry.CompoundGeometry` does not have a ``.material``
   attribute and therefore, a :class:`~sectionproperties.pre.pre.Material` cannot be assigned to a
   :class:`~sectionproperties.pre.geometry.CompoundGeometry` directly. Since a
   :class:`~sectionproperties.pre.geometry.CompoundGeometry` is simply a combination of
   :class:`~sectionproperties.pre.geometry.Geometry` objects, the Material(s) should be assigned to
   the individual :class:`~sectionproperties.pre.geometry.Geometry` objects that comprise the
   :class:`~sectionproperties.pre.geometry.CompoundGeometry`.
   :class:`~sectionproperties.pre.geometry.CompoundGeometry` objects created from a
   :class:`shapely.geometry.MultiPolygon` will have its constituent
   :class:`~sectionproperties.pre.geometry.Geometry` objects assigned the default material.


Defining Material Properties
----------------------------

Materials are defined in *sectionproperties* by creating a
:class:`~sectionproperties.pre.pre.Material` object:

..  autoclass:: sectionproperties.pre.pre.Material
    :noindex:

Each :class:`~sectionproperties.pre.geometry.Geometry` contains its own material definition,
which is stored in the ``.material`` attribute. A geometry's material may be altered at any time by
simply assigning a new :class:`~sectionproperties.pre.pre.Material` to the ``.material`` attribute.

.. warning::

   To understand how material properties affect the results reported by *sectionproperties*, see
   :ref:`label-material-results`.


Creating Section Geometries
===========================

In addition to creating geometries directly from :class:`shapely.geometry.Polygon` and/or
:class:`shapely.geometry.MultiPolygon` objects, there are other ways to create geometries for
analysis:

#. From lists of points, facets, hole regions, and control regions
#. From ``.dxf`` files
#. From ``.3dm`` files
#. From Rhino encodings
#. Using transformation methods on existing geometries and/or by applying set operations
   (e.g. `|`, `+`, `-`, `&`, `^`)
#. From *sectionproperties*'s section library

For the first two approaches, an optional ``.material`` parameter can be passed containing a
:class:`~sectionproperties.pre.pre.Material` (or list of `Material` objects) to associate with the
newly created geometry(ies). The material attribute can be altered afterward in a
:class:`~sectionproperties.pre.geometry.Geometry` object at any time by simply assigning a
different :class:`~sectionproperties.pre.pre.Material` to the ``.material`` attribute.


.. _label-from-points:

Geometry from points, facets, holes, and control points
-------------------------------------------------------

In sectionproperties ``v1.x.x``, geometries were created by specifying lists of `points`, `facets`,
`holes`, and `control_points`. This functionality has been preserved in ``v2.0.0`` by using the
:attr:`~sectionproperties.pre.geometry.Geometry.from_points()` class method.

..  autofunction:: sectionproperties.pre.geometry.Geometry.from_points
   :noindex:

For simple geometries (i.e. single-region shapes without holes), if the points are an ordered sequence of coordinates, only the `points`
argument is required (`facets`, `holes`, and `control_points` are optional). If the geometry has holes, then all arguments are required.

If the geometry has multiple regions, then the
:attr:`~sectionproperties.pre.geometry.CompoundGeometry.from_points()` class method must be used.

..  autofunction:: sectionproperties.pre.geometry.CompoundGeometry.from_points
   :noindex:

See :ref:`ref_ex_custom` for an example of this implementation.


.. _label-geometry-dxf:

Geometry from .dxf Files
------------------------

Geometries can now be created from ``.dxf`` files using the
:attr:`sectionproperties.pre.geometry.Geometry.from_dxf()` method. The returned geometry will
either be a :class:`~sectionproperties.pre.geometry.Geometry` or
:class:`~sectionproperties.pre.geometry.CompoundGeometry` object depending on the geometry in the
file (i.e. the number of contiguous regions).

..  autofunction:: sectionproperties.pre.geometry.Geometry.from_dxf
   :noindex:

..  autofunction:: sectionproperties.pre.geometry.CompoundGeometry.from_dxf
   :noindex:


.. _label-geometry-3dm:

Geometry from Rhino
-------------------

Geometries can now be created from .3dm files and BREP encodings.
Various limitations and assumptions need to be acknowledged:

* sectional analysis is based in 2d and Rhino is a 3d environment.
* the recognised Rhino geometries are limited to planer-single-surfaced BREPs.
* Rhino uses NURBS for surface boundaries and *sectionproperties* uses piecewise linear boundaries.
* a search plane is defined.

See the keyword arguments below that are used to search and simplify the Rhino geometry.

Rhino files are read via the class methods :func:`sectionproperties.pre.geometry.Geometry.from_3dm()` and
:func:`sectionproperties.pre.geometry.CompoundGeometry.from_3dm()`.
Each class method returns the respective objects.

..  autofunction:: sectionproperties.pre.geometry.Geometry.from_3dm
   :noindex:

..  autofunction:: sectionproperties.pre.geometry.CompoundGeometry.from_3dm
   :noindex:

Geometry objects can also be created from encodings of Rhino BREP.

..  autofunction:: sectionproperties.pre.geometry.Geometry.from_rhino_encoding
   :noindex:

More advanced filtering can be achieved by working with the Shapely geometries directly.
These can be accessed by :func:`~sectionproperties.pre.rhino.load_3dm()` and
:func:`~sectionproperties.pre.rhino.load_rhino_brep_encoding()`.

..  autofunction:: sectionproperties.pre.rhino.load_3dm
   :noindex:

..  autofunction:: sectionproperties.pre.rhino.load_brep_encoding
   :noindex:

.. note::

  Dependencies for importing files from rhino are not included by default. To obtain
  the required dependencies install *sectionproperties* with the rhino option:

  .. code-block:: console

    pip install sectionproperties[rhino]


.. _label-geometry-set:

Combining Geometries Using Set Operations
-----------------------------------------

Both `Geometry` and `CompoundGeometry` objects can be manipulated using Python's set operators:

- ``|``  Bitwise OR - Performs a union on the two geometries
- ``-``  Bitwise DIFFERENCE - Performs a subtraction, subtracting the second geometry from the
  first
- ``&``  Bitwise AND - Performs an intersection operation, returning the regions of geometry common
  to both
- ``^``  Bitwise XOR - Performs a symmetric difference operation, returning the regions of geometry
  that are not overlapping
- ``+``  Addition - Combines two geometries into a `CompoundGeometry`

See :ref:`label-advanced_geom` for an example using set operations.


*sectionproperties* Section Library
-----------------------------------

See :ref:`label-section-library`.


Manipulating Geometries
=======================

Each geometry instance is able to be manipulated in 2D space for the purpose of creating novel, custom section geometries
that the user may require.

.. note::
   Operations on geometries are **non-destructive**. For each operation, a new geometry object is
   returned.

   This gives *sectionproperties* geoemtries a *fluent API* meaning that transformation methods can
   be chained together. Please see :ref:`label-advanced_geom` for examples.


Aligning
--------

  ..  autofunction:: sectionproperties.pre.geometry.Geometry.align_center
      :noindex:

  ..  autofunction:: sectionproperties.pre.geometry.Geometry.align_to
      :noindex:


Mirroring
---------

  ..  autofunction:: sectionproperties.pre.geometry.Geometry.mirror_section
      :noindex:


Rotating
--------

  ..  autofunction:: sectionproperties.pre.geometry.Geometry.rotate_section
      :noindex:


Shifting
--------

  ..  autofunction:: sectionproperties.pre.geometry.Geometry.shift_points
      :noindex:

  ..  autofunction:: sectionproperties.pre.geometry.Geometry.shift_section
      :noindex:


Splitting
---------

  ..  autofunction:: sectionproperties.pre.geometry.Geometry.split_section
      :noindex:


Offsetting the Perimeter
------------------------

  ..  autofunction:: sectionproperties.pre.geometry.Geometry.offset_perimeter
      :noindex:

  ..  autofunction:: sectionproperties.pre.geometry.CompoundGeometry.offset_perimeter
      :noindex:


Visualising the Geometry
========================

Visualisation of geometry objects is best performed in the Jupyter computing environment,
however, most visualisation can also be done in any environment which supports display of
matplotlib plots.

There are generally two ways to visualise geometry objects:

#. In the Jupyter computing environment, geometry objects utilise their underlying
   ``shapely.geometry.Polygon`` object's ``_repr_svg_`` method to show the geometry
   as it's own representation.
#. By using the :func:`~sectionproperties.pre.geometry.Geometry.plot_geometry()` method

..  automethod:: sectionproperties.pre.geometry.Geometry.plot_geometry
    :noindex:

.. note::
   You can also use ``.plot_geometry()`` with ``CompoundGeometry`` objects


Generating a Mesh
=================

A finite element mesh is required to perform a cross-section analysis. After a geometry has been
created, a finite element mesh can then be created for the geometry by using the
:attr:`~sectionproperties.pre.geometry.Geometry.create_mesh()` method:

..  automethod:: sectionproperties.pre.geometry.Geometry.create_mesh
    :noindex:

..  warning:: The length of ``mesh_sizes`` must match the number of regions
  in the geometry object.

Once the mesh has been created, it is stored within the geometry object and the geometry object
can then be passed to :class:`~sectionproperties.analysis.section.Section` for analysis.

Please see :ref:`label-analysis` for further information on performing analyses.
