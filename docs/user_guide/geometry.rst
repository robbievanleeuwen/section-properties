.. _label-geometry:

Geometry
========

Before performing a cross-section analysis, the geometry of the cross-section must be
created. The geometry of a cross-section defines its shape and dimensions, and provides
a way to assign material properties for composite analyses.

There are two types of geometry objects in ``sectionproperties``:

#. The :class:`~sectionproperties.pre.geometry.Geometry` class, for simple geometries
   with a single, contiguous region.

   ..  autoclass:: sectionproperties.pre.geometry.Geometry
        :noindex:

#. The :class:`~sectionproperties.pre.geometry.CompoundGeometry` class, for complex
   geometries that comprise of two or more
   :class:`~sectionproperties.pre.geometry.Geometry` objects.

   ..  autoclass:: sectionproperties.pre.geometry.CompoundGeometry
        :noindex:

Creating Geometry Objects
-------------------------

This section will outline the many ways geometry can be created in
``sectionproperties``.

Shapely Geometry
^^^^^^^^^^^^^^^^

:class:`~sectionproperties.pre.geometry.Geometry` objects can be directly instantiated
from a shapely :class:`~shapely.Polygon`.

.. automethod:: sectionproperties.pre.geometry.Geometry.__init__
    :noindex:

.. admonition:: Example

    The following example creates a :class:`~sectionproperties.pre.geometry.Geometry`
    object from an arbitrary four-sided shapely :class:`~shapely.Polygon`.

    .. plot::
        :include-source: True
        :caption: Geometry created from points, facets and holes.

        from shapely import Polygon
        from sectionproperties.pre import Geometry

        poly = Polygon([(0, 0), (5, 2), (3, 7), (1, 6)])
        geom = Geometry(geom=poly)
        geom.plot_geometry()

:class:`~sectionproperties.pre.geometry.CompoundGeometry` objects consist of multiple
:class:`~sectionproperties.pre.geometry.Geometry` objects to form a complex region.

.. automethod:: sectionproperties.pre.geometry.CompoundGeometry.__init__
    :noindex:

.. admonition:: Example

    The following example creates a
    :class:`~sectionproperties.pre.geometry.CompoundGeometry` object from two square
    :class:`~sectionproperties.pre.geometry.Geometry` objects.

    .. plot::
        :include-source: True
        :caption: Geometry created from points, facets and holes.

        from shapely import Polygon
        from sectionproperties.pre import Geometry, CompoundGeometry

        sq1 = Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])
        sq2 = Polygon([(2, 0), (6, 0), (6, 4), (2, 4)])
        geom_sq1 = Geometry(geom=sq1)
        geom_sq2 = Geometry(geom=sq2)
        geom = CompoundGeometry(geoms=[geom_sq1, geom_sq2])
        geom.plot_geometry()

Cartesian Coordinates
^^^^^^^^^^^^^^^^^^^^^

In ``sectionproperties`` ``v1``, geometries were created by specifying lists of
``points``, ``facets``, ``holes``, and ``control_points``. This functionality has been
preserved as a legacy feature through the
:meth:`sectionproperties.pre.geometry.Geometry.from_points` and
:meth:`sectionproperties.pre.geometry.CompoundGeometry.from_points` class methods.

..  automethod:: sectionproperties.pre.geometry.Geometry.from_points
   :noindex:

..  automethod:: sectionproperties.pre.geometry.CompoundGeometry.from_points
   :noindex:

CAD Files
^^^^^^^^^

Various CAD files can be imported to creating ``sectionproperties`` geometries.
``sectionproperties`` currently supports the following formats:

#. Drawing Exchange Format - ``.dxf``
#. Rhino 3D Model Format - ``.3dm``
#. Rhino BREP Encoding

.. note::
    The dependencies used to import CAD files are not included by default in the base
    installation. To install ``sectionproperties`` with CAD import functionality, use
    the ``cad`` option:

    .. code-block:: shell

        pip install sectionproperties[cad]

``.dxf``
""""""""

.. attention::
    TODO - confirm this behaviour once PR #246 is merged,
    `link <https://github.com/robbievanleeuwen/section-properties/issues/246>`_.

:class:`~sectionproperties.pre.geometry.Geometry` objects can be created from ``.dxf``
files using the
:meth:`sectionproperties.pre.geometry.Geometry.from_dxf` method. The returned geometry
will either be a :class:`~sectionproperties.pre.geometry.Geometry` or
:class:`~sectionproperties.pre.geometry.CompoundGeometry` object depending on the
geometry in the file, i.e. the number of contiguous regions.

..  automethod:: sectionproperties.pre.geometry.Geometry.from_dxf
   :noindex:

.. admonition:: Example

    TODO

Rhino
"""""

:class:`~sectionproperties.pre.geometry.Geometry` objects can be created from ``.3dm``
files and BREP encodings. Various limitations and assumptions need to be acknowledged:

* Cross-section analysis is 2D and Rhino is a 3D environment.
* The recognised Rhino geometries are limited to planer-single-surfaced BREPs.
* Rhino uses NURBS for surface boundaries and ``sectionproperties`` uses piecewise
  linear boundaries.
* A search plane is defined.

See the keyword arguments below that are used to search and simplify the Rhino geometry.

Rhino files are read via the class methods
:func:`sectionproperties.pre.geometry.Geometry.from_3dm` and
:func:`sectionproperties.pre.geometry.CompoundGeometry.from_3dm`. Each class method
returns the respective objects.

..  automethod:: sectionproperties.pre.geometry.Geometry.from_3dm
   :noindex:

.. admonition:: Example

    TODO

..  automethod:: sectionproperties.pre.geometry.CompoundGeometry.from_3dm
   :noindex:

.. admonition:: Example

    TODO

:class:`~sectionproperties.pre.geometry.Geometry` objects can also be created from
encodings of Rhino BREP.

..  automethod:: sectionproperties.pre.geometry.Geometry.from_rhino_encoding
   :noindex:

.. admonition:: Example

    TODO

More advanced filtering can be achieved by working with the Shapely geometries directly.
These can be accessed by :func:`~sectionproperties.pre.rhino.load_3dm` and
:func:`~sectionproperties.pre.rhino.load_brep_encoding`.

Section Library
^^^^^^^^^^^^^^^

In order to make your life easier, there are a number of built-in functions that
generate typical structural cross-sections, resulting in
:class:`~sectionproperties.pre.geometry.Geometry` objects. These typical cross-sections
reside in the ``sectionproperties.pre.library`` module.

Primitive Sections
""""""""""""""""""

.. autosummary::
    :nosignatures:

    ~sectionproperties.pre.library.primitive_sections.rectangular_section
    ~sectionproperties.pre.library.primitive_sections.circular_section
    ~sectionproperties.pre.library.primitive_sections.circular_section_by_area
    ~sectionproperties.pre.library.primitive_sections.elliptical_section
    ~sectionproperties.pre.library.primitive_sections.triangular_section
    ~sectionproperties.pre.library.primitive_sections.triangular_radius_section
    ~sectionproperties.pre.library.primitive_sections.cruciform_section

Steel Sections
""""""""""""""

.. autosummary::
    :nosignatures:

    ~sectionproperties.pre.library.steel_sections.circular_hollow_section
    ~sectionproperties.pre.library.steel_sections.elliptical_hollow_section
    ~sectionproperties.pre.library.steel_sections.rectangular_hollow_section
    ~sectionproperties.pre.library.steel_sections.polygon_hollow_section
    ~sectionproperties.pre.library.steel_sections.i_section
    ~sectionproperties.pre.library.steel_sections.mono_i_section
    ~sectionproperties.pre.library.steel_sections.tapered_flange_i_section
    ~sectionproperties.pre.library.steel_sections.channel_section
    ~sectionproperties.pre.library.steel_sections.tapered_flange_channel
    ~sectionproperties.pre.library.steel_sections.tee_section
    ~sectionproperties.pre.library.steel_sections.angle_section
    ~sectionproperties.pre.library.steel_sections.cee_section
    ~sectionproperties.pre.library.steel_sections.zed_section
    ~sectionproperties.pre.library.steel_sections.box_girder_section
    ~sectionproperties.pre.library.steel_sections.bulb_section

Concrete Sections
"""""""""""""""""

.. autosummary::
    :nosignatures:

    ~sectionproperties.pre.library.concrete_sections.concrete_rectangular_section
    ~sectionproperties.pre.library.concrete_sections.concrete_column_section
    ~sectionproperties.pre.library.concrete_sections.concrete_tee_section
    ~sectionproperties.pre.library.concrete_sections.concrete_circular_section

Bridge Sections
"""""""""""""""

.. autosummary::
    :nosignatures:

    ~sectionproperties.pre.library.bridge_sections.super_t_girder_section
    ~sectionproperties.pre.library.bridge_sections.i_girder_section


Manipulating Geometry Objects
-----------------------------

kf;ldfs

Set Operations
^^^^^^^^^^^^^^

asdlkjsald

Align
^^^^^

asdjasd

Mirror
^^^^^^

asdlkjasld

Rotate
^^^^^^

adsjalsd

Shift
^^^^^

asdasd

Split
^^^^^

akjldskf

Offset
^^^^^^

sdflkjdsf


Assigning Material Properties
-----------------------------

jaskdlsad

Visualising Geometry
--------------------

asldkjasd
