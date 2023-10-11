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
        :caption: Geometry object from a shapely polygon

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
        :caption: CompoundGeometry object from two shapely polygons

        from shapely import Polygon
        from sectionproperties.pre import Geometry, CompoundGeometry

        sq1 = Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])
        sq2 = Polygon([(2, 0), (6, 0), (6, 4), (2, 4)])
        geom_sq1 = Geometry(geom=sq1)
        geom_sq2 = Geometry(geom=sq2)
        geom = CompoundGeometry(geoms=[geom_sq1, geom_sq2])
        geom.plot_geometry()

.. _label-from-points:

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
    installation. To install ``sectionproperties`` with CAD import functionality,  use
    the ``dxf`` and/or ``rhino`` options:

    .. code-block:: shell

        pip install sectionproperties[dxf]
        pip install sectionproperties[rhino]

.. _label-geometry-dxf:

``.dxf``
""""""""

:class:`~sectionproperties.pre.geometry.Geometry` objects can be created from ``.dxf``
files using the
:meth:`sectionproperties.pre.geometry.Geometry.from_dxf` method.

..  automethod:: sectionproperties.pre.geometry.Geometry.from_dxf
   :noindex:

.. admonition:: Example

    The following example loads a ``.dxf`` file and creates a
    :class:`~sectionproperties.pre.geometry.Geometry` object from its contents.

    .. plot::
        :include-source: True
        :caption: Geometry object from a ``.dxf`` file

        from sectionproperties.pre import Geometry

        # the following path is a .dxf file that describes a box section with two holes
        dxf_path = "../_static/cad_files/box_section.dxf"

        # load dxf file into a Geometry object
        geom = Geometry.from_dxf(dxf_filepath=dxf_path)
        geom.plot_geometry()

.. note::

    Loading multiple regions from a single ``.dxf`` file into a ``CompoundGeometry`` is
    not currently supported in ``sectionproperties``. A possible work around involves
    saving each region as a separate ``.dxf`` file, importing each region individually
    using ``Geometry.from_dxf()``, then combining the regions using the ``+`` operator.

.. _label-geometry-3dm:

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

    The following example loads a ``.3dm`` file and creates a
    :class:`~sectionproperties.pre.geometry.Geometry` object from its contents.

    .. plot::
        :include-source: True
        :caption: Geometry object from a ``.3dm`` file

        from sectionproperties.pre import Geometry

        # the following path is a .3dm file that describes a glazing section
        rhino_path = "../_static/cad_files/rhino.3dm"

        # load 3dm file into a Geometry object
        geom = Geometry.from_3dm(filepath=rhino_path)
        geom.plot_geometry()

..  automethod:: sectionproperties.pre.geometry.CompoundGeometry.from_3dm
   :noindex:

.. admonition:: Example

    The following example loads a ``.3dm`` file and creates a
    :class:`~sectionproperties.pre.geometry.CompoundGeometry` object from its contents.

    .. plot::
        :include-source: True
        :caption: CompoundGeometry object from a ``.3dm`` file

        from sectionproperties.pre import CompoundGeometry

        # the following path is a .3dm file that describes two distinct 2D surfaces
        rhino_path = "../_static/cad_files/rhino_compound.3dm"

        # load 3dm file into a CompoundGeometry object
        geom = CompoundGeometry.from_3dm(filepath=rhino_path)
        geom.plot_geometry()

:class:`~sectionproperties.pre.geometry.Geometry` objects can also be created from
encodings of Rhino BREP.

..  automethod:: sectionproperties.pre.geometry.Geometry.from_rhino_encoding
   :noindex:

.. admonition:: Example

    The following example loads a ``.json`` file describing a Rhino BREP and creates a
    :class:`~sectionproperties.pre.geometry.Geometry` object from its contents.

    .. plot::
        :include-source: True
        :caption: Geometry object from a Rhino BREP file

        import json
        from sectionproperties.pre import Geometry

        # the following path is a .json file that is a BREP describing a 1 x 1 square
        rhino_path = "../_static/cad_files/rhino_brep.json"

        with open(rhino_path) as rhino_file:
            brep_encoded = json.load(rhino_file)

        # load BREP file into a Geometry object
        geom = Geometry.from_rhino_encoding(r3dm_brep=brep_encoded)
        geom.plot_geometry()

More advanced filtering can be achieved by working with the Shapely geometries directly.
These can be accessed by :func:`~sectionproperties.pre.rhino.load_3dm` and
:func:`~sectionproperties.pre.rhino.load_brep_encoding`.

Section Library
^^^^^^^^^^^^^^^

In order to make your life easier, there are a number of built-in functions that
generate typical structural cross-sections, resulting in
:class:`~sectionproperties.pre.geometry.Geometry` or
:class:`~sectionproperties.pre.geometry.CompoundGeometry` objects. These typical
cross-sections reside in the ``sectionproperties.pre.library`` module.

.. _label-primitive-library:

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

.. _label-concrete-library:

Concrete Sections
"""""""""""""""""

.. autosummary::
    :nosignatures:

    ~sectionproperties.pre.library.concrete_sections.concrete_rectangular_section
    ~sectionproperties.pre.library.concrete_sections.concrete_column_section
    ~sectionproperties.pre.library.concrete_sections.concrete_tee_section
    ~sectionproperties.pre.library.concrete_sections.concrete_circular_section

.. _label-bridge-library:

Bridge Sections
"""""""""""""""

.. autosummary::
    :nosignatures:

    ~sectionproperties.pre.library.bridge_sections.super_t_girder_section
    ~sectionproperties.pre.library.bridge_sections.i_girder_section

NASTRAN Sections
""""""""""""""""

.. autosummary::
    :nosignatures:

    ~sectionproperties.pre.library.nastran_sections.nastran_bar
    ~sectionproperties.pre.library.nastran_sections.nastran_box
    ~sectionproperties.pre.library.nastran_sections.nastran_box1
    ~sectionproperties.pre.library.nastran_sections.nastran_chan
    ~sectionproperties.pre.library.nastran_sections.nastran_chan1
    ~sectionproperties.pre.library.nastran_sections.nastran_chan2
    ~sectionproperties.pre.library.nastran_sections.nastran_cross
    ~sectionproperties.pre.library.nastran_sections.nastran_fcross
    ~sectionproperties.pre.library.nastran_sections.nastran_dbox
    ~sectionproperties.pre.library.nastran_sections.nastran_gbox
    ~sectionproperties.pre.library.nastran_sections.nastran_h
    ~sectionproperties.pre.library.nastran_sections.nastran_hat
    ~sectionproperties.pre.library.nastran_sections.nastran_hat1
    ~sectionproperties.pre.library.nastran_sections.nastran_hexa
    ~sectionproperties.pre.library.nastran_sections.nastran_i
    ~sectionproperties.pre.library.nastran_sections.nastran_i1
    ~sectionproperties.pre.library.nastran_sections.nastran_l
    ~sectionproperties.pre.library.nastran_sections.nastran_rod
    ~sectionproperties.pre.library.nastran_sections.nastran_tee
    ~sectionproperties.pre.library.nastran_sections.nastran_tee1
    ~sectionproperties.pre.library.nastran_sections.nastran_tee2
    ~sectionproperties.pre.library.nastran_sections.nastran_tube
    ~sectionproperties.pre.library.nastran_sections.nastran_tube2
    ~sectionproperties.pre.library.nastran_sections.nastran_zed

Manipulating Geometry Objects
-----------------------------

Geometries in ``sectionproperties`` are able to be manipulated in 2D space for the
purpose of creating novel, custom section geometries that the user may require.

.. note::

   Operations on geometries are **non-destructive**. For each operation, a new geometry
   object is returned.

   This gives ``sectionproperties`` geometries a *fluent API*, meaning that
   transformation methods can be chained together, see
   `Advanced Geometry Creation <../examples/geometry/advanced_geometry.ipynb>`_ for
   further examples.

Align
^^^^^

There are two available align methods:

#. ``align_to()`` - aligns one geometry to another on a specified side.
#. ``align_center()`` - aligns the center of one geometry to either the center of
   another, or a specific point.

..  automethod:: sectionproperties.pre.geometry.Geometry.align_to
    :noindex:

..  automethod:: sectionproperties.pre.geometry.Geometry.align_center
    :noindex:

..  automethod:: sectionproperties.pre.geometry.CompoundGeometry.align_center
    :noindex:

Mirror
^^^^^^

Geometry can be mirrored about a specified point on either the ``x`` or ``y`` axis.

..  automethod:: sectionproperties.pre.geometry.Geometry.mirror_section
    :noindex:

Rotate
^^^^^^

Geometry can be rotated by any angle about a point.

..  automethod:: sectionproperties.pre.geometry.Geometry.rotate_section
    :noindex:

Shift
^^^^^

There are two available shift methods:

#. ``shift_section()`` - shifts the entire geometry by a vector.
#. ``shift_points()`` - shifts specific points within the geometry by either a vector,
   or to an absolute location.

..  automethod:: sectionproperties.pre.geometry.Geometry.shift_section
    :noindex:

..  automethod:: sectionproperties.pre.geometry.Geometry.shift_points
    :noindex:

Split
^^^^^

Geometry can be split either side of a straight line.

..  automethod:: sectionproperties.pre.geometry.Geometry.split_section
    :noindex:

Offset
^^^^^^

The external and/or internal perimeter of a geometry can be dilated or eroded by a set
value.

..  automethod:: sectionproperties.pre.geometry.Geometry.offset_perimeter
    :noindex:

..  automethod:: sectionproperties.pre.geometry.CompoundGeometry.offset_perimeter
    :noindex:

.. _label-geometry-set:

Set Operations
^^^^^^^^^^^^^^

Both :class:`~sectionproperties.pre.geometry.Geometry` and
:class:`~sectionproperties.pre.geometry.CompoundGeometry` objects can be manipulated
using Python's set operators. See
`Advanced Geometry Creation <../examples/geometry/advanced_geometry.ipynb>`_ for further
examples using set operations.


``|`` (Union)
"""""""""""""

.. automethod:: sectionproperties.pre.geometry.Geometry.__or__
    :noindex:

``-`` (Subtraction)
"""""""""""""""""""

.. automethod:: sectionproperties.pre.geometry.Geometry.__sub__
    :noindex:

``&`` (Intersection)
""""""""""""""""""""

.. automethod:: sectionproperties.pre.geometry.Geometry.__and__
    :noindex:

``^`` (Symmetric Difference)
""""""""""""""""""""""""""""

.. automethod:: sectionproperties.pre.geometry.Geometry.__xor__
    :noindex:

``+`` (Addition)
""""""""""""""""

.. automethod:: sectionproperties.pre.geometry.Geometry.__add__
    :noindex:


.. _label-geom-material:

Assigning Material Properties
-----------------------------

Each :class:`~sectionproperties.pre.geometry.Geometry` contains its own material
definition, which is stored in the ``.material`` attribute. The simplest way to assign
a material to a :class:`~sectionproperties.pre.geometry.Geometry` is to pass the
material as an argument to the constructor.

.. note::

    If a :class:`~sectionproperties.pre.pre.Material` is not given, then the *default
    material* is assigned to the ``Geometry.material`` attribute. The default material
    has an elastic modulus of 1, a Poisson's ratio of 0, a density of 1 and a yield
    strength of 1.

    This is equivalent to performing a purely geometric analysis of the cross-section
    and is desirable if a composite section is not being analysed.

.. warning::

    See more about how asssigning material properties affects the results reported by
    ``sectionproperties`` :ref:`here<label-material-affects-results>`.

Below are a few examples showcasing the different ways to generate geometry discussed
above:

.. admonition:: Example

    The following example assigns material properties to a number of different
    geometries:

    .. code-block:: python

        from shapely import Polygon

        from sectionproperties.pre import Material
        from sectionproperties.pre import Geometry
        from sectionproperties.pre.library import rectangular_section

        # create a steel material
        steel = Material(
            name="Steel",
            elastic_modulus=200e3,
            poissons_ratio=0.3,
            density=7.85e-6,
            yield_strength=500,
            color="grey",
        )

        # assign steel to a shapely generated geometry
        poly = Polygon([(0, 0), (5, 2), (3, 7), (1, 6)])
        geom = Geometry(geom=poly, material=steel)

        # assign steel to a geometry from points
        points = [(0, 0), (10, 5), (15, 15), (5, 10), (6, 6), (9, 7), (7, 9)]
        facets = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 4)]
        control_points = [(4, 4)]
        holes = [(7, 7)]
        geom = Geometry.from_points(
            points=points,
            facets=facets,
            control_points=control_points,
            holes=holes,
            material=steel,
        )

        # assign steel to a rectangular section
        geom = rectangular_section(d=100, b=50, material=steel)

A geometry's material may be altered at any time by simply assigning a new
:class:`~sectionproperties.pre.pre.Material` to the ``.material`` attribute. This is
also useful when creating geometry from CAD files:

.. admonition:: Example

    The following example demonstrates assigning material properties through changing
    the ``.material`` attribute.

    .. code-block:: python

        from sectionproperties.pre import Material
        from sectionproperties.pre import Geometry

        # create a steel material
        steel = Material(
            name="Steel",
            elastic_modulus=200e3,
            poissons_ratio=0.3,
            density=7.85e-6,
            yield_strength=500,
            color="grey",
        )

        # load 3dm file into a Geometry object
        geom = Geometry.from_3dm(filepath="example.3dm")

        # assign steel to the geometry
        geom.material = steel

A :class:`~sectionproperties.pre.geometry.CompoundGeometry` does not have a
``.material`` attribute and therefore, a :class:`~sectionproperties.pre.pre.Material`
cannot be directly assigned. Since a
:class:`~sectionproperties.pre.geometry.CompoundGeometry` is simply a combination of
:class:`~sectionproperties.pre.geometry.Geometry` objects, the material should be
assigned to each individual :class:`~sectionproperties.pre.geometry.Geometry` object
that make up the :class:`~sectionproperties.pre.geometry.CompoundGeometry`.

.. admonition:: Example

    The following example demonstrates assigning material properties to
    :class:`~sectionproperties.pre.geometry.CompoundGeometry` objects.

    .. plot::
        :include-source: True
        :caption: Assign materials to a ``CompoundGeometry`` object

        from shapely import Polygon

        from sectionproperties.pre import Material
        from sectionproperties.pre.library import rectangular_section
        from sectionproperties.analysis import Section

        # create steel and timber materials
        steel = Material(
            name="Steel",
            elastic_modulus=200e3,
            poissons_ratio=0.3,
            density=7.85e-6,
            yield_strength=500,
            color="grey",
        )
        timber = Material(
            name="Timber",
            elastic_modulus=8e3,
            poissons_ratio=0.35,
            density=6.5e-7,
            yield_strength=20,
            color="burlywood",
        )

        # create the individual geometries with material properties applied
        beam = rectangular_section(d=170, b=35, material=timber)
        plate = rectangular_section(d=16, b=35, material=steel)

        # combine geometries, maintaining assigned materials
        geom = beam + plate.shift_section(y_offset=-16)

        # mesh and plot
        geom.create_mesh(mesh_sizes=[20, 10])
        Section(geometry=geom).plot_mesh()

Visualising Geometry
--------------------

Visualisation of geometry objects is best performed in the Jupyter computing
environment, however, most visualisation can also be done in any environment which
supports the display of matplotlib plots.

There are generally two ways to visualise geometry objects:

#. In the Jupyter computing environment, geometry objects utilise their underlying
   :class:`shapely.geometry.Polygon` object's ``_repr_svg_`` method to show the geometry
   as it's own representation.
#. By using the :meth:`~sectionproperties.pre.geometry.Geometry.plot_geometry` method.

..  automethod:: sectionproperties.pre.geometry.Geometry.plot_geometry
    :noindex:

.. note::

   You can also use
   :meth:`~sectionproperties.pre.geometry.CompoundGeometry.plot_geometry` with
   :class:`~sectionproperties.pre.geometry.CompoundGeometry` objects
