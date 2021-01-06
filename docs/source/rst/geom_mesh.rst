.. _label-geom_mesh:

Creating a Geometry, Mesh and Material Properties
=================================================

Before performing a cross-section analysis, the geometry of the cross-section and a finite element
mesh must be created. Optionally, material properties can be applied to different regions of the
cross-section.

Cross-Section Geometry
----------------------

The geometry of a cross-section defines its dimensions and shape and involves the creation of a
:class:`~sectionproperties.pre.sections.Geometry` object. This geometry object stores all the
information needed to create a finite element mesh.

..  autoclass:: sectionproperties.pre.sections.Geometry
    :noindex:

Different regions of the geometry can be specified by defining a list of ``control_points``, which
are located within unique enclosed areas of the geometry. Different regions can be used to specify
different mesh sizes and/or different material properties within the structural cross-section. See
the :ref:`label-examples` for some example scripts in which different regions are specified through
a list of control points.


Creating Common Structural Geometries
---------------------------------------

In order to make your life easier, there are a number of built-in classes that generate typical
structural cross-sections that inherit from the :class:`~sectionproperties.pre.sections.Geometry`
class. Note that these classes automatically assign a ``control_point`` to the geometry object.

Rectangular Section
^^^^^^^^^^^^^^^^^^^
..  autoclass:: sectionproperties.pre.sections.RectangularSection
    :show-inheritance:
    :noindex:

Circular Section
^^^^^^^^^^^^^^^^
..  autoclass:: sectionproperties.pre.sections.CircularSection
    :show-inheritance:
    :noindex:

Circular Hollow Section (CHS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
..  autoclass:: sectionproperties.pre.sections.Chs
    :show-inheritance:
    :noindex:

Elliptical Section
^^^^^^^^^^^^^^^^^^
..  autoclass:: sectionproperties.pre.sections.EllipticalSection
    :show-inheritance:
    :noindex:

Elliptical Hollow Section (EHS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
..  autoclass:: sectionproperties.pre.sections.Ehs
    :show-inheritance:
    :noindex:

Rectangular Hollow Section (RHS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
..  autoclass:: sectionproperties.pre.sections.Rhs
    :show-inheritance:
    :noindex:

I-Section
^^^^^^^^^
  ..  autoclass:: sectionproperties.pre.sections.ISection
      :show-inheritance:
      :noindex:

Monosymmetric I-Section
^^^^^^^^^^^^^^^^^^^^^^^
  ..  autoclass:: sectionproperties.pre.sections.MonoISection
      :show-inheritance:
      :noindex:

Tapered Flange I-Section
^^^^^^^^^^^^^^^^^^^^^^^^
  ..  autoclass:: sectionproperties.pre.sections.TaperedFlangeISection
      :show-inheritance:
      :noindex:

Parallel Flange Channel (PFC) Section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  ..  autoclass:: sectionproperties.pre.sections.PfcSection
      :show-inheritance:
      :noindex:

Tapered Flange Channel Section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  ..  autoclass:: sectionproperties.pre.sections.TaperedFlangeChannel
      :show-inheritance:
      :noindex:

Tee Section
^^^^^^^^^^^
  ..  autoclass:: sectionproperties.pre.sections.TeeSection
      :show-inheritance:
      :noindex:

Angle Section
^^^^^^^^^^^^^
  ..  autoclass:: sectionproperties.pre.sections.AngleSection
      :show-inheritance:
      :noindex:

Cee Section
^^^^^^^^^^^
  ..  autoclass:: sectionproperties.pre.sections.CeeSection
      :show-inheritance:
      :noindex:

Zed Section
^^^^^^^^^^^
  ..  autoclass:: sectionproperties.pre.sections.ZedSection
      :show-inheritance:
      :noindex:

Cruciform Section
^^^^^^^^^^^^^^^^^
  ..  autoclass:: sectionproperties.pre.sections.CruciformSection
      :show-inheritance:
      :noindex:

Polygon Section
^^^^^^^^^^^^^^^
  ..  autoclass:: sectionproperties.pre.sections.PolygonSection
      :show-inheritance:
      :noindex:

Box Girder Section
^^^^^^^^^^^^^^^^^^
  ..  autoclass:: sectionproperties.pre.sections.BoxGirderSection
      :show-inheritance:
      :noindex:

Arbitrary Cross-Section Geometries
----------------------------------

If none of the above classes gives you what you need, you can create a
:class:`~sectionproperties.pre.sections.CustomSection` geometry object, which is defined by a list
of points (nodes), facets (node connectivities) and hole locations:

..  autoclass:: sectionproperties.pre.sections.CustomSection
    :show-inheritance:
    :noindex:

..  note:: Ensure that the ``control_points`` you choose lie within the
  :class:`~sectionproperties.pre.sections.CustomSection`. If any of the ``control_points`` are
  outside the region, on an edge or within a hole, the meshing algorithm will likely not treat
  distinct areas within the :class:`~sectionproperties.pre.sections.CustomSection` as a separate
  regions and mesh refinements may not work as anticipated.

..  note:: In order to calculate the perimeter of the cross-section be sure to enter the facet
  indices that correspond to the perimeter of your cross-section.

Merging Geometries
------------------

If you wish to merge multiple :class:`~sectionproperties.pre.sections.Geometry` objects into a
single object, you can use the :class:`~sectionproperties.pre.sections.MergedSection` class:

..  autoclass:: sectionproperties.pre.sections.MergedSection
    :show-inheritance:
    :noindex:

..  note:: There must be connectivity between the :class:`~sectionproperties.pre.sections.Geometry`
  objects that you wish to merge. It is currently not possible to analyse a cross-section that is
  composed of two or more unconnected domains.

..  note:: You may need to overwrite the perimeter facets list if predefined sections are used.
  Enabling labels while plotting the geometry is an easy way to manually identify the facet indices
  that make up the perimeter of the cross-section.

Cleaning the Geometry
---------------------

When creating a merged section often there are overlapping facets or duplicate nodes. These
geometry artefacts can cause difficulty for the meshing algorithm. It is therefore recommended to
clean the geometry after merging sections which may result in overlapping or intersecting facets,
or duplicate nodes. Cleaning the geometry can be carried out by using the
:func:`~sectionproperties.pre.sections.Geometry.clean_geometry` method:

..  automethod:: sectionproperties.pre.sections.Geometry.clean_geometry
    :noindex:

Perimeter Offset
----------------

The perimeter of a cross-section geometry can be offset by using the
:func:`~sectionproperties.pre.offset.offset_perimeter` method:

..  autofunction:: sectionproperties.pre.offset.offset_perimeter
    :noindex:

..  note:: All the built-in sections in the ``sections`` module are built using an anti-clockwise
  facet direction. As a result, side='left' will reduce the cross-section, while side='right' will
  increase the cross-section.

..  note:: The ``control_points`` may need to be manually re-assigned if reducing the cross-section
  moves the control_point outside the geometry.

..  warning:: This feature is a *beta* addition and as a result may produce some errors if the
  offsetting drastically changes the geometry.

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
