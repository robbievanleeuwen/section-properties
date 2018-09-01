.. _label-geom_mesh:

Creating a Geometry and Mesh
============================

Before performing a cross-section analysis, the geometry of the cross-section and
a finite element mesh must be created.

Cross-Section Geometry
----------------------

The geometry of a cross-section defines its dimensions, shape and material properties
and involves the creation of a :class:`~sectionproperties.pre.sections.Geometry`
object. This geometry object stores all the information needed to create a finite
element mesh.

..  autoclass:: sectionproperties.pre.sections.Geometry
    :noindex:

Different regions of the geometry can be specified by defining a list of ``control_points``,
which are located within unique enclosed areas of the geometry. Different regions
can be used to specify different mesh sizes and/or different material properties
within the structural cross-section. See the :ref:`label-examples` for some example
scripts in which different regions are specified through a list of control points.


Creating Common Structural Geometries
---------------------------------------

In order to make your life easier, there are a number of built-in classes that generate
typical structural cross-sections that also inherit from the :class:`~sectionproperties.pre.sections.Geometry`
class. Note that these classes automatically assign a ``control_point`` to the geometry
object.

# TODO: add pictures for each geometry!

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

Parallel Flange Channel (PFC) Section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  ..  autoclass:: sectionproperties.pre.sections.PfcSection
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

Arbitrary Cross-Section Geometries
----------------------------------

If none of the above classes gives you what you need, you can create a
:class:`~sectionproperties.pre.sections.CustomSection` geometry object, which is
defined by a list of points (nodes), facets (node connectivities) and hole locations:

..  autoclass:: sectionproperties.pre.sections.CustomSection
    :show-inheritance:
    :noindex:

Merging Geometries
------------------

If you wish to merge multiple :class:`~sectionproperties.pre.sections.Geometry`
objects into a single object, you can use the :class:`~sectionproperties.pre.sections.MergedSection`
class:

..  note:: There must be connectivity between the :class:`~sectionproperties.pre.sections.Geometry`
  objects that you wish to merge. It is currently not possible to analyse a cross-section
  that is composed of two or more unconnected domains.

..  autoclass:: sectionproperties.pre.sections.MergedSection
    :show-inheritance:
    :noindex:

Visualising the Geometry
------------------------

Geometry objects can be visualised by using the :func:`~sectionproperties.pre.sections.Geometry.plot_geometry`
method:

..  automethod:: sectionproperties.pre.sections.Geometry.plot_geometry
    :noindex:

Generating a Mesh
-----------------

A finite element mesh is required to perform a cross-section analysis. A finite
element mesh can be created by using the :func:`~sectionproperties.pre.sections.Geometry.create_mesh`
method:

..  automethod:: sectionproperties.pre.sections.Geometry.create_mesh
    :noindex:

..  warning:: The length of ``mesh_sizes`` must match the number of regions
  in the geometry object.
