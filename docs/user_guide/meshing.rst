Meshing
=======

A finite element mesh is required to perform a cross-section analysis. After a geometry
has been created, a finite element mesh can then be created for the geometry by using
the :meth:`sectionproperties.pre.geometry.Geometry.create_mesh` or
:meth:`sectionproperties.pre.geometry.CompoundGeometry.create_mesh` methods:

..  automethod:: sectionproperties.pre.geometry.Geometry.create_mesh
    :noindex:

..  automethod:: sectionproperties.pre.geometry.CompoundGeometry.create_mesh
    :noindex:

..  warning::

  The length of ``mesh_sizes`` must match the number of regions in the geometry object.

Once the mesh has been created, it is stored within the geometry object and the geometry
object can then be passed to :class:`~sectionproperties.analysis.section.Section` for
analysis.

Mesh quality analysis, such as plotting the mesh and displaying mesh metrics, can be
performed using the :class:`~sectionproperties.analysis.section.Section` class. Please
see :ref:`label-analysis` for further information on performing analyses.
