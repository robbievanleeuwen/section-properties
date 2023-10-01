.. _label-analysis:

Analysis
========

Section Object
--------------

The first step in running a cross-section analysis in ``sectionproperties`` involves the
creation of a :class:`~sectionproperties.analysis.section.Section` object. This object
stores the cross-section geometry and finite element mesh, providing methods to perform
various types of cross-sectional analyses.

.. automethod:: sectionproperties.analysis.section.Section.__init__
    :noindex:

Checking Mesh Quality
---------------------

Before carrying out a section analysis it is a good idea to check the quality
of the finite element mesh. Some useful methods are provided to display mesh statistics
and to plot the finite element mesh:

..  automethod:: sectionproperties.analysis.section.Section.display_mesh_info
    :noindex:

..  automethod:: sectionproperties.analysis.section.Section.plot_mesh
    :noindex:

Geometric Analysis
------------------

A geometric analysis calculates the area properties of the section.

..  automethod:: sectionproperties.analysis.section.Section.calculate_geometric_properties
    :noindex:

Plastic Analysis
----------------

A plastic analysis calculates the plastic properties of the section.

..  warning:: The plastic analysis in *sectionproperties* assumes all materials are
  able to reach their yield stress defined in the material properties. Care should be
  taken if analysing materials or cross-sections exhibiting non-linear behaviour, e.g.
  reinforced concrete or non-compact steel sections.

..  automethod:: sectionproperties.analysis.section.Section.calculate_plastic_properties
    :noindex:

Warping Analysis
----------------

A warping analysis calculates the torsion and shear properties of the section.

..  warning:: There must be connectivity between all elements of the mesh to perform a
  valid warping analysis. This is a limitation of the elastic theory that this
  implementation is based on, as there is no way to quantify the transfer of shear and
  warping between two unconnected regions.

..  automethod:: sectionproperties.analysis.section.Section.calculate_warping_properties
    :noindex:

Frame Analysis
--------------

Calculates the section properties required for a 2D or 3D frame analysis.

..  note:: This method is significantly faster than performing a geometric and
  a warping analysis and has no prerequisites.

..  automethod:: sectionproperties.analysis.section.Section.calculate_frame_properties
    :noindex:

Stress Analysis
---------------

A stress analysis calculates the section stresses arising from a set of forces
and moments. Executing this method returns a
:class:`~sectionproperties.post.stress_post.StressPost` object, which stores the
section stresses and provides stress plotting methods.

..  warning:: The stress analysis in *sectionproperties* is linear-elastic and does not
  account for the yielding of materials or any non-linearities.

..  automethod:: sectionproperties.analysis.section.Section.calculate_stress
    :noindex:
