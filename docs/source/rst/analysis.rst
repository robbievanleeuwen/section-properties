.. _label-analysis:

Running an Analysis
===================

The first step in running a section analysis is the creation of a
:class:`~sectionproperties.analysis.section.Section` object. This class
stores the structural geometry and finite element mesh and provides methods to
perform various types of sectional analyses.

..  autoclass:: sectionproperties.analysis.section.Section
    :noindex:

Checking the Mesh Quality
-------------------------

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

..  note:: A geometric analysis must be performed on the Section object before
  a plastic analysis is carried out.

..  warning:: The plastic analysis in *sectionproperties* assumes all materials are
  able to reach their yield stress defined in the material properties. Care should be
  taken if analysing materials or cross-sections exhibiting non-linear behaviour, e.g.
  reinforced concrete or non-compact steel sections.

..  automethod:: sectionproperties.analysis.section.Section.calculate_plastic_properties
    :noindex:

Warping Analysis
----------------

A warping analysis calculates the torsion and shear properties of the section.

..  note:: A geometric analysis must be performed on the Section object before
  a warping analysis is carried out.

..  warning:: There must be connectivity between all elements of the mesh to perform a
  valid warping analysis. This is a limitation of the elastic theory that this
  implementation is based on, as there is no way to quantify the transfer of shear and
  warping between two unconnected regions.

..  automethod:: sectionproperties.analysis.section.Section.calculate_warping_properties
    :noindex:

Stress Analysis
---------------

A stress analysis calculates the section stresses arising from a set of forces
and moments. Executing this method returns a :class:`~sectionproperties.analysis.section.StressResult`
object which stores the section stresses and provides stress plotting functions.

..  note:: A geometric analysis must be performed on the Section object before a stress
  analysis is carried out. Further, if the shear force or twisting moment is non-zero
  a warping analysis must also be performed.

..  warning:: The stress analysis in *sectionproperties* is linear-elastic and does not
  account for the yielding of materials or any non-linearities.

..  automethod:: sectionproperties.analysis.section.Section.calculate_stress
    :noindex:

Calculating Frame Properties
----------------------------

Calculates the section properties required for a 2D or 3D frame analysis.

..  note:: This method is significantly faster than performing a geometric and
  a warping analysis and has no prerequisites.

..  automethod:: sectionproperties.analysis.section.Section.calculate_frame_properties
    :noindex:
