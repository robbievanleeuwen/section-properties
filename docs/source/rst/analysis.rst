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

..  automethod:: sectionproperties.analysis.section.Section.calculate_plastic_properties
    :noindex:

Warping Analysis
----------------

A warping analysis calculates the torsion and shear properties of the section.

..  note:: A geometric analysis must be performed on the Section object before
  a warping analysis is carried out.

..  automethod:: sectionproperties.analysis.section.Section.calculate_warping_properties
    :noindex:

Stress Analysis
---------------

A stress analysis calculates the section stresses arising from a set of forces
and moments. Executing this method returns a :class:`~sectionproperties.analysis.section.StressResult`
object which stores the section stresses and provides stress plotting functions.

..  note:: A geometric and warping analysis must be performed on the Section
  object before a stress analysis is carried out.

..  automethod:: sectionproperties.analysis.section.Section.calculate_stress
    :noindex:

Calculating Frame Properties
----------------------------

Calculates the section properties required for a 2D or 3D frame analysis.

..  note:: This method is significantly faster than performing a geometric and
  a warping analysis and has no prerequisites.

..  automethod:: sectionproperties.analysis.section.Section.calculate_frame_properties
    :noindex:
