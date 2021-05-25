.. _label-analysis:

Running an Analysis
===================

The first step in running a cross-section analysis is the creation of a
:class:`~sectionproperties.analysis.cross_section.CrossSection` object. This class
stores the structural geometry and finite element mesh and provides methods to
perform various types of cross-section analyses.

.. autoclass:: sectionproperties.analysis.cross_section.CrossSection
   :show-inheritance:
   :noindex:

Checking the Mesh Quality
-------------------------

Before carrying out a cross-section analysis it is a good idea to check the quality
of the finite element mesh. Some useful methods are provided to display mesh statistics
and to plot the finite element mesh:

.. automethod:: sectionproperties.analysis.cross_section.CrossSection.display_mesh_info
   :noindex:

.. automethod:: sectionproperties.analysis.cross_section.CrossSection.plot_mesh
   :noindex:

Geometric Analysis
------------------

A geometric analysis calculates the area properties of the cross-section.

.. automethod:: sectionproperties.analysis.cross_section.CrossSection.calculate_geometric_properties
   :noindex:

Plastic Analysis
----------------

A plastic analysis calculates the plastic properties of the cross-section.

.. note:: A geometric analysis must be performed on the CrossSection object before
  a plastic analysis is carried out.

.. automethod:: sectionproperties.analysis.cross_section.CrossSection.calculate_plastic_properties
   :noindex:

Warping Analysis
----------------

A warping analysis calculates the torsion and shear properties of the cross-section.

.. note:: A geometric analysis must be performed on the CrossSection object before
  a warping analysis is carried out.

.. automethod:: sectionproperties.analysis.cross_section.CrossSection.calculate_warping_properties
   :noindex:

Stress Analysis
---------------

A stress analysis calculates the cross-section stresses arising from a set of forces and moments.
Executing this method returns a :class:`~sectionproperties.analysis.cross_section.StressResult`
object which stores the cross-section stresses and provides stress plotting functions.

.. note:: A geometric and warping analysis must be performed on the CrossSection
   object before a stress analysis is carried out.

.. automethod:: sectionproperties.analysis.cross_section.CrossSection.calculate_stress
   :noindex:

Calculating Frame Properties
----------------------------

Calculates the cross-section properties required for a 2D or 3D frame analysis.

.. note:: This method is significantly faster than performing a geometric and
   a warping analysis and has no prerequisites.

.. automethod:: sectionproperties.analysis.cross_section.CrossSection.calculate_frame_properties
   :noindex:
