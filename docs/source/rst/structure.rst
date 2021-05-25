Structure of an Analysis
========================

The process of performing a cross-section analysis with *sectionproperties* can
be broken down into three stages:

1. Pre-Processor: The input geometry and finite element mesh is created.
2. Solver: The cross-section properties are determined.
3. Post-Processor: The results are presented in a number of different formats.

Creating a Geometry and Mesh
----------------------------

The dimensions and shape of the cross-section to be analysed define the *geometry*
of the cross-section. The :ref:`label-sections-module` provides a number of classes
to easily generate either commonly used structural sections or an  arbitrary
cross-section, defined by a list of points, facets and holes. All of the classes
in the :ref:`label-sections-module` inherit from the
:class:`~sectionproperties.pre.sections.Geometry` class.

The final stage in the pre-processor involves generating a finite element mesh of
the *geometry* that the solver can use to calculate the cross-section properties.
This can easily be performed using the :func:`~sectionproperties.pre.sections.Geometry.create_mesh`
method that all :class:`~sectionproperties.pre.sections.Geometry` objects have
access to.

The following example creates a geometry object with a PFC cross-section. A finite element mesh is
generated with a maximum triangular area of 2.5.

.. jupyter-execute::

   import sectionproperties.pre.sections as sections
   from sectionproperties.analysis.cross_section import CrossSection
   geometry = sections.PfcSection(d=200, b=75, t_f=12, t_w=6, r=12, n_r=8)
   mesh = geometry.create_mesh(mesh_sizes=[2.5])
   section = CrossSection(geometry, mesh)
   section.plot_mesh()

If you are analysing a composite section, or would like to include material properties in your
model, material properties can be created using the :class:`~sectionproperties.pre.pre.Material`
class. The following example creates a steel material object

.. jupyter-execute::

   from sectionproperties.pre.pre import Material
   steel = Material(
       name='Steel',
       elastic_modulus=200e3,
       poissons_ratio=0.3,
       yield_strength=500,
       color='grey'
   )

Refer to :ref:`label-geom_mesh` for a more detailed explanation of the pre-processing
stage.

Running an Analysis
-------------------

The solver operates on a :class:`~sectionproperties.analysis.cross_section.CrossSection`
object and can perform four different analysis types:

- Geometric Analysis: calculates area properties.
- Plastic Analysis: calculates plastic properties.
- Warping Analysis: calculates torsion and shear properties.
- Stress Analysis: calculates cross-section stresses.

The geometric analysis can be performed individually. However in order to perform
a warping or plastic analysis, a geometric analysis must first be performed. Further,
in order to carry out a stress analysis, both a geometric and warping analysis must
have already been executed. The program will display a helpful error if you try
to run any of these analyses without first performing the prerequisite analyses.

The following example performs a geometric, plastic, and warping analysis on the
cross-section defined in the previous section with steel used as the material
property

.. jupyter-execute::

   from sectionproperties.analysis.cross_section import CrossSection
   section = CrossSection(geometry, mesh, [steel])
   section.calculate_geometric_properties()
   section.calculate_plastic_properties()
   section.calculate_warping_properties()

Refer to :ref:`label-analysis` for a more detailed explanation of the solver stage.

Viewing the Results
-------------------

Once an analysis has been performed, a number of methods belonging to the
:class:`~sectionproperties.analysis.cross_section.CrossSection` object can be called
to present the cross-section results in a number of different formats. For example
the cross-section properties can be printed to the terminal, a plot of the centroids
displayed and the cross-section stresses visualised in a contour plot.

Continuing with the example above, the cross-section properties are printed to the terminal and a
plot of the centroids is displayed.

.. jupyter-execute::

   section.display_results()

.. jupyter-execute::

   section.plot_centroids()

Refer to :ref:`label-post` for a more detailed explanation of the post-processing
stage.
