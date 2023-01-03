Changelog:
==========

v2.1.5:
-------

- Fix shapely 2.0 imports and STRtree implementation, with thanks to `@normanrichardson <https://github.com/normanrichardson>`_
- Add support for python 3.10, drop support for python 3.7

**Full changelog:** `2.1.4...2.1.5 <https://github.com/robbievanleeuwen/section-properties/compare/2.1.4...2.1.5>`_

v2.1.4:
-------

- Add side bar option to ``concrete_rectangular_section()``, thanks to `@Agent6-6-6 <https://github.com/Agent6-6-6>`_
- Fix difference operation raising an error, thanks to `@connorferster <https://github.com/connorferster>`_
- Added ``concrete_column_section()`` and ``add_bar()`` methods, thanks to `@connorferster <https://github.com/connorferster>`_

**Full changelog:** `2.1.3...2.1.4 <https://github.com/robbievanleeuwen/section-properties/compare/2.1.3...2.1.4>`_

v2.1.3:
-------

- Retrieve cross-section stresses at any point using ``get_stress_at_point()`` or
  ``get_stress_at_points()``, many thanks to `@normanrichardson <https://github.com/normanrichardson>`_
- Fix plot legend formatting, thanks to `@Agent6-6-6 <https://github.com/Agent6-6-6>`_
- Added ability for ``Geometry.align_center()`` and ``CompoundGeometry.align_center()``
  to accept an ``x``, ``y`` coordinate as a valid input, thanks to `@connorferster <https://github.com/connorferster>`_
- Only require a warping analysis to be performed for a stress analysis if the shear
  force or twisting moment is non-zero.

**Full changelog:** `2.1.2...2.1.3 <https://github.com/robbievanleeuwen/section-properties/compare/2.1.2...2.1.3>`_

v2.1.2:
-------

- Make rhino-shapley-interop an optional dependency

**Full changelog:** `2.1.1...2.1.2 <https://github.com/robbievanleeuwen/section-properties/compare/2.1.1...2.1.2>`_

v2.1.1:
-------

- Use Lagrangian multiplier for calculation of torsion properties
- Add more plotting options to ``plot_geometry()``

**Full changelog:** `2.1.0...2.1.1 <https://github.com/robbievanleeuwen/section-properties/compare/2.1.0...2.1.1>`_

v2.1.0:
-------

- Add ``bulb_section()`` to steel sections library, thanks to `@zmpulse <https://github.com/zmpulse>`_
- Add progress bar and pretty output using `rich <https://github.com/Textualize/rich>`_
- Fix logic of generating holes in CompoundGeometry using the subtraction method, thanks to `@connorferster <https://github.com/connorferster>`_
- Expand testing suite and documentation, thanks to `@czarified <https://github.com/czarified>`_
- Fix bug with plastic calculation when material properties are specified
- Add warning message for disconnected geometries when trying to calculate warping properties, thanks to `@connorferster <https://github.com/connorferster>`_
- Fix bug with material properties not being assigned when using the ``CompoundGeometry.from_points()`` method, thanks to `@connorferster <https://github.com/connorferster>`_

**Full changelog:** `2.0.3...2.1.0 <https://github.com/robbievanleeuwen/section-properties/compare/2.0.3...2.1.0>`_

v2.0.3:
-------

- Add top reinforcement to concrete section library sections
- Add option to specify concrete circle area to ``concrete_circular_section()``
- Update concrete section library to prevent overlapping geometries
- Fix implementation of ``Geometry`` and ``CompoundGeometry`` ``.__sub__()`` method
- Add method to detect overlapping geometry errors and generate warning
- Add option to create coarse mesh (no angle or area constraints)
- Update ``rhino-shapley-interop`` and ``cad_to_shapely`` requirements

**Full changelog:** `2.0.2...2.0.3 <https://github.com/robbievanleeuwen/section-properties/compare/2.0.2...2.0.3>`_

v2.0.2:
-------

- Add circular_section_by_area() in the section library
- Add option to define reinforcement by area rather than diameter for all concrete sections in the section library
- Fix bug in super_t_girder_section() which caused type 5 to be returned in all cases
- Require matplotlib >= 3.4 for CenteredNorm

**Full changelog:** `2.0.1...2.0.2 <https://github.com/robbievanleeuwen/section-properties/compare/2.0.1...2.0.2>`_

v2.0.1:
-------

- Fix issue with library module

**Full changelog:** `2.0.0...2.0.1 <https://github.com/robbievanleeuwen/section-properties/compare/2.0.0...2.0.1>`_

v2.0.0:
-------

*sectionproperties* v2 incorporates significant changes to the pre-processor, which now uses the
`Shapely <https://github.com/shapely/shapely>`_ package to power advanced geometry creation and
manipulation, and vastly improves the performance and robustness of the plastic section property
algorithm. ``v2.x.x`` introduces many breaking changes from ``v1.x.x`` when creating and manipulating
``Geometry``, refer to the `documentation <https://sectionproperties.readthedocs.io>`_ for more
information.

Pre-Processor:
^^^^^^^^^^^^^^

A special mention to `@connorferster <https://github.com/connorferster>`_ for a majority of these
fantastic additions!

- ``sections.py`` renamed to ``geometry.py``
- All ``Geometry`` objects are defined by a shapely ``Polygon``
- Addition of new geometry manipulation methods and geometry set operators
- Added .dxf import, thanks to `@aegis1980 <https://github.com/aegis1980>`_
- Added .3dm import, thanks to `@normanrichardson <https://github.com/normanrichardson>`_
- Introduction of a ``CompoundGeometry`` class for geometries with multiple regions
- ``Geometry`` objects are assigned a ``Material`` property object, ``CompoundGeometry`` objects
  can contain multiple ``Geometry`` objects (each with their own ``Material`` object)
  enabling composite analysis
- ``Geometry`` and ``CompoundGeometry`` objects contain mesh information and meshing must be
  performed before initialising a ``Section`` object
- Improved ``.offset_perimeter()`` logic
- Meshing is now performed by `triangle <https://github.com/drufat/triangle>`_, *meshpy* is no
  longer a dependency
- ``Material`` class now requires a ``.density`` parameter
- The section library (``sectionproperties.pre.library``) now contains the built-in
  *sectionproperties* geometries
- Added ``triangular_section()`` and ``triangular_radius_section()`` to the ``primitive_sections``
  library
- Added ``concrete_sections`` library - contains ``concrete_rectangular_section()``,
  ``concrete_tee_section()`` and ``concrete_circular_section()``
- Added ``bridge_section`` library, thanks to `@ccaprani <https://github.com/ccaprani>`_ - contains
  ``super_t_girder_section()`` and ``i_girder_section()``

Analysis:
^^^^^^^^^

- ``cross_section.py`` renamed to ``section.py``
- ``CrossSection`` object renamed to ``Section`` and is now initialised with only a ``Geometry`` or
  ``CompoundGeometry`` object
- Added calculation of cross-section mass
- Added calculation of weighted material properties - E_eff, G_eff, nu_eff
- The plastic algorithm is now performed by shapely, improving performance and robustness
- Added calculation of principal stresses, thanks to `@ccaprani <https://github.com/ccaprani>`_
- Shape factors are no longer calculated for composite sections (irrelevant property)

Post-Processor:
^^^^^^^^^^^^^^

- Added the ``plotting_context()`` manager, allowing easily saving files, passing kwargs to ``pyplot.subplots()``
  and much more! Many thanks to `@Spectre5 <https://github.com/Spectre5>`_
- Improved contour plotting behaviour
- Added plotting of Mohr's circle of stresses for any given point, thanks to
  `@ccaprani <https://github.com/ccaprani>`_
- ``.display_results()`` now reports E.J and E.Iw instead of G.J and G.Iw
- ``.display_results()`` now reports modulus weighted shear areas for composite sections

Misc.:
^^^^^^

- Many spelling and code style fixes, thanks to `@Spectre5 <https://github.com/Spectre5>`_
- Updated documentation to include theoretical background
- Updated examples to be performed by sphinx-gallery, thanks to
  `@normanrichardson <https://github.com/normanrichardson>`_ and
  `@Spectre5 <https://github.com/Spectre5>`_

v1.0.8:
-------

- All plots now return figure and axes objects
- Fix bug in frame_properties causing the program to crash under certain circumstances

v1.0.7:
-------

- Fix bug with geometry cleaning algorithm resulting in an infinite loop
- Added NASTRAN sections (JohnDN90)
- Added tight_layout to plots (Nils Wagner)
- Added BoxGirderSection class
- Added cross-section perimeter calculation
- Added perimeter offset method (BETA)
- Added mesh refinement example to docs

v1.0.6:
-------

- Fix bug with geometry cleaning algorithm resulting in an infinite loop
- Added NASTRAN sections (JohnDN90)
- Added tight_layout to plots (Nils Wagner)
- Added BoxGirderSection class
- Added cross-section perimeter calculation
- Added perimeter offset method (BETA)
- Added mesh refinement example to docs

v1.0.5:
-------

- Added calculation of monosymmetric constants
- Added tapered flange I Section and channel sections
- Added solid elliptical and hollow elliptical sections (BenjaminFraser)
- Added polygonal section (Agent6-6-6)
- Handle zero radius for all section classes; handle r_out < t for relevant sections
- Update Cee and Zed sections to account for short lips

v1.0.4:
-------

- Added a monosymmetric I Section class
- Extend the plastic centroid search range to the entire section
- Remove the pc_region variable from the plastic centroid calculation as it is no longer relevant
- Better verbose output for the plastic centroid calculation

v1.0.3:
-------

- Retrieve cross-section stresses using get_stress()

v1.0.2:
-------

- Fix returns for adding to geometry

v1.0.1:
-------

- Added calculate_frame_properties()
- Added methods for adding points, facets and control points to geometries
- New pypi README file

v1.0.0:
-------

- Initial release.
