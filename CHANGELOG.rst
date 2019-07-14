Change Log:
===========

v1.0.5:
-------

- Added calculation of monosymmetric constants
- Added tapered flange I-section and channel sections
- Added solid elliptical and hollow elliptcal sections (BenjaminFraser)
- Added polygonal section (Agent6-6-6)
- Handle zero radius for all section classes; handle r_out < t for relevant sections
- Update Cee and Zed sections to account for short lips

v1.0.4:
-------

- Added a monosymmetric I-section class
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
