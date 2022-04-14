r"""
.. _ref_cad_import:

Importing Geometry from CAD
---------------------------

Demonstrates loading :class:`~sectionproperties.pre.geometry.Geometry` and
:class:`~sectionproperties.pre.geometry.CompoundGeometry` objects from `.dxf` and `.3dm` (Rhino)
files.
"""

# sphinx_gallery_thumbnail_number = 8

from sectionproperties.pre.geometry import Geometry, CompoundGeometry
from sectionproperties.analysis.section import Section

# %%
# Load a geometry with a single region from a dxf file
geom = Geometry.from_dxf(dxf_filepath="files/section_holes.dxf")
geom.plot_geometry()

# %%
# Generate a mesh
geom.create_mesh([1])
sec = Section(geom)
sec.plot_mesh(materials=False)

# %%
# Conduct a geometric & plastic analysis
sec.calculate_geometric_properties()
sec.calculate_plastic_properties()
sec.plot_centroids()

# %%
# Display the geometric & plastic properties
sec.display_results()

# %%
# Load a geometry with multiple holes from a dxf file
geom = Geometry.from_dxf(dxf_filepath="files/section_holes_complex.dxf")
geom.plot_geometry()

# %%
# Generate a mesh
geom.create_mesh([1])
sec = Section(geom)
sec.plot_mesh(materials=False)

# %%
# Conduct a geometric & plastic analysis
sec.calculate_geometric_properties()
sec.calculate_plastic_properties()
sec.plot_centroids()

# %%
# Display the geometric & plastic properties
sec.display_results()

# %%
# Load a geometry from a 3dm (Rhino) file
geom = Geometry.from_3dm(filepath="files/complex_shape.3dm")
geom.plot_geometry()

# %%
# Generate a mesh
geom.create_mesh([1])
sec = Section(geom)
sec.plot_mesh(materials=False)

# %%
# Conduct a geometric & plastic analysis
sec.calculate_geometric_properties()
sec.calculate_plastic_properties()
sec.plot_centroids()

# %%
# Display the geometric & plastic properties
sec.display_results()

# %%
# Load a compound geometry with multiple regions from a 3dm (Rhino) file
geom = CompoundGeometry.from_3dm(filepath="files/compound_shape.3dm")
geom.plot_geometry()

# %%
# Generate a mesh
geom.create_mesh([1])
sec = Section(geom)
sec.plot_mesh(materials=False)

# %%
# Conduct a geometric & plastic analysis
# N.B a warping analysis would be invalid due to the lack of connectivity between the two regions
sec.calculate_geometric_properties()
sec.calculate_plastic_properties()
sec.plot_centroids()

# %%
# Display the geometric & plastic properties
sec.display_results()
