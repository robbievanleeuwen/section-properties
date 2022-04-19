r"""
.. _ref_ex_peery_beams:

Symmetric and Unsymmetric Beams in Complex Bending
--------------------------------------------------

Calculate section properties of two different beams
given in examples from 'Aircraft Structures,' by Peery. 
These cases have known results, and the output from 
SectionProperties can be compared for accuracy. These 
examples represent a more rigourous 'proof' against a 
'real' problem. Only results that have values in the 

reference material are tested here.

BibTeX Entry for reference:

    @Book{Peery,
        title = {Aircraft Structures},
        author = {David J. Peery},
        organization = {Pensylvania State University},
        publisher = {McGraw-Hill Book Company},
        year = {1950},
        edition = {First},
        ISBN = {978-0486485805}
    }

"""

# sphinz_gallery_thumbnail_number = 1

from sectionproperties.pre.library import nastran_sections
from sectionproperties.analysis.section import Section

# %%
# Example 1 in Sec. 6.2 (Symmetric Bending)
# This is a symmetric I-section with no lateral supports,
# undergoing pure unidirectional cantilever bending.
# Note that units here are **inches**, to match the text.
# 
# We'll use a very coarse mesh here, to show a conservative
# comparison for accuracy. Theoretically, with more 
# discretization, we would capture the real results more accurately.
geometry = nastran_sections.nastran_i(6,3,3,1,1,1)
geometry = geometry.shift_section(x_offset=0,y_offset=-3)
geometry = geometry.create_mesh(mesh_sizes=[0.25])
section = Section(geometry)
section.plot_mesh()

# %%
# Perform a geometric analysis on the section, and plot properties
# We don't need warping or plastic analysis for these simple checks.
section.calculate_geometric_properties()
section.plot_centroids()

