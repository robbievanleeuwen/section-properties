r"""
.. _ref-ex-peery-beams:

Symmetric and Unsymmetric Beams in Complex Bending
--------------------------------------------------

Calculate section properties of two different beams
given in examples from 'Aircraft Structures,' by Peery. 
These cases have known results, and the output from 
SectionProperties can be compared for accuracy. These 
examples represent a more rigourous 'proof' against a 
'real' problem. Only results that have values in the 
reference material are tested here.

BibTeX Entry for reference::

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

from sectionproperties.pre.library import nastran_sections
from sectionproperties.analysis.section import Section

# %%
# Example 1 in Sec. 6.2 (Symmetric Bending)
# =========================================
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
# We don't need warping analysis for these simple checks,
# but sectionproperties needs them before evaluating stress.
section.calculate_geometric_properties()
section.calculate_warping_properties()
section.plot_centroids()

# %%
# Directly from the example, we know that the 2nd moment of inertia
# resisting the bending is 43.3 in4.
section.section_props.ixx_g

# %%
# From statics, we know the max bending moment on the beam will
# be 80,000 in-lbs. We can apply this moment to the section, and 
# evaluate stress.
moment = 8e5
stress = section.calculate_stress(Mxx=moment)

# %%
# Next we can extract the max stress from the section, and let's
# go ahead and look at the calculated fringe plot. Refer to the 
# stress example for details.
numerical_result = max(stress.get_stress()[0]['sig_zz'])
stress.plot_stress_zz()

# %%
# From the book, and simple statics, we know the max stress is 
# 55,427.3 psi.
numerical_result

# %%
# This example is admittedly more simple, but it's still a nice
# check for the basics on validity of the package.
print(f'Theoretical Result = 55427 [psi]')
print(f'  Numerical Result = {numerical_result:.0f} [psi]')
acc = (55427-numerical_result)/55427
print(f'          Accuracy = {acc:%}')


# %%
# Example 1 in Sec. 7.2. (Unsymmetric Bending)
# ============================================
# Moving on to something a bit more advanced...
# This is an unsymmetric Z-section with no lateral supports.
# Note that units here are **inches**, to match the text.
base_geom = nastran_sections.nastran_zed(4,2,8,12)
base_geom = base_geom.shift_section(-5,-6)
base_geom = base_geom.create_mesh([0.25])
section = Section(base_geom)
section.calculate_geometric_properties()
section.calculate_warping_properties()
section.plot_centroids()

# %%
# Checking each property against the reference text:
props = section.section_props
print('    Property | Theoretical | Numerical')
print(f'    ixx_g    | {693.3:<12.1f}| {props.ixx_g:<.1f}')
print(f'    iyy_g    | {173.3:<12.1f}| {props.iyy_g:<.1f}')
print(f'    ixy_g    | {-240:<12.1f}| {props.ixy_g:<.1f}')
print(f'    i11_c    | {787:<12.1f}| {props.i11_c:<.1f}')
print(f'    i22_c    | {79.5:<12.1f}| {props.i22_c:<.1f}')

# %%
# Secrtion properties all look good, so we can move on to
# some stress analysis. Before we do, we will need a quick
# function to pull stress at a certain point. This is a bit 
# of a hack! Future version of sectionproperties will have 
# a much more robust system for getting stress at an 
# arbitrary location. This particular function will work for
# the locations we need, since we *know* a node will be there.
from typing import Tuple
def get_node(nodes, coord) -> Tuple[int, tuple]:
    '''
    This function will loop over the node list provided,
    finding the index of the coordinates you want.
    Returns the index in the nodes list, and the coords.
    '''
    for index,var in enumerate(nodes):
        if all(var == coord):
            return index, var
        else:
            continue
    
    raise ValueError(f'No node found with coordinates: {coord}')

# %%
# The load applied in the reference is -100,000 in-lbs about the
# x-axis, and 10,000 inb-lbs about the y-axis.
stress = section.calculate_stress(Mxx=-1e5, Myy=1e4)

# %%
# Check stress at location A (see :ref:`docs<label-testing>` page for details)
A = (-5, 4)
text_result = 1210
n, _ = get_node(section.mesh_nodes, A)
numerical_result = stress.get_stress()[0]['sig_zz'][n]
print(text_result, numerical_result)

# %%
# Check stress at location B (see :ref:`docs<label-testing>` page for details)
B = (-5, 6)
text_result = 580
n, _ = get_node(section.mesh_nodes, B)
numerical_result = stress.get_stress()[0]['sig_zz'][n]
print(text_result, numerical_result)

# %%
# Check stress at location C (see :ref:`docs<label-testing>` page for details)
C = (1, 6)
text_result = -2384
n, _ = get_node(section.mesh_nodes, C)
numerical_result = stress.get_stress()[0]['sig_zz'][n]
print(text_result, numerical_result)

# %%
# Looking at total axial stress over the section.
stress.plot_stress_zz()