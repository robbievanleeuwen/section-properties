'''
This example script shows how to perform a cross-section analysis using two
different methods for defining the geometry. The typical process for
cross-section analysis is as follows:
    1) Define geometry to be meshed through a list of points, facets and holes.
    2) Perform a cross section analysis on the geometry, given a mesh size (max.
    area) and poissons ratio (nu).
    3) Apply any combination of design actions to generate the resulting stress
    contour/vector plots. Note that this step is optional.
'''

# import main and sectionGenerator modules
import main
import sectionGenerator

# ------------------------------------------------------------------------------
# Method 1: Generate geometry using the section generator.
# ------------------------------------------------------------------------------

# Generate 200 mm x 100 mm x 9 mm RHS section, with a 15 mm external radius
# modelled by eight nodes.
(points1, facets1, holes1) = sectionGenerator.RHS(200, 100, 9, 15, 8)

# Perform a cross-section analysis with a max. mesh area of 5 mm^2 and poissons
# ratio (nu) of 0.
mesh1 = main.crossSectionAnalysis(points1, facets1, holes1, meshSize=5, nu=0)

# Plot stresses for one set of design actions
main.stressAnalysis(mesh1, Nzz=10e3, Mxx=50e6, Myy=0, M11=0, M22=0, Mzz=5e6,
    Vx=0, Vy=30e3)

# Plot stresses for another set of design actions
main.stressAnalysis(mesh1, Nzz=-10e3, Mxx=0, Myy=20e6, M11=0, M22=0, Mzz=5e6,
    Vx=10e3, Vy=0)

# ------------------------------------------------------------------------------
# Method 2: Generate geometry by defining points, facets and holes.
# ------------------------------------------------------------------------------

# Generate an asymmetric I-section with a width and depth of 100 mm.
points2 = ([(-10,0), (110,0), (100,10), (55,10), (55,90), (100,90), (110,100),
    (110,110), (-10,110), (-10,100), (0, 90), (45, 90), (45,10), (-10,10)])
facets2 = ([(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (7,8), (8,9),
    (9,10), (10,11), (11,12), (12,13), (13,0)])
holes2 = []

# Perform a cross-section analysis with a max. mesh area of 5 mm^2 and poissons
# ratio (nu) of 0.3.
mesh2 = main.crossSectionAnalysis(points2, facets2, holes2, meshSize=5, nu=0.3)

# Plot design stresses
main.stressAnalysis(mesh2, Nzz=10e3, Mxx=50e6, Myy=0, M11=0, M22=0, Mzz=5e6,
    Vx=0, Vy=30e3)
