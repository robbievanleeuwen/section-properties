import mesh2D
import femFunctions

# ------------------------------------------------------------------------------
# INPUT GEOMETRY:
# ------------------------------------------------------------------------------
# channel section
points = ([(-0.05,-0.05), (1,-0.05), (1,0.05), (0.05,0.05), (0.05, 1.95), (1, 1.95),
    (1, 2.05), (-0.05, 2.05)])
facets = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (7,0)]
holes = []
maxSize = 0.0005

# ------------------------------------------------------------------------------
# GEOMTRIC SECTION PROPERTIES:
# ------------------------------------------------------------------------------
# genereate coarse triangular mesh for geometric properties
coarseMesh = femFunctions.createMesh(points, facets, holes)

 # create mesh2D object for geometric properties
meshGeometric = mesh2D.triMesh(coarseMesh)

# plot mesh
meshGeometric.contourPlot(nodes=True, plotTitle='Mesh for Geometric Properties')

# compute geometric section properties and print
meshGeometric.computeGeometricProperties()
meshGeometric.printGeometricResults()

# compute plastic section properties and print
meshGeometric.computeGlobalPlasticProperties(points, facets, holes)
meshGeometric.printPlasticResults()

# ------------------------------------------------------------------------------
# WARPING DEPENDENT SECTION PROPERTIES:
# ------------------------------------------------------------------------------
# shift input co-ordinates
(shiftedPoints, shiftedHoles) = femFunctions.shiftGeometry(points, holes, meshGeometric.cx, meshGeometric.cy)

# genereate refined triangular mesh for warping independent properties
refinedMesh = femFunctions.createMesh(shiftedPoints, facets, shiftedHoles, maxArea=maxSize)

 # create mesh2D object for warping properties
meshWarping = mesh2D.triMesh(refinedMesh, nu=0, geometricMesh=meshGeometric)

# plot mesh
meshWarping.contourPlot(plotTitle='Mesh for Warping Properties', principalAxis=True)

# compute section properties
meshWarping.computeWarpingProperties()
meshWarping.printWarpingResults()

# ------------------------------------------------------------------------------
# CROSS-SECTION STRESSES:
# ------------------------------------------------------------------------------
# 1. Axial
meshWarping.axialStress(1e3)
meshWarping.contourPlot(z=meshWarping.sigma_zz_axial, plotTitle='Axial Stress')

# 2a. Bending (Global)
meshWarping.bendingGlobalStress(1e6, 1e6)
meshWarping.contourPlot(z=meshWarping.sigma_zz_bending, plotTitle='Bending Global Stress')

# 2a. Bending (Principal)
meshWarping.bendingPrincipalStress(1e6, 1e6)
meshWarping.contourPlot(z=meshWarping.sigma_zz_bending, plotTitle='Bending Principal Stress', principalAxis=True)

# 3. Torsion
meshWarping.torsionStress(1e6)
meshWarping.contourPlot(z=meshWarping.tau_torsion, plotTitle='Torsion Stress')
meshWarping.quiverPlot(meshWarping.tau_zx_torsion, meshWarping.tau_zy_torsion, plotTitle='Torsion Stress Vectors')

# 4. Shear
meshWarping.shearStress(1e3, 1e3)
meshWarping.contourPlot(z=meshWarping.tau_zx_shear, plotTitle='Transverse Shear (zx) Stress')
meshWarping.contourPlot(z=meshWarping.tau_zy_shear, plotTitle='Transverse Shear (zy) Stress')
meshWarping.contourPlot(z=meshWarping.tau_shear, plotTitle='Transverse Shear Stress')
meshWarping.quiverPlot(meshWarping.tau_zx_shear, meshWarping.tau_zy_shear, plotTitle='Transverse Shear Stress Vectors')

# 5a. Combined Normal Stress
meshWarping.combinedNormalStress()
meshWarping.contourPlot(z=meshWarping.sigma_zz, plotTitle='Combined Normal Stress')

# 5b. Combined Shear Stress
meshWarping.combinedShearStress()
meshWarping.contourPlot(z=meshWarping.tau, plotTitle='Combined Shear Stress')

# 6. von Mises
meshWarping.vonMisesStress()
meshWarping.contourPlot(z=meshWarping.vonMises, plotTitle='von Mises Stress')
