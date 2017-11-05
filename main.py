import mesh2D
import femFunctions
import sectionGenerator

# ------------------------------------------------------------------------------
# INPUT GEOMETRY:
# ------------------------------------------------------------------------------
# Cruciform section
(points, facets, holes) = sectionGenerator.Cruciform(200, 200, 8, 12, 8)
maxSize = 2.5

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
meshGeometric.computePrincipalPlasticProperties(points, facets, holes)
meshGeometric.printPlasticResults()

# plot mesh again
meshGeometric.contourPlot(nodes=True, plotTitle='Mesh for Geometric Properties')

# ------------------------------------------------------------------------------
# WARPING DEPENDENT SECTION PROPERTIES:
# ------------------------------------------------------------------------------
# shift input co-ordinates such that centroid lies at (0,0)
(shiftedPoints, shiftedHoles) = (femFunctions.shiftGeometry(points, holes,
    meshGeometric.cx, meshGeometric.cy))

# genereate refined triangular mesh for warping independent properties
refinedMesh = (femFunctions.createMesh(shiftedPoints, facets, shiftedHoles,
    maxArea=maxSize))

 # create mesh2D object for warping properties (load mesh with geom. properties)
meshWarping = mesh2D.triMesh(refinedMesh, nu=0, geometricMesh=meshGeometric)

# plot warping mesh
meshWarping.contourPlot(plotTitle='Mesh for Warping Properties')

# compute section properties
meshWarping.computeWarpingProperties()
meshWarping.printWarpingResults()

# ------------------------------------------------------------------------------
# CROSS-SECTION STRESSES:
# ------------------------------------------------------------------------------
# determine stresses due to unit forces/moments
meshWarping.unitStress()

# plot mesh centroids
(meshWarping.contourPlot(plotTitle='Centroids', principalAxis=True,
    centroids=True))

# apply forces/moments
(meshWarping.evaluateSectionStress(Nzz=0, Mxx=10e6, Myy=0, M11=0, M22=0,
    Mzz=1e6, Vx=0, Vy=50e3))

# 1. Axial
meshWarping.contourPlot(z=meshWarping.axialStress, plotTitle='Axial Stress')

# 2. Bending
meshWarping.contourPlot(z=meshWarping.bendingStress, plotTitle='Bending Stress')

# 3. Torsion
meshWarping.contourPlot(z=meshWarping.torsionStress, plotTitle='Torsion Stress')
(meshWarping.quiverPlot(meshWarping.torsionStress_zx,
    meshWarping.torsionStress_zy, plotTitle='Torsion Stress Vectors'))

# 4. Shear
(meshWarping.contourPlot(z=meshWarping.shearStress_zx,
    plotTitle='Transverse Shear (zx) Stress'))
(meshWarping.contourPlot(z=meshWarping.shearStress_zy,
    plotTitle='Transverse Shear (zy) Stress'))
(meshWarping.contourPlot(z=meshWarping.shearStress,
    plotTitle='Transverse Shear Stress'))
(meshWarping.quiverPlot(meshWarping.shearStress_zx, meshWarping.shearStress_zy,
    plotTitle='Transverse Shear Stress Vectors'))

# 5a. Combined Normal Stress
(meshWarping.contourPlot(z=meshWarping.sigma_zz,
    plotTitle='Combined Normal Stress'))

# 5b. Combined Shear Stress
(meshWarping.quiverPlot(meshWarping.tau_zx, meshWarping.tau_zy,
    plotTitle='Combined Shear Stress Vectors'))
meshWarping.contourPlot(z=meshWarping.tau, plotTitle='Combined Shear Stress')

# 6. von Mises
meshWarping.contourPlot(z=meshWarping.vonMises, plotTitle='von Mises Stress')
