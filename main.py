import mesh2D

# ------------------------------------------------------------------------------
# INPUT GEOMETRY:
# ------------------------------------------------------------------------------
# rectangular hollow section 1
points = [(0,0), (50,0), (50,100), (0,100), (6,6), (44, 6), (44, 94), (6, 94)]
facets = [(0,1), (1,2), (2,3), (3,0), (4,5), (5,6), (6,7), (7,4)]
holes = [(25,50)]
maxSize = 2.5

# ------------------------------------------------------------------------------
# SIMPLE SECTION PROPERTIES:
# ------------------------------------------------------------------------------
# genereate coarse triangular mesh for warping independent properties
coarseMesh = mesh2D.createMesh(points, facets, holes)

 # create mesh2D object for simple properties
meshSimple = mesh2D.triMesh(coarseMesh)

# plot mesh
meshSimple.contourPlot(principalAxis = False, nodes = True, plotTitle = 'Mesh for Simple Properties')

# compute simple section properties and print
meshSimple.computeSimpleProperties()
meshSimple.printSimpleResults()

# shift input co-ordinates
(shiftedPoints, shiftedHoles) = mesh2D.shiftGeometry(points, facets, holes, meshSimple.cx, meshSimple.cy)

# ------------------------------------------------------------------------------
# WARPING DEPENDENT SECTION PROPERTIES:
# ------------------------------------------------------------------------------
# genereate refined triangular mesh for warping independent properties
refinedMesh = mesh2D.createMesh(shiftedPoints, facets, shiftedHoles, maxArea = maxSize)

 # create mesh2D object for warping properties
meshWarping = mesh2D.triMesh(refinedMesh, 0)

# plot mesh
meshWarping.contourPlot(principalAxis = False, nodes = True, plotTitle = 'Mesh for Warping Properties')

# compute section properties
meshWarping.computeWarpingProperties(meshSimple.area, meshSimple.ixx_c, meshSimple.iyy_c, meshSimple.ixy_c)
meshWarping.printWarpingResults()

# plot results
meshWarping.contourPlot(principalAxis = False, z = meshWarping.omega, nodes = False, plotTitle = 'Warping Function')
meshWarping.contourPlot(principalAxis = False, z = meshWarping.Psi, nodes = False, plotTitle = 'Shear Function Psi')
meshWarping.contourPlot(principalAxis = False, z = meshWarping.Phi, nodes = False, plotTitle = 'Shear Function Phi')

# mesh.contourPlot(principalAxis = True, z = mesh.tau_torsion, nodes = True)
# mesh.quiverPlot(mesh.tau_zx_torsion, mesh.tau_zy_torsion)
# mesh.contourPlot(False, mesh.tau_zy_shear)
# mesh.contourPlot(False, mesh.tau_shear)
