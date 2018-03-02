'''
The main module contains a function for analysing the cross-section
(crossSectionAnalysis), and a function for applying loads and generating the
resulting stress diagrams (stressAnalysis)
'''

import mesh2D
import femFunctions


def crossSectionAnalysis(points, facets, holes, meshSize, nu, plots=True,
                         output=True):
    '''
    This function takes mesh data as input and executes a cross-section
    analaysis, printing both warping independent and warpring dependent
    properties to the console. A triMesh object is returned with unit stresses
    evaluated and stored to facilitate quick visualisation of design stresses.
    Boolean plots configures whether or not plots are drawn and boolean output
    configures whether or not output is printed to the console.
    '''
    # plot input geometry
    if plots:
        mesh2D.plotGeometry(points, facets, holes)

    # --------------------------------------------------------------------------
    # GEOMTRIC SECTION PROPERTIES:
    # --------------------------------------------------------------------------

    # genereate a coarse triangular mesh to caclulate geometric properties
    coarseMesh = femFunctions.createMesh(points, facets, holes)

    # create triMesh object from the coarse mesh
    meshGeometric = mesh2D.triMesh(coarseMesh)

    # plot the coarse mesh
    if plots:
        (meshGeometric.contourPlot(nodes=True,
                                   plotTitle='Mesh for Geometric Properties'))

    # compute geometric section properties and print to the console
    meshGeometric.computeGeometricProperties()
    if output:
        meshGeometric.printGeometricResults()

    # compute plastic section properties and print to the console
    meshGeometric.computeGlobalPlasticProperties(points, facets, holes)
    meshGeometric.computePrincipalPlasticProperties(points, facets, holes)
    if output:
        meshGeometric.printPlasticResults()

    # plot mesh again to allow the printed results to be interpreted with
    # reference to the coarse mesh
    if plots:
        (meshGeometric.contourPlot(nodes=True,
                                   plotTitle='Mesh for Geometric Properties'))

    # --------------------------------------------------------------------------
    # WARPING DEPENDENT SECTION PROPERTIES:
    # --------------------------------------------------------------------------

    # shift input co-ordinates such that centroid lies at (0,0)
    (shiftedPoints, shiftedHoles) = mesh2D.shiftGeometry(
        points, holes, meshGeometric.cx, meshGeometric.cy)

    # genereate a triangular mesh to calculate warping dependent properties
    refinedMesh = femFunctions.createMesh(
        shiftedPoints, facets, shiftedHoles, maxArea=meshSize)

    # create triMesh object from the refined mesh
    meshWarping = mesh2D.triMesh(refinedMesh, nu, geometricMesh=meshGeometric)

    # plot the refined mesh
    if plots:
        meshWarping.contourPlot(plotTitle='Mesh for Warping Properties')

    # compute warping dependent section properties and print to the console
    meshWarping.computeWarpingProperties(output=output)
    if output:
        meshWarping.printWarpingResults()

    # --------------------------------------------------------------------------
    # UNIT CROSS-SECTION STRESSES:
    # --------------------------------------------------------------------------
    # determine stresses due to unit forces/moments/shears
    meshWarping.unitStress(output=output)

    # finally plot the refined mesh with the centroids indicated
    if plots:
        (meshWarping.contourPlot(plotTitle='Centroids', principalAxis=True,
                                 centroids=True))

    # return the triMesh object to allow stress post-processing
    return meshWarping


def stressAnalysis(warpingMesh, Nzz=0, Mxx=0, Myy=0, M11=0, M22=0, Mzz=0,
                   Vx=0, Vy=0):
    '''
    This function generates plots for all stress componenets and combinations
    as a result of the input design actions.
    '''

    # apply forces/moments/shears
    warpingMesh.evaluateSectionStress(Nzz, Mxx, Myy, M11, M22, Mzz, Vx, Vy)

    # 1. Axial
    warpingMesh.contourPlot(z=warpingMesh.axialStress,
                            plotTitle='Axial Stress')

    # 2. Bending
    (warpingMesh.contourPlot(z=warpingMesh.bendingStress,
                             plotTitle='Bending Stress'))

    # 3. Torsion
    (warpingMesh.contourPlot(z=warpingMesh.torsionStress,
                             plotTitle='Torsion Stress'))
    (warpingMesh.quiverPlot(warpingMesh.torsionStress_zx,
                            warpingMesh.torsionStress_zy,
                            plotTitle='Torsion Stress Vectors'))

    # 4. Shear
    (warpingMesh.contourPlot(z=warpingMesh.shearStress_zx,
                             plotTitle='Transverse Shear (zx) Stress'))
    (warpingMesh.contourPlot(z=warpingMesh.shearStress_zy,
                             plotTitle='Transverse Shear (zy) Stress'))
    (warpingMesh.contourPlot(z=warpingMesh.shearStress,
                             plotTitle='Transverse Shear Stress'))
    (warpingMesh.quiverPlot(warpingMesh.shearStress_zx,
                            warpingMesh.shearStress_zy,
                            plotTitle='Transverse Shear Stress Vectors'))

    # 5a. Combined Normal Stress
    (warpingMesh.contourPlot(z=warpingMesh.sigma_zz,
                             plotTitle='Combined Normal Stress'))

    # 5b. Combined Shear Stress
    (warpingMesh.contourPlot(z=warpingMesh.tau,
                             plotTitle='Combined Shear Stress'))
    (warpingMesh.quiverPlot(warpingMesh.tau_zx, warpingMesh.tau_zy,
                            plotTitle='Combined Shear Stress Vectors'))

    # 6. von Mises
    (warpingMesh.contourPlot(z=warpingMesh.vonMises,
                             plotTitle='von Mises Stress'))
