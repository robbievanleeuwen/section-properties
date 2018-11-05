# import sys
# import json
# import inspect
#
# import matplotlib.pyplot as plt
#
# from sectionUtilities import sectionParse
# from otherUtilities import plotGeometry, CrossSectionSettings, LoadData
# from femUtilities import createMesh
# from CrossSectionAnalysis import CrossSectionAnalysis
#
#
# def crossSectionAnalysis(sectionTypes, sectionData, meshSizes, materials,
#                          loads, settings=None):
#     """
#     This is the main method that performs the cross-section analysis.
#     """
#
#     # check length of section inputs
#     if any(len(lst) != len(sectionTypes) for lst in [
#             sectionData, meshSizes, materials]):
#         print("Error: Mismatch in number of input sections.")
#         quit()
#
#     # load settings
#     programSettings = CrossSectionSettings(settings)
#
#     # print the loaded settings
#     if (programSettings.outputSettings):
#         programSettings.printSettings()
#
#     # generate loads from input data
#     loadData = LoadData(loads)
#
#     # generate geometry data from input data
#     (points, facets, holes, controlPoints) = sectionParse(
#         sectionTypes, sectionData, programSettings)
#
#     # check loaded geometry before proceeding
#     if (programSettings.checkGeometry):
#         # plot the input geometry
#         plotGeometry(points, facets, holes, controlPoints)
#
#         # ask the user if they want to proceed with the loaded geometry
#         while True:
#             r = str(input("Would you like to proceed with meshing? (y/n): "))
#             if r.lower().strip()[:1] == "y":
#                 plt.close("all")
#                 break
#             elif r.lower().strip()[:1] == "n":
#                 quit()
#
#     # create the finite element mesh
#     mesh = createMesh(points, facets, holes, controlPoints, meshSizes,
#                       minAngle=30, meshOrder=2, qualityMeshing=True,
#                       volumeConstraints=True, settings=programSettings)
#
#     # create the analysis object
#     analysis = CrossSectionAnalysis(mesh, materials, programSettings)
#
#     # shift the mesh origin to the centroid
#     (points, holes, controlPoints) = analysis.shiftToCentroid(
#         points, holes, controlPoints)
#
#     # check the generated mesh before proceeding
#     if (programSettings.checkMesh):
#         # plot the generated mesh
#         analysis.contourPlot(
#             nodes=False, globalAxis=True,
#             plotTitle='Mesh for Geometric Properties')
#
#         # ask the user if they want to proceed with the generated mesh
#         while True:
#             r = str(input("Would you like to proceed with analysis? (y/n): "))
#             if r.lower().strip()[:1] == "y":
#                 plt.close("all")
#                 break
#             elif r.lower().strip()[:1] == "n":
#                 quit()
#
#     # compute cross-section properties
#     analysis.computeSectionProperties(
#         points, facets, holes, controlPoints, materials)
#
#     # if the user has entered loads, evaluate stresses
#     if (loadData.containsLoads):
#         analysis.evaluateSectionStress(loadData)
#
#     # print the analysis results
#     if (programSettings.outputResults):
#         analysis.printResults(programSettings.numberFormat)
#         pass
#
#     # display plots
#     analysis.plotResults(programSettings.plots)
#
#     # ask the user if they want to quit
#     if __name__ == "__main__":
#         while True:
#             r = str(input("Enter (q) to quit: "))
#             if r.lower().strip()[:1] == "q":
#                 quit()
#
#     return analysis
#
#
# if __name__ == "__main__":
#     # read data from properly formatted json-formatted file
#     try:
#         filename = sys.argv[1]
#         with open(filename) as f:
#             print("-- Reading file '{}'...".format(filename))
#
#             data = json.load(f)  # load file
#
#             sectionTypes = data["sections"]  # load section types
#             sectionData = data["section-data"]  # load section data
#             meshSizes = data["mesh-sizes"]  # load mesh sizes
#             materials = data["materials"]  # load material data
#             loads = data["loads"]  # load load data
#             settings = data["settings"]  # load settings
#
#     except IndexError:
#         print("Usage: python efem.py <filename>")
#         quit()
#
#     except ValueError:
#         print("Error: Incorrectly formatted JSON file. Check your file at " +
#               "jsonlint.com.")
#         quit()
#
#     except KeyError as err:
#         frame = inspect.currentframe()
#         print("{}:{}: error: Key {} not found in input data file.".format(
#             __file__, frame.f_lineno, err))
#         quit()
#
#     # perform cross-section analysis
#     crossSectionAnalysis(sectionTypes, sectionData, meshSizes, materials,
#                          loads, settings)
