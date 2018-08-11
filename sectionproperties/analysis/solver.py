import time
import numpy as np
from scipy.sparse import linalg
from scipy.sparse.linalg import spsolve


def solve_cgs(k, f, m=None, tol=1e-5):
    """Solves a linear system of equations (Ku = f) using the CGS iterative
    method.

    :param k: N x N matrix of the linear system
    :type k: :class:`scipy.sparse.csc_matrix`
    :param f: N x 1 right hand side of the linear system
    :type f: :class:`numpy.ndarray`
    :param float tol: Tolerance for the solver to acheieve. The algorithm
        terminates when either the relative or the absolute residual is below
        tol.
    :param m: Preconditioner for the linear matrix approximating the inverse
        of k
    :type m: :class:`scipy.linalg.LinearOperator`

    :return: The solution vector to the linear system of equations
    :rtype: :class:`numpy.ndarray`

    :raises RuntimeError: If the CGS iterative method does not converge
    """

    (u, exit) = linalg.cgs(k, f, tol=tol, M=m)

    if (exit != 0):
        raise RuntimeError("CGS iterative method did not converge.")

    return u


def solve_cgs_lagrange(k_lg, f, tol=1e-5, m=None):
    """Solves a linear system of equations (Ku = f) using the CGS iterative
    method and the Lagrangian multiplier method.

    :param k: (N+1) x (N+1) Lagrangian multiplier matrix of the linear system
    :type k: :class:`scipy.sparse.csc_matrix`
    :param f: N x 1 right hand side of the linear system
    :type f: :class:`numpy.ndarray`
    :param float tol: Tolerance for the solver to acheieve. The algorithm
        terminates when either the relative or the absolute residual is below
        tol.
    :param m: Preconditioner for the linear matrix approximating the inverse
        of k
    :type m: :class:`scipy.linalg.LinearOperator`

    :return: The solution vector to the linear system of equations
    :rtype: :class:`numpy.ndarray`

    :raises RuntimeError: If the CGS iterative method does not converge or the
        error from the Lagrangian multiplier method exceeds the tolerance
    """

    (u, exit) = linalg.cgs(k_lg, np.append(f, 0), tol=tol, M=m)

    if (exit != 0):
        raise RuntimeError("CGS iterative method did not converge.")

    # compute error
    err = u[-1] / max(np.absolute(u))

    if err > tol:
        err = "Lagrangian multiplier method error exceeds tolerance."
        raise RuntimeError(err)

    return u[:-1]


def solve_direct(k, f):
    """Solves a linear system of equations (Ku = f) using the direct solver
    method.

    :param k: N x N matrix of the linear system
    :type k: :class:`scipy.sparse.csc_matrix`
    :param f: N x 1 right hand side of the linear system
    :type f: :class:`numpy.ndarray`

    :return: The solution vector to the linear system of equations
    :rtype: :class:`numpy.ndarray`
    """

    return spsolve(k, f)


def solve_direct_lagrange(k_lg, f):
    """Solves a linear system of equations (Ku = f) using the direct solver
    method and the Lagrangian multiplier method.

    :param k: (N+1) x (N+1) Lagrangian multiplier matrix of the linear system
    :type k: :class:`scipy.sparse.csc_matrix`
    :param f: N x 1 right hand side of the linear system
    :type f: :class:`numpy.ndarray`

    :return: The solution vector to the linear system of equations
    :rtype: :class:`numpy.ndarray`

    :raises RuntimeError: If the Lagrangian multiplier method exceeds a
        tolerance of 1e-5
    """

    u = spsolve(k_lg, np.append(f, 0))

    # compute error
    err = u[-1] / max(np.absolute(u))

    if err > 1e-5:
        err = "Lagrangian multiplier method error exceeds tolerance of 1e-5."
        raise RuntimeError(err)

    return u[:-1]


# def pcAlgorithm(tol, maxIt, u, dmin, dmax, points, facets, holes,
#                 controlPoints, nodes, elements, materials, dir):
#     """
#     Algorithm to find plastic centroid (point at which top area = bot area):
#         INPUT:
#         tol = convergence tolerance
#         maxIt = maximum iterations
#         u = unit vector in direction of axis
#         (dmin,dmax) = distance from centroid to extreme fibre of section
#         points = input points list
#         facets = input facets list
#         holes = input holes list
#         pointArray = np array containing mesh points
#         elementArray = np array containing element vertices
#         dir = 1 or 2 depending on axis direction (x or y; 11 or 22)
#
#         OUTPUT:
#         a_n = perpendicular distance from centroid to p.c.
#     """
#
#     # initialise iteration variables
#     areaConvergence_n = 0  # area convergence of step n
#     areaConvergence_n1 = 0  # area convergence of step n - 1
#     areaConvergence_n2 = 0  # area convergence of step n - 2
#     a_n = 0  # distance from centroid to pc of step n
#     a_n1 = 0  # distance from centroid to pc of step n - 1
#     a_n2 = 0  # distance from centroid to pc of step n - 2
#     iterationCount = 1
#
#     # determine vectors perpendicular to the current axis
#     if (dir == 1):
#         u_perp = np.array([u[1], -u[0]])  # u vector rotated  -90 degrees
#     elif (dir == 2):
#         u_perp = np.array([-u[1], u[0]])  # u vector rotated  90 degrees
#
#     # iterative algorithm
#     while ((abs(areaConvergence_n) > tol or iterationCount < 3) and
#            (iterationCount < maxIt)):
#         if iterationCount < 3:
#             # first two steps to setup secant method:
#             # random number between -0.5 and 0.5 multiplied by 20% of the depth
#             a_n = (np.random.rand() - 0.5) * 0.2 * (dmax - dmin)
#         else:
#             # secant method
#             a_n = (a_n2 * areaConvergence_n1 - a_n1 * areaConvergence_n2) / (
#                 areaConvergence_n1 - areaConvergence_n2)
#
#         # ensure trial axis is within section depth
#         if a_n > dmax:
#             a_n = dmax - 0.1 * abs(dmax - dmin)
#         elif a_n < dmin:
#             a_n = dmin + 0.1 * abs(dmax - dmin)
#
#         # console reporting for debugging purposes
#         # print("a_n = {}".format(a_n))
#         # print("dmin = {}".format(dmin))
#         # print("dmax = {}".format(dmax))
#
#         # determine points (p1,p2) on trial axis
#         p1 = [a_n * u_perp[0], a_n * u_perp[1]]
#         p2 = [p1[0] + u[0], p1[1] + u[1]]
#         # remesh with new trial axis included
#         (newPoints, newFacets) = divideMesh(
#             points.copy(), facets.copy(), nodes, elements,
#             p1[0], p1[1], p2[0], p2[1], abs(dmax - dmin))
#
#         ms = np.ones(len(controlPoints))  # dummy mesh sizes list
#         mesh = createMesh(
#             newPoints, newFacets, holes, controlPoints, meshSizes=ms,
#             minAngle=None, meshOrder=2, qualityMeshing=False,
#             volumeConstraints=False, settings=[])
#
#         # create section analysis object with new mesh
#         section = CrossSectionAnalysis(mesh, materials, settings=[])
#         # section.contourPlot(nodes=True)  # plot for debugging purposes
#
#         # calculate area above and below axis
#         (topA, botA, topCen, botCen) = section.computeAreaSegments(
#             u, p1[0], p1[1])
#
#         # calculate area convergence
#         areaConvergence_n = topA / botA - 1
#         # console reporting for debugging purposes
#         # print("convergence = {}".format(areaConvergence_n))
#         # print("---")
#
#         # update convergence and solution data
#         areaConvergence_n2 = areaConvergence_n1
#         areaConvergence_n1 = areaConvergence_n
#         a_n2 = a_n1
#         a_n1 = a_n
#         iterationCount += 1
#
#     if (abs(areaConvergence_n) > tol):
#         print("WARNING: Plastic centroid algorithm did not converge for the " +
#               "axis in the direction x:{0:.3f}; y:{1:.3f}\n".format(
#                   u[0], u[1]))
#
#     return (a_n, topA, botA, topCen, botCen)


def function_timer(text, function, *args):
    """Displays the message *text* and returns the time taken for a
    function, with arguments *args*, to execute. The value returned by the
    timed function is also returned.

    :param string text: Message to display
    :param function: Function to time and execute
    :type function: function
    :param args: Function arguments
    :return: Value returned from the function
    """

    start_time = time.time()

    if text != "":
        print(text)

    result = function(*args)

    if text != "":
        print("----completed in {0:.6f} seconds---".format(
            time.time() - start_time))

    return result
