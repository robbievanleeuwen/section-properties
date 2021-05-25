"""Solvers."""

import time

import numpy as np
from scipy.sparse import linalg
from scipy.sparse.linalg import spsolve


def solve_cgs(k, f, m=None, tol=1e-5):
    """Solves a linear system of equations (Ku = f) using the CGS iterative method.

    Parameters
    ----------
    k : :class:`scipy.sparse.csc_matrix`
        N x N matrix of the linear system
    f : :class:`numpy.ndarray`
        N x 1 right hand side of the linear system
    tol : float
        Tolerance for the solver to achieve. The algorithm terminates when either the relative or
        the absolute residual is below tol.
    m : :class:`scipy.sparse.linalg.LinearOperator`
        Preconditioner for the linear matrix approximating the inverse of k

    Returns
    -------
    :class:`numpy.ndarray`
        The solution vector to the linear system of equations

    Raises
    ------
    RuntimeError
        If the CGS iterative method does not converge
    """
    (u, info) = linalg.cgs(k, f, tol=tol, M=m)

    if info != 0:
        raise RuntimeError("CGS iterative method did not converge.")

    return u


def solve_cgs_lagrange(k_lg, f, tol=1e-5, m=None):
    """Solves a linear system of equations (Ku = f) using CGS and Lagrangian multipliers.

    Solves a linear system of equations (Ku = f) using the CGS iterative method and Lagrangian
    multipliers.

    Parameters
    ----------
    k_lg : :class:`scipy.sparse.csc_matrix`
        (N+1) x (N+1) Lagrangian multiplier matrix of the linear system
    f : :class:`numpy.ndarray`
        N x 1 right hand side of the linear system
    tol : float
        Tolerance for the solver to achieve. The algorithm terminates when either the relative or
        the absolute residual is below tol.
    m : :class:`scipy.sparse.linalg.LinearOperator`
        Preconditioner for the linear matrix approximating the inverse of k

    Returns
    -------
    :class:`numpy.ndarray`
        The solution vector to the linear system of equations

    Raises
    ------
    RuntimeError
        If the CGS iterative method does not converge or the error from the Lagrangian multiplier
        method exceeds the tolerance
    """
    (u, info) = linalg.cgs(k_lg, np.append(f, 0), tol=tol, M=m)

    if info != 0:
        raise RuntimeError("CGS iterative method did not converge.")

    # compute error
    err = u[-1] / max(np.absolute(u))

    if err > tol:
        err = "Lagrangian multiplier method error exceeds tolerance."
        raise RuntimeError(err)

    return u[:-1]


def solve_direct(k, f):
    """Solves a linear system of equations (Ku = f) using the direct solver method.

    Parameters
    ----------
    k : :class:`scipy.sparse.csc_matrix`
        N x N matrix of the linear system
    f : :class:`numpy.ndarray`
        N x 1 right hand side of the linear system

    Returns
    -------
    :class:`numpy.ndarray`
        The solution vector to the linear system of equations
    """
    return spsolve(k, f)


def solve_direct_lagrange(k_lg, f):
    """Solves a linear system of equations (Ku = f) using direct solver and Lagrangian multipliers.

    Solves a linear system of equations (Ku = f) using the direct solver method and Lagrangian
    multipliers.

    Parameters
    ----------
    k : :class:`scipy.sparse.csc_matrix`
        (N+1) x (N+1) Lagrangian multiplier matrix of the linear system
    f : :class:`numpy.ndarray`
        N x 1 right hand side of the linear system

    Returns
    -------
    :class:`numpy.ndarray`
        The solution vector to the linear system of equations

    Raises
    ------
    RuntimeError
        If the Lagrangian multiplier method exceeds a tolerance of 1e-5
    """
    u = spsolve(k_lg, np.append(f, 0))

    # compute error
    err = u[-1] / max(np.absolute(u))

    if err > 1e-5:
        err = "Lagrangian multiplier method error exceeds tolerance of 1e-5."
        raise RuntimeError(err)

    return u[:-1]


def function_timer(text, function, *args, **kwargs):
    r"""Displays a message and returns the time taken to run a function.

    Parameters
    ----------
    text : str
        Message to display
    function : :func:
        Function to time and execute
    \*args
        Positional arguments passed to the function
    \**kwargs
        Keyword arguments passed to the function

    Returns
    -------
    Any
        Value returned from the function
    """
    start_time = time.time()

    if text != "":
        print(text)

    result = function(*args, **kwargs)

    if text != "":
        print("----completed in {0:.6f} seconds---".format(time.time() - start_time))

    return result
