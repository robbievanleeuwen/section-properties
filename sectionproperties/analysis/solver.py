import time
import numpy as np
from rich.progress import Progress, BarColumn, ProgressColumn, TextColumn, SpinnerColumn
from rich.table import Column
from rich.text import Text
from scipy.sparse import linalg
from scipy.sparse.linalg import spsolve


def solve_cgs(k, f, m=None, tol=1e-5):
    """Solves a linear system of equations (Ku = f) using the CGS iterative method.

    :param k: N x N matrix of the linear system
    :type k: :class:`scipy.sparse.csc_matrix`
    :param f: N x 1 right hand side of the linear system
    :type f: :class:`numpy.ndarray`
    :param float tol: Tolerance for the solver to achieve. The algorithm terminates when either
        the relative or the absolute residual is below tol.
    :param m: Preconditioner for the linear matrix approximating the inverse of k
    :type m: :class:`scipy.linalg.LinearOperator`

    :return: The solution vector to the linear system of equations
    :rtype: :class:`numpy.ndarray`

    :raises RuntimeError: If the CGS iterative method does not converge
    """

    (u, info) = linalg.cgs(k, f, tol=tol, M=m)

    if info != 0:
        raise RuntimeError("CGS iterative method did not converge.")

    return u


def solve_cgs_lagrange(k_lg, f, tol=1e-5, m=None):
    """Solves a linear system of equations (Ku = f) using the CGS iterative method and the
    Lagrangian multiplier method.

    :param k: (N+1) x (N+1) Lagrangian multiplier matrix of the linear system
    :type k: :class:`scipy.sparse.csc_matrix`
    :param f: N x 1 right hand side of the linear system
    :type f: :class:`numpy.ndarray`
    :param float tol: Tolerance for the solver to achieve. The algorithm terminates when either
        the relative or the absolute residual is below tol.
    :param m: Preconditioner for the linear matrix approximating the inverse of k
    :type m: :class:`scipy.linalg.LinearOperator`

    :return: The solution vector to the linear system of equations
    :rtype: :class:`numpy.ndarray`

    :raises RuntimeError: If the CGS iterative method does not converge or the error from the
        Lagrangian multiplier method exceeds the tolerance
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

    :param k: N x N matrix of the linear system
    :type k: :class:`scipy.sparse.csc_matrix`
    :param f: N x 1 right hand side of the linear system
    :type f: :class:`numpy.ndarray`

    :return: The solution vector to the linear system of equations
    :rtype: :class:`numpy.ndarray`
    """

    return spsolve(k, f)


def solve_direct_lagrange(k_lg, f):
    """Solves a linear system of equations (Ku = f) using the direct solver method and the
    Lagrangian multiplier method.

    :param k: (N+1) x (N+1) Lagrangian multiplier matrix of the linear system
    :type k: :class:`scipy.sparse.csc_matrix`
    :param f: N x 1 right hand side of the linear system
    :type f: :class:`numpy.ndarray`

    :return: The solution vector to the linear system of equations
    :rtype: :class:`numpy.ndarray`

    :raises RuntimeError: If the Lagrangian multiplier method exceeds a tolerance of 1e-5
    """

    u = spsolve(k_lg, np.append(f, 0))

    # compute error
    err = u[-1] / max(np.absolute(u))

    if err > 1e-5:
        err = "Lagrangian multiplier method error exceeds tolerance of 1e-5."
        raise RuntimeError(err)

    return u[:-1]


class CustomTimeElapsedColumn(ProgressColumn):
    """Renders time elapsed in milliseconds."""

    def render(self, task: "Task") -> Text:
        """Show time remaining."""

        elapsed = task.finished_time if task.finished else task.elapsed

        if elapsed is None:
            return Text("-:--:--", style="progress.elapsed")

        elapsed_string = "[ {0:.4f} s ]".format(elapsed)

        return Text(elapsed_string, style="progress.elapsed")


def create_progress():
    """Returns a Rich Progress class."""

    return Progress(
        SpinnerColumn(),
        TextColumn(
            "[progress.description]{task.description}", table_column=Column(ratio=1)
        ),
        BarColumn(bar_width=None, table_column=Column(ratio=1)),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        CustomTimeElapsedColumn(),
        expand=True,
    )
