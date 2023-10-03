"""Methods used for solving linear systems and displaying info on tasks."""

from __future__ import annotations

import numpy as np
from rich.progress import (
    BarColumn,
    Progress,
    ProgressColumn,
    SpinnerColumn,
    Task,
    TextColumn,
)
from rich.table import Column
from rich.text import Text
from scipy.sparse import csc_matrix, linalg
from scipy.sparse.linalg import LinearOperator, spsolve


def solve_cgs(
    k: csc_matrix,
    f: np.ndarray,
    m: LinearOperator | None = None,
    tol: float = 1e-5,
) -> np.ndarray:
    """Solves a linear system using the CGS iterative method.

    Args:
        k: ``N x N`` matrix of the linear system
        f: ``N x 1`` right hand side of the linear system
        m: Preconditioner for the linear matrix approximating the inverse of ``k``
        tol: Tolerance for the solver to achieve. The algorithm terminates when either
            the relative or the absolute residual is below ``tol``

    Returns:
        The solution vector to the linear system of equations

    Raises:
        RuntimeError: If the CGS iterative method does not converge
    """
    u, info = linalg.cgs(A=k, b=f, tol=tol, M=m)

    if info != 0:
        raise RuntimeError("CGS iterative method did not converge.")

    return u


def solve_cgs_lagrange(
    k_lg: csc_matrix,
    f: np.ndarray,
    m: LinearOperator | None = None,
    tol: float = 1e-5,
) -> np.ndarray:
    """Solves a linear system using the CGS iterative method (Lagrangian multiplier).

    Args:
        k_lg: ``(N+1) x (N+1)`` Lagrangian multiplier matrix of the linear system
        f: ``N x 1`` right hand side of the linear system
        m: Preconditioner for the linear matrix approximating the inverse of ``k``
        tol: Tolerance for the solver to achieve. The algorithm terminates when either
            the relative or the absolute residual is below ``tol``

    Returns:
        The solution vector to the linear system of equations

    Raises:
        RuntimeError: If the CGS iterative method does not converge or the error from
            the Lagrangian multiplier method exceeds the tolerance
    """
    u, info = linalg.cgs(A=k_lg, b=np.append(f, 0), tol=tol, M=m)

    if info != 0:
        raise RuntimeError("CGS iterative method did not converge.")

    # compute error
    err = u[-1] / max(np.absolute(u))

    if err > tol:
        raise RuntimeError("Lagrangian multiplier method error exceeds tolerance.")

    return u[:-1]


def solve_direct(
    k: csc_matrix,
    f: np.ndarray,
) -> np.ndarray:
    """Solves a linear system using the direct solver method.

    Args:
        k: ``N x N`` matrix of the linear system
        f: ``N x 1`` right hand side of the linear system

    Returns:
        The solution vector to the linear system of equations
    """
    return spsolve(A=k, b=f)


def solve_direct_lagrange(
    k_lg: csc_matrix,
    f: np.ndarray,
) -> np.ndarray:
    """Solves a linear system using the direct solver method (Lagrangian multiplier).

    Args:
        k_lg: ``(N+1) x (N+1)`` Lagrangian multiplier matrix of the linear system
        f: ``N x 1`` right hand side of the linear system

    Returns:
        The solution vector to the linear system of equations

    Raises:
        RuntimeError: If the Lagrangian multiplier method exceeds a tolerance of
            ``1e-5``
    """
    u = spsolve(A=k_lg, b=np.append(f, 0))

    # compute error
    multiplier = abs(u[-1])
    rel_error = multiplier / max(np.absolute(u))

    if rel_error > 1e-7 and multiplier > 10.0 * np.finfo(float).eps:
        msg = "Lagrangian multiplier method error exceeds tolerance of 1e-5."
        raise RuntimeError(msg)

    return u[:-1]


class CustomTimeElapsedColumn(ProgressColumn):
    """Renders time elapsed in milliseconds."""

    def render(
        self,
        task: Task,
    ) -> Text:
        """Show time remaining.

        Args:
            task: Rich progress task

        Returns:
            Rich text object
        """
        elapsed = task.finished_time if task.finished else task.elapsed

        if elapsed is None:
            return Text("-:--:--", style="progress.elapsed")

        elapsed_string = f"[ {elapsed:.4f} s ]"

        return Text(elapsed_string, style="progress.elapsed")


def create_progress() -> Progress:
    """Returns a Rich Progress class.

    Returns:
        Rich Progress class containing a spinner, progress description, percentage and
        time
    """
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
