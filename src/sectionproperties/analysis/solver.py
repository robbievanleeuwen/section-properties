"""Methods used for solving linear systems and displaying info on tasks."""

from __future__ import annotations

from typing import TYPE_CHECKING

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
from scipy.sparse import coo_matrix, csc_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import cgs, splu, spsolve

if TYPE_CHECKING:
    import numpy.typing as npt
    from scipy.sparse.linalg import LinearOperator


try:
    import pypardiso

    sp_solve = pypardiso.spsolve
except ImportError:
    sp_solve = spsolve

HAS_PARDISO = sp_solve is not spsolve


def solve_cgs(
    k: csc_matrix,
    f: npt.NDArray[np.float64],
    m: LinearOperator,
    tol: float = 1e-5,
) -> npt.NDArray[np.float64]:
    """Solves a linear system using the CGS iterative method.

    Args:
        k: ``N x N`` matrix of the linear system
        f: ``N x 1`` right hand side of the linear system
        m: Preconditioner for the linear matrix approximating the inverse of ``k``
        tol: Relative tolerance for the solver to achieve. Defaults to ``1e-5``.

    Returns:
        The solution vector to the linear system of equations

    Raises:
        RuntimeError: If the CGS iterative method does not converge
    """
    u, info = cgs(A=k, b=f, rtol=tol, M=m)

    if info != 0:
        msg = "CGS iterative method did not converge."
        raise RuntimeError(msg)

    return u


def solve_cgs_lagrange(
    k_lg: csc_matrix,
    f: npt.NDArray[np.float64],
    m: LinearOperator,
    tol: float = 1e-5,
) -> npt.NDArray[np.float64]:
    """Solves a linear system using the CGS iterative method (Lagrangian multiplier).

    Args:
        k_lg: ``(N+1) x (N+1)`` Lagrangian multiplier matrix of the linear system
        f: ``N x 1`` right hand side of the linear system
        m: Preconditioner for the linear matrix approximating the inverse of ``k``
        tol: Relative tolerance for the solver to achieve. Defaults to ``1e-5``.

    Returns:
        The solution vector to the linear system of equations

    Raises:
        RuntimeError: If the CGS iterative method does not converge or the error from
            the Lagrangian multiplier method exceeds the tolerance
    """
    u, info = cgs(A=k_lg, b=np.append(f, 0), rtol=tol, M=m)

    if info != 0:
        msg = "CGS iterative method did not converge."
        raise RuntimeError(msg)

    # compute error
    err = u[-1] / max(np.absolute(u))

    if err > tol:
        msg = "Lagrangian multiplier method error exceeds tolerance."
        raise RuntimeError(msg)

    return u[:-1]


def solve_direct(
    k: csc_matrix,
    f: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Solves a linear system using the direct solver method.

    Args:
        k: ``N x N`` matrix of the linear system
        f: ``N x 1`` right hand side of the linear system

    Returns:
        The solution vector to the linear system of equations
    """
    return sp_solve(A=k, b=f)


def solve_direct_lagrange(
    k_lg: csc_matrix,
    f: npt.NDArray[np.float64],
) -> npt.NDArray[np.float64]:
    """Solves a linear system using the direct solver method (Lagrangian multiplier).

    Args:
        k_lg: ``(N+1) x (N+1)`` Lagrangian multiplier matrix of the linear system
        f: ``N x 1`` right hand side of the linear system

    Returns:
        The solution vector to the linear system of equations

    Raises:
        RuntimeError: If the Lagrangian multiplier method exceeds a relative tolerance
            of ``1e-7`` or absolute tolerance related to your machine's floating point
            precision.
    """
    u = sp_solve(A=k_lg, b=np.append(f, 0))

    # compute error
    multiplier = abs(u[-1])
    rel_error = multiplier / max(np.absolute(u))

    if rel_error > 1e-7 and multiplier > 10.0 * np.finfo(float).eps:
        msg = "Lagrangian multiplier method error exceeds the prescribed tolerance, "
        msg += "consider refining your mesh. If this error is unexpected raise an "
        msg += "issue at https://github.com/robbievanleeuwen/section-properties/issues."
        raise RuntimeError(msg)

    return u[:-1]


def lagrangian_matrix(
    row: npt.NDArray[np.int_],
    col: npt.NDArray[np.int_],
    data: npt.NDArray[np.float64],
    c: npt.NDArray[np.float64],
) -> csc_matrix:
    """Builds the Lagrangian multiplier matrix of the warping system.

    Appends the constraint vector ``c`` as an extra row and column to the stiffness
    matrix given by its COO entries.

    Args:
        row: Row indices of the stiffness matrix entries (COO format)
        col: Column indices of the stiffness matrix entries (COO format)
        data: Values of the stiffness matrix entries (COO format), duplicate entries
            are summed
        c: ``N`` constraint vector, the modulus weighted nodal areas

    Returns:
        ``(N+1) x (N+1)`` Lagrangian multiplier matrix in CSC format
    """
    n_size = len(c)
    nodes = np.arange(n_size)
    lagrange = np.full(n_size, n_size)

    k_lg = coo_matrix(
        (
            np.concatenate((data, c, c)),
            (
                np.concatenate((row, nodes, lagrange)),
                np.concatenate((col, lagrange, nodes)),
            ),
        ),
        shape=(n_size + 1, n_size + 1),
        dtype=float,
    )

    return csc_matrix(k_lg, dtype=float)


class WarpingSolver:
    r"""Direct solver for the warping system :math:`Ku = f`, with :math:`c \cdot u = 0`.

    The warping stiffness matrix ``k`` is singular: adding a constant to a warping or
    shear function leaves the equations unchanged. :func:`solve_direct_lagrange`
    removes that freedom with a Lagrangian multiplier, which appends a fully dense row
    and column to an otherwise very sparse matrix. Eliminating that row during the
    factorisation produces a large amount of fill-in.

    This solver fixes a single degree of freedom instead, which leaves a sparse,
    symmetric positive definite system, and shifts every solution by the constant that
    satisfies the constraint. The solution is the same, and the matrix is factorised
    once on construction and reused for every right hand side passed to
    :meth:`solve`.

    The null space of ``k`` contains a constant for every connected region of the
    mesh, so one degree of freedom is fixed per region. A mesh with more than one region
    is not valid for a warping analysis, but the analysis still runs, as it does with
    the Lagrangian multiplier.
    """

    def __init__(
        self,
        row: npt.NDArray[np.int_],
        col: npt.NDArray[np.int_],
        data: npt.NDArray[np.float64],
        c: npt.NDArray[np.float64],
    ) -> None:
        """Inits and factorises the WarpingSolver.

        Args:
            row: Row indices of the stiffness matrix entries (COO format)
            col: Column indices of the stiffness matrix entries (COO format)
            data: Values of the stiffness matrix entries (COO format), duplicate entries
                are summed
            c: ``N`` constraint vector, the modulus weighted nodal areas
        """
        n_size = len(c)
        self.c = c

        k = csc_matrix(
            coo_matrix((data, (row, col)), shape=(n_size, n_size), dtype=float),
            dtype=float,
        )

        # fix the best conditioned degree of freedom of every connected region
        on_diagonal = row == col
        diagonal = np.bincount(
            row[on_diagonal], weights=data[on_diagonal], minlength=n_size
        )
        n_regions, region = connected_components(csgraph=k, directed=False)
        self.pinned = np.array(
            [
                int(np.flatnonzero(region == idx)[np.argmax(diagonal[region == idx])])
                for idx in range(n_regions)
            ]
        )

        # replace the equations of the fixed degrees of freedom with u[dof] = 0
        is_pinned = np.zeros(n_size, dtype=bool)
        is_pinned[self.pinned] = True
        keep = ~(is_pinned[row] | is_pinned[col])

        k_pinned = csc_matrix(
            coo_matrix(
                (
                    np.concatenate((data[keep], np.ones(n_regions))),
                    (
                        np.concatenate((row[keep], self.pinned)),
                        np.concatenate((col[keep], self.pinned)),
                    ),
                ),
                shape=(n_size, n_size),
                dtype=float,
            ),
            dtype=float,
        )

        # PARDISO keeps the factorisation of the last matrix it was given, SuperLU is
        # told the matrix is symmetric positive definite: symmetric ordering and no
        # partial pivoting
        self.k_pinned = k_pinned
        self.lu = (
            None
            if HAS_PARDISO
            else splu(
                k_pinned,
                permc_spec="MMD_AT_PLUS_A",
                diag_pivot_thresh=0.0,
                options={"SymmetricMode": True},
            )
        )

    def solve(
        self,
        f: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """Solves the warping system for one or more right hand sides.

        Args:
            f: ``N`` right hand side of the linear system, or ``N x M`` right hand
                sides solved together

        Returns:
            The solution, of the same shape as ``f``, that satisfies the constraint
        """
        rhs = np.array(f, dtype=float)
        rhs[self.pinned, ...] = 0.0

        if self.lu is None:
            u = np.asarray(sp_solve(A=self.k_pinned, b=rhs), dtype=float)
        else:
            u = np.asarray(self.lu.solve(rhs), dtype=float)

        return u - (self.c @ u) / float(self.c.sum())


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
