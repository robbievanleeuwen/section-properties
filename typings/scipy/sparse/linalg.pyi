from collections.abc import Callable
from typing import Any

import numpy as np
import numpy.typing as npt
from scipy.sparse import coo_matrix, csc_matrix

class LinearOperator:
    def __init__(
        self,
        shape: tuple,
        matvec: Callable,
    ) -> None: ...

def spsolve(
    A: npt.NDArray[np.float64] | csc_matrix | coo_matrix,
    b: npt.NDArray[np.float64] | csc_matrix | coo_matrix,
    permc_spec: str = ...,
    use_umfpack: bool = ...,
) -> npt.NDArray[np.float64]: ...
def cgs(
    A: npt.NDArray[np.float64] | csc_matrix | coo_matrix,
    b: npt.NDArray[np.float64],
    x0: npt.NDArray[np.float64] = ...,
    rtol: float = ...,
    atol: float = ...,
    maxiter: int = ...,
    M: LinearOperator = ...,
    callback: Callable[[Any], Any] = ...,
) -> tuple[npt.NDArray[np.float64], int]: ...
def spilu(A: csc_matrix | coo_matrix) -> SuperLU: ...

class SuperLU:
    def solve(
        self,
        rhs: npt.NDArray[np.float64] | tuple[Any],
    ) -> npt.NDArray[np.float64] | tuple[Any]: ...
