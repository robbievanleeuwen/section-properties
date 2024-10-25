from collections.abc import Callable

import numpy as np
from scipy.sparse import coo_matrix, csc_matrix

class LinearOperator:
    def __init__(
        self,
        shape: tuple,
        matvec: Callable,
    ) -> None: ...

def spsolve(
    A: np.ndarray | csc_matrix | coo_matrix,
    b: np.ndarray | csc_matrix | coo_matrix,
    permc_spec: str = ...,
    use_umfpack: bool = ...,
) -> np.ndarray: ...
def cgs(
    A: np.ndarray | csc_matrix | coo_matrix,
    b: np.ndarray,
    x0: np.ndarray = ...,
    rtol: float = ...,
    atol: float = ...,
    maxiter: int = ...,
    M: LinearOperator = ...,
    callback: Callable = ...,
) -> tuple[np.ndarray, int]: ...
def spilu(A: csc_matrix | coo_matrix) -> SuperLU: ...

class SuperLU:
    def solve(
        self,
        rhs: np.ndarray | tuple,
    ) -> np.ndarray | tuple: ...
