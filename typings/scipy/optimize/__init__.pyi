from collections.abc import Callable
from typing import Any

def brentq(
    f: Callable[[float, Any, Any, Any], float],
    a: float,
    b: float,
    args: tuple[Any, ...] = ...,
    xtol: float = ...,
    rtol: float = ...,
    maxiter: int = ...,
    full_output: bool = ...,
    disp: bool = ...,
) -> tuple[float, RootResults]: ...

class RootResults:
    root: float
    iterations: int
    function_calls: int
    converged: bool
    flag: str
    method: str
