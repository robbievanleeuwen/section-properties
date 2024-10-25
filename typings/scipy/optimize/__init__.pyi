from collections.abc import Callable

def brentq(
    f: Callable,
    a: float,
    b: float,
    args: tuple = ...,
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
