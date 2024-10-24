from typing import Any

import numpy as np

def spsolve(
    A,
    b,
    factorize: bool = ...,
    squeeze: bool = ...,
    solver: Any = ...,
    *args,
    **kwargs,
) -> np.ndarray: ...
