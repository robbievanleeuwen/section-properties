from collections.abc import Callable
from typing import Any

def njit(
    cache: bool,
    nogil: bool,
) -> Callable[[Any], Any]: ...
