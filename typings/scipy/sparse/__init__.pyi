import numpy as np

class csc_matrix:
    def __init__(
        self,
        data: coo_matrix,
        dtype: np.dtype | type,
    ) -> None: ...

class coo_matrix:
    def __init__(
        self,
        data: tuple,
        shape: tuple[int, int],
        dtype: np.dtype | type,
    ) -> None: ...
