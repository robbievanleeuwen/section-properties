from typing import Any

import numpy as np
import numpy.typing as npt
from scipy.sparse import coo_matrix, csc_matrix

def spsolve(
    A: npt.NDArray[np.float64] | csc_matrix | coo_matrix,
    b: npt.NDArray[np.float64] | csc_matrix | coo_matrix,
) -> npt.NDArray[np.float64]: ...
