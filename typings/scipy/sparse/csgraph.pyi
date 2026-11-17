import numpy as np
import numpy.typing as npt
from scipy.sparse import csc_matrix

def connected_components(
    csgraph: csc_matrix,
    directed: bool = ...,
) -> tuple[int, npt.NDArray[np.int32]]: ...
