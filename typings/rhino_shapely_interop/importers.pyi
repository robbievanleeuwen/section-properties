from collections.abc import Iterator

import numpy as np
from shapely.geometry import Polygon

class RhImporter:
    @classmethod
    def from_file(
        cls,
        file_name: str,
    ) -> RhImporter: ...
    @classmethod
    def from_serialzed_brep(
        cls,
        s_brep: str,
    ) -> RhImporter: ...
    def get_planer_brep(
        self,
        refine_num: int = ...,
        tol: float = ...,
        vec1: np.ndarray = ...,
        vec2: np.ndarray = ...,
        plane_distance: float = ...,
        project: bool = ...,
        parallel: bool = ...,
    ) -> Iterator[Polygon]: ...
