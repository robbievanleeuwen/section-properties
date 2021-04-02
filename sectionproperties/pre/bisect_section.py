from typing import Tuple, Union

import numpy as np
from shapely.geometry import Point, LineString

def create_line_segment(
    point_on_line: Union[Tuple[float, float], np.ndarray],
    vector: np.ndarray, 
    bounds: tuple,
    axis: str,
    ):
    """
    Return a LineString of a line that contains 'point_on_line' in the direction of 'unit_vector' 
    bounded by 'bounds'.
    'bounds' is a tuple of float containing a max ordinate and min ordinate.
    'axis' is a str either "x" or "y" and corresponds to the axis of which 'bounds' are given
    """
    c_x, c_y = point_on_line
    b_2 = max(bounds)
    b_1 = min(bounds)
    if axis == "x":
        scale_factor_2 = (b_2 - c_x) / vector[0]
        y_2 = scale_factor_2 * vector[1] + c_y
        
        scale_factor_1 = (b_1 - c_x) / vector[0]
        y_1 = scale_factor_1 * vector[1] + c_y
        return LineString([(b_1, y_1), (b_2, y_2)])
    else:
        scale_factor_2 = (b_2 - c_y) / vector[1]
        x_2 = scale_factor_2 * vector[0] + c_x
        
        scale_factor_1 = (b_1 - c_y) / vector[1]
        x_1 = scale_factor_1 * vector[0] + c_x
        return LineString([(x_1, b_1), (x_2, b_2)])