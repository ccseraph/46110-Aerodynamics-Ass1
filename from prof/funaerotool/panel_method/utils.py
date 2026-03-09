
from typing import Union

import numpy as np

ScalarOrArray = Union[float, np.ndarray]


def _maybe_scalar(value: np.ndarray) -> ScalarOrArray:
    return float(value) if value.shape == () else value


def broadcast_float_arrays(*values: ScalarOrArray) -> list[np.ndarray]:
    """Broadcast scalars/arrays to a common shape as float ndarrays."""

    arrays = [np.asarray(v, dtype=float) for v in values]
    return list(np.broadcast_arrays(*arrays))
