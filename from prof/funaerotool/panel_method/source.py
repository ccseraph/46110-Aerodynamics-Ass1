
import numpy as np

from .utils import ScalarOrArray, _maybe_scalar


def point_source_induced_velocity(
    x: ScalarOrArray,
    y: ScalarOrArray,
    x_source: ScalarOrArray,
    y_source: ScalarOrArray,
    strength: ScalarOrArray = 1.0,
) -> tuple[ScalarOrArray, ScalarOrArray]:
    """Velocity induced by a point source at arbitrary evaluation points."""

    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    x0_arr = np.asarray(x_source, dtype=float)
    y0_arr = np.asarray(y_source, dtype=float)
    strength_arr = np.asarray(strength, dtype=float)

    x_b, y_b, x0_b, y0_b, strength_b = np.broadcast_arrays(
        x_arr, y_arr, x0_arr, y0_arr, strength_arr
    )

    dx = x_b - x0_b
    dy = y_b - y0_b
    r2 = dx**2 + dy**2
    if np.any(r2 == 0.0):
        raise ValueError("Evaluation point coincides with source location.")

    coeff = strength_b / (2.0 * np.pi * r2)
    ux = coeff * dx
    uy = coeff * dy

    return _maybe_scalar(ux), _maybe_scalar(uy)


def source_panel_induced_velocity_local(
    sx: ScalarOrArray, sy: ScalarOrArray, panel_length: ScalarOrArray
) -> tuple[ScalarOrArray, ScalarOrArray]:
    """Induced velocity from a unit-strength source panel in its local frame.

    The equations are those used in the provided MATLAB script (Eqs. 16/18 in the note).

    Self-induction is returned when the evaluation point is at the panel
    control point (sx and sy both < 1e-6 in magnitude).
    """
    sx_arr = np.asarray(sx, dtype=float)
    sy_arr = np.asarray(sy, dtype=float)
    panel_len_arr = np.asarray(panel_length, dtype=float)

    if np.any(panel_len_arr <= 0.0):
        raise ValueError("'panel_length' must be positive.")

    sx_b, sy_b, panel_len_b = np.broadcast_arrays(sx_arr, sy_arr, panel_len_arr)

    ut = np.empty_like(sx_b, dtype=float)
    un = np.empty_like(sx_b, dtype=float)

    self_mask = (np.abs(sx_b) < 1e-6) & (np.abs(sy_b) < 1e-6)
    if np.any(self_mask):
        ut[self_mask] = 0.0
        un[self_mask] = 0.5

    other_mask = ~self_mask
    if np.any(other_mask):
        half = 0.5 * panel_len_b[other_mask]
        sx_o = sx_b[other_mask]
        sy_o = sy_b[other_mask]
        num = (sx_o + half) ** 2 + sy_o**2
        den = (sx_o - half) ** 2 + sy_o**2

        ut[other_mask] = np.log(num / den) / (4.0 * np.pi)
        un[other_mask] = (
            np.arctan2(sy_o, sx_o - half) - np.arctan2(sy_o, sx_o + half)
        ) / (2.0 * np.pi)

    return _maybe_scalar(ut), _maybe_scalar(un)
