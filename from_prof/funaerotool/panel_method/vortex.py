
import numpy as np

from .source import point_source_induced_velocity, source_panel_induced_velocity_local
from .utils import ScalarOrArray, _maybe_scalar


def point_vortex_induced_velocity(
    x: ScalarOrArray,
    y: ScalarOrArray,
    x_vortex: ScalarOrArray,
    y_vortex: ScalarOrArray,
    circulation: ScalarOrArray = 1.0,
) -> tuple[ScalarOrArray, ScalarOrArray]:
    """Velocity induced by a point vortex at arbitrary evaluation points."""

    # Vortex velocity is a 90° rotation of the radial source velocity
    ux_src, uy_src = point_source_induced_velocity(
        x=x, y=y, x_source=x_vortex, y_source=y_vortex, strength=circulation
    )

    ux_src_arr = np.asarray(ux_src, dtype=float)
    uy_src_arr = np.asarray(uy_src, dtype=float)

    ux = -uy_src_arr
    uy = ux_src_arr

    return _maybe_scalar(ux), _maybe_scalar(uy)


def constant_vortex_distribution(panel_lengths: np.ndarray) -> np.ndarray:
    """Constant vortex distribution normalized by panel lengths.

    Returns a per-panel array ``gamma_hat`` such that
    ``sum(gamma_hat * panel_lengths) == 1``.
    """
    lengths = np.asarray(panel_lengths, dtype=float)
    if lengths.ndim != 1:
        raise ValueError("'panel_lengths' must be a 1D array.")
    if lengths.size < 1:
        raise ValueError("At least one panel is required.")
    if np.any(lengths <= 0.0):
        raise ValueError("All panel lengths must be positive.")

    total_length = float(np.sum(lengths))
    if np.isclose(total_length, 0.0):
        raise ValueError("Degenerate panel lengths encountered.")
    return np.full(lengths.shape, 1.0 / total_length, dtype=float)


def parabolic_vortex_distribution(panel_lengths: np.ndarray) -> np.ndarray:
    """Parabolic vortex distribution used in the provided MATLAB reference."""
    lengths = np.asarray(panel_lengths, dtype=float)
    if lengths.ndim != 1:
        raise ValueError("'panel_lengths' must be a 1D array.")
    n = lengths.size
    if n < 2:
        raise ValueError("At least two panels are required.")
    if np.any(lengths <= 0.0):
        raise ValueError("All panel lengths must be positive.")

    idx = np.arange(n, dtype=float)
    numer = idx * (n - 1.0 - idx)
    denom = np.sum(numer * lengths)
    if np.isclose(denom, 0.0):
        raise ValueError("Degenerate panel weighting encountered.")
    return numer / denom


def vortex_panel_induced_velocity_local(
    sx: ScalarOrArray, sy: ScalarOrArray, panel_length: ScalarOrArray
) -> tuple[ScalarOrArray, ScalarOrArray]:
    """Induced velocity from a unit-strength constant vortex panel in local frame.

    Self-induction is returned when the evaluation point is at the panel
    control point (sx and sy both < 1e-6 in magnitude).
    """

    ut_source, un_source = source_panel_induced_velocity_local(
        sx=sx,
        sy=sy,
        panel_length=panel_length,
    )

    return -un_source, ut_source
