
import numpy as np


def _validate_contour(x: np.ndarray, y: np.ndarray) -> None:
    """Validate a 2D contour array pair used for panel discretization."""
    if x.ndim != 1 or y.ndim != 1:
        raise ValueError("Contour coordinates must be 1D arrays.")
    if x.size != y.size:
        raise ValueError("Contour arrays 'x' and 'y' must have equal length.")
    if x.size < 3:
        raise ValueError("At least 3 contour points are required.")


def flip_contour(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return contour arrays reversed end-to-end while preserving shape.

    This is useful to fix airfoil contours that are ordered opposite the
    expected TE->upper->LE->lower->TE convention enforced by _validate_contour.
    """
    x_arr = np.asarray(x, dtype=float).ravel()
    y_arr = np.asarray(y, dtype=float).ravel()
    return x_arr[::-1], y_arr[::-1]


def panel_geometry(x: np.ndarray, y: np.ndarray) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Compute panel geometry vectors from an ordered contour.

    Returns
    -------
    tuple
        (plength, xp, yp, Tx, Ty, Nx, Ny)
    """
    _validate_contour(x, y)

    dx = x[1:] - x[:-1]
    dy = y[1:] - y[:-1]
    plength = np.hypot(dx, dy)

    if np.any(plength <= 0.0):
        raise ValueError("Panels must have positive length.")

    xp = 0.5 * (x[1:] + x[:-1])
    yp = 0.5 * (y[1:] + y[:-1])

    Tx = -dx / plength
    Ty = -dy / plength
    Nx = -Ty
    Ny = Tx

    return plength, xp, yp, Tx, Ty, Nx, Ny
