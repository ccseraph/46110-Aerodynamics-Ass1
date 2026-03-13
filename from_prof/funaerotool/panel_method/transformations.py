
import numpy as np

from .utils import ScalarOrArray, _maybe_scalar


def local_to_global(
    ut: ScalarOrArray,
    un: ScalarOrArray,
    tx: ScalarOrArray,
    ty: ScalarOrArray,
    nx: ScalarOrArray,
    ny: ScalarOrArray,
) -> tuple[ScalarOrArray, ScalarOrArray]:
    """Transform velocity from panel-local to global frame (scalar or array)."""
    ut_a, un_a, tx_a, ty_a, nx_a, ny_a = (
        np.asarray(arr, dtype=float) for arr in (ut, un, tx, ty, nx, ny)
    )
    ux = ut_a * tx_a + un_a * nx_a
    uy = ut_a * ty_a + un_a * ny_a
    return _maybe_scalar(ux), _maybe_scalar(uy)


def global_to_local(
    ux: ScalarOrArray,
    uy: ScalarOrArray,
    tx: ScalarOrArray,
    ty: ScalarOrArray,
    nx: ScalarOrArray,
    ny: ScalarOrArray,
) -> tuple[ScalarOrArray, ScalarOrArray]:
    """Project global velocity into panel-local tangential/normal components."""
    ux_a, uy_a, tx_a, ty_a, nx_a, ny_a = (
        np.asarray(arr, dtype=float) for arr in (ux, uy, tx, ty, nx, ny)
    )
    ut = ux_a * tx_a + uy_a * ty_a
    un = ux_a * nx_a + uy_a * ny_a
    return _maybe_scalar(ut), _maybe_scalar(un)
