
from typing import Literal

import numpy as np

from .source import source_panel_induced_velocity_local
from .transformations import global_to_local, local_to_global
from .vortex import vortex_panel_induced_velocity_local


def global_panel_induced_velocity_matrices(
    x_eval: np.ndarray,
    y_eval: np.ndarray,
    xp: np.ndarray,
    yp: np.ndarray,
    Tx: np.ndarray,
    Ty: np.ndarray,
    Nx: np.ndarray,
    Ny: np.ndarray,
    panel_lengths: np.ndarray,
    panel_type: Literal["source", "vortex"] = "source",
) -> tuple[np.ndarray, np.ndarray]:
    """Build induced velocity matrices for arbitrary evaluation points."""
    x_eval = np.asarray(x_eval, dtype=float).ravel()
    y_eval = np.asarray(y_eval, dtype=float).ravel()
    if x_eval.size != y_eval.size:
        raise ValueError("'x_eval' and 'y_eval' must have the same length.")

    xp = np.asarray(xp, dtype=float).ravel()
    yp = np.asarray(yp, dtype=float).ravel()
    Tx = np.asarray(Tx, dtype=float).ravel()
    Ty = np.asarray(Ty, dtype=float).ravel()
    Nx = np.asarray(Nx, dtype=float).ravel()
    Ny = np.asarray(Ny, dtype=float).ravel()
    panel_lengths = np.asarray(panel_lengths, dtype=float).ravel()

    n_panels = panel_lengths.size
    if not (xp.size == yp.size == Tx.size == Ty.size == Nx.size == Ny.size == n_panels):
        raise ValueError("All panel arrays must have the same length.")

    if panel_type not in ("source", "vortex"):
        raise ValueError(f"panel_type must be 'source' or 'vortex', got {panel_type!r}")
    # Relative positions of evaluation points to panel control points (broadcasted)
    dx = x_eval[:, None] - xp[None, :]
    dy = y_eval[:, None] - yp[None, :]

    dt, dn = global_to_local(
        ux=dx,
        uy=dy,
        tx=Tx[None, :],
        ty=Ty[None, :],
        nx=Nx[None, :],
        ny=Ny[None, :],
    )

    panel_lengths_b = panel_lengths[None, :]
    if panel_type == "vortex":
        ut, un = vortex_panel_induced_velocity_local(
            sx=dt,
            sy=dn,
            panel_length=panel_lengths_b,
        )
    else:
        ut, un = source_panel_induced_velocity_local(
            sx=dt,
            sy=dn,
            panel_length=panel_lengths_b,
        )

    A_global, B_global = local_to_global(
        ut=ut,
        un=un,
        tx=Tx[None, :],
        ty=Ty[None, :],
        nx=Nx[None, :],
        ny=Ny[None, :],
    )

    return np.asarray(A_global, dtype=float), np.asarray(B_global, dtype=float)
