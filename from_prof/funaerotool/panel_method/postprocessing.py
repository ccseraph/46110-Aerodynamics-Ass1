
import numpy as np

from .freestream import freestream_components
from .induction_matrix import global_panel_induced_velocity_matrices
from .preprocessing import panel_geometry
from .source import point_source_induced_velocity
from .utils import broadcast_float_arrays
from .vortex import point_vortex_induced_velocity


def compute_pressure_coefficient(
    ux: np.ndarray, uy: np.ndarray, U_inf: float = 1.0
) -> np.ndarray:
    """Compute pressure coefficient from velocity components."""
    if U_inf <= 0.0:
        raise ValueError("'U_inf' must be positive.")
    speed = np.hypot(ux, uy)
    return 1.0 - (speed / U_inf) ** 2


def _points_inside_polygon(
    x: np.ndarray, y: np.ndarray, x_poly: np.ndarray, y_poly: np.ndarray
) -> np.ndarray:
    """Return boolean mask of points inside a polygon (ray-casting)."""
    x_flat = x.ravel()
    y_flat = y.ravel()

    poly_x = np.asarray(x_poly, dtype=float).ravel()
    poly_y = np.asarray(y_poly, dtype=float).ravel()

    if poly_x[0] != poly_x[-1] or poly_y[0] != poly_y[-1]:
        poly_x = np.concatenate([poly_x, poly_x[:1]])
        poly_y = np.concatenate([poly_y, poly_y[:1]])

    n_vert = poly_x.size - 1
    inside = np.zeros_like(x_flat, dtype=bool)

    for i in range(n_vert):
        x0, y0 = poly_x[i], poly_y[i]
        x1, y1 = poly_x[i + 1], poly_y[i + 1]

        intersects = (y0 > y_flat) != (y1 > y_flat)
        denom = (y1 - y0) if (y1 - y0) != 0.0 else 1e-300
        x_int = x0 + (x1 - x0) * (y_flat - y0) / denom
        inside ^= intersects & (x_flat < x_int)

    return inside.reshape(x.shape)


def compute_panel_flow_field(
    x_field: np.ndarray,
    y_field: np.ndarray,
    x_contour: np.ndarray,
    y_contour: np.ndarray,
    sigma: np.ndarray | None = None,
    gamma: np.ndarray | None = None,
    U_inf: float = 1.0,
    aoa_deg: float | None = None,
    mask_inside: bool = True,
) -> dict:
    """Compute panel-induced flow field velocities and pressure coefficient on a grid."""
    x_field = np.asarray(x_field, dtype=float)
    y_field = np.asarray(y_field, dtype=float)

    if x_field.shape != y_field.shape:
        raise ValueError("x_field and y_field must share shape.")

    x_contour = np.asarray(x_contour, dtype=float).ravel()
    y_contour = np.asarray(y_contour, dtype=float).ravel()

    x_eval = x_field.ravel()
    y_eval = y_field.ravel()

    ux_total = np.zeros_like(x_eval, dtype=float)
    uy_total = np.zeros_like(y_eval, dtype=float)

    plength, xp, yp, Tx, Ty, Nx, Ny = panel_geometry(x_contour, y_contour)

    if sigma is not None:
        A_global, B_global = global_panel_induced_velocity_matrices(
            x_eval=x_eval,
            y_eval=y_eval,
            xp=xp,
            yp=yp,
            Tx=Tx,
            Ty=Ty,
            Nx=Nx,
            Ny=Ny,
            panel_lengths=plength,
            panel_type="source",
        )

        ux_total = ux_total + A_global @ sigma
        uy_total = uy_total + B_global @ sigma

    if gamma is not None:
        A_global_vort, B_global_vort = global_panel_induced_velocity_matrices(
            x_eval=x_eval,
            y_eval=y_eval,
            xp=xp,
            yp=yp,
            Tx=Tx,
            Ty=Ty,
            Nx=Nx,
            Ny=Ny,
            panel_lengths=plength,
            panel_type="vortex",
        )
        ux_total = ux_total + A_global_vort @ gamma
        uy_total = uy_total + B_global_vort @ gamma

    if aoa_deg is not None:
        u_fs_x, u_fs_y = freestream_components(aoa_deg=aoa_deg, U_inf=U_inf)
        ux_total = ux_total + u_fs_x
        uy_total = uy_total + u_fs_y

    ux_grid = ux_total.reshape(x_field.shape)
    uy_grid = uy_total.reshape(y_field.shape)

    Cp = compute_pressure_coefficient(ux=ux_grid, uy=uy_grid, U_inf=U_inf)

    if mask_inside:
        inside = _points_inside_polygon(
            x_field, y_field, x_poly=x_contour, y_poly=y_contour
        )
        ux_grid = ux_grid.copy()
        uy_grid = uy_grid.copy()
        Cp = Cp.copy()
        ux_grid[inside] = np.nan
        uy_grid[inside] = np.nan
        Cp[inside] = np.nan
    u = np.hypot(ux_grid, uy_grid)
    return {"ux": ux_grid, "uy": uy_grid, "Cp": Cp, "u": u}


def compute_point_flow_field(
    x_field: np.ndarray,
    y_field: np.ndarray,
    Sigma: np.ndarray | None = None,
    x_sigma: np.ndarray | None = None,
    y_sigma: np.ndarray | None = None,
    Gamma: np.ndarray | None = None,
    x_gamma: np.ndarray | None = None,
    y_gamma: np.ndarray | None = None,
    U_inf: float = 1.0,
    aoa_deg: float | None = None,
) -> dict:
    """Compute flow field induced by point sources/vortices plus optional freestream."""

    x_field = np.asarray(x_field, dtype=float)
    y_field = np.asarray(y_field, dtype=float)

    if x_field.shape != y_field.shape:
        raise ValueError("x_field and y_field must share shape.")

    x_eval = x_field.ravel()
    y_eval = y_field.ravel()

    ux_total = np.zeros_like(x_eval, dtype=float)
    uy_total = np.zeros_like(y_eval, dtype=float)

    if Sigma is not None:
        if x_sigma is None or y_sigma is None:
            raise ValueError("x_sigma and y_sigma are required when Sigma is provided.")

        s_b, xs_b, ys_b = broadcast_float_arrays(Sigma, x_sigma, y_sigma)

        for s_val, xs_val, ys_val in zip(s_b.ravel(), xs_b.ravel(), ys_b.ravel()):
            ux_p, uy_p = point_source_induced_velocity(
                x=x_eval, y=y_eval, x_source=xs_val, y_source=ys_val, strength=s_val
            )
            ux_total = ux_total + ux_p
            uy_total = uy_total + uy_p

    if Gamma is not None:
        if x_gamma is None or y_gamma is None:
            raise ValueError("x_gamma and y_gamma are required when Gamma is provided.")

        g_b, xg_b, yg_b = broadcast_float_arrays(Gamma, x_gamma, y_gamma)

        for g_val, xg_val, yg_val in zip(g_b.ravel(), xg_b.ravel(), yg_b.ravel()):
            ux_p, uy_p = point_vortex_induced_velocity(
                x=x_eval,
                y=y_eval,
                x_vortex=xg_val,
                y_vortex=yg_val,
                circulation=g_val,
            )
            ux_total = ux_total + ux_p
            uy_total = uy_total + uy_p

    if aoa_deg is not None:
        u_fs_x, u_fs_y = freestream_components(aoa_deg=aoa_deg, U_inf=U_inf)
        ux_total = ux_total + u_fs_x
        uy_total = uy_total + u_fs_y

    ux_grid = ux_total.reshape(x_field.shape)
    uy_grid = uy_total.reshape(y_field.shape)

    Cp = compute_pressure_coefficient(ux=ux_grid, uy=uy_grid, U_inf=U_inf)

    u = np.hypot(ux_grid, uy_grid)
    return {"ux": ux_grid, "uy": uy_grid, "Cp": Cp, "u": u}
