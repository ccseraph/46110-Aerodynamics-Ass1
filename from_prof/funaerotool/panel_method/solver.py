import numpy as np

from funaerotool.panel_method.freestream import freestream_components
from funaerotool.panel_method.induction_matrix import (
    global_panel_induced_velocity_matrices,
)
from funaerotool.panel_method.postprocessing import compute_pressure_coefficient
from funaerotool.panel_method.preprocessing import panel_geometry
from funaerotool.panel_method.transformations import global_to_local
from funaerotool.panel_method.vortex import parabolic_vortex_distribution


def solve_closed_contour_panel_method(
    x: np.ndarray,
    y: np.ndarray,
    aoa_deg: float = 0.0,
    U_inf: float = 1.0,
    kutta_condition: bool = True,
) -> dict:
    """Solve a closed-contour panel system with optional Kutta circulation.

    Parameters
    ----------
    x, y : np.ndarray
        Closed contour coordinates (first point == last point).
    aoa_deg : float
        Angle of attack in degrees.
    U_inf : float
        Freestream speed magnitude.
    kutta_condition : bool
        When True, solve a coupled system with a parabolic vortex distribution
        to enforce a Kutta condition.
    """
    # Flatten inputs to 1D arrays so geometry math stays vectorized
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()

    # Panel geometry: midpoints (xp, yp), tangents (Tx, Ty), normals (Nx, Ny), lengths
    plength, xp, yp, Tx, Ty, Nx, Ny = panel_geometry(x, y)
    n_panels = plength.size

    # Influence matrices for source panels evaluated at panel control points
    # A_g_source[i,j] = x-wise velocity at panel i due to panel j
    # B_g_source[i,j] = y-wise velocity at panel i due to panel j
    A_g_source, B_g_source = global_panel_induced_velocity_matrices(
        x_eval=xp,
        y_eval=yp,
        xp=xp,
        yp=yp,
        Tx=Tx,
        Ty=Ty,
        Nx=Nx,
        Ny=Ny,
        panel_lengths=plength,
        panel_type="source",
    )

    # Transform induced velocities into each panel's local frame
    # A_l_source[i,j] = normal velocity at panel i due to panel j
    # B_l_source[i,j] = tangential velocity at panel i due to panel j
    B_l_source, A_l_source = global_to_local(
        ux=A_g_source,
        uy=B_g_source,
        tx=Tx[:, None],
        ty=Ty[:, None],
        nx=Nx[:, None],
        ny=Ny[:, None],
    )

    # Freestream components in global coordinates
    u_inf_x, u_inf_y = freestream_components(aoa_deg=aoa_deg, U_inf=U_inf)

    # Freestream components in each panel's local frame
    u_inf_t, u_inf_n = global_to_local(
        ux=u_inf_x,
        uy=u_inf_y,
        tx=Tx,
        ty=Ty,
        nx=Nx,
        ny=Ny,
    )

    if kutta_condition is False:
        # Source-only solve: A * σ = F  →  σ = A⁻¹ F
        sigma = np.linalg.solve(A_l_source, -u_inf_n)
        circulation = 0.0
        gamma = None
    else:
        # Parabolic vortex distribution used to satisfy the Kutta condition
        gamma_hat = parabolic_vortex_distribution(plength)

        # Influence matrices for vortex panels evaluated at the same control points
        A_g_vortex, B_g_vortex = global_panel_induced_velocity_matrices(
            x_eval=xp,
            y_eval=yp,
            xp=xp,
            yp=yp,
            Tx=Tx,
            Ty=Ty,
            Nx=Nx,
            Ny=Ny,
            panel_lengths=plength,
            panel_type="vortex",
        )

        # Transform vortex-induced velocities into local panel frames
        B_l_vortex, A_l_vortex = global_to_local(
            ux=A_g_vortex,
            uy=B_g_vortex,
            tx=Tx[:, None],
            ty=Ty[:, None],
            nx=Nx[:, None],
            ny=Ny[:, None],
        )

        # Assemble coupled source-vortex linear system enforcing Kutta condition
        lhs = np.zeros((n_panels + 1, n_panels + 1), dtype=float)
        lhs[:n_panels, :n_panels] = A_l_source
        lhs[:n_panels, n_panels] = A_l_vortex @ gamma_hat
        lhs[n_panels, :n_panels] = B_l_source[0, :] + B_l_source[-1, :]
        lhs[n_panels, n_panels] = (B_l_vortex[0, :] + B_l_vortex[-1, :]) @ gamma_hat

        # Right-hand side: no-penetration plus Kutta tangency constraint
        rhs = np.zeros(n_panels + 1, dtype=float)
        rhs[:n_panels] = -u_inf_n
        rhs[-1] = -(u_inf_t[0] + u_inf_t[-1])

        # Solve for source strengths and circulation, then recover distributed gamma
        solution = np.linalg.solve(lhs, rhs)
        sigma = solution[:n_panels]
        circulation = solution[-1]
        gamma = circulation * gamma_hat

    # Total tangential and normal velocities including optional circulation
    Vt = B_l_source @ sigma + u_inf_t
    Vn = A_l_source @ sigma + u_inf_n
    if gamma is not None:
        Vt = Vt + B_l_vortex @ gamma
        Vn = Vn + A_l_vortex @ gamma

    # Pressure coefficient from local velocities
    Cp = compute_pressure_coefficient(ux=Vt, uy=Vn, U_inf=U_inf)

    # Lift coefficient from Kutta circulation; chord inferred from x-span
    chord = np.max(x) - np.min(x)
    Cl = -2.0 * circulation / (chord * U_inf)

    # Bundle outputs for downstream analysis and plotting
    return {
        "sigma": sigma,
        "strengths": sigma,
        "circulation": circulation,
        "gamma": gamma,
        "Vt": Vt,
        "Vn": Vn,
        "Cp": Cp,
        "Cl": Cl,
        "xp": xp,
        "yp": yp,
        "panel_lengths": plength,
        "Tx": Tx,
        "Ty": Ty,
        "Nx": Nx,
        "Ny": Ny,
    }
