
import numpy as np

from .panel_method.postprocessing import compute_pressure_coefficient


def _validate_inputs(R: float, U_inf: float) -> None:
    """Validate core cylinder-flow parameters.

    Parameters
    ----------
    R
        Cylinder radius.
    U_inf
        Freestream velocity magnitude.
    """
    if R <= 0.0:
        raise ValueError("Cylinder radius 'R' must be positive.")
    if U_inf <= 0.0:
        raise ValueError("Freestream velocity 'U_inf' must be positive.")


def _resolve_circulation(
    circulation: float | None, R: float, U_inf: float, aoa_deg: float
) -> float:
    """Return circulation or apply Kutta condition when None."""
    if circulation is None:
        return cylinder_circulation_for_kutta_condition(
            R=R, U_inf=U_inf, aoa_deg=aoa_deg
        )
    return circulation


def cylinder_complex_potential(
    z: np.ndarray,
    R: float,
    U_inf: float,
    circulation: float | None = None,
    aoa_deg: float = 0.0,
) -> np.ndarray:
    """Compute the complex potential W(z) for flow around a cylinder.

    Parameters
    ----------
    z
        Complex array of points where the potential is evaluated.
    R
        Cylinder radius.
    U_inf
        Freestream velocity magnitude.
    circulation
        Circulation strength (positive for counterclockwise).
    aoa_deg
        Angle of attack in degrees (freestream flow direction).
    """
    _validate_inputs(R, U_inf)
    circulation = _resolve_circulation(circulation, R=R, U_inf=U_inf, aoa_deg=aoa_deg)
    aoa = np.deg2rad(aoa_deg)
    with np.errstate(divide="ignore", invalid="ignore"):
        W = U_inf * (z * np.exp(-1j * aoa) + (R**2) / z * np.exp(1j * aoa)) + 1j * (
            circulation / (2.0 * np.pi)
        ) * np.log(z / R)
    return W


def cylinder_complex_velocity(
    z: np.ndarray,
    R: float,
    U_inf: float,
    circulation: float | None = None,
    aoa_deg: float = 0.0,
) -> np.ndarray:
    """Compute the complex velocity dW/dz for flow around a cylinder."""
    _validate_inputs(R, U_inf)
    circulation = _resolve_circulation(circulation, R=R, U_inf=U_inf, aoa_deg=aoa_deg)
    aoa = np.deg2rad(aoa_deg)
    with np.errstate(divide="ignore", invalid="ignore"):
        dWdz = (
            U_inf * (np.exp(-1j * aoa) - (R**2) / (z**2) * np.exp(1j * aoa))
            + 1j * (circulation / (2.0 * np.pi)) / z
        )
    return dWdz


def cylinder_flow_field(
    x: np.ndarray,
    y: np.ndarray,
    R: float,
    U_inf: float,
    circulation: float | None = None,
    aoa_deg: float = 0.0,
    mask_inside: bool = True,
    eps: float = 5e-3,
) -> dict:
    """Compute cylinder flow field (ux, uy, speed, Cp) on a grid.

    Returns a dictionary with keys ``ux``, ``uy``, ``u`` (speed), and ``Cp``.
    Points inside the cylinder are masked as NaN when ``mask_inside`` is True.
    """
    z = x + 1j * y
    dWdz = cylinder_complex_velocity(
        z, R=R, U_inf=U_inf, circulation=circulation, aoa_deg=aoa_deg
    )
    if mask_inside:
        inside = np.abs(z) <= R - eps
        dWdz[inside] = np.nan + 1j * np.nan
    ux = np.real(dWdz)
    uy = -np.imag(dWdz)
    speed = np.hypot(ux, uy)
    Cp = 1.0 - (speed / U_inf) ** 2
    return {"ux": ux, "uy": uy, "u": speed, "Cp": Cp}


def cylinder_surface_velocity(
    theta: np.ndarray,
    R: float,
    U_inf: float,
    circulation: float | None = None,
    aoa_deg: float = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute polar surface velocity components (u_r, u_theta) at r=a."""
    u_complex = cylinder_complex_velocity(
        R * np.exp(1j * theta),
        R=R,
        U_inf=U_inf,
        circulation=circulation,
        aoa_deg=aoa_deg,
    )
    u_x = np.real(u_complex)
    u_y = -np.imag(u_complex)
    u_r = u_x * np.cos(theta) + u_y * np.sin(theta)
    u_theta = -u_x * np.sin(theta) + u_y * np.cos(theta)
    return u_r, u_theta


def cylinder_pressure_coefficient_surface(
    theta: np.ndarray,
    R: float,
    U_inf: float,
    circulation: float | None = None,
    aoa_deg: float = 0.0,
) -> np.ndarray:
    """Compute pressure coefficient C_p on the cylinder surface."""
    ux, uy = cylinder_surface_velocity(
        theta, R=R, U_inf=U_inf, circulation=circulation, aoa_deg=aoa_deg
    )
    return compute_pressure_coefficient(ux=ux, uy=uy, U_inf=U_inf)


def cylinder_circulation_for_kutta_condition(
    R: float, U_inf: float, aoa_deg: float
) -> float:
    """Compute circulation that satisfies the Joukowsky/Kutta condition.

    For the classical circle-to-airfoil (Joukowsky) setup, this selects the
    circulation such that a stagnation point is located at the trailing-edge
    pre-image on the circle.

    The relation is
        Gamma = 4 * pi * a * U_inf * sin(aoa)
    where ``aoa`` is the angle of attack in radians.
    """
    _validate_inputs(R, U_inf)
    aoa = np.deg2rad(aoa_deg)
    return float(4.0 * np.pi * R * U_inf * np.sin(aoa))


def cylinder_lift_coefficient(
    R: float, U_inf: float, circulation: float | None = None, aoa_deg: float = 0.0
) -> float:
    """Compute lift coefficient for the cylinder from the circulation."""
    _validate_inputs(R, U_inf)
    circulation = _resolve_circulation(
        circulation=circulation, R=R, U_inf=U_inf, aoa_deg=aoa_deg
    )
    chord = 2.0 * R
    return 2.0 * circulation / (chord * U_inf)
