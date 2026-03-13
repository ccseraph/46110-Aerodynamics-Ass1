
import numpy as np


def freestream_components(aoa_deg: float, U_inf: float = 1.0) -> tuple[float, float]:
    """Return freestream components (u_x, u_y) for given angle of attack."""
    if U_inf <= 0.0:
        raise ValueError("'U_inf' must be positive.")
    alpha = np.deg2rad(aoa_deg)
    return float(U_inf * np.cos(alpha)), float(U_inf * np.sin(alpha))
