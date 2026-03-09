
import numpy as np


def generate_circle_contour(
    n_points: int = 51, radius: float = 1.0
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a closed circular contour with points from 0 to 2*pi.

    Parameters
    ----------
    n_points
        Number of contour points (includes repeated start/end point).
    radius
        Circle radius.
    """
    if n_points < 3:
        raise ValueError("'n_points' must be at least 3.")
    if radius <= 0.0:
        raise ValueError("'radius' must be positive.")

    t = np.linspace(0.0, 2.0 * np.pi, n_points)
    x = radius * np.cos(t)
    y = radius * np.sin(t)
    return x, y


def naca4_parameters_from_code(code: str) -> tuple[float, float, float]:
    """Parse a 4-digit NACA code into ``(m, p, t)`` fractions.

    Examples
    --------
    - ``"2412"`` -> ``m=0.02``, ``p=0.4``, ``t=0.12``
    - ``"0012"`` -> ``m=0.00``, ``p=0.0``, ``t=0.12``
    """
    if len(code) != 4 or not code.isdigit():
        raise ValueError("'code' must be a 4-digit string, e.g. '2412'.")

    m = int(code[0]) / 100.0
    p = int(code[1]) / 10.0
    t = int(code[2:]) / 100.0

    if t <= 0.0:
        raise ValueError("NACA thickness must be positive.")
    if m > 0.0 and not (0.0 < p < 1.0):
        raise ValueError(
            "For cambered airfoils, camber location must satisfy 0 < p < 1."
        )

    return m, p, t


def naca4_surfaces(
    m: float,
    p: float,
    t: float,
    n_points: int = 201,
    closed_te: bool = True,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Generate NACA 4-series upper/lower surfaces over ``x in [0, 1]``.

    Returns ``(x_upper, y_upper, x_lower, y_lower)`` where both surfaces are
    ordered from leading edge to trailing edge.
    """
    if n_points < 3:
        raise ValueError("'n_points' must be at least 3.")
    if not (0.0 <= m <= 0.09):
        raise ValueError("'m' must be in [0, 0.09].")
    if m > 0.0 and not (0.0 < p < 1.0):
        raise ValueError("For cambered airfoils, 'p' must satisfy 0 < p < 1.")
    if not (0.0 < t <= 0.4):
        raise ValueError("'t' must be in (0, 0.4].")

    beta = np.linspace(0.0, np.pi, n_points)
    x = 0.5 * (1.0 - np.cos(beta))

    te_coeff = -0.1036 if closed_te else -0.1015
    y_t = (
        5.0
        * t
        * (
            0.2969 * np.sqrt(x)
            - 0.1260 * x
            - 0.3516 * x**2
            + 0.2843 * x**3
            + te_coeff * x**4
        )
    )

    y_c = np.zeros_like(x)
    dyc_dx = np.zeros_like(x)

    if m > 0.0:
        left = x < p
        right = ~left

        y_c[left] = m / p**2 * (2.0 * p * x[left] - x[left] ** 2)
        dyc_dx[left] = 2.0 * m / p**2 * (p - x[left])

        y_c[right] = (
            m / (1.0 - p) ** 2 * ((1.0 - 2.0 * p) + 2.0 * p * x[right] - x[right] ** 2)
        )
        dyc_dx[right] = 2.0 * m / (1.0 - p) ** 2 * (p - x[right])

    theta = np.arctan(dyc_dx)

    x_upper = x - y_t * np.sin(theta)
    y_upper = y_c + y_t * np.cos(theta)
    x_lower = x + y_t * np.sin(theta)
    y_lower = y_c - y_t * np.cos(theta)

    return x_upper, y_upper, x_lower, y_lower


def generate_naca4_contour(
    naca_code: str = "2412",
    n_points: int = 401,
    closed_te: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a closed NACA 4-series contour from a 4-digit code.

    The contour is ordered from trailing edge (upper surface) to leading edge,
    then back to trailing edge along the lower surface.
    """
    if n_points < 5 or n_points % 2 == 0:
        raise ValueError("'n_points' must be an odd integer >= 5.")

    m, p, t = naca4_parameters_from_code(naca_code)
    n_surface = (n_points + 1) // 2
    x_u, y_u, x_l, y_l = naca4_surfaces(
        m=m,
        p=p,
        t=t,
        n_points=n_surface,
        closed_te=closed_te,
    )

    x = np.concatenate((x_u[::-1], x_l[1:]))
    y = np.concatenate((y_u[::-1], y_l[1:]))

    # Normalize x so LE sits at 0 and TE at 1 exactly (within floating error)
    x_min = float(np.min(x))
    x_max = float(np.max(x))
    span = x_max - x_min
    if span <= 0.0:
        raise ValueError("Generated contour has zero x-span; check inputs.")
    x = (x - x_min) / span
    return x, y
