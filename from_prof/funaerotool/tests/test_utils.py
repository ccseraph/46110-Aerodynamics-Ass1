
import numpy as np

from funaerotool.utils import (
    generate_circle_contour,
    generate_naca4_contour,
    naca4_parameters_from_code,
    naca4_surfaces,
)


def test_generate_circle_contour_is_closed_and_radius_constant() -> None:
    x, y = generate_circle_contour(n_points=41, radius=1.7)

    assert x.shape == (41,)
    assert y.shape == (41,)
    assert np.isclose(x[0], x[-1])
    assert np.isclose(y[0], y[-1])
    assert np.allclose(np.hypot(x, y), 1.7)


def test_naca4_parameters_from_code_parses_2412() -> None:
    m, p, t = naca4_parameters_from_code("2412")

    assert np.isclose(m, 0.02)
    assert np.isclose(p, 0.4)
    assert np.isclose(t, 0.12)


def test_naca4_parameters_from_code_rejects_invalid_input() -> None:
    import pytest

    with pytest.raises(ValueError):
        naca4_parameters_from_code("241")

    with pytest.raises(ValueError):
        naca4_parameters_from_code("24A2")

    with pytest.raises(ValueError):
        naca4_parameters_from_code("2400")


def test_naca4_surfaces_are_symmetric_for_0012() -> None:
    x_u, y_u, x_l, y_l = naca4_surfaces(m=0.0, p=0.0, t=0.12, n_points=101)

    assert np.allclose(x_u, x_l)
    assert np.allclose(y_u, -y_l)


def test_generate_naca4_contour_is_closed_and_ordered() -> None:
    x, y = generate_naca4_contour(naca_code="2412", n_points=401)

    assert x.shape == (401,)
    assert y.shape == (401,)
    assert np.isfinite(x).all()
    assert np.isfinite(y).all()

    assert np.isclose(x[0], 1.0, atol=5e-3)
    assert np.isclose(x[-1], 1.0, atol=5e-3)
    assert np.isclose(y[0], 0.0, atol=5e-3)
    assert np.isclose(y[-1], 0.0, atol=5e-3)

    le_distance = np.min(np.hypot(x, y))
    assert le_distance < 2e-3


def test_generate_naca4_contour_requires_odd_n_points() -> None:
    import pytest

    with pytest.raises(ValueError):
        generate_naca4_contour(naca_code="0012", n_points=400)
