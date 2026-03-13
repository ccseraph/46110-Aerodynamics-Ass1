"""Tests for cylinder potential-flow utilities (complex form).

This suite focuses on the remaining complex-potential API:
- `_validate_inputs`
- `cylinder_complex_potential` / `cylinder_complex_velocity`
- `cylinder_flow_field`
- `surface_velocity`
- `pressure_coefficient_surface` and `_kutta`
- `gamma_for_joukowsky_kutta`
"""


import numpy as np
import pytest

from funaerotool.cylinder_potential_flow import (
    _validate_inputs,
    cylinder_circulation_for_kutta_condition,
    cylinder_complex_potential,
    cylinder_complex_velocity,
    cylinder_flow_field,
    cylinder_lift_coefficient,
    cylinder_pressure_coefficient_surface,
    cylinder_surface_velocity,
)


def test_validate_inputs_accepts_positive_values() -> None:
    _validate_inputs(R=1.0, U_inf=2.0)


@pytest.mark.parametrize(
    ("a", "u_inf"),
    [
        (0.0, 1.0),
        (-1.0, 1.0),
        (1.0, 0.0),
        (1.0, -1.0),
    ],
)
def test_validate_inputs_raises_for_non_positive(a: float, u_inf: float) -> None:
    with pytest.raises(ValueError):
        _validate_inputs(R=a, U_inf=u_inf)


def test_cylinder_complex_velocity_matches_manual_expression() -> None:
    z = np.array([2.0 + 1.0j])
    R = 1.0
    U_inf = 1.3
    circulation = -0.7
    aoa_deg = 12.5
    aoa = np.deg2rad(aoa_deg)

    vel = cylinder_complex_velocity(
        z, R=R, U_inf=U_inf, circulation=circulation, aoa_deg=aoa_deg
    )
    expected = (
        U_inf * (np.exp(-1j * aoa) - (R**2) / (z**2) * np.exp(1j * aoa))
        + 1j * (circulation / (2.0 * np.pi)) / z
    )

    assert np.allclose(vel, expected)


def test_cylinder_flow_field_matches_complex_velocity_and_masks_inside() -> None:
    x = np.array([[0.0, 2.0]])
    y = np.array([[0.0, 0.5]])
    R = 1.0
    U_inf = 1.0
    circulation = 0.3
    aoa_deg = -5.0

    flow = cylinder_flow_field(
        x=x,
        y=y,
        R=R,
        U_inf=U_inf,
        circulation=circulation,
        aoa_deg=aoa_deg,
        mask_inside=True,
    )

    # Inside point masked
    assert np.isnan(flow["ux"][0, 0])
    assert np.isnan(flow["uy"][0, 0])
    assert np.isnan(flow["Cp"][0, 0])

    # Outside point matches direct complex-velocity evaluation
    z_out = x[0, 1] + 1j * y[0, 1]
    vel_out = cylinder_complex_velocity(
        z_out, R=R, U_inf=U_inf, circulation=circulation, aoa_deg=aoa_deg
    )
    ux_expected = np.real(vel_out)
    uy_expected = -np.imag(vel_out)
    speed_expected = np.hypot(ux_expected, uy_expected)

    assert np.isclose(flow["ux"][0, 1], ux_expected)
    assert np.isclose(flow["uy"][0, 1], uy_expected)
    assert np.isclose(flow["u"][0, 1], speed_expected)
    assert np.isclose(flow["Cp"][0, 1], 1.0 - (speed_expected / U_inf) ** 2)


def test_cylinder_complex_potential_returns_expected_value() -> None:
    z = np.array([2.0 + 0.5j])
    R = 1.0
    U_inf = 1.5
    circulation = 0.9
    aoa_deg = 7.5
    aoa = np.deg2rad(aoa_deg)

    W = cylinder_complex_potential(
        z, R=R, U_inf=U_inf, circulation=circulation, aoa_deg=aoa_deg
    )
    expected = U_inf * (
        z * np.exp(-1j * aoa) + (R**2) / z * np.exp(1j * aoa)
    ) + 1j * (circulation / (2.0 * np.pi)) * np.log(z / R)

    assert np.allclose(W, expected)


def test_surface_velocity_matches_boundary_expression() -> None:
    theta = np.array([0.0, np.pi / 2.0])
    R = 1.0
    U_inf = 3.0
    circulation = 2.0
    aoa_deg = 8.5
    aoa = np.deg2rad(aoa_deg)

    u_r, u_theta = cylinder_surface_velocity(
        theta, R=R, U_inf=U_inf, circulation=circulation, aoa_deg=aoa_deg
    )
    u_theta_expected = -2.0 * U_inf * np.sin(theta - aoa) - circulation / (
        2.0 * np.pi * R
    )

    assert np.allclose(u_r, 0.0)
    assert np.allclose(u_theta, u_theta_expected)


def test_pressure_coefficient_surface_matches_zero_circulation_relation() -> None:
    theta = np.array([0.0, np.pi / 2.0, np.pi])
    cp = cylinder_pressure_coefficient_surface(
        theta, R=1.0, U_inf=2.0, circulation=0.0, aoa_deg=0.0
    )
    cp_expected = 1.0 - 4.0 * np.sin(theta) ** 2
    assert np.allclose(cp, cp_expected)


def test_pressure_coefficient_surface_kutta_matches_surface_with_kutta_gamma() -> None:
    theta = np.array([0.0, np.pi / 3.0, np.pi])
    R = 1.0
    U_inf = 2.5
    aoa_deg = 9.0

    cp_kutta = cylinder_pressure_coefficient_surface(
        theta, R=R, U_inf=U_inf, circulation=None, aoa_deg=aoa_deg
    )

    gamma = cylinder_circulation_for_kutta_condition(R=R, U_inf=U_inf, aoa_deg=aoa_deg)
    cp_expected = cylinder_pressure_coefficient_surface(
        theta, R=R, U_inf=U_inf, circulation=gamma, aoa_deg=aoa_deg
    )

    assert np.allclose(cp_kutta, cp_expected)


def test_cylinder_lift_coefficient_defaults_to_kutta() -> None:
    R = 1.1
    U_inf = 2.3
    aoa_deg = 6.0
    cl = cylinder_lift_coefficient(R=R, U_inf=U_inf, circulation=None, aoa_deg=aoa_deg)
    gamma = cylinder_circulation_for_kutta_condition(R=R, U_inf=U_inf, aoa_deg=aoa_deg)
    assert np.isclose(cl, 2.0 * gamma / (2.0 * R * U_inf))


def test_gamma_for_joukowsky_kutta_scales_with_sin_aoa() -> None:
    R = 1.2
    U_inf = 3.4
    aoa_deg = 15.0
    gamma = cylinder_circulation_for_kutta_condition(R=R, U_inf=U_inf, aoa_deg=aoa_deg)
    expected = 4.0 * np.pi * R * U_inf * np.sin(np.deg2rad(aoa_deg))
    assert np.isclose(gamma, expected)
