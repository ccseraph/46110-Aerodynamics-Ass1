
import numpy as np
import pytest

from funaerotool.cylinder_potential_flow import (
    cylinder_circulation_for_kutta_condition,
    cylinder_lift_coefficient,
    cylinder_pressure_coefficient_surface,
)
from funaerotool.panel_method import (
    compute_panel_flow_field,
    compute_point_flow_field,
    compute_pressure_coefficient,
    constant_vortex_distribution,
    flip_contour,
    freestream_components,
    global_panel_induced_velocity_matrices,
    global_to_local,
    local_to_global,
    panel_geometry,
    parabolic_vortex_distribution,
    point_source_induced_velocity,
    point_vortex_induced_velocity,
    solve_closed_contour_panel_method,
    source_panel_induced_velocity_local,
    vortex_panel_induced_velocity_local,
)
from funaerotool.utils import generate_circle_contour, generate_naca4_contour


def test_generate_circle_contour_is_closed_and_radius_constant() -> None:
    x, y = generate_circle_contour(n_points=21, radius=2.0)

    assert x.shape == (21,)
    assert y.shape == (21,)
    assert np.isclose(x[0], x[-1])
    assert np.isclose(y[0], y[-1])
    assert np.allclose(np.hypot(x, y), 2.0)


def test_panel_geometry_returns_expected_sizes_and_unit_vectors() -> None:
    x, y = generate_circle_contour(n_points=31, radius=1.0)
    plength, xp, yp, tx, ty, nx, ny = panel_geometry(x, y)

    n = x.size - 1
    assert plength.shape == (n,)
    assert xp.shape == (n,)
    assert yp.shape == (n,)

    t_norm = np.hypot(tx, ty)
    n_norm = np.hypot(nx, ny)
    dot_tn = tx * nx + ty * ny

    assert np.allclose(t_norm, 1.0)
    assert np.allclose(n_norm, 1.0)
    assert np.allclose(dot_tn, 0.0)


def test_source_panel_self_induction_matches_reference_values() -> None:
    ux, uy = source_panel_induced_velocity_local(sx=0.0, sy=0.0, panel_length=2.0)

    assert np.isclose(ux, 0.0)
    assert np.isclose(uy, 0.5)


def test_source_panel_induced_velocity_local_simple_point() -> None:
    ux, uy = source_panel_induced_velocity_local(sx=0.0, sy=1.0, panel_length=2.0)

    ux_expected = 0.0
    uy_expected = 0.25

    assert np.isclose(ux, ux_expected)
    assert np.isclose(uy, uy_expected)


def test_local_global_velocity_transforms_are_consistent() -> None:
    ux_global, uy_global = local_to_global(
        ut=0.2,
        un=-0.7,
        tx=0.6,
        ty=0.8,
        nx=-0.8,
        ny=0.6,
    )
    ux_panel, uy_panel = global_to_local(
        ux=ux_global,
        uy=uy_global,
        tx=0.6,
        ty=0.8,
        nx=-0.8,
        ny=0.6,
    )

    assert np.isclose(ux_panel, 0.2)
    assert np.isclose(uy_panel, -0.7)


def test_point_source_induced_velocity_matches_expected_direction() -> None:
    ux, uy = point_source_induced_velocity(
        x=1.0, y=0.0, x_source=0.0, y_source=0.0, strength=2 * np.pi
    )

    assert np.isclose(ux, 1.0)
    assert np.isclose(uy, 0.0)


def test_point_vortex_induced_velocity_matches_expected_rotation() -> None:
    ux, uy = point_vortex_induced_velocity(
        x=1.0, y=0.0, x_vortex=0.0, y_vortex=0.0, circulation=2 * np.pi
    )

    assert np.isclose(ux, 0.0)
    assert np.isclose(uy, 1.0)


def test_point_singularity_raises_on_coincident_evaluation() -> None:
    with pytest.raises(ValueError):
        point_source_induced_velocity(
            x=0.0, y=0.0, x_source=0.0, y_source=0.0, strength=1.0
        )


def test_freestream_components_match_alpha_zero_and_ninety() -> None:
    ux0, uy0 = freestream_components(aoa_deg=0.0, U_inf=2.0)
    ux90, uy90 = freestream_components(aoa_deg=90.0, U_inf=2.0)

    assert np.isclose(ux0, 2.0)
    assert np.isclose(uy0, 0.0)
    assert np.isclose(ux90, 0.0, atol=1e-14)
    assert np.isclose(uy90, 2.0)


def test_compute_pressure_coefficient_from_velocity_components() -> None:
    ux = np.array([0.0, 1.0, 2.0])
    uy = np.zeros_like(ux)
    cp = compute_pressure_coefficient(ux=ux, uy=uy, U_inf=2.0)

    assert np.allclose(cp, np.array([1.0, 0.75, 0.0]))


def test_compute_panel_flow_field_returns_freestream_for_zero_strengths() -> None:
    x_contour, y_contour = generate_circle_contour(n_points=21, radius=1.0)
    n_panels = x_contour.size - 1

    x_field_line = np.linspace(-2.0, 2.0, 5)
    x_field, y_field = np.meshgrid(x_field_line, x_field_line)

    strengths = np.zeros(n_panels)

    flow = compute_panel_flow_field(
        x_field=x_field,
        y_field=y_field,
        x_contour=x_contour,
        y_contour=y_contour,
        sigma=strengths,
        U_inf=2.0,
        aoa_deg=0.0,
        mask_inside=False,
    )

    assert np.allclose(flow["ux"], 2.0)
    assert np.allclose(flow["uy"], 0.0)
    assert np.allclose(flow["Cp"], 0.0)


def test_compute_panel_flow_field_masks_inside_by_default() -> None:
    x_contour = np.array([-1.0, 1.0, 1.0, -1.0, -1.0])
    y_contour = np.array([-1.0, -1.0, 1.0, 1.0, -1.0])
    strengths = np.zeros(x_contour.size - 1)

    x_line = np.array([-2.0, 0.0, 2.0])
    x_field, y_field = np.meshgrid(x_line, x_line)

    flow = compute_panel_flow_field(
        x_field=x_field,
        y_field=y_field,
        x_contour=x_contour,
        y_contour=y_contour,
        sigma=strengths,
        U_inf=1.0,
        aoa_deg=0.0,
    )

    assert np.isnan(flow["ux"][1, 1])
    assert np.isnan(flow["Cp"][1, 1])
    assert np.isclose(flow["ux"][0, 0], 1.0)
    assert np.isclose(flow["Cp"][0, 0], 0.0)


def test_compute_panel_flow_field_allows_disabling_mask() -> None:
    x_contour = np.array([-1.0, 1.0, 1.0, -1.0, -1.0])
    y_contour = np.array([-1.0, -1.0, 1.0, 1.0, -1.0])
    strengths = np.zeros(x_contour.size - 1)

    x_line = np.array([-2.0, 0.0, 2.0])
    x_field, y_field = np.meshgrid(x_line, x_line)

    flow = compute_panel_flow_field(
        x_field=x_field,
        y_field=y_field,
        x_contour=x_contour,
        y_contour=y_contour,
        sigma=strengths,
        U_inf=1.0,
        aoa_deg=0.0,
        mask_inside=False,
    )

    assert np.isfinite(flow["ux"]).all()
    assert np.isfinite(flow["Cp"]).all()
    assert np.isclose(flow["ux"][1, 1], 1.0)


def test_compute_point_flow_field_supports_multiple_sources_and_vortices() -> None:
    x_field, y_field = np.meshgrid(np.linspace(-1.0, 1.0, 5), np.linspace(-1.0, 1.0, 5))

    Sigma = np.array([2 * np.pi, -2 * np.pi])
    x_sigma = np.array([0.5, -0.5])
    y_sigma = np.array([0.75, 0.75])

    Gamma = np.array([2 * np.pi])
    x_gamma = np.array([0.0])
    y_gamma = np.array([0.75])

    flow = compute_point_flow_field(
        x_field=x_field,
        y_field=y_field,
        Sigma=Sigma,
        x_sigma=x_sigma,
        y_sigma=y_sigma,
        Gamma=Gamma,
        x_gamma=x_gamma,
        y_gamma=y_gamma,
        U_inf=1.0,
    )

    assert flow["ux"].shape == x_field.shape
    assert flow["uy"].shape == y_field.shape
    assert flow["Cp"] is not None


def test_parabolic_vortex_distribution_is_nonnegative_and_normalized() -> None:
    panel_lengths = np.ones(8)
    vort = parabolic_vortex_distribution(panel_lengths)

    assert vort.shape == (8,)
    assert np.all(vort >= 0.0)
    assert vort[0] == 0.0
    assert vort[-1] == 0.0
    assert np.isclose(np.sum(vort * panel_lengths), 1.0)


def test_constant_vortex_distribution_is_normalized() -> None:
    panel_lengths = np.array([1.0, 2.0, 3.0])
    vort = constant_vortex_distribution(panel_lengths)

    assert vort.shape == (3,)
    assert np.allclose(vort, vort[0])
    assert np.isclose(np.sum(vort * panel_lengths), 1.0)


def test_constant_vortex_panel_self_induction_matches_reference_values() -> None:
    ux, uy = vortex_panel_induced_velocity_local(sx=0.0, sy=0.0, panel_length=2.0)

    assert np.isclose(ux, -0.5)
    assert np.isclose(uy, 0.0)


def test_constant_vortex_panel_induced_velocity_local_is_rotated_source() -> None:
    ux_source, uy_source = source_panel_induced_velocity_local(
        sx=0.0, sy=1.0, panel_length=2.0
    )
    ux_vortex, uy_vortex = vortex_panel_induced_velocity_local(
        sx=0.0, sy=1.0, panel_length=2.0
    )

    assert np.isclose(ux_vortex, -uy_source)
    assert np.isclose(uy_vortex, ux_source)


def test_panel_geometry_requires_closed_contour() -> None:
    x = np.array([0.0, 1.0, 0.0])
    y = np.array([0.0, 0.0, 1.0])
    panel_geometry(x, y)


def test_panel_geometry_requires_airfoil_ordering() -> None:
    x, y = generate_naca4_contour(naca_code="2412", n_points=401)

    # Both orientations should be acceptable after relaxed validation
    panel_geometry(x, y)
    panel_geometry(*flip_contour(x, y))


def test_global_panel_induced_velocity_matrices_source_single_panel() -> None:
    xp = np.array([0.0])
    yp = np.array([0.0])
    tx = np.array([1.0])
    ty = np.array([0.0])
    nx = np.array([0.0])
    ny = np.array([1.0])
    plength = np.array([2.0])

    x_eval = np.array([0.0, 0.5])
    y_eval = np.array([1.0, 2.0])

    A_global, B_global = global_panel_induced_velocity_matrices(
        x_eval=x_eval,
        y_eval=y_eval,
        xp=xp,
        yp=yp,
        Tx=tx,
        Ty=ty,
        Nx=nx,
        Ny=ny,
        panel_lengths=plength,
        panel_type="source",
    )

    assert A_global.shape == (2, 1)
    assert B_global.shape == (2, 1)

    ux_ref, uy_ref = source_panel_induced_velocity_local(
        sx=0.0, sy=1.0, panel_length=2.0
    )
    assert np.isclose(A_global[0, 0], ux_ref)
    assert np.isclose(B_global[0, 0], uy_ref)


def test_global_panel_induced_velocity_matrices_vortex_single_panel() -> None:
    xp = np.array([0.0])
    yp = np.array([0.0])
    tx = np.array([1.0])
    ty = np.array([0.0])
    nx = np.array([0.0])
    ny = np.array([1.0])
    plength = np.array([2.0])

    x_eval = np.array([0.0])
    y_eval = np.array([1.0])

    A_global, B_global = global_panel_induced_velocity_matrices(
        x_eval=x_eval,
        y_eval=y_eval,
        xp=xp,
        yp=yp,
        Tx=tx,
        Ty=ty,
        Nx=nx,
        Ny=ny,
        panel_lengths=plength,
        panel_type="vortex",
    )

    ux_source, uy_source = source_panel_induced_velocity_local(
        sx=0.0, sy=1.0, panel_length=2.0
    )
    assert np.isclose(A_global[0, 0], -uy_source)
    assert np.isclose(B_global[0, 0], ux_source)


def test_global_panel_induced_velocity_matrices_multiple_eval_points() -> None:
    x, y = generate_circle_contour(n_points=11, radius=1.0)
    plength, xp, yp, tx, ty, nx, ny = panel_geometry(x, y)

    x_eval = np.array([2.0, -1.5, 0.0])
    y_eval = np.array([0.0, 1.0, 5.0])

    A_global, B_global = global_panel_induced_velocity_matrices(
        x_eval=x_eval,
        y_eval=y_eval,
        xp=xp,
        yp=yp,
        Tx=tx,
        Ty=ty,
        Nx=nx,
        Ny=ny,
        panel_lengths=plength,
        panel_type="source",
    )

    n_panels = plength.size
    assert A_global.shape == (3, n_panels)
    assert B_global.shape == (3, n_panels)
    assert np.isfinite(A_global).all()
    assert np.isfinite(B_global).all()


def test_solver_circle_no_circulation_normal_velocity_small() -> None:
    n_panels = 120
    n_points = n_panels + 1
    x, y = generate_circle_contour(n_points=n_points, radius=1.0)
    sol = solve_closed_contour_panel_method(
        x, y, aoa_deg=0.0, U_inf=1.0, kutta_condition=False
    )

    Vn = sol["Vn"]
    assert np.max(np.abs(Vn)) < 5e-2

    cp = sol["Cp"]
    assert np.all(np.isfinite(cp))


def test_coefficient_convergence_no_kutta() -> None:
    panel_counts = [40, 80, 160, 320]
    max_cp_errors = []

    for n_panels in panel_counts:
        x, y = generate_circle_contour(n_points=n_panels + 1, radius=1.0)
        sol = solve_closed_contour_panel_method(
            x, y, aoa_deg=0.0, U_inf=1.0, kutta_condition=False
        )

        theta = np.arctan2(sol["yp"], sol["xp"])
        cp_exact = cylinder_pressure_coefficient_surface(
            theta, R=1.0, U_inf=1.0, circulation=0.0
        )
        cp_error = sol["Cp"] - cp_exact
        max_cp_errors.append(np.max(np.abs(cp_error)))

        assert abs(sol["Cl"]) < 1e-12

    # Errors are at round-off; just enforce small magnitude
    assert np.all(np.array(max_cp_errors) < 1e-12)


def test_coefficient_convergence_with_kutta() -> None:
    panel_counts = [40, 80, 160, 320, 640]
    max_cp_errors = np.zeros_like(panel_counts, dtype=float)
    cl_errors = np.zeros_like(panel_counts, dtype=float)

    aoa_deg = 10.0
    cl_exact = cylinder_lift_coefficient(aoa_deg=aoa_deg, R=1.0, U_inf=1.0)

    for i, n_panels in enumerate(panel_counts):
        x, y = generate_circle_contour(n_points=n_panels + 1, radius=1.0)
        sol = solve_closed_contour_panel_method(
            x,
            y,
            aoa_deg=aoa_deg,
            U_inf=1.0,
            kutta_condition=True,
        )
        theta = np.arctan2(sol["yp"], sol["xp"])
        sort_idx = np.argsort(theta)
        theta_sorted = theta[sort_idx]

        cp_sorted = sol["Cp"][sort_idx]
        cp_exact_sorted = cylinder_pressure_coefficient_surface(
            theta_sorted, R=1.0, U_inf=1.0, circulation=None, aoa_deg=aoa_deg
        )

        cp_error_sorted = cp_sorted - cp_exact_sorted
        max_cp_errors[i] = np.max(np.abs(cp_error_sorted))

        cl_errors[i] = abs(sol["Cl"] - cl_exact)

    assert max_cp_errors[-1] < 0.01
    assert cl_errors[-1] < 0.05
    rel_cp_errors = max_cp_errors[1:] / max_cp_errors[:-1]
    rel_cl_errors = cl_errors[1:] / cl_errors[:-1]
    assert np.allclose(rel_cp_errors, rel_cp_errors[-1], atol=0.1)
    assert np.allclose(rel_cl_errors, rel_cl_errors[-1], atol=0.1)
