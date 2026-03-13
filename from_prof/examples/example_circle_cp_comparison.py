import matplotlib.pyplot as plt
import numpy as np

from funaerotool import (
    cylinder_pressure_coefficient_surface,
    solve_closed_contour_panel_method,
)
from funaerotool.utils import generate_circle_contour


# Configuration
radius = 1.0
U_inf = 1.0
aoa_deg = 10.0
panel_counts = [40, 80, 160, 320]

fig, (ax_cp, ax_err) = plt.subplots(2, 1, figsize=(9, 8), constrained_layout=True)
max_errors = []

for n_panels in panel_counts:
    n_points = n_panels + 1  # closed contour
    x, y = generate_circle_contour(n_points=n_points, radius=radius)

    sol = solve_closed_contour_panel_method(
        x, y, aoa_deg=aoa_deg, U_inf=U_inf, kutta_condition=True
    )

    theta = np.arctan2(sol["yp"], sol["xp"])
    sort_idx = np.argsort(theta)
    theta_sorted = theta[sort_idx]

    cp_num = sol["Cp"][sort_idx]
    cp_exact_sorted = cylinder_pressure_coefficient_surface(
        theta_sorted, R=radius, U_inf=U_inf, circulation=None, aoa_deg=aoa_deg
    )
    cp_error_sorted = cp_num - cp_exact_sorted

    # Map back to original order for plotting raw theta if desired
    inv_idx = np.argsort(sort_idx)
    cp_exact = cp_exact_sorted[inv_idx]
    cp_error = cp_error_sorted[inv_idx]
    max_errors.append(np.max(np.abs(cp_error)))

    ax_cp.plot(theta_sorted, cp_num, label=f"Panel (n={n_panels})")
    ax_err.plot(theta_sorted, cp_error_sorted, label=f"Error (n={n_panels})")

ax_cp.plot(theta_sorted, cp_exact_sorted, "--k", linewidth=1.5, label="Analytical")
ax_cp.set_xlabel("θ [rad]")
ax_cp.set_ylabel("C_p")
ax_cp.set_title("Pressure coefficient on a circle: panel vs analytical")
ax_cp.legend()
ax_cp.grid(True)

ax_err.set_xlabel("θ [rad]")
ax_err.set_ylabel("C_p error (num - exact)")
ax_err.set_title("C_p error for varying panel counts")
ax_err.legend()
ax_err.grid(True)

fig_err, ax_log = plt.subplots(figsize=(6, 4), constrained_layout=True)
ax_log.loglog(panel_counts, max_errors, "-o", label="Max |Cp error|")

log_panels = np.log10(panel_counts)
log_errors = np.log10(max_errors)
slope, intercept = np.polyfit(log_panels, log_errors, 1)
fit_errors = 10 ** (intercept + slope * log_panels)

ax_log.loglog(panel_counts, fit_errors, "--", label=f"Fit slope ≈ {slope:.2f}")
ax_log.set_xlabel("Number of panels")
ax_log.set_ylabel("Max |C_p error|")
ax_log.set_title("Error convergence vs panel count")
ax_log.grid(True, which="both")
ax_log.legend()

plt.show()
