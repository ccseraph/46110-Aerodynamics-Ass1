import matplotlib.pyplot as plt
import numpy as np

from funaerotool.cylinder_potential_flow import (
    cylinder_flow_field,
    cylinder_lift_coefficient,
    cylinder_pressure_coefficient_surface,
)
from funaerotool.plotting import plot_cp_distribution, plot_flow_field
from funaerotool.utils import generate_circle_contour

# Parse command-line options

# Define physical setup parameters
R = 1.0
U_inf = 1.0
aoa_deg = 10.0
circulation = None  # None = Kutta condition

# If you want to lock a specific circulation value, set it explicitly:
# circulation = cylinder_circulation_for_kutta_condition(R=R, U_inf=U_inf, aoa_deg=aoa_deg)

# Compute nondimensional lift coefficient from the chosen circulation
cl = cylinder_lift_coefficient(
    R=R, U_inf=U_inf, aoa_deg=aoa_deg, circulation=circulation
)

# Build Cartesian grid for field evaluation
x_vec = np.linspace(-5.0, 5.0, 200)
y_vec = np.linspace(-5.0, 5.0, 198)
x_grid, y_grid = np.meshgrid(x_vec, y_vec)

# Evaluate flow field on the grid (includes Cp)
flow = cylinder_flow_field(
    x=x_grid,
    y=y_grid,
    R=R,
    U_inf=U_inf,
    circulation=circulation,
    aoa_deg=aoa_deg,
    mask_inside=True,
)

# Build contour and evaluate surface Cp along the cylinder
x_contour, y_contour = generate_circle_contour(n_points=361, radius=R)
theta = np.arctan2(y_contour, x_contour)
c_p = cylinder_pressure_coefficient_surface(
    theta, R=R, U_inf=U_inf, circulation=circulation, aoa_deg=aoa_deg
)

# Build 1x2 figure similar to panel-method example
fig, (ax_cp, ax_flow) = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)

# Cp subplot (-C_p vs x/R)
plot_cp_distribution(xp=x_contour / R, cp=c_p, ax=ax_cp)
ax_cp.set_xlabel(r"$x/R$")
ax_cp.set_title("Pressure coefficient vs x/R")
ax_cp.set_xlim(-1.1, 1.1)

# Flow-field subplot: streamlines + Cp heatmap (computed from |u|)
Cp_field = flow["Cp"]
Cp_field = np.ma.filled(np.ma.masked_invalid(Cp_field), np.nan)
plot_flow_field(
    X=x_grid,
    Y=y_grid,
    ux=flow["ux"],
    uy=flow["uy"],
    Cp=flow["Cp"],
    x_contour=x_contour,
    y_contour=y_contour,
    ax=ax_flow,
    stream_density=1.2,
)
ax_flow.set_xlim(-5, 5)
ax_flow.set_ylim(-3, 3)
ax_flow.set_title(f"Flow field around cylinder (Cl={cl:.3f})")

# Print key scalar outputs to terminal
print(f"Cl = {cl:.6f}")
print(f"aoa = {aoa_deg:.3f} deg")

# Save figure or show interactively
plt.show()
