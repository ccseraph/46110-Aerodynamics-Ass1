from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from funaerotool.panel_method import (
    compute_panel_flow_field,
    solve_closed_contour_panel_method,
)
from funaerotool.plotting import plot_cp_distribution, plot_flow_field
from funaerotool.utils import generate_naca4_contour


aoa_deg = 10.0
U_inf = 1.0
n_points = 161

# Generate NACA 4412 airfoil contour
x, y = generate_naca4_contour("4412", n_points=n_points, closed_te=False)

# Solve panel method
sol = solve_closed_contour_panel_method(
    x, y, aoa_deg=aoa_deg, U_inf=U_inf, kutta_condition=True
)
fig, (ax_cp, ax_flow) = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)

# Cp subplot
plot_cp_distribution(xp=sol["xp"], cp=sol["Cp"], ax=ax_cp)

# Flow-field subplot: streamlines + Cp heatmap
grid_x = np.linspace(-1.0, 2.0, 200)
grid_y = np.linspace(-1.0, 1.0, 198)
X, Y = np.meshgrid(grid_x, grid_y)

flow = compute_panel_flow_field(
    x_field=X,
    y_field=Y,
    x_contour=x,
    y_contour=y,
    sigma=sol["sigma"],
    gamma=sol.get("gamma"),
    U_inf=U_inf,
    aoa_deg=aoa_deg,
)
plot_flow_field(
    X=X,
    Y=Y,
    ux=flow["ux"],
    uy=flow["uy"],
    Cp=flow["Cp"],
    x_contour=x,
    y_contour=y,
    ax=ax_flow,
    stream_density=0.8,
)

figures_dir = Path("figures")

