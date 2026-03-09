import matplotlib.pyplot as plt
import numpy as np
from solver import solve_closed_contour_panel_method

from funaerotool.panel_method import compute_panel_flow_field, compute_point_flow_field
from funaerotool.plotting import plot_cp_distribution, plot_flow_field
from funaerotool.utils import generate_circle_contour, generate_naca4_contour

aoa_deg = 10.0
U_inf = 1.0
n_points = 161
kutta_condition = False  # Set to False to see effect of not enforcing Kutta condition
use_cylinder = True  # Set to True to run the circle Cp comparison instead of the airfoil flow field example

if use_cylinder:
    x, y = generate_circle_contour(n_points=n_points, radius=0.5)
    x += 0.5  # Shift circle to the right so the point vortex is not inside the contour
else:
    # Generate NACA 4412 airfoil contour
    x, y = generate_naca4_contour("0018", n_points=n_points, closed_te=False)

# Solve panel method
sol = solve_closed_contour_panel_method(
    x, y, aoa_deg=aoa_deg, U_inf=U_inf, kutta_condition=kutta_condition
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
    gamma=sol["gamma"],
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

# Flow-field visualization: streamlines and Cp heatmap
plt.show()
