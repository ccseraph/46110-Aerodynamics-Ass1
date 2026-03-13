import matplotlib.pyplot as plt
import numpy as np

from funaerotool.panel_method import compute_panel_flow_field, compute_point_flow_field
from funaerotool.plotting import plot_flow_field
from funaerotool.utils import generate_circle_contour


# Grid for all cases
x_line = np.linspace(-2.0, 2.0, 220)
y_line = np.linspace(-2.0, 2.0, 220)
X, Y = np.meshgrid(x_line, y_line)

# Freestream only (aoa=0, aoa=15)
freestream_0 = compute_point_flow_field(
    x_field=X,
    y_field=Y,
    Sigma=None,
    Gamma=None,
    U_inf=1.0,
    aoa_deg=0.0,
)
freestream_15 = compute_point_flow_field(
    x_field=X,
    y_field=Y,
    Sigma=None,
    Gamma=None,
    U_inf=1.0,
    aoa_deg=15.0,
)

# Point singularities: source and vortex at origin
point_source = compute_point_flow_field(
    x_field=X,
    y_field=Y,
    Sigma=1,
    x_sigma=0.0,
    y_sigma=0.0,
    Gamma=None,
    U_inf=1.0,
    aoa_deg=None,
)
point_vortex = compute_point_flow_field(
    x_field=X,
    y_field=Y,
    Sigma=None,
    Gamma=1,
    x_gamma=0.0,
    y_gamma=0.0,
    U_inf=1.0,
    aoa_deg=None,
)

# Panel singularities on a circular contour
panel_x = np.array([-0.5, 0.0, 0.5])
panel_y = np.array([0.0, 0.0, 0.0])
n_panels = panel_x.size - 1
sigma_uniform = np.ones(n_panels)
gamma_uniform = np.ones(n_panels)

panel_source = compute_panel_flow_field(
    x_field=X,
    y_field=Y,
    x_contour=panel_x,
    y_contour=panel_y,
    sigma=sigma_uniform,
    gamma=None,
    U_inf=1.0,
    aoa_deg=None,
)
panel_vortex = compute_panel_flow_field(
    x_field=X,
    y_field=Y,
    x_contour=panel_x,
    y_contour=panel_y,
    sigma=None,
    gamma=gamma_uniform,
    U_inf=1.0,
    aoa_deg=None,
)

fig, axs = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)

# Column 1: freestream
plot_flow_field(
    X=X,
    Y=Y,
    ux=freestream_0["ux"],
    uy=freestream_0["uy"],
    Cp=freestream_0["Cp"],
    ax=axs[0, 0],
    stream_density=1.0,
)
axs[0, 0].set_title("Freestream: aoa=0°")

plot_flow_field(
    X=X,
    Y=Y,
    ux=freestream_15["ux"],
    uy=freestream_15["uy"],
    Cp=freestream_15["Cp"],
    ax=axs[1, 0],
    stream_density=1.0,
)
axs[1, 0].set_title("Freestream: aoa=15°")

# Column 2: point singularities
plot_flow_field(
    X=X,
    Y=Y,
    ux=point_source["ux"],
    uy=point_source["uy"],
    Cp=point_source["Cp"],
    ax=axs[0, 1],
)
axs[0, 1].set_title("Point source (Γ=0)")

plot_flow_field(
    X=X,
    Y=Y,
    ux=point_vortex["ux"],
    uy=point_vortex["uy"],
    Cp=point_vortex["Cp"],
    ax=axs[1, 1],
)
axs[1, 1].set_title("Point vortex (Γ=2π)")

# Column 3: panel singularities on a circle
plot_flow_field(
    X=X,
    Y=Y,
    ux=panel_source["ux"],
    uy=panel_source["uy"],
    Cp=panel_source["Cp"],
    x_contour=panel_x,
    y_contour=panel_y,
    ax=axs[0, 2],
)
axs[0, 2].set_title("Panel source (uniform σ)")

plot_flow_field(
    X=X,
    Y=Y,
    ux=panel_vortex["ux"],
    uy=panel_vortex["uy"],
    Cp=panel_vortex["Cp"],
    x_contour=panel_x,
    y_contour=panel_y,
    ax=axs[1, 2],
)
axs[1, 2].set_title("Panel vortex (uniform γ)")

for ax in axs.ravel():
    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(-2.0, 2.0)

plt.show()
