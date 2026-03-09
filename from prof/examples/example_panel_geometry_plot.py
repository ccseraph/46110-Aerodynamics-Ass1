import matplotlib.pyplot as plt

from funaerotool.panel_method import panel_geometry
from funaerotool.plotting import plot_panel_geometry
from funaerotool.utils import generate_naca4_contour


# Generate NACA 4412 airfoil contour
x, y = generate_naca4_contour("4412", n_points=31, closed_te=False)

# Compute panel geometry
plength, xp, yp, Tx, Ty, Nx, Ny = panel_geometry(x, y)

fig, ax = plot_panel_geometry(
    x=x,
    y=y,
    xp=xp,
    yp=yp,
    Tx=Tx,
    Ty=Ty,
    Nx=Nx,
    Ny=Ny,
    scale=0.05,
)

plt.show()
