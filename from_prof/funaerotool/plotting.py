
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np


def plot_cp_distribution(
    xp: np.ndarray,
    cp: np.ndarray,
    ax: Optional[plt.Axes] = None,
    **kwargs_plot,
) -> Tuple[plt.Figure, plt.Axes]:
    """Plot pressure coefficient distribution along the contour control points."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 4))
    else:
        fig = ax.figure

    kwargs_plot.setdefault("color", "k")
    kwargs_plot.setdefault("linewidth", 2.0)
    kwargs_plot.setdefault("ls", "-")

    ax.plot(xp, cp, **kwargs_plot)
    ax.set_xlabel("x (control points)")
    ax.set_ylabel("$C_p$")
    ax.set_title("Pressure coefficient (panel method)")
    ax.grid(True)
    ax.invert_yaxis()
    if "label" in kwargs_plot:
        ax.legend()
    return fig, ax


def plot_flow_field(
    X: np.ndarray,
    Y: np.ndarray,
    ux: np.ndarray,
    uy: np.ndarray,
    Cp: np.ndarray,
    x_contour: np.ndarray | None = None,
    y_contour: np.ndarray | None = None,
    ax: Optional[plt.Axes] = None,
    stream_density: float = 0.8,
    cmap: str = "coolwarm",
    linewidth_min: float = 0.4,
    linewidth_max: float = 1.6,
) -> Tuple[plt.Figure, plt.Axes]:
    """Plot flow field with streamlines and C_p colormap; contour overlay optional."""
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))
    else:
        fig = ax.figure

    speed = np.hypot(ux, uy)
    speed_masked = np.ma.masked_invalid(speed)
    speed_norm = speed_masked / speed_masked.max()
    stream_linewidth = linewidth_min + (
        linewidth_max - linewidth_min
    ) * speed_norm.filled(0.0)

    cp_plot = ax.pcolormesh(X, Y, Cp, shading="auto", cmap=cmap)
    ax.streamplot(
        X,
        Y,
        ux,
        uy,
        color="k",
        linewidth=stream_linewidth,
        density=stream_density,
        arrowsize=1.1,
    )
    if x_contour is not None and y_contour is not None:
        ax.fill(x_contour, y_contour, color="k", alpha=1.0, zorder=4)
        ax.plot(x_contour, y_contour, color="k", linewidth=1.0, zorder=5)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Flow field: streamlines + $C_p$ heatmap")
    fig.colorbar(cp_plot, ax=ax, orientation="vertical", shrink=0.9, label="$C_p$")
    return fig, ax


def plot_panel_geometry(
    x: np.ndarray,
    y: np.ndarray,
    xp: np.ndarray,
    yp: np.ndarray,
    Tx: np.ndarray,
    Ty: np.ndarray,
    Nx: np.ndarray,
    Ny: np.ndarray,
    ax: Optional[plt.Axes] = None,
    scale: float = 0.2,
    color_contour: str = "k",
    color_tangent: str = "tab:blue",
    color_normal: str = "tab:red",
    color_edges: str = "0.6",
) -> Tuple[plt.Figure, plt.Axes]:
    """Plot panel contour with control points and tangent/normal vectors."""

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 5))
    else:
        fig = ax.figure

    ax.plot(x, y, color=color_contour, linewidth=1.2, label="Contour")
    ax.scatter(xp, yp, s=12, color=color_contour, zorder=3, label="Control points")

    # Panel edges: short orthogonal line segments centered on each control poin
    half_span = 1e-2

    x0 = np.zeros_like(x)
    y0 = np.zeros_like(y)
    x1 = np.zeros_like(x)
    y1 = np.zeros_like(y)
    x0[:-1] = x[:-1] - half_span * Nx
    y0[:-1] = y[:-1] - half_span * Ny
    x1[:-1] = x[:-1] + half_span * Nx
    y1[:-1] = y[:-1] + half_span * Ny
    x0[-1] = x[-1] - half_span * Nx[-1]
    y0[-1] = y[-1] - half_span * Ny[-1]
    x1[-1] = x[-1] + half_span * Nx[-1]
    y1[-1] = y[-1] + half_span * Ny[-1]

    xs = np.empty(x.size * 3)
    ys = np.empty(x.size * 3)
    xs[0::3] = x0
    xs[1::3] = x1
    xs[2::3] = np.nan
    ys[0::3] = y0
    ys[1::3] = y1
    ys[2::3] = np.nan

    ax.plot(
        xs,
        ys,
        color=color_edges,
        linewidth=2.0,
        alpha=0.9,
        zorder=2,
        label="Panel edges",
    )

    ax.quiver(
        xp,
        yp,
        Tx,
        Ty,
        angles="xy",
        scale_units="xy",
        scale=1.0 / scale,
        color=color_tangent,
        width=0.005,
        label="Tangent",
        zorder=4,
    )
    ax.quiver(
        xp,
        yp,
        Nx,
        Ny,
        angles="xy",
        scale_units="xy",
        scale=1.0 / scale,
        color=color_normal,
        width=0.005,
        label="Normal",
        zorder=3,
    )

    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_ylim(-0.5, 0.5)
    ax.set_title("Panel geometry: contour, control points, tangents, normals")
    ax.legend()
    ax.grid(True, linestyle=":", linewidth=0.8)
    return fig, ax
