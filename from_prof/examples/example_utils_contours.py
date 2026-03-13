import matplotlib.pyplot as plt

from funaerotool.utils import generate_circle_contour, generate_naca4_contour


# Example 1: circle contour
x_circle, y_circle = generate_circle_contour(n_points=201, radius=1.0)

# Example 2: NACA 4-series contour
x_naca, y_naca = generate_naca4_contour(naca_code="2412", n_points=401)

fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
ax1, ax2 = axes

ax1.plot(x_circle, y_circle, "-k", linewidth=2)
ax1.set_aspect("equal", adjustable="box")
ax1.grid(True)
ax1.set_title("Circle contour")
ax1.set_xlabel("x")
ax1.set_ylabel("y")

ax2.plot(x_naca, y_naca, "-k", linewidth=2)
ax2.scatter([0.0, 1.0], [0.0, 0.0], color=["tab:red", "tab:blue"], s=35, zorder=5)
ax2.annotate(
    "LE", (0.0, 0.0), xytext=(5, 6), textcoords="offset points", color="tab:red"
)
ax2.annotate(
    "TE mid",
    (1.0, 0.0),
    xytext=(5, 6),
    textcoords="offset points",
    color="tab:blue",
)
ax2.set_aspect("equal", adjustable="box")
ax2.grid(True)
ax2.set_title("NACA 2412 contour")
ax2.set_xlabel("x")
ax2.set_ylabel("y")

plt.show()
