import matplotlib.pyplot as plt
import numpy as np

from funaerotool.utils import generate_naca4_contour


fig, axes = plt.subplots(1, 3, figsize=(14, 4.8), constrained_layout=True)
colors = ["tab:blue", "tab:orange", "tab:green"]


def naca_code_from_params(m: float, p: float, t: float) -> str:
    return f"{int(round(100 * m)):1d}{int(round(10 * p)):1d}{int(round(100 * t)):02d}"


base = {"m": 0.04, "p": 0.4, "t": 0.12}

variation_definitions = [
    ("m", np.array([0.00, 0.02, 0.04]), "max camber m"),
    ("p", np.array([0.20, 0.40, 0.60]), "camber location p"),
    ("t", np.array([0.09, 0.12, 0.18]), "thickness t"),
]

for ax, (parameter_name, values, title) in zip(axes, variation_definitions):
    for idx, value in enumerate(values):
        params = dict(base)
        params[parameter_name] = float(value)
        code = naca_code_from_params(params["m"], params["p"], params["t"])
        x, y = generate_naca4_contour(naca_code=code, n_points=401)
        ax.plot(
            x,
            y,
            linewidth=2,
            color=colors[idx % len(colors)],
            label=f"NACA {code}",
        )

    ax.scatter([0.0, 1.0], [0.0, 0.0], color=["tab:red", "k"], s=25, zorder=5)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"Vary {title}")
    ax.legend(fontsize=8)

plt.show()
