"""
Author: Lorenzo Ceresole
Class: 46110 Fundamentals of Aerodynamics
----------------------
---- Solve with XFOIL ----
----------------------
"""

import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# path to xfoil
xfoil_path = r"C:\Users\Lorenzo Ceresole\Documents\Wind Turby\Xfoil\xfoil.exe"

#naca Airfoil codes
airfoils = ["2312", "2324", "4412", "4424"]

# cases for xfoil
cases = ["free", "fixed"]

# flow settings
re = 1.5e6
mach = 0.0
ncrit = 9


alpha_start = -10
alpha_end = 15
alpha_step = 1

alpha_cp = 10

# output folder
results_folder = "xfoil_results"
os.makedirs(results_folder, exist_ok=True)


# --------------------------
# ---- Run XFOIL case   ----
# --------------------------
def run_case(airfoil, case):
    polar_file = f"{airfoil}_{case}_polar.txt"
    cp_file = f"{airfoil}_{case}_cp.txt"

    polar_path = os.path.join(results_folder, polar_file)
    cp_path = os.path.join(results_folder, cp_file)

    # delete old files so xfoil does not append to them
    if os.path.exists(polar_path):
        os.remove(polar_path)

    if os.path.exists(cp_path):
        os.remove(cp_path)

    commands = []
    commands.append(f"NACA {airfoil}")
    commands.append("OPER")
    commands.append(f"VISC {re}")
    commands.append(f"MACH {mach}")
    commands.append("ITER 150")
    commands.append("VPAR")
    commands.append(f"N {ncrit}")

    # free transition = do nothing
    # fixed transition = upper side forced at x/c = 0.1
    if case == "fixed":
        commands.append("XTR 0.1 1")

    commands.append("")
    commands.append("PACC")
    commands.append(polar_file)
    commands.append("")
    commands.append(f"ASEQ {alpha_start} {alpha_end} {alpha_step}")
    commands.append("PACC")
    commands.append(f"ALFA {alpha_cp}")
    commands.append(f"CPWR {cp_file}")
    commands.append("")
    commands.append("QUIT")

    text = "\n".join(commands) + "\n"

    subprocess.run(xfoil_path,input=text,text=True,cwd=results_folder)
    return polar_path, cp_path


# ------------------------------
# ---- Read XFOIL polar file ----
# ------------------------------
def read_polar_file(file_path):
    alpha = []
    cl = []
    cd = []

    with open(file_path, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            try:
                a = float(parts[0])
                cl_val = float(parts[1])
                cd_val = float(parts[2])
                alpha.append(a)
                cl.append(cl_val)
                cd.append(cd_val)
            except:
                pass

    alpha = np.array(alpha)
    cl = np.array(cl)
    cd = np.array(cd)

    if len(alpha) == 0:
        return alpha, cl, cd

    # remove duplicates
    alpha_unique, idx = np.unique(alpha, return_index=True)
    cl = cl[idx]
    cd = cd[idx]
    alpha = alpha_unique

    # sort by alpha
    sort_idx = np.argsort(alpha)
    alpha = alpha[sort_idx]
    cl = cl[sort_idx]
    cd = cd[sort_idx]

    return alpha, cl, cd


# ---------------------------
# ---- Read XFOIL Cp file ---
# ---------------------------
def read_cp_file(file_path):
    x = []
    cp = []

    with open(file_path, "r") as f:
        for line in f:
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                x_val = float(parts[0])
                cp_val = float(parts[1])
                x.append(x_val)
                cp.append(cp_val)
            except:
                pass

    return np.array(x), np.array(cp)


# --------------------------
# ---- Solve all cases  ----
# --------------------------
all_data = {}

for airfoil in airfoils:
    all_data[airfoil] = {}

    for case in cases:
        polar_path, cp_path = run_case(airfoil, case)
        alpha, cl, cd = read_polar_file(polar_path)
        x_cp, cp = read_cp_file(cp_path)

        all_data[airfoil][case] = {
            "alpha": alpha,
            "cl": cl,
            "cd": cd,
            "x_cp": x_cp,
            "cp": cp
        }


# ----------------------------
# ---- Part 2 Lift plots  ----
# ----------------------------
for airfoil in airfoils:
    plt.figure()

    for case in cases:
        alpha = all_data[airfoil][case]["alpha"]
        cl = all_data[airfoil][case]["cl"]
        plt.plot(alpha, cl, marker="o", label=case)

    plt.xlabel("AoA (deg)")
    plt.ylabel("$C_l$")
    plt.title(f"XFOIL - $C_l$ vs AoA [NACA {airfoil}]")
    plt.grid()
    plt.legend()
    plt.savefig(os.path.join(results_folder, f"{airfoil}_cl_vs_aoa.png"), dpi=300, bbox_inches="tight")
    plt.show()


# ----------------------------
# ---- Part 2 Drag plots  ----
# ----------------------------
for airfoil in airfoils:
    plt.figure()

    for case in cases:
        alpha = all_data[airfoil][case]["alpha"]
        cd = all_data[airfoil][case]["cd"]
        plt.plot(alpha, cd, marker="o", label=case)

    plt.xlabel("AoA (deg)")
    plt.ylabel("$C_d$")
    plt.title(f"XFOIL - $C_d$ vs AoA [NACA {airfoil}]")
    plt.grid()
    plt.legend()
    plt.savefig(os.path.join(results_folder, f"{airfoil}_cd_vs_aoa.png"), dpi=300, bbox_inches="tight")
    plt.show()


# --------------------------------------
# ---- Part 4 Pressure distribution ----
# --------------------------------------
for airfoil in airfoils:
    plt.figure()

    for case in cases:
        x_cp = all_data[airfoil][case]["x_cp"]
        cp = all_data[airfoil][case]["cp"]
        plt.plot(x_cp, cp, marker=".", label=case)

    plt.gca().invert_yaxis()
    plt.xlabel("x/c")
    plt.ylabel("$C_p$")
    plt.title(f"XFOIL - $C_p$ distribution at $\\alpha$ = {alpha_cp}° [NACA {airfoil}]")
    plt.grid()
    plt.legend()
    plt.savefig(os.path.join(results_folder, f"{airfoil}_cp_distribution.png"), dpi=300, bbox_inches="tight")
    plt.show()
    