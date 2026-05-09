
import sys
sys.path.insert(0, 'from_prof')         # fixes our path problem
from airfoil_toolbox import run_case, read_polar_file
import numpy as np
import matplotlib.pyplot as plt
import os

xfoil_path = r"C:\Users\mcook\OneDrive - Danmarks Tekniske Universitet\Documents\DTU\Spring 2026\46110 Aerodynamics\xfoil.exe"          # or full path to your xfoil executable
results_folder = "xfoil_results"  # folder where polar files get written
os.makedirs(results_folder, exist_ok=True) 
os.makedirs("Figures", exist_ok=True)

# Run free transition (no turbulator)
polar_free = run_case(
    airfoil="2402",
    xfoil_path=xfoil_path,
    results_folder=results_folder,
    re=10640,
    mach=0.0,
    ncrit=9,
    alpha_start=-5,
    alpha_end=15,
    alpha_step=0.5
)

alpha, cl, cd = read_polar_file(polar_free)

# Cl vs AoA
plt.figure()
plt.plot(alpha, cl)
plt.xlabel("AoA (deg)")
plt.ylabel("$C_l$")
plt.grid(True)
plt.savefig("Figures/task4_cl_aoa.png", dpi=150)

# Cl vs Cd polar
plt.figure()
plt.plot(cd, cl)
plt.xlabel("$C_d$")
plt.ylabel("$C_l$")
plt.grid(True)
plt.savefig("Figures/task4_cl_cd.png", dpi=150)