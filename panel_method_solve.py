"""
Author: (Cc) Matthew Cook
Group: 17
Class: 46110 Fundamentals of Aerodynamics
-----------------------------
---- Solve Panel Method ----
-----------------------------
"""

# Magic function to solve all our panel method problems
from from_prof.funaerotool.panel_method.solver import solve_closed_contour_panel_method  

# our functions
from airfoil_toolbox import solve_panel_method
import numpy as np
import matplotlib.pyplot as plt

# set your NACA airfoil code
# maybe:  2312, 2324, 4412, 4424
NACA_code = 2312

# set number of panels 
# (this should probably be odd and on the order of 100)
N = 201

# ----------------------------
# ---- Part 2 Ranging AoA ---- 
# ----------------------------
AoA_range = np.linspace(-10, 15, 26)    # range of alpha -10 to 25
Cl_list = []                            # saving lift coeficients

for AoA in AoA_range:
    results = solve_panel_method(NACA_code, AoA, N)     # halelujah chorus
    Cl_list.append(results["Cl"])                       # collecting our Cl's

plt.plot(AoA_range, Cl_list, 'salmon', label=f"NACA {NACA_code}")
plt.xlabel("AoA (°)")
plt.ylabel("$C_L$")
plt.title("Panel Method - $C_l$ vs AoA")
plt.legend()
plt.grid()
plt.show()


# ------------------------------------
# ---- Part 3 Pressure difference ---- 
# ------------------------------------
# set angle of attack in degrees
AoA = 10
# halelujah chorus
results = solve_panel_method(NACA_code, AoA, N)

Cp = results["Cp"]
xp = results["xp"]

# There's probably an easier way than to split xp's but here we are
split = len(xp) // 2

# upper surface reverse (LE → TE)
xp_upper = xp[:split][::-1]
Cp_upper = Cp[:split][::-1]

# lower surface (LE→TE)
xp_lower = xp[split:]
Cp_lower = Cp[split:]

# interpolate lower onto upper x points for subtraction
Cp_lower_interp = np.interp(xp_upper, xp_lower, Cp_lower)

# delta Cp = upper - lower (I think?) (What was that other equation?)
delta_Cp = Cp_upper - Cp_lower_interp

fig, ax = plt.subplots()
ax.plot(xp_upper, delta_Cp, 'salmon', label=f"$\\alpha$ = {AoA}°")
ax.set_title(f'Panel Method - $\\Delta C_p$ distribution [NACA {NACA_code}]')
ax.set_xlabel("x/c")
ax.set_ylabel("$\\Delta C_p$")
ax.legend(loc='best')
plt.grid()
plt.show()


# --------------------------------------
# ---- Part 4 Pressure distribution ---- 
# --------------------------------------

# Same angle of attack, Cp, and xp as part 3

plt.close('all') 
fig, ax = plt.subplots()
ax.plot(xp, Cp, 'salmon', label=f"$\\alpha$ = {AoA}°")
ax.set_title(f'Panel Method - $C_P$ distribution [NACA {NACA_code}]')
ax.set_xlabel("x/c")
ax.set_ylabel("$C_P$")
ax.legend(loc='best')
plt.grid()
plt.show()



