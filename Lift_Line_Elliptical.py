"""
Author: (Cc) Matthew Cook
Group: 17
Class: 46110 Fundamentals of Aerodynamics
-----------------------------
---- Lifting Line -----------
---- Elliptical wing --------
-----------------------------
"""

import sys
sys.path.insert(0, 'from_prof')         # fixes our path problem
import numpy as np
import matplotlib.pyplot as plt

# our function
from airfoil_toolbox import solve_panel_method

# set your NACA airfoil code
NACA_code = 2410
# set Reynold's number
Re = 5*10**6
# Aspect ratio values
AR = np.array([4, 6, 8, 10, 10**5])         # I just did a big number instead of infinity
# Angles of attack
alpha_geo = np.arange(-4, 9, 1)

# get zero lift angle using panel method from last assingment
Cl_list = []                                                              
for AoA in alpha_geo:                               # loop using solver for Cl's
    results = solve_panel_method(NACA_code, AoA)    # panel method
    Cl_list.append(results["Cl"])                   # saving lift coeficients     
Cl_list = np.array(Cl_list) 
alpha_zero = np.interp(0, Cl_list, alpha_geo)       # finding zero lift angle with panel method
print(f"Zero lift AoA: {alpha_zero:.3f} degrees")   

# lifting line C_L vs AoA
C_L_alpha = (2 * np.pi) / (1 + (2 / AR))
C_L_alpha_deg = C_L_alpha / (180 / np.pi) 

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))

for i, ar in enumerate(AR):
    label = f"AR = {ar}" if ar < 20 else "AR = ∞"
    C_L = C_L_alpha_deg[i] * (alpha_geo - alpha_zero)
    alpha_i = C_L / (np.pi * AR[i])             # induced AoA. Equations from lecture
    alpha_i_deg = np.degrees(alpha_i)           # convert to degrees for plotting
    C_di = C_L * alpha_i                        # our induced drag
    
    ax1.plot(alpha_geo, C_L, marker='x', label=f"{label}")
    ax2.plot(alpha_geo, C_di, marker='x', label=f"{label}")
    ax3.plot(alpha_geo, alpha_i_deg, marker='x', label=label)

# now we plot everything       
ax1.set_xlabel("Angle of Attack α (degrees)")
ax1.set_ylabel("$C_L$")
ax1.set_title("Lift Coefficient $C_L$")
ax1.legend()
ax1.grid()

ax2.set_xlabel("Angle of Attack α (degrees)")
ax2.set_ylabel("$C_{D,i}$")
ax2.set_title("Induced Drag Coefficient $C_{D,i}$")
ax2.legend()
ax2.grid()

ax3.set_xlabel("Angle of Attack α (degrees)")
ax3.set_ylabel("$α_i$ (degrees)")
ax3.set_title("Induced Angle of Attack $α_i$")
ax3.legend()
ax3.grid()

plt.suptitle("Elliptic Wing — Lifting Line Theory (NACA 2410)")
plt.tight_layout()
plt.show()

# I should include something to use xfoil and aquire C_d.
# Then I can add a plot of C_D for the whole wingspan