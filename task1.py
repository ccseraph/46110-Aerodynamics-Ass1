import numpy as np
import matplotlib.pyplot as plt

# ──────────────────────────────────────────────
# Known parameters
# ──────────────────────────────────────────────
m_total   = 1.8          # kg
g         = 3.72         # m/s^2
rho       = 0.0174       # kg/m^3

D         = 1.2          # m
R         = D/2
A         = np.pi * R**2

rpm       = 2800
Omega     = rpm * 2*np.pi/60

N_blades  = 2
Cd0       = 0.02
gamma     = 1.15

# Fraction of lower disk area inside upper FAR wake
f_overlap = 0.5  # your "half area sees 2*v_upper"

W = m_total * g  # hover force balance: W = Tu + Tl

# ──────────────────────────────────────────────
# Solve equal-power condition Pu = Pl
# With partial overlap: r = v_l / v_u solves:
# r^3 + 4f r^2 + 4f r - 1 = 0
# ──────────────────────────────────────────────
coeff = [1, 4*f_overlap, 4*f_overlap, -1]
roots = np.roots(coeff)

# pick the positive real root
r_candidates = np.real(roots[np.isreal(roots)])
r = r_candidates[r_candidates > 0][0]

# Thrust ratio: Tl/Tu = r(2f + r)
tl_over_tu = r * (2*f_overlap + r)

# Force balance gives thrusts
T_upper = W / (1 + tl_over_tu)
T_lower = W - T_upper

# Upper induced velocity from Tu = 2 rho A v_u^2
v_upper = np.sqrt(T_upper / (2 * rho * A))

# Lower induced increment v_l = r * v_u
v_lower = r * v_upper

# Consistency check for lower thrust:
# Tl = 2 rho A v_l (2f v_u + v_l)
T_lower_check = 2 * rho * A * v_lower * (2*f_overlap*v_upper + v_lower)

# Induced powers
P_i_upper = gamma * 2 * rho * A * v_upper**3
P_i_lower = gamma * 2 * rho * A * v_lower * (
    f_overlap * (2*v_upper + v_lower)**2 + (1 - f_overlap) * (v_lower**2)
)

# ──────────────────────────────────────────────
# Profile power vs chord sweep
# P0 = (1/8) rho c Nb Cd0 Omega^3 R^4
# ──────────────────────────────────────────────
def profile_power(c):
    return (1/8) * rho * c * N_blades * Cd0 * Omega**3 * R**4

c_over_R = np.linspace(1/50, 1/5, 300)
c_vals   = c_over_R * R
P0_vals  = profile_power(c_vals)

# Total power per rotor
P_upper = P_i_upper + P0_vals
P_lower = P_i_lower + P0_vals
P_total = P_upper + P_lower

# ──────────────────────────────────────────────
# Print key results
# ──────────────────────────────────────────────
print("="*75)
print("Coaxial Hover Model — Partial Far-Wake Overlap + Equal Power Split")
print("="*75)
print(f"Overlap fraction f = {f_overlap:.2f}")
print(f"Weight W = {W:.4f} N\n")

print("--- Solved velocity ratio ---")
print(f"r = v_lower / v_upper = {r:.6f}\n")

print("--- Thrust split (W = Tu + Tl) ---")
print(f"T_upper = {T_upper:.4f} N  ({T_upper/W*100:.2f}% of W)")
print(f"T_lower = {T_lower:.4f} N  ({T_lower/W*100:.2f}% of W)")
print(f"T_lower check (model) = {T_lower_check:.4f} N\n")

print("--- Induced velocities ---")
print(f"v_upper = {v_upper:.4f} m/s")
print(f"v_lower = {v_lower:.4f} m/s")
print(f"Upper far wake velocity = 2*v_upper = {2*v_upper:.4f} m/s\n")

print("--- Induced powers (should be equal) ---")
print(f"P_i_upper = {P_i_upper:.6f} W")
print(f"P_i_lower = {P_i_lower:.6f} W")

# ──────────────────────────────────────────────
# Plot power vs c/R
# ──────────────────────────────────────────────
plt.figure(figsize=(8.6, 5.2))
plt.plot(c_over_R, P_upper, label="Upper rotor power")
plt.plot(c_over_R, P_lower, "--", label="Lower rotor power")
plt.plot(c_over_R, P_total, "k", linewidth=2, label="Total system power")

plt.xlabel("Chord ratio c/R [-]")
plt.ylabel("Power [W]")
plt.title("Power vs chord ratio c/R (1/50 to 1/5)\nPartial far-wake overlap, equal power split")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()