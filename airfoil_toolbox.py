"""
Author: (Cc) Matthew Cook
Group: 17
Class: 46110 Fundamentals of Aerodynamics
------------------------------
---- Aerodynamcs Toolbox ----
------------------------------

This houses all the functions we create to run our code.
Scripts can import functions from here. 
"""
import numpy as np 
from from_prof.funaerotool.panel_method.solver import solve_closed_contour_panel_method  
# Magic function to solve all our panel method problems ^^

def parse_naca(code):
    # give this function the 4 digit code of the NACA airfoil
    m = code // 1000
    p = (code // 100) % 10
    xx = code % 100
    return m, p, xx


def shape_naca(m, p, xx, c=1, N=200):
    """
    This function takes inputs as the NACA airfoil codes (2312, 4412, etc.)
    and plots the airfoil surface. I think it will break if camber does not equal 1.  
    """
    m = m/100           # 1st char in NACA code, max camber as a percentage of the chord
    p = p/10            # 2nd char in NACA code, location of the max camber as a percentage of the chord
    xx = xx/100         # 3rd and 4th char in NACA code, thickness from the last 2 characters as a percentage of chord length
    t = c*xx            # max thickness
    x_vals = np.linspace(0, c, N)   # linspace is better than arange since we need c inclusive

    pos_camber = np.zeros(shape=(N,2))          # arrays for storing data to plot
    pos_upper = np.zeros(shape=(N,2))    
    pos_lower = np.zeros(shape=(N,2))    

    for i, x in enumerate(x_vals):            # enumerate to track iteration counter
        
        xi = x/c                      # xi is the position between 0 and c along chord
        if xi >= 0 and xi <= p:                     # camber line based on class equations
            y_camber = (m/p**2)*(2*p*xi - xi**2)    # this include derivatives I calulated so could be wrong
            dy_dx = 2*m/(c*p**2)*(p - xi)
        elif xi > p and xi <= 1:
            y_camber = (m/(1-p)**2)*(1 - 2*p + 2*p*xi - xi**2)
            dy_dx = 2*m/(c*(1-p)**2)*(p - xi)
        else:
            break
        
        pos_camber[i] = [x, y_camber]       # storing camber data    
        
        y_thick = 5*t*(0.2969*np.sqrt(xi) - 0.1260*xi - 0.3516*(xi)**2 + 0.2843*(xi)**3 - 0.1015*(xi)**4)        # Equation for thickness
        
        theta = np.arctan(dy_dx)                  # upper and lower surface positions calculated 
        x_u = x - y_thick*np.sin(theta)         # using theta
        y_u = y_camber + y_thick*np.cos(theta)
        x_l = x + y_thick*np.sin(theta)
        y_l = y_camber - y_thick*np.cos(theta)

        pos_upper[i] = [x_u, y_u]               # upper surface position array
        pos_lower[i] = [x_l, y_l]               # lower

    # force exact trailing edge
    # might not need to force it
    # pos_upper[-1] = [1.0, 0.0]
    # pos_lower[-1] = [1.0, 0.0]
    return pos_camber, pos_upper, pos_lower


def solve_panel_method(NACA_code, AoA, N=201):
    """
    This function solves your panel method
    """
    # defining airfoil shape
    m, p, xx = parse_naca(NACA_code)

    # getting our surfaces
    pos_camber, pos_upper, pos_lower = shape_naca(m, p, xx, c=1, N=N)

    # taping them together into x's and y's
    airfoil_surface = np.concatenate([pos_upper[::-1], pos_lower[1:]])
    xs = airfoil_surface[:, 0]
    ys = airfoil_surface[:, 1]

    # magic solving function from class
    results = solve_closed_contour_panel_method(xs, ys, aoa_deg=AoA, kutta_condition=True)
    """
    results has these arguments:
        "sigma": sigma,
        "strengths": sigma,
        "circulation": circulation,
        "gamma": gamma,
        "Vt": Vt,
        "Vn": Vn,
        "Cp": Cp,
        "Cl": Cl,
        "xp": xp,
        "yp": yp,
        "panel_lengths": plength,
        "Tx": Tx,
        "Ty": Ty,
        "Nx": Nx,
        "Ny": Ny,
    use them by saying something like 
    Cl = results["Cl"]
    """
    return results

# consider adding a function here which can mamke all our plots
# look identical and nice looking for the overleaf


