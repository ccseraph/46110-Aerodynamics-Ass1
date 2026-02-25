"""

Created by Cc for Aerodynamics Group 17
24-02-2026
46110 Fundamentals of Aerodynamics

This code should plot the models for 
NACA 2312 and 2324 airfoils,
and the NACA 4412 and 4424 airfoils
including camber line

"""

import numpy as np # trying to import numpy functions
import matplotlib.pyplot as plt # I think I will try to use this to plot


def shape_naca(m, p, xx, c=1, N=200):
    """
    This function takes inputs as the NACA airfoil codes (2312, 4412, etc.)
    and plots the airfoil surface. I think it will break if camber does not equal 1.  
    """
    m = m/100           # 1stchar in NACA code, max camber as a percentage of the chord
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
        
        y_thick = 5*t*(0.2969*np.sqrt(x) - 0.1260*x - 0.3516*(x)**2 + 0.2843*(x)**3 - 0.1015*(x)**4)        # Equation for thickness
        
        theta = np.atan(dy_dx)                  # upper and lower surface positions calculated 
        x_u = x - y_thick*np.sin(theta)         # using theta
        y_u = y_camber + y_thick*np.cos(theta)
        x_l = x + y_thick*np.sin(theta)
        y_l = y_camber - y_thick*np.cos(theta)

        pos_upper[i] = [x_u, y_u]               # upper surface position array
        pos_lower[i] = [x_l, y_l]               # lower

    return pos_camber, pos_upper, pos_lower


# NACA 2312

camber, upper, lower = shape_naca(2, 3, 12)             # NACA code is the input, not the percentages/decimals
plt.close('all')                                    
fig, ax = plt.subplots()
ax.plot(upper[:, 0], upper[:, 1], 'r', label="$NACA 2312$")
ax.plot(lower[:, 0], lower[:, 1], 'r')
ax.plot(camber[:,0], camber[:,1], 'k')

# NACA 2324

camber, upper, lower = shape_naca(2, 3, 24)
ax.plot(upper[:, 0], upper[:, 1], 'b', label="$NACA 2324$")
ax.plot(lower[:, 0], lower[:, 1], 'b')
ax.plot(camber[:,0], camber[:,1], 'k')
ax.axis('equal')
ax.set_title('Compare two NACA airfoils')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.legend(loc='best')
plt.grid()
plt.show()

# NACA 4412

plt.close('all') 
camber, upper, lower = shape_naca(4, 4, 12)
fig, ax = plt.subplots()
ax.plot(upper[:, 0], upper[:, 1], 'r', label=f"$NACA 4412$")
ax.plot(lower[:, 0], lower[:, 1], 'r')
ax.plot(camber[:,0], camber[:,1], 'k')

# NACA 4424

camber, upper, lower = shape_naca(4, 4, 24)
ax.plot(upper[:, 0], upper[:, 1], 'b', label=f"$NACA 4424$")
ax.plot(lower[:, 0], lower[:, 1], 'b')
ax.plot(camber[:,0], camber[:,1], 'k')
ax.axis('equal')
ax.set_title('Compare two NACA airfoils')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.legend(loc='best')
plt.grid()
plt.show()

"""
:D
"""