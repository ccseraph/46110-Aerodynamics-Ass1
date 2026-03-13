"""
Author: (Cc) Matthew Cook
Group: 17
Class: 46110 Fundamentals of Aerodynamics

------------------------
---- Plot Airfoils ----
------------------------

This code should plot the models for 
NACA 2312 and 2324 airfoils,
and the NACA 4412 and 4424 airfoils
including camber line

"""

import matplotlib.pyplot as plt # I think I will try to use this to plot
from airfoil_toolbox import shape_naca

# NACA 2312

camber, upper, lower = shape_naca(2, 3, 12)             # NACA code is the input, not the percentages/decimals
plt.close('all')                                    
fig, ax = plt.subplots()
ax.plot(upper[:, 0], upper[:, 1], 'r', label="NACA 2312")
ax.plot(lower[:, 0], lower[:, 1], 'r')
ax.plot(camber[:,0], camber[:,1], 'k', linestyle='--')

# NACA 2324

camber, upper, lower = shape_naca(2, 3, 24)
ax.plot(upper[:, 0], upper[:, 1], 'b', label="NACA 2324")
ax.plot(lower[:, 0], lower[:, 1], 'b')
ax.plot(camber[:,0], camber[:,1], 'k', linestyle='--')
ax.axis('equal')
ax.set_title('Compare NACA 2312/24 airfoils')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.legend(loc='best')
plt.grid()
plt.show()

# NACA 4412

plt.close('all') 
camber, upper, lower = shape_naca(4, 4, 12)
fig, ax = plt.subplots()
ax.plot(upper[:, 0], upper[:, 1], 'r', label="NACA 4412")
ax.plot(lower[:, 0], lower[:, 1], 'r')
ax.plot(camber[:,0], camber[:,1], 'k', linestyle='--')

# NACA 4424

camber, upper, lower = shape_naca(4, 4, 24)
ax.plot(upper[:, 0], upper[:, 1], 'b', label="NACA 4424")
ax.plot(lower[:, 0], lower[:, 1], 'b')
ax.plot(camber[:,0], camber[:,1], 'k', linestyle='--')
ax.axis('equal')
ax.set_title('Compare NACA 4412/24 airfoils')
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.legend(loc='best')
plt.grid()
plt.show()

"""
:D
"""