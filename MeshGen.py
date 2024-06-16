# ===========================================================
# Comp Fluid Dyn - Final Project
# This script generates a structured mesh for the Wedge
# ===========================================================
import numpy as np
import matplotlib.pyplot as plt
from math import exp
import os

# gemoetry configuration
AB = 1
AE = 2.4
DE = 4
angle = 15

mu = -1 # stretching factor

# corner points
x_l = y_b = 0
x_B = AB
x_r = DE
y_t = AE
y_C = (DE - AB) * np.tan(np.deg2rad(angle)) + y_b

# grid resolution
nx = 121 # Number of grid points in x-direction
ny = 61 # Number of grid points in y-direction
nx1 = int(nx * (AB / DE))
nx2 = nx + 1 - nx1
nx = nx1 + nx2 - 1

# generate domain
X_b1 = np.linspace(x_l , x_B , nx1)
X_b2 = np.linspace(x_B , x_r , nx2)
Y_b1 = np.ones(nx1) * y_b
Y_b2 = y_b + (X_b2 - x_B) * (np.deg2rad(angle))
X_t1 = np.linspace(x_l , x_B , nx1)
X_t2 = np.linspace(x_B , x_r , nx2)
X_t = np.linspace(x_l , x_r , nx)
X_bottom = np.concatenate ((X_b1 , X_b2[1:]))
Y_bottom = np.concatenate ((Y_b1 , Y_b2[1:]))
X_top = np.concatenate ((X_t1 , X_t2[1:]))
Y_top = np.ones(nx) * y_t

# generate grid
x = np.zeros ((nx , ny)) # x coordinate matrix
y = np.zeros ((nx , ny)) # y coordinate matrix
for i in range(0,nx):
    for j in range(0,ny):
        x[i, j] = X_bottom[i] + (X_top[i] - X_bottom[i]) * (j) /(ny-1)
        rescaling_j = (j) / (ny - 1)
        rescaling_j = (exp(mu * rescaling_j) - 1) / (exp(mu) - 1)
        y[i, j] = Y_bottom[i] + (Y_top[i] - Y_bottom[i]) * rescaling_j
        
# visualize
plt.figure ()
plt.plot(x,y,'k',linewidth=0.5)
plt.plot(np.transpose(x),np.transpose(y),'k',linewidth=0.5)
plt.axis('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

# save mesh
np.save('Mesh/x.npy', x)
np.save('Mesh/y.npy', y)

# parameters for calculations
param = open('Mesh/param.ini', 'w')
param.write(str(nx) + '\n')
param.write(str(ny) + '\n')
