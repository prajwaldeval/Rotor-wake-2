import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pandas as pd
import numpy as np
from defs import *

number_of_blades = 3
R = 50  # m
r_hub = 0.2 * R
U0 = 10
a = 0.25
nr_blade_elements = 7
rho = 1.225
# Constant or cosine element spacing on the blade
r, dr = geometry_constant(r_hub, R, nr_blade_elements)
# r, dr = geometry_cosine(r_hub, R, nr_blade_elements)
iterations = 200
gamma_convergence_weight = 0.3
error_limit = 0.01

wakelength = 0.3 # how many diameters long the wake shall be prescribed for
nt = 50
tip_speed_ratio = 6


cps, fils = wake_system_generation(r, dr, U0, a, wakelength, number_of_blades, tip_speed_ratio)

plane_dist = 1#m

res = 20
y,z = np.meshgrid(np.linspace(-50,50,res),np.linspace(-50,50,res))
v_ind = np.zeros((res,res))
w_ind = np.zeros((res,res))
print('generated mesh')
for i in range(res):
    print(i)

    for j in range(res):
        print(j)
        for k in range(len(fils['x1'])):
            print(k)
            u,v,w = biot_savart_function(fils["x1"][k],fils["x2"][k],fils["z1"][k],fils["x2"][k],fils["y2"][k],fils["z2"][k],plane_dist,y[i,j],z[i,j],0,0.1)
            v_ind[i,j] =+ v
            w_ind[i, j] = + w
plt.quiver(y,z,v,w)
plt.show()












