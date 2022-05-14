import numpy as np
import matplotlib.pyplot as plt

from defs import *

number_of_blades = 3
R = 50  # m
r_hub = 0.2 * R
U0 = 10
a = 0.25

# Constant or cosine element spacing on the blade
r, dr = geometry_constant(r_hub, R, 20)
# r,dr = geometry_cosine(r_hub,R,30)

w =
wakelength = 4  # how many dimaters long the wake shall be prescribed for
nt = 50


rotor_vortex_system = wake_system_generation(r, dr, U0, w, a, wakelength, nt, 10, number_of_blades)  # don't know where the 10 comes from