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


wakelength = 4  # how many dimaters long the wake shall be prescribed for
nt = 50
tip_speed_ratio = 8


rotor_vortex_system = wake_system_generation(r, dr, U0, a, wakelength, number_of_blades, tip_speed_ratio)

