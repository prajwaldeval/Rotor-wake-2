import numpy as np
import matplotlib.pyplot as plt
from defs import *
from old_defs import *

number_of_blades = 3
R = 50  # m
r_hub = 0.2 * R
U0 = 10
a = 0.25

# Constant or cosine element spacing on the blade
r, dr = geometry_constant(r_hub, R, 20)
# r,dr = geometry_cosine(r_hub,R,30)


wakelength = 2  # how many dimaters long the wake shall be prescribed for
nt = 50
tip_speed_ratio = 8


rotor_vortex_system = wake_system_generation(r, dr, U0, a, wakelength, number_of_blades, tip_speed_ratio)
Ua, Va, Wa = unit_strength_induction_matrix(cps, fils, nr_blade_elements, number_of_blades)

a_new, gamma, Fax_ll, Faz_ll = iteration(iterations, Ua, Va, Wa,
                                         cps, tip_speed_ratio,
                                         gamma_convergence_weight, error_limit, R)


# alphas = np.zeros(n)
# phis = np.zeros(n)
# a = np.zeros(n)
# ap = np.zeros(n)
# Fax = np.zeros(n)
# Faz = np.zeros(n)
# gamma = np.zeros(n)
# Cax = np.zeros(n)
# Caz = np.zeros(n)
# cls = np.zeros(n)
# cds = np.zeros(n)
#
#
#
# CT = 0
# CP = 0
# CQ = 0
#
# for i in range(len(r_arr)):
#     a[i],ap[i],Fax[i],Faz[i],gamma[i], phi, alpha, cls[i], cds[i] = solve_anul(U0,r_arr[i],dr_arr[i],c[i],b[i],Prandtl=True, Glauert=True)
#     alphas[i] = alpha
#     phis[i] = phi
#     Cax[i] = Fax[i] * nb/(0.5*rho*U0**2*np.pi*R**2*c[i])
#     Caz[i] = Faz[i] * nb/(0.5*rho*U0**2*np.pi*R**2*c[i])
#     CT += dr_arr[i] * Fax[i] * nb/(0.5*rho*U0**2*np.pi*R**2)
#     CP += dr_arr[i] * Faz[i] * nb * R * omega/(0.5*rho*U0**3*np.pi*R**2)
#     CQ += r_arr[i] * dr_arr[i] * Faz[i] * nb/(0.5*rho*U0**2*np.pi*R**2)
#     print(i,'finished')
#
# print('Finished with CT =',CT,'and CP =',CP)
# r_R = r_arr/R
# clcd = cls/cds
# plt.plot(r_R[1:],clcd[1:]*c[1:])
# plt.title('TSR= '+str(TSR))
# plt.grid()
# plt.xlabel(r'$r/R$')
# plt.ylabel(r'($C_{l}c$')
# plt.show()
