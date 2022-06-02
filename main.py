import numpy as np
import matplotlib.pyplot as plt
import defs as ll
import bem_defs as bem

number_of_blades = 3
R = 50  # m
r_hub = 0.2 * R
U0 = 10
a = 0.25
tip_speed_ratio = 8
rho = 1.225
U0 = 10

# values
number_elements = 25
wake_length = 2.5
dtheta = np.pi / 30
iterations = 200
gamma_convergence_weight = 0.3
error_limit = 0.001

# Constant or cosine element spacing on the blade
r, dr = ll.geometry_constant(r_hub, R, number_elements)

## LL solution
cps, fils = ll.wake_system_generation(r, dr, U0, a, wake_length, number_of_blades, tip_speed_ratio, dtheta)
print('wake generated')
Ua, Va, Wa = ll.unit_strength_induction_matrix(cps, fils, number_elements, number_of_blades)
print('induction matrices created')
a_new, gamma, Fax_ll, Faz_ll, phi_ll, alpha_ll = ll.iteration(iterations, Ua, Va, Wa,
                                            cps, tip_speed_ratio,
                                            gamma_convergence_weight, error_limit, R, U0, rho)
print('converged LL')
CT_ll, CP_ll = ll.coefficients(Fax_ll, Faz_ll, r, dr, number_of_blades, rho, U0, tip_speed_ratio)
## BEM Solution
c = bem.chord(r / R)
b = bem.twist(r / R)
omega = U0*tip_speed_ratio/R

alphas_b = np.zeros(number_elements)
phis_b = np.zeros(number_elements)
a_arr_b = np.zeros(number_elements)
ap_b = np.zeros(number_elements)
Fax_b = np.zeros(number_elements)
Faz_b = np.zeros(number_elements)
gamma_b = np.zeros(number_elements)
Cax_b = np.zeros(number_elements)
Caz_b = np.zeros(number_elements)
cls_b = np.zeros(number_elements)
cds_b = np.zeros(number_elements)

CT_b = 0
CP_b = 0

for i in range(len(r)):
    a_arr_b[i], ap_b[i], Fax_b[i], Faz_b[i], gamma_b[i], phis_b[i], alphas_b[i], cls_b[i], cds_b[i] = bem.solve_anul(U0, r[i], dr[i], c[i], b[i], True, True, omega)

    Cax_b[i] = Fax_b[i] * number_of_blades / (0.5 * rho * U0 ** 2 * np.pi * R ** 2 * c[i])
    Caz_b[i] = Faz_b[i] * number_of_blades / (0.5 * rho * U0 ** 2 * np.pi * R ** 2 * c[i])
    CT_b += dr[i] * Fax_b[i] * number_of_blades / (0.5 * rho * U0 ** 2 * np.pi * R ** 2)
    CP_b += dr[i] * Faz_b[i] * r[i] * omega * number_of_blades / (0.5 * rho * U0 ** 3 * np.pi * R ** 2)
    print(r[i], 'Converged BEM')

##plotting


#alpha
