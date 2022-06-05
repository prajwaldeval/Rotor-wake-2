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

# values
number_elements = 50
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
a_new, gamma_ll, Fax_ll, Faz_ll, alpha_ll, phi_ll = ll.iteration(iterations, Ua, Va, Wa,
                                            cps, tip_speed_ratio,
                                            gamma_convergence_weight, error_limit, R, U0, rho)
Cax_ll = Fax_ll * number_of_blades/ (0.5 * rho * U0 ** 2 * R)
Caz_ll = Faz_ll * number_of_blades/ (0.5 * rho * U0 ** 2 * R)

print('converged LL')
CT_ll, CP_ll = ll.coefficients(Fax_ll, Faz_ll, r, dr, number_of_blades, rho, U0, tip_speed_ratio)
## BEM Solution
c = bem.chord(r / R)
b = bem.twist(r / R) + 2
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

    Cax_b[i] = Fax_b[i] * number_of_blades / (0.5 * rho * U0 ** 2 * R)
    Caz_b[i] = Faz_b[i] * number_of_blades / (0.5 * rho * U0 ** 2 * R)
    CT_b += dr[i] * Fax_b[i] * number_of_blades / (0.5 * rho * U0 ** 2 * np.pi * R ** 2)
    CP_b += dr[i] * Faz_b[i] * r[i] * omega * number_of_blades / (0.5 * rho * U0 ** 3 * np.pi * R ** 2)
    print(r[i], 'Converged BEM')




##plotting ll = darkolivegreen bem= indianred
'''
#α&φ
alphi = plt.figure()
plt.plot(r/R, alphas_b, c='indianred', label=r'$\alpha$ BEM')
plt.plot(r[:number_elements-1]/R, alpha_ll[:number_elements-1], c='darkolivegreen', label=r'$\alpha$ LL')

plt.plot(r/R, phis_b, c='indianred', label=r'$\phi$ BEM', linestyle='dashed')
plt.plot(r[:number_elements-1]/R, phi_ll[:number_elements-1], c='darkolivegreen', label=r'$\phi$ LL', linestyle = 'dashed')

plt.ylabel(r'angle [$^\circ$]')
plt.xlabel(r'$r/R$')
plt.legend()
plt.show()

#alpha
alphafig = plt.figure()
plt.plot(r/R, alphas_b, c='indianred', label='BEM')
plt.plot(r[:number_elements-1]/R, alpha_ll[:number_elements-1], c='darkolivegreen', label='LL')

plt.ylabel(r'$\alpha$ [$^\circ$]')
plt.xlabel(r'$r/R$')
plt.legend()
plt.show()
#phi
phifig = plt.figure()
plt.plot(r/R, phis_b, c='indianred', label='BEM')
plt.plot(r[:number_elements-1]/R, phi_ll[:number_elements-1], c='darkolivegreen', label='LL')

plt.ylabel(r'$\phi$ [$^\circ$]')
plt.xlabel(r'$r/R$')
plt.legend()
plt.show()

#vorticity
gammafig = plt.figure()
gnorm = np.pi*U0*U0/omega/number_of_blades

plt.plot(r/R, gamma_b/gnorm, c='indianred', label='BEM')
plt.plot(r[:-1]/R, gamma_ll[:number_elements-1]/gnorm, c='darkolivegreen', label='LL')
plt.ylabel(r'$\Gamma \cdot \frac{\omega n_{blades}}{\pi * U_{\inf}^{2}}$')
plt.xlabel(r'$r/R$')
plt.legend()
plt.show()




#Cax
caxfig = plt.figure()
plt.plot(r/R, Cax_b, c='indianred', label='BEM')
plt.plot(r[:number_elements-1]/R, Cax_ll[:number_elements-1], c='darkolivegreen', label='LL')
plt.ylabel(r'$C_{ax}$')
plt.xlabel(r'$r/R$')
plt.legend()
plt.show()

#Caz
cazfig = plt.figure()
plt.plot(r/R, Caz_b, c='indianred', label='BEM')
plt.plot(r[:number_elements-1]/R, Caz_ll[:number_elements-1], c='darkolivegreen', label='LL')
plt.ylabel(r'$C_{az}$')
plt.xlabel(r'$r/R$')
plt.legend()
plt.show()
'''