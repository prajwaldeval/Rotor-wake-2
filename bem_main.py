import matplotlib.pyplot as plt
import numpy as np
from defs import *


constant = True # True is we want constant, False if we want cosine

r_cos,dr_cos = geometry_cosine(n)
r_con,dr_con = geometry_constant(n)

if constant:
    c = chord(r_con / R)
    b = twist(r_con / R)
    r_arr = r_con
    dr_arr = dr_con
else:
    c = chord(r_cos / R)
    b = twist(r_cos / R)
    r_arr = r_cos
    dr_arr = dr_cos

'''Plotting for constant TSR'''
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


'''Plotting annulus expansion & entropy jump'''
# r1,r23,r4,a1,a23,a4 = areas(a,r_arr,dr_arr)
# jump = enthalpy_jump(a23,Fax)
# p = 101325 #Pa, sea level pressure
# E1 = (p/rho + U0**2/2)*np.ones(n) #entropy at 1
# E2 = E1 #entropy at 2
# E3 = E2 + jump #entropy at 3
# E4 = E3 #entropy at 4
#
# r1c = (r1[:-1]+r1[1:])/(2*max(r1))
# r23c = (r23[:-1]+r23[1:])/(2*max(r23))
# r4c = (r4[:-1]+r4[1:])/(2*max(r4))
#
# plt.plot(r1c,E1,label='Start of stream tube',linestyle='--',marker='o')
# plt.plot(r23c,E2,label='Before rotor')
# plt.plot(r23c,E3,label='After rotor')
# plt.plot(r4c[:-1],E4[:-1],label='End of stream tube')
# plt.xlabel(r'$r/R$')
# plt.ylabel(r'$h_{s}$')
# plt.legend(loc='upper left')
# plt.grid()
# plt.show()

# fig = plt.figure()
# ax = fig.add_subplot(111, polar=True)
# theta = np.linspace(0,2*np.pi,1000)
# r = np.ones(len(theta))
# for i in range(len(r1)-1):
#     ax.plot(theta,r*r23[i])
# plt.show()





# plt.figure()
# # plt.plot(r_R[1:],Cax[1:],'r',label='Axial loading coefficient')
# clcd = cls/cds
# plt.plot(r_R[1:],clcd[1:]/117.50503018108652,'b')#,label='Azimuthal loading coefficient')
# plt.title('TSR = ' + str(TSR))
# plt.xlabel('Blade location r/R')
# plt.ylabel(r'$C_{l}/C_{d}$')
# # plt.legend(loc='upper right')
# plt.grid()
#
# plt.show()

'''PLOTS FOR VARYING TSR'''
alphas = np.zeros(n)
phis = np.zeros(n)
a = np.zeros(n)
ap = np.zeros(n)
Fax = np.zeros(n)
Faz = np.zeros(n)
gamma = np.zeros(n)
Cax = np.zeros(n)
Caz = np.zeros(n)
cls = np.zeros(n)
cds = np.zeros(n)

TSRar = np.arange(6,11)
CTar = np.zeros(len(TSRar))
CPar = np.zeros(len(TSRar))
CQar = np.zeros(len(TSRar))
for j in range(len(TSRar)):

    CT=0
    CP=0
    CQ=0
    TSR = TSRar[j]
    omega = U0*TSR/R
    for i in range(len(r_arr)):
        a[i],ap[i],Fax[i],Faz[i],gamma[i], phi, alpha, cl ,cd = solve_anul(U0,r_arr[i],dr_arr[i],c[i],b[i],omega=omega,Prandtl=True, Glauert=True)

        CT += dr_arr[i] * Fax[i] * nb/(0.5*rho*U0**2*np.pi*R**2)
        CP += dr_arr[i] * Faz[i] * nb * r_arr[i] * omega/(0.5*rho*U0**3*np.pi*R**2)
        CQ += dr_arr[i] * Faz[i] * nb /(0.5*rho*U0**2*np.pi*R**2)

        # print(i,'finished')
    CPar[j] = CP
    CTar[j] = CT
    CQar[j] = CQ
    print('Finished with CT =',CT,'and CP =',CP)
plt.plot(TSRar,CQar)
plt.xlabel(r'$\lambda$ [-]')
plt.ylabel(r'$C_{Q}$')
plt.show()
# plt.plot(TSRar,CPar)
# plt.xlabel(r'$\lambda$ [-]')
# plt.ylabel(r'$C_{P}$')
# plt.show()




