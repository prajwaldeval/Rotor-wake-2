import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# TSR = 10
n = 100 #number of discretised blade elements

R = 50 #m
r_hub = 0.2*R
pitch = -2 #deg
U0 = 10 #m/s
# omega= U0*TSR/R #rad/s
nb = 3 #number of blades
rho = 1.225

#constant for Glauert
CT1 = 1.816
CT2 = 2*np.sqrt(CT1) - CT1

def chord(mu):
    return 3*(1-mu) +1

def twist(mu):
    return 14*(1-mu) + pitch

def geometry_cosine(n):
    r = np.linspace(0,np.pi,n+1)
    r = (1-np.cos(r))/2
    blen = R-r_hub
    r = r_hub + blen*r
    dr = r[1:]-r[:-1]
    r = (r[1:]+r[:-1])/2
    return r, dr

def geometry_constant(n):
    r = np.linspace(r_hub, R, n+1)
    dr = r[1:] - r[:-1]
    r = (r[1:] + r[:-1]) / 2
    return  r,dr


def aero_coeffs(alpha,filename = 'polar DU95W180.xlsx'):
    data = np.array(pd.read_excel(filename))[3:,:]
    data = np.array(data,dtype='float64')
    cl = np.interp(alpha,data[:,0],data[:,1])
    cd = np.interp(alpha,data[:,0],data[:,2])
    return cl, cd

def BE_loads(a,ap,r,dr,b,c,omega):
    Vtan = omega*r*(1+ap)
    Vax = U0*(1-a) #1-a
    Vps= Vtan**2+Vax**2
    phi = np.arctan2(Vax,Vtan)
    alpha = phi*180/np.pi - b
    cl, cd = aero_coeffs(alpha)
    L = 0.5*c*rho*Vps*cl
    D = 0.5*c*rho*Vps*cd
    Faz = L*np.sin(phi) - D*np.cos(phi)
    Fax = L*np.cos(phi) + D*np.sin(phi)
    gamma = 0.5 * np.sqrt(Vps) * cl * c  # vorticity
    phi = phi*180/np.pi
    return Vax,Vtan,Fax,Faz,gamma,phi,alpha,cl,cd

def MT_induction(Fax,Faz,r,dr,b,c,Prandtl,Glauert,omega):
    CT = (Fax*nb*dr)/(0.5*rho*(U0**2)*2*np.pi*r*dr)

    a = 0.5 - 0.5*np.sqrt(1-CT)
    # print(a)
    if Glauert:
        a1 = 1 - np.sqrt(CT1)/2
        if a < a1:
            CT = CT1 - 4*(np.sqrt(CT1)-1)*(1-a)
            if CT<CT2:
                a = 0.5 - 0.5 * np.sqrt(1 - CT)
            else:
                print('Using Glauert', a, '!!!!')
                a = 1 + (CT - CT1)/(4*np.sqrt(CT1)-4)
                print(a)


    ap = (Faz * nb) / (2 * rho * (2 * np.pi * r) * U0**2 * (1 - a) * r * omega)  # 1-a
    if Prandtl:
        mu = r/R
        l = omega*R/U0
        p1 = -nb/2*(1-mu)/(mu)*np.sqrt(1+(l**2*mu**2)/(1-a)**2)
        ftip = 2/np.pi * np.arccos(np.exp(p1))
        p2 = -nb/2*(mu-r_hub/R)/(mu)*np.sqrt(1+(l**2*mu**2)/(1-a)**2)
        froot = 2/np.pi * np.arccos(np.exp(p2))
        ftot = ftip * froot
        if (ftot < 0.0001):
            ftot = 0.0001
        a = a/ftot
        ap = ap/ftot

    return a,ap

def solve_anul(Uinf, r, dr, c, b,Prandtl ,Glauert,omega):

    #initialization of variables
    a = 0.2   # axial induction factor
    ap = 0  # tangential induction factor

    convergence_error = 0.0001
    max_nr_iterations = 100

    for i in range(max_nr_iterations):

        Vax, Vtan, Fax, Faz, gamma,phi,alpha,cl,cd = BE_loads(a,ap,r,dr,b,c,omega)
        an, apn = MT_induction(Fax,Faz,r,dr,b,c, Prandtl, Glauert,omega)

        anext = 0.25*a + 0.75*an
        apnext = 0.25*ap + 0.75*apn

        # anext = an
        # apnext = apn

        if (np.abs(a-anext)<convergence_error):
            # print('N =', i,'converged with a =', a, 'and ap =',ap, 'and r=', r/R)
            break
        if i==max_nr_iterations-1:
            # print('Not converged and r=', r/R)
            a = anext
            ap = apnext

    return a,ap,Fax,Faz,gamma, phi, alpha, cl, cd

def areas(a,r,dr):
    r23 = np.zeros(len(r)+1)
    r23[:-1]=r-0.5*dr
    r23[-1]=r[-1]+0.5*dr[-1]
    r1 = np.zeros(len(r23))
    r4 = np.zeros(len(r23))
    r1[0] = r23[0]
    r4[0] = r23[0]
    a23 = np.zeros(len(r))
    a1 = np.zeros(len(r))
    a4 = np.zeros(len(r))

    for i in range(1,len(r)+1):
        r1[i] = np.sqrt( (1 - a[i - 1]) * (r23[i] ** 2 - r23[i - 1] ** 2) + r1[i - 1] ** 2 )
        r4[i] = np.sqrt( ((1 - a[i - 1]) * (r23[i] ** 2 - r23[i - 1] ** 2))/(1 - 2 * a[i-1]) + r4[i - 1] ** 2)
        print('r1=',r1[i],', r23=', r23[i], ', r4=', r4[i])
        a1[i-1] = np.pi*(r1[i]**2-r1[i-1]**2)
        a23[i-1] = np.pi*(r23[i]**2-r23[i-1]**2)
        a4[i - 1] = np.pi * (r4[i] ** 2 - r4[i - 1] ** 2)
    return r1,r23,r4,a1,a23,a4

def enthalpy_jump(a23,Fax):
    f = Fax/a23
    jump = f/rho
    return -1*jump






