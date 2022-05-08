import matplotlib.pyplot as plt
import numpy as np

def geometry_cosine(r_hub,R,n):
    r = np.linspace(0,np.pi,n+1)
    r = (1-np.cos(r))/2
    blen = R-r_hub
    r = r_hub + blen*r
    dr = r[1:]-r[:-1]
    r = (r[1:]+r[:-1])/2
    return r, dr

def geometry_constant(r_hub,R,n):
    r = np.linspace(r_hub, R, n+1)
    dr = r[1:] - r[:-1]
    r = (r[1:] + r[:-1]) / 2
    return  r,dr

def wake_gen(r,dr,U0,w,a,wakelength,nt):
    Ur = U0*(1-a) #axial velocity just after rotor
    UD = U0*(1-2*a) #axial velocity 'infinity' downwind of rotor
    Uw = np.linspace(Ur,UD,nt) #axial velocity at every discretised location in the wake
    D = (r[-1]+0.5*dr[-1])*2 #rotor diameter
    zw = np.linspace(0,wakelength*D,nt) #z(axial) coordinate of wake discretised points
    tw = np.zeros(nt)
    for i in range(1,nt):
        tw[i] = tw[i-1] + wakelength*D/Uw[i]/nt

    n = len(r) #number of horseshoe vortices per blade
    rw = np.zeros((nt,n)) #each row is one timestep, each column spanwise station
    thw = np.zeros((nt,n))
    print(rw)
    # for i in range(n):

    return 0

if __name__ == '__main__':
    R = 50  # m
    r_hub = 0.2 * R
    U0 = 10
    a = 0.25
    r,dr = geometry_constant(r_hub,R,2)
    # r,dr = geometry_cosine(r_hub,R,30)

    wake = wake_gen(r,dr,U0,a,4,50,10)