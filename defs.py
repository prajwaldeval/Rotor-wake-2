import matplotlib.pyplot as plt
import numpy as np


def geometry_cosine(r_hub, R, n):
    r = np.linspace(0, np.pi, n + 1)
    r = (1 - np.cos(r)) / 2
    blen = R - r_hub
    r = r_hub + blen * r
    dr = r[1:] - r[:-1]
    r = (r[1:] + r[:-1]) / 2
    return r, dr


def geometry_constant(r_hub, R, n):
    r = np.linspace(r_hub, R, n + 1)
    dr = r[1:] - r[:-1]
    r = (r[1:] + r[:-1]) / 2
    return r, dr

def blade_geometry(radius):
    radial = [0, 0.3, .5, .8, 1]
    chorddist = [.05, .04, .03, .02, .015]
    twistdist = [-12, -8, -5, -4, 0.]
    pitch = 2
    chord = 3 * (1 - radius) + 1
    twist = -14 * (1 - radius)
    result = [chord, twist + pitch]
    return result

def wake_system_generation(r_array, dr, U0, w, a, wakelength, nt, number_of_blades):

    controlpoints = []
    filaments = []

    Ur = U0 * (1 - a)  # axial velocity just after rotor
    UD = U0 * (1 - 2 * a)  # axial velocity 'infinity' downwind of rotor

    # U_wake = np.linspace(Ur, UD, nt)  # axial velocity at every discretised location in the wake
    U_wake = Ur  # at least initially, as indicated by the tutorial

    D = (r[-1] + 0.5 * dr[-1]) * 2  # rotor diameter
    z_wake = np.linspace(0, wakelength * D, nt)  # z(axial) coordinate of wake discretised points
    t_wake = np.zeros(nt)

    for blade_nr in range(number_of_blades):

        angle_rotation = (2 * np.pi) / (blade_nr * number_of_blades)
        cos_rotation = np.cos(angle_rotation)
        sin_rotation = np.sin(angle_rotation)

        for r in r_array:

            geodef = blade_geometry(r)
            angle = geodef[1] * np.pi/180

            # define control points
            temp1 = {"coordinates": [0, r, 0],
                     "chord": geodef[0],
                     "normal": [np.cos(angle), 0, -1 * np.sin(angle)],
                     "tangential": [-1 * np.sin(angle), 0, -1 * np.cos(angle)]}

            # rotate blade to position
            temp1.coordinates = [0, temp1.coordinates[1] * cos_rotation - temp1.coordinates[2] * sin_rotation,
                                 temp1.coordinates[1] * sin_rotation + temp1.coordinates[2] * cos_rotation]

            temp1.normal = [temp1.normal[0], temp1.normal[1] * cos_rotation - temp1.normal[2] * sin_rotation,
                            temp1.normal[1] * sin_rotation + temp1.normal[2] * cos_rotation]

            temp1.tangential = [temp1.tangential[0], temp1.tangential[1] * cos_rotation - temp1.tangential[2] * sin_rotation,
                                temp1.tangential[1] * sin_rotation + temp1.tangential[2] * cos_rotation]

            controlpoints.append(temp1)

            # define bound vortex filament
            temp1 = {"x1": 0,
                     "y1": r_array[r],
                     "z1": 0,
                     "x2": 0,
                     "y2": r_array[r + 1],
                     "z2": 0,
                     "Gamma": 0}

            # rotate filament to position
            filaments.append(temp1)

            # create trailing filaments, at x1 of bound filament
            geodef = blade_geometry(r_array[r] / r)
            angle = geodef[1] * np.pi / 180

            temp1 = {"x1": geodef[0] * np.sin(-1 * angle),
                     "y1": r_array[r],
                     "z1": -1 * geodef[0] * np.cos(angle),
                     "x2": 0,
                     "y2": r_array[r],
                     "z2": 0,
                     "Gamma": 0}

            filaments.append(temp1)

            # define the trailing vortex filaments
            for j in range(1, nt):
                t_wake[j] = t_wake[j - 1] + wakelength * D / U_wake[j] / nt

            n = len(r)  # number of horseshoe vortices per blade
            rw = np.zeros((nt, n))  # each row is one timestep, each column spanwise station
            thw = np.zeros((nt, n))
            print(rw)
            # for i in range(n):

    return 0
