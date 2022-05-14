import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3Dimport
from mpl_toolkits import mplot3d
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
    # radial = [0, 0.3, .5, .8, 1]
    # chorddist = [.05, .04, .03, .02, .015]
    # twistdist = [-12, -8, -5, -4, 0.]
    pitch = 2
    chord = 3 * (1 - radius) + 1
    twist = -14 * (1 - radius)
    result = [chord, twist + pitch]
    return result


def wake_system_generation(r_array, dr, U0, a, wakelength, number_of_blades, tip_speed_ratio):
    controlpoints = []
    filaments = []

    Ur = U0 * (1 - a)  # axial velocity just after rotor
    UD = U0 * (1 - 2 * a)  # axial velocity 'infinity' downwind of rotor

    # U_wake = np.linspace(Ur, UD, nt)  # axial velocity at every discretised location in the wake
    U_wake = Ur  # at least initially, as indicated by the tutorial

    D = (r_array[-1] + 0.5 * dr[-1]) * 2  # rotor diameter
    # z_wake = np.linspace(0, wakelength * D, nt)  # z(axial) coordinate of wake discretised points
    # t_wake = np.zeros(nt)

    n_rotations = (tip_speed_ratio * np.pi * wakelength) / (1 - a)
    theta_array = np.arange(0, n_rotations * 2 * np.pi, np.pi / 10)

    for blade_nr in range(number_of_blades):

        angle_rotation = ((2 * np.pi) / number_of_blades) * blade_nr
        print('rot=', angle_rotation)
        cos_rotation = np.cos(angle_rotation)
        sin_rotation = np.sin(angle_rotation)

        # apply rotation to whole theta_array
        # theta_array = theta_array - angle_rotation

        for ri in range(len(r_array)-1):
            geodef = blade_geometry(r_array[ri])  # chord, twist+pitch=phi
            angle = geodef[1] * np.pi / 180

            # define control points
            temp1 = {"coordinates": [0, r_array[ri], 0],
                     "chord": geodef[0],
                     "normal": [np.cos(angle), 0, -1 * np.sin(angle)],
                     "tangential": [-1 * np.sin(angle), 0, -1 * np.cos(angle)],
                     "blade": blade_nr}

            # rotate blade to position
            temp1["coordinates"] = [0,
                                    temp1["coordinates"][1] * cos_rotation - temp1["coordinates"][2] * sin_rotation,
                                    temp1["coordinates"][1] * sin_rotation + temp1["coordinates"][2] * cos_rotation]

            temp1["normal"] = [temp1["normal"][0],
                               temp1["normal"][1] * cos_rotation - temp1["normal"][2] * sin_rotation,
                               temp1["normal"][1] * sin_rotation + temp1["normal"][2] * cos_rotation]

            temp1["tangential"] = [temp1["tangential"][0],
                                   temp1["tangential"][1] * cos_rotation - temp1["tangential"][2] * sin_rotation,
                                   temp1["tangential"][1] * sin_rotation + temp1["tangential"][2] * cos_rotation]

            controlpoints.append(temp1)

            # define bound vortex filament & rotate filament to position
            temp1 = {"x1": 0,
                     "y1": r_array[ri] * cos_rotation,
                     "z1": r_array[ri] * sin_rotation,
                     "x2": 0,
                     "y2": r_array[ri + 1] * cos_rotation,
                     "z2": r_array[ri + 1] * sin_rotation,
                     "Gamma": 0,
                     "blade": blade_nr}

            filaments.append(temp1)

            # create trailing filaments, at x1 (1st point) of bound filament
            geodef = blade_geometry(r_array[ri] / R)
            angle = geodef[1] * np.pi / 180

            temp1 = {"x1": geodef[0] * np.sin(-1 * angle),
                     "y1": r_array[ri],
                     "z1": -1 * geodef[0] * np.cos(angle),
                     "x2": 0,
                     "y2": r_array[ri],
                     "z2": 0,
                     "Gamma": 0,
                     "blade": blade_nr}

            filaments.append(temp1)

            for j in range(len(theta_array) - 1):
                xt = filaments[len(filaments) - 1]["x1"]
                yt = filaments[len(filaments) - 1]["y1"]
                zt = filaments[len(filaments) - 1]["z1"]

                dx = (theta_array[j + 1] - theta_array[j]) / tip_speed_ratio * D / 2
                dy = (np.cos(theta_array[j + 1]) - np.cos(theta_array[j])) * r_array[ri]
                dz = (np.sin(-theta_array[j + 1]) - np.sin(-theta_array[j])) * r_array[ri]


                temp1 = {"x1": xt + dx,
                         "y1": yt + dy,
                         "z1": zt + dz,
                         "x2": xt,
                         "y2": yt,
                         "z2": zt,
                         "Gamma": 0,
                         "blade": blade_nr}

                filaments.append(temp1)

            # create trailing filaments, at x2 (2nd point) of bound filament
            geodef = blade_geometry(r_array[ri + 1] / R)
            angle = geodef[1] * np.pi / 180

            temp1 = {"x1": geodef[0] * np.sin(-1 * angle),
                     "y1": r_array[ri + 1],
                     "z1": -1 * geodef[0] * np.cos(angle),
                     "x2": 0,
                     "y2": r_array[ri + 1],
                     "z2": 0,
                     "Gamma": 0,
                     "blade": blade_nr}

            filaments.append(temp1)
            for j in range(len(theta_array) - 1):
                xt = filaments[len(filaments) - 1]["x2"]
                yt = filaments[len(filaments) - 1]["y2"]
                zt = filaments[len(filaments) - 1]["z2"]

                dx = (theta_array[j + 1] - theta_array[j]) / tip_speed_ratio * D / 2
                dy = (np.cos(theta_array[j + 1]) - np.cos(theta_array[j])) * r_array[ri + 1]
                dz = (np.sin(-theta_array[j + 1]) - np.sin(-theta_array[j])) * r_array[ri + 1]


                temp1 = {"x1": xt,
                         "y1": yt,
                         "z1": zt,
                         "x2": xt + dx,
                         "y2": yt + dy,
                         "z2": zt + dz,
                         "Gamma": 0,
                         "blade": blade_nr}

                filaments.append(temp1)

    return controlpoints, filaments


if __name__ == '__main__':
    number_of_blades = 3
    R = 50  # m
    r_hub = 0.2 * R
    U0 = 10
    a = 0.25

    # Constant or cosine element spacing on the blade
    r, dr = geometry_constant(r_hub, R, 20)
    # r,dr = geometry_cosine(r_hub,R,30)

    wakelength = 0.05  # how many diameters long the wake shall be prescribed for
    nt = 50
    tip_speed_ratio = 8

    cps, fils = wake_system_generation(r, dr, U0, a, wakelength, number_of_blades, tip_speed_ratio)

    fig = plt.figure()
    ax = plt.axes(projection="3d")

    for i in range(len(fils)):
        x, y, z = [fils[i]["x1"], fils[i]["x2"]], [fils[i]["y1"], fils[i]["y2"]], [fils[i]["z1"], fils[i]["z2"]]
        ax.plot(x, y, z, color='black')

    for i in range(len(cps)):

        x, y, z = cps[i]["coordinates"]
        ax.scatter(x, y, z, c='red', s=100)

    plt.show()
