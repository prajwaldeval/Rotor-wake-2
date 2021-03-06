import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import pandas as pd
import numpy as np
import copy
from time import process_time


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


def aero_coeffs(alpha, filename='polar DU95W180.xlsx'):
    data = np.array(pd.read_excel(filename))[3:, :]
    data = np.array(data, dtype='float64')
    # if alpha > 30.06 or alpha < -16.0:
    #     print("AoA out of bounds", alpha)
    cl = np.interp(alpha, data[:, 0], data[:, 1])
    cd = np.interp(alpha, data[:, 0], data[:, 2])
    return cl, cd


def blade_geometry(r_R):
    # radial = [0, 0.3, .5, .8, 1]
    # chorddist = [.05, .04, .03, .02, .015]
    # twistdist = [-12, -8, -5, -4, 0.]
    pitch = 2
    chord = 3 * (1 - r_R) + 1
    twist = 14 * (1 - r_R)
    result = [chord, twist + pitch]
    return result


def wake_system_generation(r_array, dr, U0, a, wakelength, number_of_blades, tip_speed_ratio, increment):
    controlpoints = []
    filaments = []

    Ur = U0 * (1 - a)  # axial velocity just after rotor
    UD = U0 * (1 - 2 * a)  # axial velocity 'infinity' downwind of rotor

    D = (r_array[-1] + 0.5 * dr[-1]) * 2  # rotor diameter
    R = D/2

    n_rotations = (tip_speed_ratio * wakelength) / ((1 - a) * np.pi)
    theta_array = np.arange(0, n_rotations * 2 * np.pi, increment)

    angle_rotation_first_blade = 0
    cos_rotation = np.cos(angle_rotation_first_blade)
    sin_rotation = np.sin(angle_rotation_first_blade)

    for ri in range(len(r_array) - 1):
        geodef = blade_geometry(r_array[ri] / R)  # chord, twist+pitch=phi
        angle = geodef[1] * np.pi / 180

        """" Defining Control Points """
        temp1 = {"coordinates": [0, r_array[ri], 0],
                 "chord": geodef[0],
                 "normal": [np.cos(angle), 0, -1 * np.sin(angle)],
                 "tangential": [-1 * np.sin(angle), 0, -1 * np.cos(angle)],
                 "angle": 0,
                 "radius": r_array[ri],
                 "dr": dr[ri],
                 "y_offset": 0}

        controlpoints.append(temp1)

        """" Define Bound Vortex Filaments """
        temp1 = {"x1": 0,
                 "y1": r_array[ri] - 0.5 * dr[ri],
                 "z1": 0,
                 "x2": 0,
                 "y2": (r_array[ri + 1] - 0.5 * dr[ri + 1]),
                 "z2": 0,
                 "Gamma": 0,
                 "Blade": 0,
                 "Horse": ri}

        filaments.append(temp1)

        """" create trailing filaments, at x1 (1st point) of bound filament """
        geodef = blade_geometry(r_array[ri] / R)
        angle = geodef[1] * np.pi / 180

        temp1 = {"x1": geodef[0] * np.sin(-angle),
                 "y1": (r_array[ri] - 0.5 * dr[ri]) * cos_rotation,
                 "z1": (r_array[ri] - 0.5 * dr[ri]) * sin_rotation,
                 "x2": 0,
                 "y2": (r_array[ri] - 0.5 * dr[ri]) * cos_rotation,
                 "z2": (r_array[ri] - 0.5 * dr[ri]) * sin_rotation,
                 "Gamma": 0,
                 "Blade": 0,
                 "Horse": ri}

        filaments.append(temp1)

        for j in range(len(theta_array) - 1):
            xt = filaments[len(filaments) - 1]["x1"]
            yt = filaments[len(filaments) - 1]["y1"]
            zt = filaments[len(filaments) - 1]["z1"]

            dx = (theta_array[j + 1] - theta_array[j]) / tip_speed_ratio * D / 2
            dy = (np.cos(-theta_array[j + 1]) - np.cos(-theta_array[j])) * (r_array[ri] - 0.5 * dr[ri])
            dz = (np.sin(-theta_array[j + 1]) - np.sin(-theta_array[j])) * (r_array[ri] - 0.5 * dr[ri])

            temp1 = {"x1": xt + dx,
                     "y1": yt + dy,
                     "z1": zt + dz,
                     "x2": xt,
                     "y2": yt,
                     "z2": zt,
                     "Gamma": 0,
                     "Blade": 0,
                     "Horse": ri}

            filaments.append(temp1)

        """"" create trailing filaments, at x2 (1st point) of bound filament """""
        geodef = blade_geometry(r_array[ri] / R)
        angle = geodef[1] * np.pi / 180

        temp2 = {"x1": geodef[0] * np.sin(angle),
                 "y1": (r_array[ri] + 0.5 * dr[ri]) * cos_rotation,
                 "z1": (r_array[ri] + 0.5 * dr[ri]) * sin_rotation,
                 "x2": 0,
                 "y2": (r_array[ri] + 0.5 * dr[ri]) * cos_rotation,
                 "z2": (r_array[ri] + 0.5 * dr[ri]) * sin_rotation,
                 "Gamma": 0,
                 "Blade": 0,
                 "Horse": ri}

        filaments.append(temp2)

        for j in range(len(theta_array) - 1):
            xt = filaments[len(filaments) - 1]["x2"]
            yt = filaments[len(filaments) - 1]["y2"]
            zt = filaments[len(filaments) - 1]["z2"]

            dx = (theta_array[j + 1] - theta_array[j]) / tip_speed_ratio * D / 2
            dy = (np.cos(-theta_array[j + 1]) - np.cos(-theta_array[j])) * (r_array[ri + 1] - 0.5 * dr[ri + 1])
            dz = (np.sin(-theta_array[j + 1]) - np.sin(-theta_array[j])) * (r_array[ri + 1] - 0.5 * dr[ri + 1])

            temp1 = {"x1": xt,
                     "y1": yt,
                     "z1": zt,
                     "x2": xt + dx,
                     "y2": yt + dy,
                     "z2": zt + dz,
                     "Gamma": 0,
                     "Blade": 0,
                     "Horse": ri
                     }

            filaments.append(temp1)

    angle_rotation = ((2 * np.pi) / number_of_blades)

    fils_new_blades = []
    cps_new_blades = []

    for blade_nr in range(1, number_of_blades):
        theta = angle_rotation * blade_nr

        for i in range(len(controlpoints)):
            # rotation around X-axis
            temp1 = {"coordinates": [0,
                                     controlpoints[i]["coordinates"][1] * np.cos(theta) -
                                     controlpoints[i]["coordinates"][2] * np.sin(
                                         theta),
                                     controlpoints[i]["coordinates"][1] * np.sin(theta) +
                                     controlpoints[i]["coordinates"][2] * np.cos(
                                         theta)],

                     "chord": 0,

                     "normal": [controlpoints[i]["normal"][0],
                                controlpoints[i]["normal"][1] * np.cos(theta) - controlpoints[i]["normal"][2] * np.sin(
                                    theta),
                                controlpoints[i]["normal"][1] * np.sin(theta) + controlpoints[i]["normal"][2] * np.cos(
                                    theta)],

                     "tangential": [controlpoints[i]["tangential"][0],
                                    controlpoints[i]["tangential"][1] * np.cos(theta) - controlpoints[i]["tangential"][
                                        2] * np.sin(theta),
                                    controlpoints[i]["tangential"][1] * np.sin(theta) + controlpoints[i]["tangential"][
                                        2] * np.cos(theta)],
                     "angle": theta,

                     "radius": controlpoints[i]["radius"],

                     "dr": controlpoints[i]["dr"],

                     "y_offset": 0}

            cps_new_blades.append(temp1)

        for i in range(len(filaments)):
            temp_dict = {"x1": filaments[i]["x1"],
                         "y1": filaments[i]["y1"] * np.cos(theta) - filaments[i]["z1"] * np.sin(theta),
                         "z1": filaments[i]["y1"] * np.sin(theta) + filaments[i]["z1"] * np.cos(theta),
                         "x2": filaments[i]["x2"],
                         "y2": filaments[i]["y2"] * np.cos(theta) - filaments[i]["z2"] * np.sin(theta),
                         "z2": filaments[i]["y2"] * np.sin(theta) + filaments[i]["z2"] * np.cos(theta),
                         "Gamma": 0,
                         "Blade": blade_nr,
                         "Horse": filaments[i]["Horse"]
                         }

            fils_new_blades.append(temp_dict)

    cps_all = controlpoints + cps_new_blades
    fils_all = filaments + fils_new_blades

    fils_dict = {"x1": [],
                 "y1": [],
                 "z1": [],
                 "x2": [],
                 "y2": [],
                 "z2": [],
                 "Gamma": [],
                 "Blade": [],
                 "Horse": []}

    for i in range(len(fils_all)):
        fils_dict["x1"].append(fils_all[i]["x1"])
        fils_dict["y1"].append(fils_all[i]["y1"])
        fils_dict["z1"].append(fils_all[i]["z1"])
        fils_dict["x2"].append(fils_all[i]["x2"])
        fils_dict["y2"].append(fils_all[i]["y2"])
        fils_dict["z2"].append(fils_all[i]["z2"])
        fils_dict["Gamma"].append(fils_all[i]["Gamma"])
        fils_dict["Blade"].append(fils_all[i]["Blade"])
        fils_dict["Horse"].append(fils_all[i]["Horse"])
    fils_dict["x1"] = np.asarray(fils_dict["x1"])
    fils_dict["y1"] = np.asarray(fils_dict["y1"])
    fils_dict["z1"] = np.asarray(fils_dict["z1"])
    fils_dict["x2"] = np.asarray(fils_dict["x2"])
    fils_dict["y2"] = np.asarray(fils_dict["y2"])
    fils_dict["z2"] = np.asarray(fils_dict["z2"])
    fils_dict["Gamma"] = np.asarray(fils_dict["Gamma"])
    fils_dict["Blade"] = np.asarray(fils_dict["Blade"])
    fils_dict["Horse"] = np.asarray(fils_dict["Horse"])

    return cps_all, fils_dict


def biot_savart_function(fil_x1, fil_y1, fil_z1,
                         fil_x2, fil_y2, fil_z2,
                         cp_x, cp_y, cp_z,
                         gamma, core):
    R1 = np.sqrt((cp_x - fil_x1) ** 2 + (cp_y - fil_y1) ** 2 + (cp_z - fil_z1) ** 2)
    R2 = np.sqrt((cp_x - fil_x2) ** 2 + (cp_y - fil_y2) ** 2 + (cp_z - fil_z2) ** 2)

    R1_2X = (cp_y - fil_y1) * (cp_z - fil_z2) - (cp_z - fil_z1) * (cp_y - fil_y2)
    R1_2Y = -(cp_x - fil_x1) * (cp_z - fil_z2) + (cp_z - fil_z1) * (cp_x - fil_x2)
    R1_2Z = (cp_x - fil_x1) * (cp_y - fil_y2) - (cp_y - fil_y1) * (cp_x - fil_x2)

    R12_sqr = R1_2X ** 2 + R1_2Y ** 2 + R1_2Z ** 2

    R01 = (fil_x2 - fil_x1) * (cp_x - fil_x1) + (fil_y2 - fil_y1) * (cp_y - fil_y1) + (fil_z2 - fil_z1) * (
            cp_z - fil_z1)
    R02 = (fil_x2 - fil_x1) * (cp_x - fil_x2) + (fil_y2 - fil_y1) * (cp_y - fil_y2) + (fil_z2 - fil_z1) * (
            cp_z - fil_z2)

    if R12_sqr < core ** 2:
        R12_sqr = core ** 2

    if R1 < core:
        R1 = core

    if R2 < core:
        R2 = core

    K = gamma / (4 * np.pi * R12_sqr) * (R01 / R1 - R02 / R2)

    u = K * R1_2X
    v = K * R1_2Y
    w = K * R1_2Z

    return u, v, w


def unit_strength_induction_matrix(cps, fils, n, number_of_blades):
    """""
     set up unit strength induction matrix
     & initalize and calculate matrices for velocity induced by horseshoe vortex rings
    """""

    # setting up unit induction matrix with the biot-savart function

    unitU_ind = np.zeros((len(cps), (n - 1) * number_of_blades))
    unitV_ind = np.zeros((len(cps), (n - 1) * number_of_blades))
    unitW_ind = np.zeros((len(cps), (n - 1) * number_of_blades))



    for i_cp in range(len(cps)):
        x_cp, y_cp, z_cp = cps[i_cp]["coordinates"]

        j = 0
        for i_blade in range(number_of_blades):
            for i_horse in range(n - 1):

                x1s = fils['x1'][np.logical_and(fils['Blade'] == i_blade, fils['Horse'] == i_horse)]
                x2s = fils['x2'][np.logical_and(fils['Blade'] == i_blade, fils['Horse'] == i_horse)]
                y1s = fils['y1'][np.logical_and(fils['Blade'] == i_blade, fils['Horse'] == i_horse)]
                y2s = fils['y2'][np.logical_and(fils['Blade'] == i_blade, fils['Horse'] == i_horse)]
                z1s = fils['z1'][np.logical_and(fils['Blade'] == i_blade, fils['Horse'] == i_horse)]
                z2s = fils['z2'][np.logical_and(fils['Blade'] == i_blade, fils['Horse'] == i_horse)]

                for i_fil in range(len(x1s)):
                    u, v, w = biot_savart_function(x1s[i_fil], y1s[i_fil], z1s[i_fil],
                                                   x2s[i_fil], y2s[i_fil], z2s[i_fil],
                                                   x_cp, y_cp, z_cp, gamma=1, core=0.5)

                    unitU_ind[i_cp, j] = unitU_ind[i_cp, j] + u
                    unitV_ind[i_cp, j] = unitV_ind[i_cp, j] + v
                    unitW_ind[i_cp, j] = unitW_ind[i_cp, j] + w

                j = j + 1

    return unitU_ind, unitV_ind, unitW_ind


def BE_loads(Vax, Vtan, beta, c, rho):
    Vps = Vtan ** 2 + Vax ** 2
    phi = np.arctan(Vax / Vtan)
    alpha = phi * 180 / np.pi - beta  # could be minus beta, to check
    cl, cd = aero_coeffs(alpha)
    L = 0.5 * c * rho * Vps * cl
    D = 0.5 * c * rho * Vps * cd
    Faz = L * np.sin(phi) - D * np.cos(phi)
    Fax = L * np.cos(phi) + D * np.sin(phi)
    gamma = 0.5 * np.sqrt(Vps) * cl * c  # vorticity
    phi = phi * 180 / np.pi

    return Fax, Faz, gamma, phi, alpha, cl, cd


def iteration(iterations, Ua, Va, Wa, cps, tsr, gamma_convergence_weight, error_limit, R, U0, rho):
    gamma = np.ones(len(cps))
    gamma_new = np.ones(len(cps))
    a_new = np.empty(len(cps))
    Fax_ll = np.empty(len(cps))
    Faz_ll = np.empty(len(cps))
    alpha = np.empty(len(cps))
    phi = np.empty(len(cps))
    for i in range(iterations):

        u_ind = Ua.dot(gamma)
        v_ind = Va.dot(gamma)
        w_ind = Wa.dot(gamma)

        # print(u_ind)

        for i_cp in range(len(cps)):
            # position and geometry at current control point
            # radius = np.sqrt(cps[i_cp]["coordinates"][1] ** 2 + cps[i_cp]["coordinates"][2] ** 2)
            radius = cps[i_cp]['radius']
            yoffset = cps[i_cp]['y_offset']
            geodef = blade_geometry(radius / R)
            c = geodef[0]
            twist_and_pitch = geodef[1]
            #
            # print("y_offset", yoffset)

            # velocities at current control point
            omega = tsr * U0 * 1 / R

            coord_vect = [cps[i_cp]["coordinates"][0], cps[i_cp]["coordinates"][1] - yoffset, cps[i_cp]["coordinates"][2]]

            Vrot = np.cross([-omega, 0, 0], coord_vect)

            # axial and tangential
            vel_perceived = [U0 + u_ind[i_cp] + Vrot[0], v_ind[i_cp] + Vrot[1], w_ind[i_cp] + Vrot[2]]
            Vax = vel_perceived[0]

            azimdir = np.cross([-1 / radius, 0, 0], coord_vect)
            Vtan = np.dot(azimdir, vel_perceived)

            # angle = cps[i_cp]["angle"]
            # Vtan = -v_ind[i_cp] * np.sin(angle) - w_ind[i_cp] * np.cos(angle) - omega * radius

            # send to BEM model function
            # print("beta is ", twist_and_pitch, "deg")
            # print("axial velocity", Vax, "Vtan", Vtan)
            Fax_ll[i_cp], Faz_ll[i_cp], gamma_new[i_cp], phi[i_cp], alpha[i_cp], cl, cd = BE_loads(Vax, Vtan, twist_and_pitch, c, rho)

            # update axial induction factor
            a_new[i_cp] = 1 - Vax / U0

        # update circulation
        gamma = (1 - gamma_convergence_weight) * gamma + gamma_convergence_weight * gamma_new

        # # check convergence
        # ref_error = max(np.abs(gamma_new))
        # ref_error = max(ref_error, 0.001)  # define scale of bound circulation

        error = (np.absolute(gamma_new - gamma)).max()  # difference betweeen iterations
        # error = error / ref_error  # relative error
        # print("we are at iteration", i)
        if error < error_limit:
            print("convergence threshold met at iteration", i)
            break

    return a_new, gamma, Fax_ll, Faz_ll, alpha, phi


def create_second_rotor_wake(cps, fils, y_offset, phase_difference):
    theta = phase_difference * np.pi / 180
    cps_second = copy.deepcopy(cps)
    for i in range(len(cps_second)):
        # rotation around X-axis
        cps_second[i] = {"coordinates": [0,
                                         cps_second[i]["coordinates"][1] * np.cos(theta) -
                                         cps_second[i]["coordinates"][2] * np.sin(theta),
                                         cps_second[i]["coordinates"][1] * np.sin(theta) +
                                         cps_second[i]["coordinates"][2] * np.cos(theta)],

                         "chord": 0,

                         "normal": [cps_second[i]["normal"][0],
                                    cps_second[i]["normal"][1] * np.cos(theta) - cps_second[i]["normal"][2] * np.sin(
                                        theta),
                                    cps_second[i]["normal"][1] * np.sin(theta) + cps[i]["normal"][2] * np.cos(theta)],

                         "tangential": [cps_second[i]["tangential"][0],
                                        cps_second[i]["tangential"][1] * np.cos(theta) - cps_second[i]["tangential"][
                                            2] * np.sin(theta),
                                        cps_second[i]["tangential"][1] * np.sin(theta) + cps_second[i]["tangential"][
                                            2] * np.cos(theta)],
                         "angle": theta,

                         "radius": cps_second[i]["radius"],

                         "y_offset": y_offset}

        # translation in y-direction with the offset value
        cps_second[i]["coordinates"][1] = cps_second[i]["coordinates"][1] + y_offset

        global fils_second
        fils_second = copy.deepcopy(fils)

        for j in range(len(fils_second["x1"])):
            fils_second["y1"][j] = fils["y1"][j] * np.cos(theta) - fils["z1"][j] * np.sin(theta) + y_offset
            fils_second["z1"][j] = fils["y1"][j] * np.sin(theta) + fils["z1"][j] * np.cos(theta)
            fils_second["y2"][j] = fils["y2"][j] * np.cos(theta) - fils["z2"][j] * np.sin(theta) + y_offset
            fils_second["z2"][j] = fils["y2"][j] * np.sin(theta) + fils["z2"][j] * np.cos(theta)
            fils_second["Blade"][j] = fils_second["Blade"][j] #+ 3



    return cps_second, fils_second

def coefficients(Fax, Faz, r, dr, nb, rho, U0, TSR):
    CT = 0
    CP = 0
    R = r[-1] + 0.5 * dr[-1]
    omega = TSR * U0 / R
    for i in range(len(r) - 1):
        CT = CT + (dr[i] * Fax[i] * nb) / (0.5 * rho * U0 ** 2 * np.pi * R **2)
        CP = CP + (dr[i] * Faz[i] * r[i] * omega * nb) / (0.5 * rho * U0 ** 3 * np.pi * R ** 2)
    return CT, CP

def coefficients1(Fax, Faz, R, rho, U0, TSR, cps):
    CT = 0
    CP = 0
    omega = TSR * U0 / R
    for i in range(len(cps)):
        CT = CT + (cps[i]["dr"] * Fax[i]) / (0.5 * rho * U0 ** 2 * np.pi * R **2)
        CP = CP + (cps[i]["dr"] * Faz[i] * cps[i]["radius"] * omega) / (0.5 * rho * U0 ** 3 * np.pi * R ** 2)
    return CT, CP

if __name__ == '__main__':
    number_of_blades = 3
    R = 50  # m
    r_hub = 0.2 * R
    U0 = 10
    a = 0.25
    nr_blade_elements = 25
    rho = 1.225
    tip_speed_ratio = 8
    iterations = 200
    gamma_convergence_weight = 0.3
    error_limit = 0.001
    wake_length = 2.5  # how many diameters long the wake shall be prescribed for
    nt = 50
    increment = np.pi / 30

    # Constant or cosine element spacing on the blade
    r, dr = geometry_constant(r_hub, R, nr_blade_elements)
    # r, dr = geometry_cosine(r_hub, R, nr_blade_elements)

    cps, fils = wake_system_generation(r, dr, U0, a, wake_length, number_of_blades, tip_speed_ratio, increment)

    Ua, Va, Wa = unit_strength_induction_matrix(cps, fils, nr_blade_elements, number_of_blades)

    a_new_single, gamma_single, Fax_ll_single, Faz_ll_single, alpha_ll_single, phi_ll_single = iteration(iterations, Ua,
                                                                                                         Va, Wa,
                                                                                                         cps,
                                                                                                         tip_speed_ratio,
                                                                                                         gamma_convergence_weight,
                                                                                                         error_limit, R,
                                                                                                         U0, rho)

    CT_single, CP_single = coefficients(Fax_ll_single, Faz_ll_single, r, dr, number_of_blades, rho, U0, tip_speed_ratio)
    Cax_ll_single = Fax_ll_single * number_of_blades / (0.5 * rho * U0 ** 2 * R)
    Caz_ll_single = Faz_ll_single * number_of_blades / (0.5 * rho * U0 ** 2 * R)

    # 2nd rotor generation
    offset = 1
    y_offset = offset * 2 * R
    phase_diff = 60

    cps_second, fils_second = create_second_rotor_wake(cps, fils, y_offset, phase_diff)

    U00, V00, W00 = unit_strength_induction_matrix(cps, fils, nr_blade_elements, number_of_blades)
    U01, V01, W01 = unit_strength_induction_matrix(cps, fils_second, nr_blade_elements, number_of_blades)
    U10, V10, W10 = unit_strength_induction_matrix(cps_second, fils, nr_blade_elements, number_of_blades)
    U11, V11, W11 = unit_strength_induction_matrix(cps_second, fils_second, nr_blade_elements, number_of_blades)
    Ua = np.block([
        [U00, U01],
        [U10, U11]
    ])

    Va = np.block([
        [V00, V01],
        [V10, V11]
    ])

    Wa = np.block([
        [W00, W01],
        [W10, W11]
    ])

    cp_total = cps + cps_second
    a_new, gamma, Fax_ll, Faz_ll, alpha_ll, phi_ll = iteration(iterations, Ua, Va, Wa, cp_total, tip_speed_ratio,
                                                               gamma_convergence_weight, error_limit, R, U0, rho)
    CT, CP = coefficients1(Fax_ll, Faz_ll, R, rho, U0, tip_speed_ratio, cps)
    Cax_ll = Fax_ll * number_of_blades / (0.5 * rho * U0 ** 2 * R)
    Caz_ll = Faz_ll * number_of_blades / (0.5 * rho * U0 ** 2 * R)

    # # Plotting of wake system
    # fig = plt.figure()
    # ax = plt.axes(projection="3d")
    # for i in range(len(fils["x1"])):
    #     x, y, z = [fils["x1"][i], fils["x2"][i]], [fils["y1"][i], fils["y2"][i]], [fils["z1"][i], fils["z2"][i]]
    #     ax.plot(x, y, z, color='black')
    #
    # for i in range(len(fils_second["x1"])):
    #     x, y, z = [fils_second["x1"][i], fils_second["x2"][i]], \
    #               [fils_second["y1"][i], fils_second["y2"][i]], \
    #               [fils_second["z1"][i], fils_second["z2"][i]]
    #     ax.plot(x, y, z, color='black')
    #
    # for i in range(len(cps_second)):
    #     x, y, z = cps_second[i]["coordinates"]
    #     ax.scatter(x, y, z, c='red', s=50)
    #
    # for i in range(len(cps)):
    #     x, y, z = cps[i]["coordinates"]
    #     ax.scatter(x, y, z, c='blue', s=50)

    # plt.show()

    # Plotting of results
    r_flip = r / R

    fig1 = plt.figure()
    plt.plot(r_flip[:-1], alpha_ll_single, label='single rotor')
    plt.plot(r_flip[:-1], alpha_ll[0:len(r_flip)-1], label='2 rotors, blade 1')
    plt.plot(r_flip[:-1], alpha_ll[len(r_flip):2*len(r_flip)-1], label='2 rotors, blade 2')
    plt.plot(r_flip[:-1], alpha_ll[2*len(r_flip):3*len(r_flip)-1], label='2 rotors, blade 3')
    plt.xlabel(r'Angle of attack')
    plt.ylabel(r'r/R')
    plt.legend()
    plt.show()
    #
    # fig2 = plt.figure()
    # plt.plot(offsets, CPs, label='Two rotor cases')
    # plt.xlabel(r'offsets in multiple of diameter')
    # plt.ylabel(r'$C_{P}$')
    # plt.axhline(y=0.333,label='Single rotor case')
    # plt.legend()
    # plt.show()





    # time_lst = []
    # wake_lst = []
    # elem_lst = []
    #
    # for tip_speed_ratio in range(6, 10, 2):
    #
    #     t1_start = process_time()
    #
    #     r, dr = geometry_constant(r_hub, R, nr_blade_elements)
    #
    #     cps, fils = wake_system_generation(r, dr, U0, a, wake_length, number_of_blades, tip_speed_ratio, increment)
    #
    #     Ua, Va, Wa = unit_strength_induction_matrix(cps, fils, nr_blade_elements, number_of_blades)
    #
    #     a_new, gamma, Fax_ll, Faz_ll, alpha_ll, phi_ll = iteration(iterations, Ua, Va, Wa,
    #                                              cps, tip_speed_ratio,
    #                                              gamma_convergence_weight, error_limit, R, U0, rho)
    #
    #     CT, CP = coefficients(Fax_ll, Faz_ll, r, dr, number_of_blades, rho, U0, tip_speed_ratio)
    #
    #     t1_stop = process_time()
    #
    #     print("Elapsed time for wakelength of", wake_length, "blade sections", nr_blade_elements, "was",
    #           t1_stop - t1_start)
    #
    #     r_flip = r/R
    #     norm = 0.5 * rho * U0 ** 2 * R
    #
    #     fig1 = plt.figure()
    #     plt.plot(r_flip[0:(nr_blade_elements - 1)], a_new[0:(nr_blade_elements - 1)])
    #     plt.title("Axial Induction Factor 'a'")
    #     plt.show()
    #
    #     fig2 = plt.figure()
    #     plt.plot(r_flip[0:(nr_blade_elements - 1)], Fax_ll[0:(nr_blade_elements - 1)]/norm)
    #     plt.title("Axial Normalised Force 'F_ax'")
    #     plt.show()
    #
    #     fig3 = plt.figure()
    #     plt.plot(r_flip[0:(nr_blade_elements - 1)], Faz_ll[0:(nr_blade_elements - 1)]/norm)
    #     plt.title("Azimuthal Normalised Force 'F_az'")
    #     plt.show()