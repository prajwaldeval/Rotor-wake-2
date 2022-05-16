import numpy as np


def biot_savart_function(fil_x1, fil_y1, fil_z1,
                         fil_x2, fil_y2, fil_z2,
                         cp_x, cp_y, cp_z,
                         gamma):

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

    K = gamma / (4 * np.pi * R12_sqr) * (R01 / R1 - R02 / R2)

    u = K * R1_2X
    v = K * R1_2Y
    w = K * R1_2Z

    return u, v, w


def unit_strength_induction_matrix():
    """""
     set up unit strength induction matrix
     & initalize and calculate matrices for velocity induced by horseshoe vortex rings       
    """""

    # setting up unit induction matrix with the biot-savart function

    unitU_ind = []
    unitV_ind = []
    unitW_ind = []

    for ii in range(len(controlpoints)):
        for nn in range(len(controlpoints)):
            for jj = 1:2 * N_trail:
                for nb in range(nr_of_blades)

                    fil_x1 =
                    fil_y1 =
                    fil_z1 =
                    fil_x2 =
                    fil_y2 =
                    fil_z2 =
                    cp_x =
                    cp_y =
                    cp_z =
                    gamma = 1
                    uu, vv, ww = biot_savart_function(fil_x1, fil_y1, fil_z1,
                                                        fil_x2, fil_y2, fil_z2,
                                                        cp_x, cp_y, cp_z,
                                                        gamma)

    unitU_ind(ii, nn) = unitU_ind(ii, nn) + uu
    unitV_ind(ii, nn) = unitV_ind(ii, nn) + vv
    unitW_ind(ii, nn) = unitW_ind(ii, nn) + ww

    return [unitU_ind, unitV_ind, unitW_ind]