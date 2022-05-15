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

    return [u, v, w]
