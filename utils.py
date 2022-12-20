import numpy as np


#  +------------------------------------------------+
#  |           Matrix Library Utils                 |
#  +------------------------------------------------+

def get_point(mat, p, check = True):
    if not point_exists(mat, p):
        return None
    return mat[p[0]][p[1]][p[2]]

def set_point(mat, p, new_val, check = True):
    if not point_exists(mat, p):
        return None
    mat[p[0]][p[1]][p[2]] = new_val

def make_point(dx=0.5, dy=0.5, dz=0.5, k=1, cp=1, p=1):
    return [dx, dy, dz, k, cp, p]

def get_dx(p):
    return p[0]

def get_dy(p):
    return p[1]

def get_dz(p):
    return p[2]

def get_k(p):
    return p[3]

def get_cp(p):
    return p[4]

def get_p(p):
    return p[5]


#  +------------------------------------------------+
#  |             Helper Functions                   |
#  +------------------------------------------------+


def get_z_temp(mat_t, z=0):
    return mat_t[:, :, z]

def point_exists(mat_d, p):
    x = p[0]
    y = p[1]
    z = p[2]

    x_check = not (z < 0 or z >= mat_d.shape[2])
    y_check = not (y < 0 or y >= mat_d.shape[1])
    z_check = not (x < 0 or x >= mat_d.shape[0])

    return x_check and y_check and z_check


