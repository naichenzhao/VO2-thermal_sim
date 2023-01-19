import numpy as np


#  +------------------------------------------------------------------+
#  |                                                                  |
#  |                        FDM Functions                             |
#  |                                                                  |
#  +------------------------------------------------------------------+


#  +------------------------------------------------+
#  |               Set Point Functions              |
#  +------------------------------------------------+

def get_min_timestep(mat):
    min_t = 100000
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            for k in range(mat.shape[2]):
                curr_p = get_point(mat, (i, j, k))
                a = get_cp(curr_p)*get_p(curr_p)
                b = 3*get_dx(curr_p)*get_dy(curr_p)*get_dz(curr_p)
                curr_t = 0.9 * a/b
                min_t = min(curr_t, min_t)
    return min_t



#  +------------------------------------------------+
#  |               Set Point Functions              |
#  +------------------------------------------------+

def set_points(t_mat, temps, coordinates):
    mask = np.zeros(t_mat.shape)

    for i in range(len(temps)):
        curr_loc = coordinates[i]
        curr_temp = temps[i]
        set_point(t_mat, curr_loc, curr_temp)
        set_point(mask, curr_loc, 1)

    return mask


def set_mat(t_mat, temp, coordinates):
    mask = np.zeros(t_mat.shape)

    for i in range(len(coordinates)):
        curr_loc = coordinates[i]
        set_point(t_mat, curr_loc, temp)
        set_point(mask, curr_loc, 1)

    return mask


def set_heat_mat(h_mat, temp, coordinates):
    for curr in coordinates:
        h_mat[curr[0], curr[1]] = temp







#  +------------------------------------------------------------------+
#  |                                                                  |
#  |                      Library Functions                           |
#  |                                                                  |
#  +------------------------------------------------------------------+

#  +------------------------------------------------+
#  |           Matrix Library Utils                 |
#  +------------------------------------------------+

def get_point(mat, p, check = True):
    if not point_exists(mat, p):
        return None
    return mat[p[0]][p[1]][p[2]]


def set_point(mat, p, new_val, check=True):
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
#  |           Heat Transfer Equations              |
#  +------------------------------------------------+


def apply_heat(mat_t, h_state, mat_h):
    mat_t[:, :, 0] += mat_h * h_state



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


