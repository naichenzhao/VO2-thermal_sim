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
#  |             Calculation Utils                  |
#  +------------------------------------------------+

def calculate_dt(mat_d, mat_t, p, dt, r=1):
    curr_point = get_point(mat_d, p)
    adj_points = get_adj(mat_d, p)

    curr_t = get_point(mat_t, p)
    cp = get_cp(curr_point)
    p = get_p(curr_point)

    a = 1 - (dt/(p*cp))*(len(adj_points)/r**3)
    b = (dt/(p*cp))*(1/r**3)
    sum_t = sum([get_point(mat_t, i) for i in adj_points])
    new_t = curr_t*a + b*sum_t

    return new_t, set(adj_points)


def calculate_dt_l(mat_d, mat_t, p, dt):
    a_sum = 0
    b_sum = 0

    x = p[0]
    y = p[1]
    z = p[2]

    adj_points = get_adj(mat_d, p)
    curr_p = get_point(mat_d, p)
    curr_t = get_point(mat_t, p)
    
    # ---------- Calculate X points ----------
    A = curr_p[2] * curr_p[1] # dz * dy

    if x+1 < mat_d.shape[0]:  # Get point for (1, 0, 0)
        target_p = get_point(mat_d, (x+1, y, z), check=False)
        target_t = get_point(mat_t, (x+1, y, z), check=False)

        a_sum += 1/((curr_p[0] + target_p[0]) * A)
        b_sum += target_t/((curr_p[0] + target_p[0]) * A)
    
    if x-1 >= 0:  # Get point for (-1, 0, 0)
        target_p = get_point(mat_d, (x-1, y, z), check=False)
        target_t = get_point(mat_t, (x-1, y, z), check=False)

        a_sum += 1/((curr_p[0] + target_p[0]) * A)
        b_sum += target_t/((curr_p[0] + target_p[0]) * A)

    # ---------- Calculate Y points ----------
    A = curr_p[2] * curr_p[0] # dz * dx
    if y+1 < mat_d.shape[1]:  # Get point for (0, 1, 0)
        target_p = get_point(mat_d, (x, y+1, z), check=False)
        target_t = get_point(mat_t, (x, y+1, z), check=False)

        a_sum += 1/((curr_p[1] + target_p[1]) * A)
        b_sum += target_t/((curr_p[1] + target_p[1]) * A)

    if y-1 >= 0:  # Get point for (0, -1, 0)
        target_p = get_point(mat_d, (x, y-1, z), check=False)
        target_t = get_point(mat_t, (x, y-1, z), check=False)

        a_sum += 1/((curr_p[1] + target_p[1]) * A)
        b_sum += target_t/((curr_p[1] + target_p[1]) * A)

    # ---------- Calculate Z points ----------
    A = curr_p[0] * curr_p[1]  # dx * dy
    if z+1 < mat_d.shape[2]:  # Get point for (0, 0, 1)
        target_p = get_point(mat_d, (x, y, z+1), check=False)
        target_t = get_point(mat_t, (x, y, z+1), check=False)

        a_sum += 1/((curr_p[2] + target_p[2]) * A)
        b_sum += target_t/((curr_p[2] + target_p[2]) * A)
    
    if z-1 >= 0:  # Get point for (0, 0, -1)
        target_p = get_point(mat_d, (x, y, z-1), check=False)
        target_t = get_point(mat_t, (x, y, z-1), check=False)

        a_sum += 1/((curr_p[2] + target_p[2]))/A
        b_sum += target_t/((curr_p[2] + target_p[2]))/A

    a = curr_t * (1 - (dt/(get_cp(curr_p) * get_p(curr_p))) * a_sum/4)
    b = (dt/(get_cp(curr_p) * get_p(curr_p))) * b_sum/4
    return a + b, adj_points


#  +------------------------------------------------+
#  |             Helper Functions                   |
#  +------------------------------------------------+


def get_z_temp(mat_t, z=0):
    return mat_t[:, :, z]


def get_adj(mat_d, p):
    ret = set()

    x = p[0]
    y = p[1]
    z = p[2]

    if (x + 1 < mat_d.shape[0]):
        ret.add((x+1, y, z))
    if (y + 1 < mat_d.shape[1]):
        ret.add((x, y+1, z))
    if (z + 1 < mat_d.shape[2]):
        ret.add((x, y, z+1))
    
    if (x - 1 >= 0):
        ret.add((x-1, y, z))
    if (y - 1 >= 0):
        ret.add((x, y-1, z))
    if (z - 1 >= 0):
        ret.add((x, y, z-1))

    return ret


def point_exists(mat_d, p):
    x = p[0]
    y = p[1]
    z = p[2]

    x_check = not (z < 0 or z >= mat_d.shape[2])
    y_check = not (y < 0 or y >= mat_d.shape[1])
    z_check = not (x < 0 or x >= mat_d.shape[0])

    return x_check and y_check and z_check


def p_add(p1, p2):
    return (p1[0] + p2[0], p1[1] + p2[1], p1[2] + p2[2])


def p_sub(p1, p2):
    return (p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2])
