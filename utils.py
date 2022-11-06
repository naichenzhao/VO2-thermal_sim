import numpy as np
from itertools import product

def get_point(mat, p):
    if not point_exists(mat, p):
        return None
    return mat[p[0]][p[1]][p[2]]


def set_point(mat, p, new_p):
    if not point_exists(mat, p):
        return
    mat[p[0]][p[1]][p[2]] = new_p


def make_point(t=0, dx=0.5, dy=0.5, dz=0.5, k=1, cp=1, p=4, s=False):
    return [t, dx, dy, dz, s, k, cp, p]


def set_t(p, t):
    if p[4] == 0:
        p[0] = t


def get_t(p):
    return p[0]

def get_dx(p):
    return p[1]

def get_dy(p):
    return p[2]

def get_dz(p):
    return p[3]

def get_s(p):
    return p[4]

def get_k(p):
    return p[5]

def get_cp(p):
    return p[6]

def get_p(p):
    return p[7]


def calculate_dt(mat, p, dt, r):
    curr_point = get_point(mat, p)
    adj_points = get_adj(mat, p)

    curr_t = get_t(curr_point)
    cp = get_cp(curr_point)
    p=get_p(curr_point)

    a = 1 - (dt/(p*cp))*(len(adj_points)/r**3)
    b = (dt/(p*cp))*(1/r**3)
    sum_t = sum([get_t(get_point(mat, i)) for i in adj_points])
    new_t = curr_t*a + b*sum_t

    return new_t, adj_points


def get_mat_temp(mat):
    x_dim = mat.shape[0]
    y_dim = mat.shape[1]
    z_dim = mat.shape[2]

    t_mat = np.zeros((x_dim, y_dim, z_dim))
    space = np.array([*product(range(x_dim), range(y_dim), range(z_dim))])

    for i in space:
        set_point(t_mat, i, get_t(get_point(mat, i)))
    
    return t_mat


def get_adj(mat, p):
    ret = []
    x = p[0]
    y = p[1]
    z = p[2]

    if point_exists(mat, (x+1, y, z)):
        ret.append((x+1, y, z))
    if point_exists(mat, (x-1, y, z)):
        ret.append((x-1, y, z))
    if point_exists(mat, (x, y+1, z)):
        ret.append((x, y+1, z))
    if point_exists(mat, (x, y-1, z)):
        ret.append((x, y-1, z))
    if point_exists(mat, (x, y, z+1)):
        ret.append((x, y, z+1))
    if point_exists(mat, (x, y, z-1)):
        ret.append((x, y, z-1))
    return ret


def point_exists(mat, p):
    x = p[0]
    y = p[1]
    z = p[2]

    x_check = not (z < 0 or z >= mat.shape[2])
    y_check = not (y < 0 or y >= mat.shape[1])
    z_check = not (x < 0 or x >= mat.shape[0])

    return x_check and y_check and z_check



