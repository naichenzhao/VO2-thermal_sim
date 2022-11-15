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


def calculate_dt_l(mat, p, dt):
    a_sum = 0
    b_sum = 0
    adj_points = get_adj(mat, p)
    curr_p = get_point(mat, p)
    
    # ---------- Calculate X points ----------
    print(p_add(p, (1, 0, 0)))
    if point_exists(mat, p_add(p, (1, 0, 0))):  # Get point for (1, 0, 0)
        target = get_point(mat, p_add(p, (1, 0, 0)))
        A = get_dz(curr_p) * get_dy(curr_p) * 4
        a_sum += (1/A) * (1/(get_dx(curr_p) + get_dx(target)))
        b_sum += (1/A) * (get_t(target)/(get_dx(curr_p) + get_dx(target)))
        print(get_dx(curr_p))
    
    if point_exists(mat, p_add(p, (-1, 0, 0))):  # Get point for (-1, 0, 0)
        target = get_point(mat, p_add(p, (-1, 0, 0)))
        A = get_dz(curr_p) * get_dy(curr_p) * 4
        a_sum += (1/A) * (1/(get_dx(curr_p) + get_dx(target)))
        b_sum += (1/A) * (get_t(target)/(get_dx(curr_p) + get_dx(target)))

    # ---------- Calculate Y points ----------
    if point_exists(mat, p_add(p, (0, 1, 0))):  # Get point for (0, 1, 0)
        target = get_point(mat, p_add(p, (0, 1, 0)))
        A = get_dz(curr_p) * get_dx(curr_p) * 4
        a_sum += (1/A) * (1/(get_dy(curr_p) + get_dy(target)))
        b_sum += (1/A) * (get_t(target)/(get_dy(curr_p) + get_dy(target)))

    if point_exists(mat, p_add(p, (0, -1, 0))):  # Get point for (0, -1, 0)
        target = get_point(mat, p_add(p, (0, -1, 0)))
        A = get_dz(curr_p) * get_dx(curr_p) * 4
        a_sum += (1/A) * (1/(get_dy(curr_p) + get_dy(target)))
        b_sum += (1/A) * (get_t(target)/(get_dy(curr_p) + get_dy(target)))

    # ---------- Calculate Z points ----------
    if point_exists(mat, p_add(p, (0, 0, 1))):  # Get point for (0, 0, 1)
        target = get_point(mat, p_add(p, (0, 0, 1)))
        A = get_dx(curr_p) * get_dy(curr_p) * 4
        a_sum += (1/A) * (1/(get_dz(curr_p) + get_dz(target)))
        b_sum += (1/A) * (get_t(target)/(get_dz(curr_p) + get_dz(target)))
    
    if point_exists(mat, p_add(p, (0, 0, -1))):  # Get point for (0, 0, -1)
        target = get_point(mat, p_add(p, (0, 0, -1)))
        A = get_dx(curr_p) * get_dy(curr_p) * 4
        a_sum += (1/A) * (1/(get_dz(curr_p) + get_dz(target)))
        b_sum += (1/A) * (get_t(target)/(get_dz(curr_p) + get_dz(target)))

    a = get_t(curr_p) * (1 - (dt/(get_cp(curr_p) * get_p(curr_p))) * a_sum)
    b = (dt/(get_cp(curr_p) * get_p(curr_p))) * b_sum

    new_t =  a + b

    return new_t, adj_points

    


def get_mat_temp(mat):
    return mat[:,:,:,0]


def get_adj(mat, p):
    ret = []
    mask = [(1,0,0), (0,1,0), (0,0,1)]

    for i in mask:
        if point_exists(mat, p_add(p, i)):
            ret.append(p_add(p, i))
        if point_exists(mat, p_sub(p, i)):
            ret.append(p_sub(p, i))
    return ret


def point_exists(mat, p):
    x = p[0]
    y = p[1]
    z = p[2]

    x_check = not (z < 0 or z >= mat.shape[2])
    y_check = not (y < 0 or y >= mat.shape[1])
    z_check = not (x < 0 or x >= mat.shape[0])

    return x_check and y_check and z_check


def p_add(p1, p2):
    return (p1[0] + p2[0], p1[1] + p2[1], p1[2] + p2[2])


def p_sub(p1, p2):
    return (p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2])
