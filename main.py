import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt
from utils import *

    
#  +-------------------------------------------+
#  |           Constant values                 |
#  +-------------------------------------------+
 
X_GRID = 11
Y_GRID = 1
Z_GRID = 1
r = 1
dt = 1


#  +-------------------------------------------+
#  |        Data point Organization            |
#  +-------------------------------------------+
'''
[t, dx, dy, dz, s, k, cp, p]

    t (0) = temperature
    dx (1) = d-distance
    dy (2) = y-distance
    dz (3) = z-distance

    s (4) = if node is static (1 or 0)

    k (5) = thermal conductivity
    cp (6) = heat capacity
    p (7) = density (rho)
'''

dp_length = 8




def main():
    #  +-------------------------------------------+
    #  |           Setup Matrix                    |
    #  +-------------------------------------------+

    mat = np.empty((X_GRID, Y_GRID, Z_GRID, dp_length))

    for i in range(X_GRID):
        for j in range(Y_GRID):
            for k in range(Z_GRID):
                set_point(mat, (i, j, k), make_point())
    
    next_points = set_mat(mat, (100, 100), ((0, 0, 0), (10, 0, 0)))


    #  +-------------------------------------------+
    #  |              Run Loop                     |
    #  +-------------------------------------------+
    for i in tqdm(range(3)):
        next_points = sim(mat, next_points, dt)

    print(get_mat_temp(mat))



def set_mat(mat, temps, coordinates):
    r_set = set()

    for i in range(len(temps)):
        curr_loc = coordinates[i]
        curr_temp = temps[i]

        set_point(mat, curr_loc, make_point(t=curr_temp, s = 1))

        r_set.add(curr_loc)
        r_set.update(get_adj(mat, curr_loc))
    
    return r_set


def sim(mat, vals, dt):
    ret_set = set()
    temp_mat = np.zeros((mat.shape[0], mat.shape[1], mat.shape[2]))

    for p in vals:
        t, points = calculate_dt_l(mat, p, dt)
        set_point(temp_mat, p, t)
        ret_set.update(points)

    for p in vals:
        set_t(get_point(mat, p), get_point(temp_mat, p))

    return ret_set




if __name__ == '__main__':
    main()

    

