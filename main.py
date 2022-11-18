import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt
from utils import *

    
#  +-------------------------------------------+
#  |           Constant values                 |
#  +-------------------------------------------+
 
X_GRID = 51
Y_GRID = 51
Z_GRID = 1
r = 1


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
    print("Setting up matrix... ")

    mat = np.empty((X_GRID, Y_GRID, Z_GRID, dp_length))

    for i in range(X_GRID):
        for j in range(Y_GRID):
            for k in range(Z_GRID):
                set_point(mat, (i, j, k), make_point())

    min_time = get_min_timestep(mat)
    
    TEMPS = (100, 100, 100, 100)
    POINTS = ((0, 0, 0), (50, 50, 0), (3, 0, 0), (25, 25, 0))

    next_points = set_mat(mat, TEMPS, POINTS)


    #  +-------------------------------------------+
    #  |              Run Loop                     |
    #  +-------------------------------------------+
    print("Running loop... ")
    dt = 0.25
    NUM_CYCLES = 400

    if dt > min_time:
        print("defined time is too large for stability")
        print("new timestep of", min_time, "has been selected")
        dt = min_time



    for i in tqdm(range(NUM_CYCLES)):
        next_points = sim(mat, next_points, dt)
    


    #  +-------------------------------------------+
    #  |           Print Values                    |
    #  +-------------------------------------------+
    print("Printing... ")

    xy_grid = np.transpose(get_z_temp(mat))

    # print(xy_grid)

    p = plt.imshow(xy_grid, extent=[0,X_GRID,0,Y_GRID], cmap='gist_heat', vmax = 120)
    plt.colorbar(p)
    plt.show()



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
    r_set = set()
    new_temps = {}

    for p in vals:
        t, points = calculate_dt_l(mat, p, dt)
        new_temps[p] = t
        r_set.update(points)

    for p in new_temps.keys():
        set_t(get_point(mat, p), new_temps[p])

    return r_set



def get_min_timestep(mat):
    min_t = 100000;
    for i in range(X_GRID):
        for j in range(Y_GRID):
            for k in range(Z_GRID):
                curr_p = get_point(mat, (i, j, k))
                a = get_cp(curr_p)*get_p(curr_p)
                b = 3*get_dx(curr_p)*get_dy(curr_p)*get_dz(curr_p)
                curr_t = 0.9 * a/b
                min_t = min(curr_t, min_t)
    return min_t


if __name__ == '__main__':
    main()

    

