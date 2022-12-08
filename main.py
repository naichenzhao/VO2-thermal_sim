import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt
from utils import *
from state_matrix import *


# TODO:
#   - try to get faster computation (maybe using np.pad)
#   - make a way of just calculating all points - make use of caching 



X_GRID = 100
Y_GRID = 100
Z_GRID = 4



#  +-------------------------------------------+
#  |        Data point Organization            |
#  +-------------------------------------------+
''' [dx, dy, dz, k, cp, p]
    dx (0) = d-distance
    dy (1) = y-distance
    dz (2) = z-distance

    k (3) = thermal conductivity
    cp (4) = heat capacity
    p (5) = density (rho)
'''


def main():
    #  +-------------------------------------------+
    #  |             Setup Defaults                |
    #  +-------------------------------------------+
    print("Setting up Constants... ")
    dt = 0.2
    NUM_CYCLES = 15000

    # Matrix for storing the point data
    mat_d = np.empty((X_GRID, Y_GRID, Z_GRID, 6))

    # Set default values for data matrix
    for i in range(X_GRID):
        for j in range(Y_GRID):
            for k in range(Z_GRID):
                set_point(mat_d, (i, j, k), make_point())


    # Matrix for storing the point temperatures
    mat_t = np.empty((X_GRID, Y_GRID, Z_GRID))

    # Set default values for temperature matrix
    for i in range(X_GRID):
        for j in range(Y_GRID):
            for k in range(Z_GRID):
                set_point(mat_t, (i, j, k), 0)




    #  +-------------------------------------------+
    #  |           Setup Matrix                    |
    #  +-------------------------------------------+
    print("Setting up Matrix... ")

    # Check for min timestep/ stab ility
    min_time = get_min_timestep(mat_d)
    if dt > min_time:
        print("defined time is too large for stability")
        print("new timestep of", min_time, "has been selected")
        dt = min_time

    # Calculate state transition matrix
    a_state, b_state = gen_state_matrix(mat_d, dt)

    # set heated elements
    y_n_plane = [(x, Y_GRID-1, 0) for x in range(X_GRID)]
    x_0_plane = [(0, y, 0) for y in range(Y_GRID)]
    x_n_plane = [(X_GRID-1, y, 0) for y in range(Y_GRID)]

    POINTS = [*y_n_plane, *x_0_plane, *x_n_plane]
    set_mat(mat_t, 100, POINTS)


    print("Printing initial conditions")
    xy_grid = np.transpose(get_z_temp(mat_t))
    print_plane(xy_grid)
    



    #  +-------------------------------------------+
    #  |              Run Loop                     |
    #  +-------------------------------------------+
    print("Running loop... ")

    # Create loop helpers
    new_temps = np.empty((X_GRID, Y_GRID, Z_GRID))
    comp_mat = np.empty((X_GRID, Y_GRID, Z_GRID, 6))
    initial_temps = mat_t.copy()
    
    

    # Run loop
    p_bar = tqdm(range(NUM_CYCLES), desc="Running Sim")
    for i in p_bar:
        new_temps = transition_state(mat_t, comp_mat, a_state, b_state)
        mat_t = set_mat_temps(mat_t, new_temps, initial_temps)





    #  +-------------------------------------------+
    #  |           Print Values                    |
    #  +-------------------------------------------+
    print("Printing final result... ")

    print("-----------------------------")

    xy_grid = np.transpose(get_z_temp(mat_t))
    print_plane(xy_grid)







    
#  +------------------------------------------------+
#  |             Helper Functions                   |
#  +------------------------------------------------+

def print_plane(plane):
    p = plt.imshow(plane, extent=[0, X_GRID,
                   0, Y_GRID], cmap='gist_heat', vmax=115)
    plt.colorbar(p)
    plt.show()


def get_min_timestep(mat):
    min_t = 100000
    for i in range(X_GRID):
        for j in range(Y_GRID):
            for k in range(Z_GRID):
                curr_p = get_point(mat, (i, j, k))
                a = get_cp(curr_p)*get_p(curr_p)
                b = 3*get_dx(curr_p)*get_dy(curr_p)*get_dz(curr_p)
                curr_t = 0.9 * a/b
                min_t = min(curr_t, min_t)
    return min_t


#  +------------------------------------------------+
#  |            Set Point Functions                 |
#  +------------------------------------------------+

def set_points(mat_t, temps, coordinates):
    r_set = set()

    for i in range(len(temps)):
        curr_loc = coordinates[i]
        curr_temp = temps[i]
        set_point(mat_t, curr_loc, curr_temp)


def set_mat(t_mat, temp, coordinates):
    
    for i in range(len(coordinates)):
        curr_loc = coordinates[i]
        set_point(t_mat, curr_loc, temp)


def set_points_r(mat_t, temps, coordinates):
    r_set = set()

    for i in range(len(temps)):
        curr_loc = coordinates[i]
        curr_temp = temps[i]
        set_point(mat_t, curr_loc, curr_temp)

        r_set.add(curr_loc)
        r_set.update(get_adj(mat_t, curr_loc))
    return r_set


def set_mat_r(t_mat, temp, coordinates):
    r_set = set()
    for i in range(len(coordinates)):
        # Set the temperature values
        curr_loc = coordinates[i]
        set_point(t_mat, curr_loc, temp)

        r_set.add(curr_loc)
    return r_set



#  +------------------------------------------------+
#  |            Set Point Functions                 |
#  +------------------------------------------------+

def sim_priority(mat, mat_t, vals, dt):
    r_set = set()

    for p in vals:
        t, points = calculate_dt_l(mat, p, dt)
        set_point(mat_t, p, t)
        r_set.update(points)

    for p in vals:
        set_point(get_point(mat, p), get_point(mat_t, p))

    return r_set



if __name__ == '__main__':
    main()

    

