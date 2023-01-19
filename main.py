import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt
from utils import *
from state_matrix import *


# TODO:
#   - Test sample transient/steady state solutions


#  +------------------------------------------------------------------+
#  |                                                                  |
#  |                       Constant Values                            |
#  |                                                                  |
#  +------------------------------------------------------------------+

X_GRID = 300
Y_GRID = 100
Z_GRID = 1



#  +------------------------------------------------------------------+
#  |                                                                  |
#  |                    Set Point Functions                           |
#  |                                                                  |
#  +------------------------------------------------------------------+
''' 
Finite Difference for 3D Heat Transfer

Temperature Matrix: mat_t 
    3D matrix of point temperatures
    

Data matrix: mat_d - [dx, dy, dz, k, cp, p]
    dx (0) = x-distance
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
    dt = 0.1
    NUM_CYCLES = 100000

    # Create primary matrices to use
    mat_d = np.zeros((X_GRID, Y_GRID, Z_GRID, 6))
    mat_t = np.zeros((X_GRID, Y_GRID, Z_GRID))
    mat_h = np.zeros((X_GRID, Y_GRID))

    # Set default values for data matrix
    for i in range(X_GRID):
        for j in range(Y_GRID):
            for k in range(Z_GRID):
                set_point(mat_d, (i, j, k), make_point())
    
    
    
    
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
    h_state = gen_h_state(mat_d, dt)


    # Set values for heat transfer matrix
    # set_heat_mat(mat_h, 100000, ((50, 40), (50, 60)))


    # set heated matrix elements
    z_heat = 0
    
    x_0_plane = [(0, y, z_heat) for y in range(Y_GRID)]
    x_n_plane = [(X_GRID-1, y, z_heat) for y in range(Y_GRID)]
    y_n_plane = [(x, Y_GRID-1, z_heat) for x in range(1, X_GRID-1)]
    POINTS = [*x_0_plane, *x_n_plane, *y_n_plane]
    mask = set_mat(mat_t, 100, POINTS)

    y_0_plane = [(x, 0, z_heat) for x in range(1, X_GRID-1)]
    POINTS2 = [*y_0_plane]
    mask2 = set_mat(mat_t, 0, POINTS2)

    mask = mask + mask2


    # Print the initial Conditions
    print("Printing initial conditions")
    print_mat(mat_t)
    
    
    
    
    
    #  +-------------------------------------------+
    #  |              Run Loop                     |
    #  +-------------------------------------------+
    print("Running loop... ")

    # Create loop helpers
    new_temps = np.zeros((X_GRID, Y_GRID, Z_GRID))
    comp_mat = np.zeros((X_GRID, Y_GRID, Z_GRID, 6))
    initial_temps = mat_t.copy()
    
    
    
    
    # Run loop
    p_bar = tqdm(range(NUM_CYCLES), desc="Running Sim")
    for i in p_bar:
        new_temps = transition_state(mat_t, comp_mat, a_state, b_state).copy()
        mat_t = set_mat_temps(mask, initial_temps, new_temps)
        apply_heat(mat_t, h_state, mat_h)



    #  +-------------------------------------------+
    #  |           Print Values                    |
    #  +-------------------------------------------+
    print("Printing final result... ")

    print("-----------------------------")

    print_mat(mat_t)





#  +------------------------------------------------------------------+
#  |                                                                  |
#  |                       Helper Functions                           |
#  |                                                                  |
#  +------------------------------------------------------------------+

def print_plane(planes):
    C = 'gist_heat'
    MAXT = 115

    num_args = len(planes)
    x = int(np.sqrt(num_args))
    y = int(np.ceil(num_args/x))

    if x ==1 and y == 1:
        p = plt.imshow(np.transpose(planes[0]), extent=[0, X_GRID,
                                      0, Y_GRID], cmap='gist_heat', vmax=115)
        plt.colorbar(p)
        plt.show()
    else:
        fig, axes = plt.subplots(nrows=x, ncols=y)
        i = 0;
        for ax in axes.flat:
            pl = np.transpose(planes[i])
            i += 1
            im = ax.imshow(pl, extent=[0, X_GRID, 0, Y_GRID], cmap=C, vmax=MAXT)
            ax.set_title(f'Layer {i}', fontsize=8)
        fig.colorbar(im, ax=axes.ravel().tolist())

    plt.show()

def print_mat(mat_t):
    grid1 = get_z_temp(mat_t)
    # grid2 = get_z_temp(mat_t, 1)
    # grid3 = get_z_temp(mat_t, 2)
    # grid4 = get_z_temp(mat_t, 3)
    print_plane(np.array([grid1]))


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


#  +------------------------------------------------------------------+
#  |                                                                  |
#  |                    Set Point Functions                           |
#  |                                                                  |
#  +------------------------------------------------------------------+


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


if __name__ == '__main__':
    main()
