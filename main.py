from PySpice import *
from PySpice.Spice.Netlist import Circuit
import PySpice.Logging.Logging as Logging

import numpy as np
from tqdm import tqdm
import pandas as pd
import os
from matplotlib import pyplot as plt
from thermo_utils import *
from electro_utils import *
from state_matrix import *



#  +------------------------------------------------------------------+
#  |                                                                  |
#  |                       Constant Values                            |
#  |                                                                  |
#  +------------------------------------------------------------------+


# Number of runs
NUM_CYCLES = 1000  # 1 second of time

# Grid dimensionns
X_GRID = 250
Y_GRID = 250
Z_GRID = 50

# Electrostatic Dimensions
STARTX = 84
STARTY = 117

X_ESIM = 82
Y_ESIM = 16
Z_ESIM = 2

CONTACT_LENGTH = 16
SCALE = 2
VOLTAGE = 30


LASER_POWER = 0.05/50





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
    #  |              Setup Matrix                 |
    #  +-------------------------------------------+
    print("Setting up Matrix... ")
    dt = 1 # Set constant of 1 and allow program to calculate bounds

    # Create primary matrices to use
    mat_d = np.zeros((X_GRID, Y_GRID, Z_GRID, 6))
    mat_t = np.ones((X_GRID, Y_GRID, Z_GRID)) * (273.15 + 20) # Set to 20C
    h_add = np.zeros((X_GRID, Y_GRID))

    # Set default values for data matrix
    mat_d[:, :, :] = make_point()
    # for i in range(X_GRID):
    #     for j in range(Y_GRID):
    #         for k in range(Z_GRID):
    #             set_point(mat_d, (i, j, k), make_point())

    
    
    #  +-------------------------------------------+
    #  |           Setup Thermo                    |
    #  +-------------------------------------------+
    print("Setting up Thermo... ")
    
    '''Check for min timestep/ stability'''
    min_time = get_min_timestep(mat_d)
    if dt > min_time:
        print("\n=====================================")
        print("defined time:", dt, "is too large for stability")
        print("new timestep of", min_time, "has been selected")
        print("=====================================\n")
        dt = min_time


    '''Calculate state transition matrix'''
    a_state, b_state = gen_state_matrix(mat_d, dt)
    hstate_t = get_hstate_thermo(mat_d, dt)
    mask = None

    ''' Set values for heat transfer matrix '''
    laser_points = [[i, 125] for i in range(100, 151)]
    set_added_heat(h_add, LASER_POWER, laser_points)

    
    '''Set up boundry Temperatures'''
    mask = np.zeros((X_GRID, Y_GRID, Z_GRID))

    # set edges to 0
    mask[0, :, :] = 1
    mask[X_GRID-1, :, :] = 1
    mask[:, 0, :] = 1
    mask[:, Y_GRID-1, :] = 1
    mask[:, :, Z_GRID-1] = 1

    if mask is None:
        mask = np.zeros((X_GRID, Y_GRID, Z_GRID))


    #  +-------------------------------------------+
    #  |           Setup Electrostatics            |
    #  +-------------------------------------------+
    print("Setting up Electrostatics... ")

    # check if scaling works
    if X_ESIM % SCALE != 0 or Y_ESIM % SCALE != 0 or Z_ESIM % SCALE != 0 or CONTACT_LENGTH % SCALE != 0:
        print("\n=====================================")
        print("[ERROR]: INVALID SCALE VAUE")
        print("Thermo dimensions dont work, try changing simulation area")
        print("=====================================\n")
        quit()

    # Make node matrix
    nodes = np.zeros((X_ESIM//SCALE, Y_ESIM//SCALE,
                     Z_ESIM//SCALE, 3), dtype=np.int8)
    for i in range(X_ESIM//SCALE):
        for j in range(Y_ESIM//SCALE):
            for k in range(Z_ESIM//SCALE):
                nodes[i, j, k, 0] = i
                nodes[i, j, k, 1] = j
                nodes[i, j, k, 2] = k
    hstate_elec = get_hstate_elec(get_selected_area(mat_d, STARTX, STARTY, X_ESIM, Y_ESIM, Z_ESIM), dt, SCALE)
    
    # Make reistor matrix
    mat_r = get_res_matrix(mat_t, mat_d[0, 0, 0, 0], STARTX, STARTY, X_ESIM, Y_ESIM, Z_ESIM, SCALE)
    
    # setup resistors
    circuit = Circuit('sim')  # Remake the circuit
    circuit.V('input', 'vin', circuit.gnd, VOLTAGE)
    circuit.Vinput.minus.add_current_probe
    setup_resistors(circuit, mat_r, nodes, CONTACT_LENGTH//SCALE)


    #Print the initial Conditions
    print("Printing initial conditions")
    print("    Printing Thermal Conditions")
    print_mat_t(mat_t)
    print("    Printing Electrostatic Conditions")
    print_mat_e(mat_t, STARTX, STARTY, X_ESIM, Y_ESIM, CONTACT_LENGTH)





    #  +-----------------------------------------------------+
    #  |                                                     |
    #  |                       Run Loop                      |
    #  |                                                     |
    #  +-----------------------------------------------------+
    print("Running loop... ")

    # Create loop helpers
    new_temps = np.zeros((X_GRID, Y_GRID, Z_GRID))
    initial_temps = mat_t.copy()

    comp_mat = np.zeros((X_GRID, Y_GRID, Z_GRID, 6))
    res_heat = np.zeros((X_ESIM//SCALE, Y_ESIM//SCALE, Z_ESIM//SCALE))
    dv_mat = np.zeros((X_ESIM//SCALE, Y_ESIM//SCALE, Z_ESIM//SCALE, 6))
    
    power_draw = 0

    mat_v = np.zeros((X_ESIM//SCALE, Y_ESIM//SCALE, Z_ESIM//SCALE))
    
    # Run loop
    bar = tqdm(range(NUM_CYCLES), desc="Running Sim")
    for i in bar:
        # -------------------------------
        #   Thermal Sim
        # -------------------------------
        transition_state(mat_t, comp_mat, a_state, b_state, new_temps)
        mat_t = set_mat_temps(mask, initial_temps, new_temps)
        apply_heat(mat_t, hstate_t, h_add)


        
        # # -------------------------------
        # #   Electrostatic Sim
        # # -------------------------------
        # mat_r = get_res_matrix(mat_t, mat_d[0,0,0,0], STARTX, STARTY, X_ESIM, Y_ESIM, Z_ESIM, SCALE)

        # '''For every N cycles, reset spice to make sure we dotn use too much memory'''
        # N = 1000
        # if i>0 and i%N == 0 :
        #     # Gotta ngspice or else theres a memory leak
        #     ngspice = simulator.factory(circuit).ngspice
        #     ngspice.remove_circuit()
        #     ngspice.destroy()

        #     circuit = Circuit('sim')  # Remake the circuit
        #     circuit.V('input', 'vin', circuit.gnd, VOLTAGE)

        #     # Re-create resistor array
        #     setup_resistors(circuit, mat_r, nodes, CONTACT_LENGTH//SCALE)
        #     res_heat, simulator, mat_v = get_heat(circuit, mat_r, dv_mat, SCALE)

        # else:
        #     # Keep using the same simulator
        #     update_resistors(circuit, mat_r, nodes, CONTACT_LENGTH//SCALE)
        #     res_heat, simulator, mat_v = get_heat(circuit, mat_r, dv_mat, SCALE)

        # add_head(mat_t, res_heat, hstate_elec, STARTX, STARTY, X_ESIM, Y_ESIM, Z_ESIM)


        # # Get probe temperature
        # probe_t = mat_t[X_GRID//2, Y_GRID//2, 0] - (273.15 + 20)
        # bar.set_postfix({'probe temp: ': probe_t})

        # for node in simulator.operating_point().branches.values():
        #     power_draw += dt * float(node) * VOLTAGE
        

    #  +-------------------------------------------+
    #  |           Print Values                    |
    #  +-------------------------------------------+
    

    print("Saving temperature data... ")
    print("-----------------------------")
    save_data(mat_t, mat_r, mat_v)

    



    print("Printing final result... ")
    print("-----------------------------")

    total_time = dt * NUM_CYCLES
    print("Circuit power draw:", -power_draw)
    print("total time that has passed:", total_time)
    print_matrix(mat_v)
    print_matrix(mat_r)
    print_mat_t(mat_t)








#  +------------------------------------------------------------------+
#  |                                                                  |
#  |                     Printing Functions                           |
#  |                                                                  |
#  +------------------------------------------------------------------+


def print_plane(planes):
    C = 'gist_heat'
    MAXT = 320

    num_args = len(planes)
    x = int(np.sqrt(num_args))
    y = int(np.ceil(num_args/x))

    if x ==1 and y == 1:
        p = plt.imshow(np.transpose(planes[0]), extent=[0, X_GRID,
                                                        0, Y_GRID], cmap='gist_heat', vmin=293.15)
        plt.colorbar(p)
        plt.show()
    else:
        fig, axes = plt.subplots(nrows=x, ncols=y)
        i = 0
        for ax in axes.flat:
            pl = np.transpose(planes[i])
            i += 1
            im = ax.imshow(pl, extent=[0, X_GRID, 0, Y_GRID], cmap=C)
            ax.set_title(f'Layer {i}', fontsize=8)
        fig.colorbar(im, ax=axes.ravel().tolist())

    plt.show()


def save_data(t, r, v):
    filepath = "sim_output"
    # Save temperature
    for i in range(t.shape[2]):
        DF = pd.DataFrame(np.transpose(t[:, :, i]))
        file_name = os.path.join( filepath, "final_temps", "temps_" + str(i) + ".csv")
        DF.to_csv(file_name, header=False, index=False)
        
    # Save resistance
    for i in range(r.shape[2]):
        DF = pd.DataFrame(np.transpose(r[:, :, i]))
        file_name = os.path.join(
            filepath, "final_resistance", "res_" + str(i) + ".csv")
        DF.to_csv(file_name, header=False, index=False)
    
    # Save voltage
    for i in range(v.shape[2]):
        DF = pd.DataFrame(np.transpose(v[:, :, i]))
        file_name = os.path.join(
            filepath, "final_voltage", "v_" + str(i) + ".csv")
        DF.to_csv(file_name, header=False, index=False)

def print_mat_t(mat_t):
    grid1 = get_z_temp(mat_t)
    # grid2 = get_z_temp(mat_t, 1)
    # grid3 = get_z_temp(mat_t, 2)
    # grid4 = get_z_temp(mat_t, 3)
    print_plane(np.array([grid1]))


def print_mat_e(mat_t, startx, starty, x, y, contact_length):
    endx = startx + x
    endy = starty + y

    printmat = np.zeros((mat_t.shape[0], mat_t.shape[1]))
    printmat[startx:endx, starty:endy] = 1
    printmat[startx:startx+contact_length, starty:endy] = 2
    printmat[endx-contact_length:endx, starty:endy] = 2


    plt.imshow(np.transpose(printmat), extent=[0,X_GRID,0,Y_GRID], cmap='plasma')
    plt.show()

def print_matrix(mat):
    p = plt.imshow(np.transpose(mat[:,:,0]), extent=[
        0, mat.shape[0], 0, mat.shape[1]], cmap='GnBu')
    plt.colorbar(p)
    plt.show()


if __name__ == '__main__':
    main()
