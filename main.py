from PySpice import *
from PySpice.Spice.Netlist import Circuit
import PySpice.Logging.Logging as Logging

import torch
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

# Set type of TF acceleration
'''Setup device for PyTorch calculation for Thermo
    cpu: Uses the computer's CPU. 
        - This is required if you want to also have electrostatic sim
        - Can support 64-bit floating point
        - Generally slower than the other two ptions

    cuda: Uses Nvidia's CUDA GPU processors
        - This is only available on GPUs with Nvidia graphics
        - This supports 64-bit floating point
        - Probably the fatsets option

    mps: Uses apple's metal GPU acceleration
        - This is what is used to test GPU acceleration
        - Only supports 32-bit floating point (sadge)
        - Similar to cuda but 32-bit is not accurate enough for our puroses
    
    For testing purposes, I am usually using gpu. For running the actual trial, its probably
        best to use CUDA if available for ~7.5x speed increase from my testing

'''
device = "cuda"

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
    dt = 1.15e-12  # Set timestep, should be within bounds

    # Create primary matrices to use
    mat_d_np = np.zeros((X_GRID, Y_GRID, Z_GRID, 6))
    mat_t = torch.ones((X_GRID, Y_GRID, Z_GRID)) * (273.15 + 20)  # Set to 20C
    h_add = torch.zeros((X_GRID, Y_GRID))
    print("     Finished creating empty matrices ")

    # Set default values for data matrix
    mat_d_np[:, :, :] = make_point()
    mat_d = torch.tensor(mat_d_np)
    print("     Finished setting points ")

    #  +-------------------------------------------+
    #  |           Setup Thermo                    |
    #  +-------------------------------------------+
    print("Setting up Thermo... ")
    
    '''Check for min timestep/ stability'''
    min_time = get_min_timestep(mat_d_np)
    if dt > min_time:
        print("\n=====================================")
        print("defined time:", dt, "is too large for stability")
        print("new timestep of", min_time, "has been selected")
        print("=====================================\n")
        dt = min_time
    print("     - Finished checking timestep ")

    '''Calculate state transition matrix'''
    a_state_np, b_state_np = gen_state_matrix(mat_d_np, dt)
    hstate_t_np = get_hstate_thermo(mat_d_np, dt)

    a_state = torch.tensor(a_state_np).type(torch.float64)
    b_state = torch.tensor(b_state_np).type(torch.float64)
    hstate_t = torch.tensor(hstate_t_np).type(torch.float64)
    print("     - Finished creating a, b, and hstate matrices ")
    
    mask = None

    ''' Set values for heat transfer matrix '''
    laser_points = [[i, 125] for i in range(100, 151)]
    set_added_heat(h_add, LASER_POWER, laser_points)
    print("     - Finished setting laser points ")

    
    '''Set up boundry Temperatures'''
    mask = torch.zeros((X_GRID, Y_GRID, Z_GRID))

    # set edges to 0
    mask[0, :, :] = 1
    mask[X_GRID-1, :, :] = 1
    mask[:, 0, :] = 1
    mask[:, Y_GRID-1, :] = 1
    mask[:, :, Z_GRID-1] = 1

    if mask is None:
        mask = torch.zeros((X_GRID, Y_GRID, Z_GRID))

    print("     - Finished setting boundary temps ")


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
    nodes = torch.zeros((X_ESIM//SCALE, Y_ESIM//SCALE,
                     Z_ESIM//SCALE, 3), dtype=torch.int8)
    for i in range(X_ESIM//SCALE):
        for j in range(Y_ESIM//SCALE):
            for k in range(Z_ESIM//SCALE):
                nodes[i, j, k, 0] = i
                nodes[i, j, k, 1] = j
                nodes[i, j, k, 2] = k
    hstate_elec_np = get_hstate_elec(get_selected_area(mat_d, STARTX, STARTY, X_ESIM, Y_ESIM, Z_ESIM), dt, SCALE)
    hstate_elec = hstate_elec_np.type(torch.float32)
    print("     - Finished setting electro matrices ")
    
    # Make reistor matrix
    mat_r = get_res_matrix(mat_t, mat_d[0, 0, 0, 0], STARTX, STARTY, X_ESIM, Y_ESIM, Z_ESIM, SCALE)
    print("     - Finished getting initial resistances ")
    
    # setup resistors
    circuit = Circuit('sim')  # Remake the circuit
    circuit.V('input', 'vin', circuit.gnd, VOLTAGE)
    circuit.Vinput.minus.add_current_probe
    setup_resistors(circuit, mat_r, nodes, CONTACT_LENGTH//SCALE)
    print("     - Finished creating initial circuit ")


    # #Print the initial Conditions
    # print("Printing initial conditions")
    # print("    Printing Thermal Conditions")
    # print_mat_t(mat_t)
    # print("    Printing Electrostatic Conditions")
    # print_mat_e(mat_t, STARTX, STARTY, X_ESIM, Y_ESIM, CONTACT_LENGTH)





    #  +-----------------------------------------------------+
    #  |                                                     |
    #  |                       Run Loop                      |
    #  |                                                     |
    #  +-----------------------------------------------------+
    print("Running loop... ")

    # Create loop helpers
    new_temps = torch.zeros((X_GRID, Y_GRID, Z_GRID))
    initial_temps = torch.tensor(np.array(mat_t).copy())

    # pre-allocate matrices used for calculation
    comp_mat = torch.zeros((X_GRID, Y_GRID, Z_GRID, 6))
    res_heat = torch.zeros((X_ESIM//SCALE, Y_ESIM//SCALE, Z_ESIM//SCALE))
    dv_mat = torch.zeros((X_ESIM//SCALE, Y_ESIM//SCALE, Z_ESIM//SCALE, 6))
    
    power_draw = 0
    ref_length = mat_d[0, 0, 0, 0]

    # Create matrix of voltages
    mat_v = torch.zeros((X_ESIM//SCALE, Y_ESIM//SCALE, Z_ESIM//SCALE))
    print("     - Finished Loop setup")

    # Move all matrices to device 
    new_temps = new_temps.to(device=device)
    initial_temps = initial_temps.to(device=device)
    comp_mat = comp_mat.to(device=device)
    mat_t = mat_t.to(device=device)
    hstate_t = hstate_t.to(device=device)
    h_add = h_add.to(device=device)
    mask = mask.to(device=device)

    a_state = a_state.to(device=device)
    b_state = b_state.to(device=device)

    print("     - Finished Moving to pyTorch device")
    print("Running Loop using device:", device)

    # Run loop
    bar = tqdm(range(NUM_CYCLES), desc="Running Sim")
    for i in bar:
        # -------------------------------
        #   Thermal Sim
        # -------------------------------
        transition_state(mat_t, comp_mat, a_state, b_state, new_temps)
        mat_t = set_mat_temps(mask, initial_temps, new_temps)
        apply_heat(mat_t, hstate_t, h_add)


        '''Commenting out the Electro Sim so we can use CUDA acceleration for purely thermo   '''

        '''=============== START COMMENT ==============='''
        # # -------------------------------
        # #   Electrostatic Sim
        # # -------------------------------
        # mat_r = get_res_matrix(mat_t, ref_length, STARTX, STARTY, X_ESIM, Y_ESIM, Z_ESIM, SCALE)

        # '''For every N cycles, reset spice to make sure we dotn use too much memory'''
        # N = 1000
        # if i%N == 0 :
        #     # Gotta ngspice or else theres a memory leak
        #     if (i != 0) :
        #         ngspice = simulator.factory(circuit).ngspice
        #         ngspice.remove_circuit()
        #         ngspice.destroy()

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

        '''=============== END COMMENT ==============='''
        
    new_temps = new_temps.to(device="cpu")
    initial_temps = initial_temps.to(device="cpu")
    comp_mat = comp_mat.to(device="cpu")
    mat_t = mat_t.to(device="cpu")
    hstate_t = hstate_t.to(device="cpu")
    h_add = h_add.to(device="cpu")
    mask = mask.to(device="cpu")

    a_state = a_state.to(device="cpu")
    b_state = b_state.to(device="cpu")

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
        p = plt.imshow(torch.t(planes[0]), extent=[0, X_GRID,
                                                        0, Y_GRID], cmap='gist_heat', vmin=293.15)
        plt.colorbar(p)
        plt.show()
    else:
        fig, axes = plt.subplots(nrows=x, ncols=y)
        i = 0
        for ax in axes.flat:
            pl = torch.t(planes[i])
            i += 1
            im = ax.imshow(pl, extent=[0, X_GRID, 0, Y_GRID], cmap=C)
            ax.set_title(f'Layer {i}', fontsize=8)
        fig.colorbar(im, ax=axes.ravel().tolist())

    plt.show()


def save_data(t, r, v):
    # t = np.array(t_tf)
    # r = np.array(r_tf)
    # v = np.array(v_tf)




    filepath = "sim_output"
    # Save temperature
    for i in range(t.shape[2]):
        DF = pd.DataFrame(torch.t(t[:, :, i]))
        file_name = os.path.join( filepath, "final_temps", "temps_" + str(i) + ".csv")
        DF.to_csv(file_name, header=False, index=False)
        
    # Save resistance
    for i in range(r.shape[2]):
        DF = pd.DataFrame(torch.t(r[:, :, i]))
        file_name = os.path.join(
            filepath, "final_resistance", "res_" + str(i) + ".csv")
        DF.to_csv(file_name, header=False, index=False)
    
    # Save voltage
    for i in range(v.shape[2]):
        DF = pd.DataFrame(torch.t(v[:, :, i]))
        file_name = os.path.join(
            filepath, "final_voltage", "v_" + str(i) + ".csv")
        DF.to_csv(file_name, header=False, index=False)

def print_mat_t(mat_t):
    # mat_t = np.array(mat_t_tf)
    grid1 = get_z_temp(mat_t)
    # grid2 = get_z_temp(mat_t, 1)
    # grid3 = get_z_temp(mat_t, 2)
    # grid4 = get_z_temp(mat_t, 3)
    print_plane([grid1])


def print_mat_e(mat_t, startx, starty, x, y, contact_length):
    # mat_t = np.array(mat_t_tf)
    endx = startx + x
    endy = starty + y

    printmat = torch.zeros((mat_t.shape[0], mat_t.shape[1]))
    printmat[startx:endx, starty:endy] = 1
    printmat[startx:startx+contact_length, starty:endy] = 2
    printmat[endx-contact_length:endx, starty:endy] = 2


    plt.imshow(torch.t(printmat), extent=[0,X_GRID,0,Y_GRID], cmap='plasma')
    plt.show()

def print_matrix(mat):
    # mat = np.array(mat_tf)
    p = plt.imshow(torch.t(mat[:,:,0]), extent=[
        0, mat.shape[0], 0, mat.shape[1]], cmap='GnBu')
    plt.colorbar(p)
    plt.show()


if __name__ == '__main__':
    main()
