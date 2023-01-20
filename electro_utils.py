from PySpice import *
from PySpice.Spice.Netlist import Circuit
import PySpice.Logging.Logging as Logging
import numpy as np


#  +------------------------------------------------+
#  |                 Main Functions                 |
#  +------------------------------------------------+

def get_voltage_matrix():
    # set up variables:
    return None








#  +------------------------------------------------+
#  |               PySpice Functions                |
#  +------------------------------------------------+


def update_resistors(circuit, r_mat, nodes, contact_length):
    X_LEN = nodes.shape[0]
    Y_LEN = nodes.shape[1]
    Z_LEN = nodes.shape[2]
    c = 0

    # Set resistors for main array IJK
    for i in range(X_LEN-1):
        for j in range(Y_LEN-1):
            for k in range(Z_LEN-1):
                circuit['R' + str(c)].resistance = r_mat[i, j, k] + r_mat[i+1, j, k]
                c += 1
                circuit['R' + str(c)].resistance = r_mat[i, j, k] + r_mat[i, j+1, k]
                c += 1
                circuit['R' + str(c)].resistance = r_mat[i,j, k] + r_mat[i, j, k+1]
                c += 1

    # set final JK
    for j in range(Y_LEN-1):
        for k in range(Z_LEN-1):
            circuit['R' + str(c)].resistance = r_mat[X_LEN-1, j, k] + r_mat[X_LEN-1, j+1, k]
            c += 1
            circuit['R' + str(c)].resistance = r_mat[X_LEN-1, j, k] + r_mat[X_LEN-1, j, k+1]
            c += 1

    # set final IK
    for i in range(X_LEN-1):
        for k in range(Z_LEN-1):
            circuit['R' + str(c)].resistance = r_mat[i, Y_LEN-1, k] + r_mat[i+1, Y_LEN-1, k]
            c += 1
            circuit['R' + str(c)].resistance = r_mat[i, Y_LEN-1, k] + r_mat[i, Y_LEN-1, k+1]
            c += 1

    # set final IJ
    for i in range(X_LEN-1):
        for j in range(Y_LEN-1):
            circuit['R' + str(c)].resistance = r_mat[i, j, Z_LEN-1] + r_mat[i+1, j, Z_LEN-1]
            c += 1
            circuit['R' + str(c)].resistance = r_mat[i, j, Z_LEN-1] + r_mat[i, j+1, Z_LEN-1]
            c += 1

    # set final I
    for i in range(X_LEN-1):
        circuit['R' + str(c)].resistance = r_mat[i, Y_LEN-1, Z_LEN-1] + r_mat[i+1, Y_LEN-1, Z_LEN-1]
        c += 1

    # set final J
    for j in range(Y_LEN-1):
        circuit['R' + str(c)].resistance = r_mat[X_LEN-1, j, Z_LEN-1] + r_mat[X_LEN-1, j+1, Z_LEN-1]
        c += 1

    # set final K
    for k in range(Z_LEN-1):
        circuit['R' + str(c)].resistance = r_mat[X_LEN-1, Y_LEN-1, k] + r_mat[X_LEN-1, Y_LEN-1, k+1]
        c += 1

    # Setup the voltage sources and contacts
    for j in range(Y_LEN):
        for i in range(contact_length):
            circuit['R' + str(c)].resistance = r_mat[i, j, 0]
            c += 1
            circuit['R' + str(c)].resistance = r_mat[i, j, 0]
            c += 1





def setup_resistors(circuit, r_mat, nodes, contact_length):
    X_LEN = nodes.shape[0]
    Y_LEN = nodes.shape[1]
    Z_LEN = nodes.shape[2]
    c = 0

    # Set resistors for main array IJK
    for i in range(X_LEN-1):
        for j in range(Y_LEN-1):
            for k in range(Z_LEN-1):
                n = get_node(nodes, i, j, k)
                circuit.R(c, n, get_node(nodes, i+1, j, k),
                        r_mat[i, j, k] + r_mat[i+1, j, k])
                c += 1
                circuit.R(c, n, get_node(nodes, i, j+1, k),
                        r_mat[i, j, k] + r_mat[i, j+1, k])
                c += 1
                circuit.R(c, n, get_node(nodes, i, j, k+1),
                        r_mat[i, j, k] + r_mat[i, j, k+1])
                c += 1

    # set final JK
    for j in range(Y_LEN-1):
        for k in range(Z_LEN-1):
            circuit.R(c, get_node(nodes, X_LEN-1, j, k), get_node(nodes, X_LEN-1, j+1, k),
                      r_mat[X_LEN-1, j, k] + r_mat[X_LEN-1, j+1, k])
            c += 1
            circuit.R(c, get_node(nodes, X_LEN-1, j, k), get_node(nodes, X_LEN-1, j, k+1),
                      r_mat[X_LEN-1, j, k] + r_mat[X_LEN-1, j, k+1])
            c += 1

    # set final IK
    for i in range(X_LEN-1):
        for k in range(Z_LEN-1):
            circuit.R(c, get_node(nodes, i, Y_LEN-1, k), get_node(nodes, i+1, Y_LEN-1, k),
                      r_mat[i, Y_LEN-1, k] + r_mat[i+1, Y_LEN-1, k])
            c += 1
            circuit.R(c, get_node(nodes, i, Y_LEN-1, k), get_node(nodes, i, Y_LEN-1, k+1),
                      r_mat[i, Y_LEN-1, k] + r_mat[i, Y_LEN-1, k+1])
            c += 1

    # set final IJ
    for i in range(X_LEN-1):
        for j in range(Y_LEN-1):
            circuit.R(c, get_node(nodes, i, j, Z_LEN-1), get_node(nodes, i+1, j, Z_LEN-1),
                      r_mat[i, j, Z_LEN-1] + r_mat[i+1, j, Z_LEN-1])
            c += 1
            circuit.R(c, get_node(nodes, i, j, Z_LEN-1), get_node(nodes, i, j+1, Z_LEN-1),
                      r_mat[i, j, Z_LEN-1] + r_mat[i, j+1, Z_LEN-1])
            c += 1

    # set final I
    for i in range(X_LEN-1):
        circuit.R(c, get_node(nodes, i, Y_LEN-1, Z_LEN-1), get_node(nodes, i+1, Y_LEN-1, Z_LEN-1),
                  r_mat[i, Y_LEN-1, Z_LEN-1] + r_mat[i+1, Y_LEN-1, Z_LEN-1])
        c += 1

    # set final J
    for j in range(Y_LEN-1):
        circuit.R(c, get_node(nodes, X_LEN-1, j, Z_LEN-1), get_node(nodes, X_LEN-1, j+1, Z_LEN-1),
                  r_mat[X_LEN-1, j, Z_LEN-1] + r_mat[X_LEN-1, j+1, Z_LEN-1])
        c += 1

    # set final K
    for k in range(Z_LEN-1):
        circuit.R(c, get_node(nodes, X_LEN-1, Y_LEN-1, k), get_node(nodes, X_LEN-1, Y_LEN-1, k+1),
                  r_mat[X_LEN-1, Y_LEN-1, k] + r_mat[X_LEN-1, Y_LEN-1, k+1])
        c += 1

    # Setup the voltage sources and contacts
    for j in range(Y_LEN):
        for i in range(contact_length):
            circuit.R(c, 'vin', get_node(nodes, i, j, 0), r_mat[i, j, 0])
            c += 1
            circuit.R(c, circuit.gnd, get_node(nodes, X_LEN - i- 1, j, 0), r_mat[i, j, 0])
            c += 1



def get_node(nodes, i, j, k):
    return str(nodes[i, j, k, 0]) + "/" + str(nodes[i, j, k, 1]) + "/" + str(nodes[i, j, k, 2])

#  +------------------------------------------------+
#  |                Helper Functions                |
#  +------------------------------------------------+

def get_selected_area(mat, startx, starty, x, y, z):
    endx = startx + x
    endy = starty + y
    return mat[startx:endx, starty:endy, 0:z]



#  +------------------------------------------------+
#  |             Resistance Functions               |
#  +------------------------------------------------+

def get_res_matrix(mat_t, mat_d, STARTX, STARTY, X_ESIM, Y_ESIM, Z_ESIM):
    emat_t = get_selected_area(mat_t, STARTX, STARTY, X_ESIM, Y_ESIM, Z_ESIM)
    emat_d = get_selected_area(mat_d[:,:,:,0], STARTX, STARTY, X_ESIM, Y_ESIM, Z_ESIM)

    r_mat = get_resistivity(emat_t, alpha = -0.001)/(emat_d * 2)
    return r_mat


def get_resistivity(M, r_0 = 1.68e-8, alpha = 0.00386, T0 = 293.15):
    T0_mat = np.ones(M.shape) * T0
    return r_0 * (1 + alpha*(M-T0_mat))




