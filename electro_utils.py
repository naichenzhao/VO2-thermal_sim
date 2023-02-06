from PySpice import *
from PySpice.Spice.Netlist import Circuit
import PySpice.Logging.Logging as Logging
import numpy as np


#  +------------------------------------------------+
#  |                 Main Functions                 |
#  +------------------------------------------------+

def get_heat(circuit, r_mat, dv_mat):
    # set up variables:
    simulator = circuit.simulator()
    analysis = simulator.operating_point()

    volt_matrix = np.zeros(r_mat.shape)


    for node in analysis.nodes.values():
        if str(node) == 'vin':
            continue
        currP = str(node).split('/')
        volt_matrix[int(currP[0]), int(currP[1]), int(currP[2])] = float(node)
    
    X = r_mat.shape[0]
    Y = r_mat.shape[1]
    Z = r_mat.shape[2]

    # ---------- Calculate X points ----------
    dv_mat[:X-1, :, :, 0] = volt_matrix[:X-1, :, :] - volt_matrix[1:, :, :]
    dv_mat[1:, :, :, 1] = volt_matrix[1:, :, :] - volt_matrix[:X-1, :, :]

    # ---------- Calculate Y points ----------
    dv_mat[:, :Y-1, :, 2] = volt_matrix[:, :Y-1, :] - volt_matrix[:, 1:, :]
    dv_mat[:, 1:, :, 3] = volt_matrix[:, 1:, :] - volt_matrix[:, :Y-1, :]

    # ---------- Calculate Z points ----------
    dv_mat[:, :, :Z-1, 4] = volt_matrix[:, :, :Z-1] - volt_matrix[:, :, 1:]
    dv_mat[:, :, 1:, 5] = volt_matrix[:, :, 1:] - volt_matrix[:, :, :Z-1]
    
    return np.sum( dv_mat * dv_mat, axis=3)/r_mat, simulator, volt_matrix


def add_head(mat_t, r_heat, dt, startx, starty, x, y, z, S):

    heat_added = np.zeros((x, y, z))
    for i in range (r_heat.shape[0]):
        for j in range (r_heat.shape[1]):
            for k in range (r_heat.shape[2]):
                heat_added[S*i:S*(i+1), S*j:S*(j+1), S*k:S*(k+1)] = r_heat[i, j, k]
    endx = startx + x
    endy = starty + y
    mat_t[startx:endx, starty:endy, 0:z] = mat_t[startx:endx, starty:endy, 0:z] + heat_added * dt






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

def get_res_matrix(mat_t, L, STARTX, STARTY, X_ESIM, Y_ESIM, Z_ESIM, S):
    em = get_selected_area(mat_t, STARTX, STARTY, X_ESIM, Y_ESIM, Z_ESIM)

    t_scaled = np.zeros((em.shape[0]//S, em.shape[1]//S, em.shape[2]//S))

    for i in range(t_scaled.shape[0]):
        for j in range(t_scaled.shape[1]):
            for k in range(t_scaled.shape[2]):
                t_scaled[i, j, k] = np.mean(em[S*i:S*(i+1), S*j:S*(j+1), S*k:S*(k+1)])

    r_mat = get_resistivity(t_scaled, alpha=0, r_0=0.1)/(L * 2)
    return r_mat


def get_resistivity(M, r_0 = 1.68e-8, alpha = 0.00386, T0 = 293.15):
    T0_mat = np.ones(M.shape) * T0
    return np.maximum(r_0 * (1 + alpha*(M-T0_mat)), np.zeros(M.shape))




