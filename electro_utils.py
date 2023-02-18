from PySpice import *
from PySpice.Spice.Netlist import Circuit
import PySpice.Logging.Logging as Logging
import numpy as np
import torch



#  +------------------------------------------------+
#  |                 Main Functions                 |
#  +------------------------------------------------+

def get_heat(circuit, r_mat, dv_mat, scale):
    # set up variables:
    simulator = circuit.simulator()
    analysis = simulator.operating_point()

    volt_matrix = torch.zeros(r_mat.shape)
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
    
    #print(dv_mat)
    #print(volt_matrix)
    return size_up(torch.sum(dv_mat * dv_mat, axis=3)/(r_mat * 4), scale), simulator, volt_matrix


def add_head(mat_t, r_heat, hstate, startx, starty, x, y, z):
    endx = startx + x
    endy = starty + y

    #print((r_heat*hstate))
    mat_t[startx:endx, starty:endy, 0:z] = mat_t[startx:endx,starty:endy, 0:z] + (r_heat*hstate)


def get_hstate_elec(mat_d, dt, s):
    X_GRID = mat_d.shape[0]
    Y_GRID = mat_d.shape[1]
    Z_GRID = mat_d.shape[2]

    h_state = np.zeros((X_GRID, Y_GRID, Z_GRID))
    mass = (8*mat_d[:,:,:,0]*mat_d[:,:,:,1]*mat_d[:,:,:,2])*mat_d[:,:,:,5]*(s**3)
    h_state = dt/(mass * mat_d[:,:,:,4])

    return h_state




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

def size_up(m, s):
    new_m = torch.zeros( (m.shape[0]*s, m.shape[1]*s, m.shape[2]*s) )
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            for k in range(m.shape[2]):
                new_m[s*i:s*(i+1), s*j:s*(j+1), s*k:s*(k+1)] = m[i, j, k]
    return new_m

def size_down(m, s):
    new_m = torch.zeros((m.shape[0]//s, m.shape[1]//s, m.shape[2]//s))
    for i in range(m.shape[0]//s):
        for j in range(m.shape[1]//s):
            for k in range(m.shape[2]//s):
                new_m[i, j, k] = torch.mean(
                    m[s*i:s*(i+1), s*j:s*(j+1), s*k:s*(k+1)])
    return new_m



#  +------------------------------------------------+
#  |             Resistance Functions               |
#  +------------------------------------------------+

def get_res_matrix(mat_t, L, STARTX, STARTY, X_ESIM, Y_ESIM, Z_ESIM, S):
    em = get_selected_area(mat_t, STARTX, STARTY, X_ESIM, Y_ESIM, Z_ESIM)
    r_mat = get_resistivity(size_down(em, S))/(L * S * 4)
    return r_mat


def get_resistivity(M):
    r_mat = (9.225e-11) * M**5 + (-1.665e-7) * M**4 + (0.0001994) * M**3 + (-0.04253) * M**2 + (7.511) * M + (-525.9)
    return torch.maximum(r_mat, torch.zeros(M.shape))




