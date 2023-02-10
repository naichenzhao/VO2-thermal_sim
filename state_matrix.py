import numpy as np
from thermo_utils import *


def gen_state_matrix(mat_d, dt):
    ''' state translation: for b-state [A, B, C, D, E, F]
    A: x + 1
    B: x - 1
    C: y + 1
    D: y - 1
    E: z + 1
    F: z - 1
    '''
    X_GRID = mat_d.shape[0]
    Y_GRID = mat_d.shape[1]
    Z_GRID = mat_d.shape[2]

    a_state = np.zeros((X_GRID, Y_GRID, Z_GRID), dtype=np.float64)
    b_state = np.zeros((X_GRID, Y_GRID, Z_GRID, 6), dtype=np.float64)

    for i in range(X_GRID):
        for j in range(Y_GRID):
            for k in range(Z_GRID):
                curr_p = mat_d[i, j, k]
                a_sum = 0
                b_const = (dt/(curr_p[4] * curr_p[5]))


                # ---------- Calculate X points ----------
                A = curr_p[2] * curr_p[1] * 4  # dz * dy

                if i + 1 < X_GRID:  # Get point for (1, 0, 0)
                    target = mat_d[i+1, j, k]
                    a_sum += 1/((curr_p[0] + target[0]) * A)
                    b_state[i, j, k, 0] = b_const/((curr_p[0] + target[0]) * A)

                if i - 1 >= 0:  # Get point for (-1, 0, 0)
                    target = mat_d[i-1, j, k]
                    a_sum += 1/((curr_p[0] + target[0]) * A)
                    b_state[i, j, k, 1] = b_const/((curr_p[0] + target[0]) * A)

                # ---------- Calculate Y points ----------
                A = curr_p[2] * curr_p[0] * 4  # dz * dx

                if j + 1 < Y_GRID:  # Get point for (0, 1, 0)
                    target = mat_d[i, j+1, k]
                    a_sum += 1/((curr_p[1] + target[1]) * A)
                    b_state[i, j, k, 2] = b_const/((curr_p[1] + target[1]) * A)

                if j - 1 >= 0:  # Get point for (0, -1, 0)
                    target = mat_d[i, j-1, k]
                    a_sum += 1/((curr_p[1] + target[1]) * A)
                    b_state[i, j, k, 3] = b_const/((curr_p[1] + target[1]) * A)

                # ---------- Calculate Z points ----------
                A = curr_p[0] * curr_p[1] * 4  # dx * dy

                if k + 1 < Z_GRID:  # Get point for (0, 0, 1)
                    target = mat_d[i, j, k+1]
                    a_sum += 1/((curr_p[2] + target[2]) * A)
                    b_state[i, j, k, 4] = b_const/((curr_p[2] + target[2]) * A)

                if k - 1 >= 0:  # Get point for (0, 0, -1)
                    target = mat_d[i, j, k-1]
                    a_sum += 1/((curr_p[2] + target[2]) * A)
                    b_state[i, j, k, 5] = b_const/((curr_p[2] + target[2]))/A

                a_state[i, j, k] = (1 - (dt/(curr_p[4] * curr_p[5])) * a_sum)

    return a_state, b_state


def get_hstate_thermo(mat_d, dt):
    X_GRID = mat_d.shape[0]
    Y_GRID = mat_d.shape[1]

    h_state = np.zeros((X_GRID, Y_GRID), dtype=np.float64)
    mass = (8*mat_d[:, :, 0, 0]*mat_d[:, :, 0, 1]*mat_d[:, :, 0, 2])*mat_d[:,:,0,5]
    h_state = dt/(mass * mat_d[:, :, 0, 4])

    return h_state


def transition_state(mat_t, comp_mat, a_state, b_state):
    X_GRID = mat_t.shape[0]
    Y_GRID = mat_t.shape[1]
    Z_GRID = mat_t.shape[2]

    # ---------- Calculate X points ----------
    comp_mat[:X_GRID-1, :, :, 0] = mat_t[1:, :, :]
    comp_mat[1:, :, :, 1] = mat_t[:X_GRID-1, :, :]
    

    # ---------- Calculate Y points ----------
    comp_mat[:, :Y_GRID-1, :, 2] = mat_t[:, 1:, :]
    comp_mat[:, 1:, :, 3] = mat_t[:, :Y_GRID-1, :]
    

    # ---------- Calculate Z points ----------
    comp_mat[:, :, :Z_GRID-1, 4] = mat_t[:, :, 1:]
    comp_mat[:, :, 1:, 5] = mat_t[:, :, :Z_GRID-1]
    

    return mat_t * a_state + np.sum(comp_mat*b_state, axis=3)


def set_mat_temps(mask, initial, new):
    return np.where(mask > 0, initial, new)
