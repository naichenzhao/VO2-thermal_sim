import numpy as np
from utils import *

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

    a_state = np.empty((X_GRID, Y_GRID, Z_GRID))
    b_state = np.empty((X_GRID, Y_GRID, Z_GRID, 6))

    for i in range(X_GRID):
        for j in range(Y_GRID):
            for k in range(Z_GRID):
                curr_p = mat_d[i, j, k]
                a_sum = 0
                b_const = (dt/(curr_p[4] * curr_p[5])/4)

                A = curr_p[2] * curr_p[1]  # dz * dy
                if i + 1 < X_GRID:  # Get point for (1, 0, 0)
                    target = mat_d[i+1, j, k]
                    a_sum += 1/((curr_p[0] + target[0]) * A)
                    b_state[i, j, k, 0] = b_const/((curr_p[0] + target[0]) * A)

                if i - 1 >= 0:  # Get point for (-1, 0, 0)
                    target = mat_d[i-1, j, k]
                    a_sum += 1/((curr_p[0] + target[0]) * A)
                    b_state[i, j, k, 1] = b_const/((curr_p[0] + target[0]) * A)

                # ---------- Calculate Y points ----------
                A = curr_p[2] * curr_p[0]  # dz * dx
                if j + 1 < Y_GRID:  # Get point for (0, 1, 0)
                    target = mat_d[i, j+1, k]
                    a_sum += 1/((curr_p[1] + target[1]) * A)
                    b_state[i, j, k, 2] = b_const/((curr_p[1] + target[1]) * A)

                if j - 1 >= 0:  # Get point for (0, -1, 0)
                    target = mat_d[i, j-1, k]
                    a_sum += 1/((curr_p[1] + target[1]) * A)
                    b_state[i, j, k, 3] = b_const/((curr_p[1] + target[1]) * A)

                # ---------- Calculate Z points ----------
                A = curr_p[0] * curr_p[1]  # dx * dy
                if k + 1 < Z_GRID:  # Get point for (0, 0, 1)
                    target = mat_d[i, j, k+1]
                    a_sum += 1/((curr_p[2] + target[2]) * A)
                    b_state[i, j, k, 4] = b_const/((curr_p[2] + target[2]) * A)

                if k - 1 >= 0:  # Get point for (0, 0, -1)
                    target = mat_d[i, j, k-1]
                    a_sum += 1/((curr_p[2] + target[2]) * A)
                    b_state[i, j, k, 5] = b_const/((curr_p[2] + target[2]))/A

                a_state[i, j, k] = (1 - (dt/(curr_p[4] * curr_p[5])) * a_sum/4)

    return a_state, b_state



def transition_state(mat_t, comp_mat, a_state, b_state):
    X_GRID = mat_t.shape[0]
    Y_GRID = mat_t.shape[1]
    Z_GRID = mat_t.shape[2]

    comp_mat[:, :, :, 0] = np.roll(mat_t, 1, axis= 0)
    comp_mat[0, :, :, 0] = np.zeros((Y_GRID, Z_GRID))

    comp_mat[:, :, :, 1] = np.roll(mat_t, -1, axis=0)
    comp_mat[X_GRID-1, :, :, 1] = np.zeros((Y_GRID, Z_GRID))

    comp_mat[:, :, :, 2] = np.roll(mat_t, 1, axis=1)
    comp_mat[:, 0, :, 2] = np.zeros((X_GRID, Z_GRID))

    comp_mat[:, :, :, 3] = np.roll(mat_t, -1, axis = 1)
    comp_mat[:, Y_GRID-1, :, 3] = np.zeros((X_GRID, Z_GRID))

    comp_mat[:, :, :, 4] = np.roll(mat_t, 1, axis=2)
    comp_mat[:, :, 0, 4] = np.zeros((X_GRID, Y_GRID))

    comp_mat[:, :, :, 5] = np.roll(mat_t, -1, axis=2)
    comp_mat[:, :, Z_GRID-1, 5] = np.zeros((X_GRID, Y_GRID))

    return mat_t * a_state + np.sum(comp_mat*b_state, axis=3)


def set_mat_temps(mat_t, new, initial):
    X_GRID = mat_t.shape[0]
    Y_GRID = mat_t.shape[1]
    Z_GRID = mat_t.shape[2]
    return np.where(initial > 0, initial, new)


