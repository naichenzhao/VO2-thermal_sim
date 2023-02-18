import numpy as np
from thermo_utils import *
import threading
import torch


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

    a_state = np.zeros((X_GRID, Y_GRID, Z_GRID))
    b_state = np.zeros((X_GRID, Y_GRID, Z_GRID, 6))

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

    h_state = torch.zeros((X_GRID, Y_GRID))
    mass = (8*mat_d[:, :, 0, 0]*mat_d[:, :, 0, 1]*mat_d[:, :, 0, 2])*mat_d[:,:,0,5]
    h_state = dt/(mass * mat_d[:, :, 0, 4])

    return h_state


def transition_state(mat_t, comp_mat, a_state, b_state, out):
    X_GRID = mat_t.shape[0]
    Y_GRID = mat_t.shape[1]
    Z_GRID = mat_t.shape[2]

    ''' Previous way of calculating the comp_mat
    '''
    # ---------- Calculate X points ----------
    comp_mat[:X_GRID-1, :, :, 0] = mat_t[1:, :, :]
    comp_mat[1:, :, :, 1] = mat_t[:X_GRID-1, :, :]
    

    # ---------- Calculate Y points ----------
    comp_mat[:, :Y_GRID-1, :, 2] = mat_t[:, 1:, :]
    comp_mat[:, 1:, :, 3] = mat_t[:, :Y_GRID-1, :]
    

    # ---------- Calculate Z points ----------
    comp_mat[:, :, :Z_GRID-1, 4] = mat_t[:, :, 1:]
    comp_mat[:, :, 1:, 5] = mat_t[:, :, :Z_GRID-1]


    ''' Previous way(s) of getting the resut
    '''
    # sum_mat = np.sum(comp_mat*b_state, axis=3)
    # return pyfma.fma(mat_t, a_state, sum_mat)
    # return ne.evaluate('a*b+c', local_dict={'a': mat_t, 'b': a_state, 'c': np.sum(comp_mat*b_state, axis=3)})
    out[:] =  mat_t * a_state + torch.sum(comp_mat*b_state, axis=3)

    # '''New approach to calculation uses multi-threaded approach to speed up calculations'''
    # par_comp_mat(comp_mat, mat_t, func=eq)
    # par_fmadd(mat_t, a_state, comp_mat, b_state, 10, out)


def set_mat_temps(mask, initial, new):
    return torch.where(mask > 0, initial, new)


def fmadd(a, b, c, d, out):
    out[:] = a*b+torch.sum(c*d, axis=3)


def eq(a, b):
    a[:] = b[:]


def par_fmadd(a, b, c, d, numCores, out, func=fmadd):

    threads = []
    pl = a.shape[0]//numCores
    for i in range(numCores - 1):
        th = threading.Thread(target=func, args=(
            a[i*pl:(i+1)*pl, :, :], 
            b[i*pl:(i+1)*pl, :, :], 
            c[i*pl:(i+1)*pl, :, :, :], 
            d[i*pl:(i+1)*pl, :, :, :], 
            out[i*pl:(i+1)*pl, :, :]))
        th.start()
        threads.append(th)

    for th in threads:
        th.join()


def par_comp_mat(comp_mat, mat_t, func=eq):

    X_GRID = mat_t.shape[0]
    Y_GRID = mat_t.shape[1]
    Z_GRID = mat_t.shape[2]

    threads = []

    # ---------- Calculate X points ----------
    th0 = threading.Thread(target=func, args=(
        comp_mat[:X_GRID-1, :, :, 0], mat_t[1:, :, :] ))
    th1 = threading.Thread(target=func, args=(
        comp_mat[1:, :, :, 1], mat_t[:X_GRID-1, :, :]))
    th0.start()
    th1.start()
    threads.append(th0)
    threads.append(th1)

    # ---------- Calculate X points ----------
    th2 = threading.Thread(target=func, args=(
        comp_mat[:, :Y_GRID-1, :, 2], mat_t[:, 1:, :]))
    th3 = threading.Thread(target=func, args=(
        comp_mat[:, 1:, :, 3], mat_t[:, :Y_GRID-1, :]))
    th2.start()
    th3.start()
    threads.append(th2)
    threads.append(th3)

    # ---------- Calculate X points ----------
    th4 = threading.Thread(target=func, args=(
        comp_mat[:, :, :Z_GRID-1, 4], mat_t[:, :, 1:]))
    th5 = threading.Thread(target=func, args=(
        comp_mat[:, :, 1:, 5], mat_t[:, :, :Z_GRID-1]))
    th4.start()
    th5.start()
    threads.append(th4)
    threads.append(th5)

    for th in threads:
        th.join()


