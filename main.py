import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt

from datapoint import *
from utils import *

    
#  +-------------------------------------------+
#  |           Constant values                 |
#  +-------------------------------------------+
 
X_GRID = 1001
Y_GRID = 20
Z_GRID = 20
r = 1                 
dt = 1                




def main():
    #  +-------------------------------------------+
    #  |           Setup Matrix                    |
    #  +-------------------------------------------+

    mat = np.empty((X_GRID, Y_GRID, Z_GRID), dtype=object)

    for i in range(X_GRID):
        for j in range(Y_GRID):
            for k in range(Z_GRID):
                set_point(mat, (i, j, k), DataPoint())

    next_points = set_mat(mat, (100, 100), ((0, 0, 0), (1000, 0, 0)))


    #  +-------------------------------------------+
    #  |              Run Loop                     |
    #  +-------------------------------------------+
    for i in tqdm(range(5)):
        next_points = sim(mat, next_points, dt)


    print(mat)






def set_mat(mat, temps, coordinates):
    r_set = set()

    for i in range(len(temps)):
        curr_loc = coordinates[i]
        curr_temp = temps[i]

        set_point(mat, curr_loc, DataPoint(t = curr_temp, s=True))

        r_set.add(curr_loc)
        r_set.update(get_adj(mat, curr_loc))
    
    return r_set


def sim(mat, vals, dt):
    ret_set = set()
    temp_mat = np.zeros((mat.shape[0], mat.shape[1], mat.shape[2]))

    for p in vals:
        t, points = calculate_dt(mat, p, dt, r)
        temp_mat[p[0]][p[1]][p[2]] = t 
        ret_set.update(points)

    for p in vals:
        set_temp(mat, p, temp_mat[p[0]][p[1]][p[2]])              

    return ret_set




if __name__ == '__main__':
    main()

    

