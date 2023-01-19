from PySpice import *
from PySpice.Spice.Netlist import Circuit

import numpy as np

import PySpice.Logging.Logging as Logging
logger = Logging.setup_logging()


circuit = Circuit('voltage bridge')
circuit.V('input', 'vin', circuit.gnd, '10V')

# Make node coordinates
linex = range(0, 100)  # x
liney = range(0, 3)  # y
linez = range(0, 100)  # z
nodes = np.array([[[str(i) + str(j) + str(k) for k in linez]
                 for j in liney]for i in linex])

x = nodes.shape[0]
y = nodes.shape[1]
z = nodes.shape[2]

# Make reistor matrix
r_mat = np.ones((x, y, z))*500

print(nodes)
print(r_mat)

c = 0

# Set resistors for main array IJK
for i in range(x-1):
    for j in range(y-1):
        for k in range(z-1):
            n = nodes[i, j, k]
            circuit.R(c, n, nodes[i+1, j, k],
                      r_mat[i, j, k] + r_mat[i+1, j, k])
            c += 1

            circuit.R(c, n, nodes[i, j+1, k],
                      r_mat[i, j, k] + r_mat[i, j+1, k])
            c += 1

            circuit.R(c, n, nodes[i, j, k+1],
                      r_mat[i, j, k] + r_mat[i, j, k+1])
            c += 1


ex = x-1
ey = y-1
ez = z-1

# set final JK
for j in range(y-1):
    for k in range(z-1):
        circuit.R(c, nodes[ex, j, k], nodes[ex, j+1, k],
                  r_mat[ex, j, k] + r_mat[ex, j+1, k])
        c += 1
        circuit.R(c, nodes[ex, j, k], nodes[ex, j, k+1],
                  r_mat[ex, j, k] + r_mat[ex, j, k+1])
        c += 1

# set final IK
for i in range(x-1):
    for k in range(z-1):
        circuit.R(c, nodes[i, ey, k], nodes[i+1, ey, k],
                  r_mat[i, ey, k] + r_mat[i+1, ey, k])
        c += 1
        circuit.R(c, nodes[i, ey, k], nodes[i, ey, k+1],
                  r_mat[i, ey, k] + r_mat[i, ey, k+1])
        c += 1

# set final IJ
for i in range(x-1):
    for j in range(y-1):
        circuit.R(c, nodes[i, j, ez], nodes[i+1, j, ez],
                  r_mat[i, j, ez] + r_mat[i+1, j, ez])
        c += 1
        circuit.R(c, nodes[i, j, ez], nodes[i, j+1, ez],
                  r_mat[i, j, ez] + r_mat[i, j+1, ez])
        c += 1


# set final I
for i in range(x-1):
    circuit.R(c, nodes[i, ey, ez], nodes[i+1, ey, ez],
              r_mat[i, ey, ez] + r_mat[i+1, ey, ez])
    c += 1

# set final J
for j in range(y-1):
    circuit.R(c, nodes[ex, j, ez], nodes[ex, j+1, ez],
              r_mat[ex, j, ez] + r_mat[ex, j+1, ez])
    c += 1


# set final K
for k in range(z-1):
    circuit.R(c, nodes[ex, ey, k], nodes[ex, ey, k+1],
              r_mat[ex, ey, k] + r_mat[ex, ey, k+1])
    c += 1


# set vin/grounds
for i in range(x):
    for j in range(y):
        circuit.R('vin' + str(i) + "-" +  str(j), 'vin',
                  nodes[i, j, 0], r_mat[i, j, 0])
        circuit.R('gnd' + str(i) + "-" + str(j), circuit.gnd,
                  nodes[i, j, ez], r_mat[i, j, ez])


print("finished setup")
# print(circuit)

simulator = circuit.simulator()
analysis = simulator.operating_point()


print("finished analysis")

final = np.zeros((x, y, z))

for node in analysis.nodes.values():
    if str(node) == 'vin':
        continue
    currP = str(node)
    final[int(currP[0]), int(currP[1]), int(
        currP[2])] = np.round(float(node), 2)

print(final)
