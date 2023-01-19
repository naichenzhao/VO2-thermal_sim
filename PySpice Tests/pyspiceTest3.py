from PySpice import *
from PySpice.Spice.Netlist import Circuit

import numpy as np

import PySpice.Logging.Logging as Logging
logger = Logging.setup_logging()


circuit = Circuit('voltage bridge')
circuit.V('input', 'vin', circuit.gnd, '10V')

x = 50
y = 15
z = 10

# Make node coordinates
nodes = np.zeros((x, y, z, 3), dtype=np.int8)

for i in range(x):
    for j in range(y):
        for k in range(z):
            nodes[i, j, k, 0] = i
            nodes[i, j, k, 1] = j
            nodes[i, j, k, 2] = k

# Make reistor matrix
r_mat = np.ones((x, y, z))*500


def get_node(nodes, i, j, k):
    return str(nodes[i,j,k,0]) + "/" + str(nodes[i,j,k,1]) + "/" + str(nodes[i,j,k,2])



print(nodes[:,:,:,0])
c = 0

# Set resistors for main array IJK
for i in range(x-1):
    for j in range(y-1):
        for k in range(z-1):
            n = get_node(nodes, i, j, k)
            circuit.R(c, n, get_node(nodes, i+1, j, k), r_mat[i, j, k] + r_mat[i+1, j, k])
            c += 1

            circuit.R(c, n, get_node(nodes, i, j+1, k),
                      r_mat[i, j, k] + r_mat[i, j+1, k])
            c += 1

            circuit.R(c, n, get_node(nodes, i, j, k+1),
                      r_mat[i, j, k] + r_mat[i, j, k+1])
            c += 1


ex = x-1
ey = y-1
ez = z-1

# set final JK
for j in range(y-1):
    for k in range(z-1):
        circuit.R(c, get_node(nodes, ex, j, k), get_node(nodes, ex, j+1, k),
                  r_mat[ex, j, k] + r_mat[ex, j+1, k])
        c += 1
        circuit.R(c, get_node(nodes, ex, j, k), get_node(nodes, ex, j, k+1),
                  r_mat[ex, j, k] + r_mat[ex, j, k+1])
        c += 1

# set final IK
for i in range(x-1):
    for k in range(z-1):
        circuit.R(c, get_node(nodes, i, ey, k), get_node(nodes, i+1, ey, k),
                  r_mat[i, ey, k] + r_mat[i+1, ey, k])
        c += 1
        circuit.R(c, get_node(nodes, i, ey, k), get_node(nodes, i, ey, k+1),
                  r_mat[i, ey, k] + r_mat[i, ey, k+1])
        c += 1

# set final IJ
for i in range(x-1):
    for j in range(y-1):
        circuit.R(c, get_node(nodes, i, j, ez), get_node(nodes, i+1, j, ez),
                  r_mat[i, j, ez] + r_mat[i+1, j, ez])
        c += 1
        circuit.R(c, get_node(nodes, i, j, ez), get_node(nodes, i, j+1, ez),
                  r_mat[i, j, ez] + r_mat[i, j+1, ez])
        c += 1


# set final I
for i in range(x-1):
    circuit.R(c, get_node(nodes, i, ey, ez), get_node(nodes, i+1, ey, ez),
              r_mat[i, ey, ez] + r_mat[i+1, ey, ez])
    c += 1

# set final J
for j in range(y-1):
    circuit.R(c, get_node(nodes, ex, j, ez), get_node(nodes, ex, j+1, ez),
              r_mat[ex, j, ez] + r_mat[ex, j+1, ez])
    c += 1


# set final K
for k in range(z-1):
    circuit.R(c, get_node(nodes, ex, ey, k), get_node(nodes, ex, ey, k+1),
              r_mat[ex, ey, k] + r_mat[ex, ey, k+1])
    c += 1


# set vin/grounds
for i in range(x):
    for j in range(y):
        circuit.R('vin' + str(i) + "-" +  str(j), 'vin',
                  get_node(nodes, i, j, 0), r_mat[i, j, 0])
        circuit.R('gnd' + str(i) + "-" + str(j), circuit.gnd,
                  get_node(nodes, i, j, ez), r_mat[i, j, ez])


print("finished setup")
# print(circuit)

simulator = circuit.simulator()
analysis = simulator.operating_point()


print("finished analysis")

final = np.zeros((x, y, z))

for node in analysis.nodes.values():
    if str(node) == 'vin':
        continue
    currP = str(node).split('/')
    final[int(currP[0]), int(currP[1]), int(
        currP[2])] = np.round(float(node), 2)

print(final)
