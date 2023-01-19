from PySpice import *
from PySpice.Spice.Netlist import Circuit

import numpy as np

import PySpice.Logging.Logging as Logging
logger = Logging.setup_logging()



circuit = Circuit('voltage bridge')
circuit.V('input', 'vin', circuit.gnd, '10V')

# Make node coordinates
linex = range(0, 4) # x
liney = range(0, 7)  # y
nodes = np.array([[str(i) + str(j) for j in linex] for i in liney])

x = nodes.shape[0]
y = nodes.shape[1]

# Make reistor matrix
r_mat = np.ones((x, y))*500

print(nodes)
print(r_mat)

c = 0;

print(x, y)

# Set resistors for main array
for i in range(x-1):
    for j in range(y-1):
        n = nodes[i, j]
        circuit.R(c, n, nodes[i, j+1], r_mat[i, j] + r_mat[i, j+1])
        c += 1

        circuit.R(c, n, nodes[i+1, j], r_mat[i, j] + r_mat[i+1, j])
        c += 1
        
# set final row
for j in range(y-1):
    circuit.R(c, nodes[x-1, j], nodes[x-1, j+1], r_mat[x-1, j] + r_mat[x-1, j+1])
    c += 1

# set final column
for i in range(x-1):
    circuit.R(c, nodes[i, y-1], nodes[i+1, y-1], r_mat[i, y-1] + r_mat[i+1, y-1])
    c += 1

# set vin/grounds
for i in range(x):
    circuit.R('vin' + str(i), 'vin', nodes[i, 0], r_mat[i, 1])
    circuit.R('gnd' + str(i), circuit.gnd, nodes[i, y-1], r_mat[i, y-1])


print("finished setup")
# print(circuit)

simulator = circuit.simulator()
analysis = simulator.operating_point()


print("finished analysis")

final = np.zeros((x, y))

for node in analysis.nodes.values():
    if str(node) == 'vin':
        continue
    currP = str(node)
    final[int(currP[0]), int(currP[1])] = np.round(float(node), 2)

print(final)


