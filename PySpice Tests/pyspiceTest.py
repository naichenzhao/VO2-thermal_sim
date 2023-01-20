from PySpice.Unit import *
from PySpice.Spice.Netlist import Circuit
import PySpice.Logging.Logging as Logging
logger = Logging.setup_logging()


MAX_LEN = 5

circuit = Circuit('voltage bridge')
circuit.V('input', circuit.gnd, 1, '10V')


for i in range(MAX_LEN):
    circuit.R(i, 1, circuit.gnd, kilo(1))
    circuit["R" + str(i)].resistance = 100



for i in range(MAX_LEN):
    curr = circuit["R" + str(i)]
    curr.minus.add_current_probe(circuit)  # to get positive value

print("finished setup")

simulator = circuit.simulator()
analysis = simulator.operating_point()

print("finished analysis")

for node in analysis.branches.values():
    print('Node {}: {:5.3f} A'.format(str(node), float(node)))
