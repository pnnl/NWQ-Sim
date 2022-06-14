import numpy as np
from math import pi
from qiskit import QuantumCircuit, transpile, execute
from qiskit import ClassicalRegister, QuantumRegister
from qiskit.tools.visualization import plot_histogram
from qiskit_nwqsim_provider import NWQSimProvider

nwqsim = NWQSimProvider('NWQSim')
print (nwqsim.backends)
sim = nwqsim.backends['dmsim_cpu']

circ = QuantumCircuit(3,3)
circ.h(0)
circ.cx(0,1)
circ.cx(1,2)
circ.measure([0,1,2], [0,1,2])

result = sim.run(transpile(circ,sim)).result()
counts = result.get_counts(0)
print (counts)
plot_histogram(counts, title="ideal counts for 3-qubit ghz state")
