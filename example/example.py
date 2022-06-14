import numpy as np
from math import pi
from qiskit import QuantumCircuit, transpile, execute
from qiskit import ClassicalRegister, QuantumRegister
from qiskit_nwqsim_provider import NWQSimProvider


nwqsim = NWQSimProvider('NWQSim')
print (nwqsim.backends)
backend = nwqsim.backends['svsim_cpu']
backend.set_n_cpus(4)

#q_a = QuantumRegister(3, name='qa')
#q_b = QuantumRegister(13, name='qb')
#c_a = ClassicalRegister(3)
#c_b = ClassicalRegister(13)

q_a = QuantumRegister(3, name='qa')
q_b = QuantumRegister(5, name='qb')
c_a = ClassicalRegister(3)
c_b = ClassicalRegister(5)

circuit = QuantumCircuit(q_a, q_b, c_a, c_b)

circuit.x(q_a[1])
circuit.x(q_b[1])
circuit.x(q_b[2])
circuit.x(q_b[4])
circuit.barrier()
circuit.h(q_a)
circuit.barrier(q_a)
circuit.h(q_b)
circuit.cswap(q_b[0], q_b[1], q_b[2])
circuit.cswap(q_b[2], q_b[3], q_b[4])
circuit.cswap(q_b[3], q_b[4], q_b[0])
circuit.barrier(q_b)
circuit.measure(q_a, c_a)
circuit.measure(q_b, c_b);

result = execute(circuit, backend, seed_transpiler=111).result()
print(result.get_counts(circuit))


