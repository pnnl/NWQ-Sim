import numpy as np
from math import pi
from qiskit import QuantumCircuit, transpile, execute
from qiskit import ClassicalRegister, QuantumRegister
from qiskit_nwqsim_provider import NWQSimProvider
from qiskit.quantum_info import SuperOp
from qiskit.providers.aer.noise import NoiseModel
from qiskit.providers.aer.noise import depolarizing_error
from qiskit.providers.aer.noise import amplitude_damping_error
from qiskit.providers.aer.noise import phase_damping_error

nwqsim = NWQSimProvider('NWQSim')
print (nwqsim.backends)
backend = nwqsim.backends['dmsim_cpu']
backend.set_n_cpus(4)

#Bell State
circuit = QuantumCircuit(2,2)
circuit.h(0)
circuit.cx(0,1)
circuit.measure([0,1], [0,1])

noise_model = NoiseModel()
#Add amplitude and phase damping error for 1-qubit
param_amplitude_damping = 0.05 
param_phase_damping = 0.03
amplitude_error = amplitude_damping_error(param_amplitude_damping)
phase_error = phase_damping_error(param_phase_damping)
qerror_q1 = phase_error.compose(amplitude_error)

#Add depolarizing error for 2-qubit
param_q2 = 0.08
qerror_q2 = depolarizing_error(param_q2, 2)
q1_superop = SuperOp(qerror_q1)
q2_superop = SuperOp(qerror_q2)
backend.set_noise_1q_sop(q1_superop)
backend.set_noise_2q_sop(q2_superop)

result = execute(circuit, backend, seed_transpiler=111).result()
print(result.get_counts(circuit))
print(result.get_statevector())
