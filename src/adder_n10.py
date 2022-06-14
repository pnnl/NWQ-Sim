# ---------------------------------------------------------------------------
# NWQsim: Northwest Quantum Circuit Simulation Environment
# ---------------------------------------------------------------------------
# Ang Li, Senior Computer Scientist
# Pacific Northwest National Laboratory(PNNL), U.S.
# Homepage: http://www.angliphd.com
# GitHub repo: http://www.github.com/pnnl/DM-Sim
# PNNL-IPID: 31919-E, ECCN: EAR99, IR: PNNL-SA-143160
# BSD Lincese.
# ---------------------------------------------------------------------------
# File: adder_n10_omp.py
# A 10-qubit adder example using Python API.
# Single-node multiple NVIDIA GPUs; No inter-GPU communication required.
# Requires: PyBind11 (https://github.com/pybind/pybind11)
#           CUDA-10.0 or newer (required by pybind11 for Python API)
# ---------------------------------------------------------------------------
import numpy as np
import sys
from qiskit_nwqsim_provider import libdmsim as dmsim

## Call via: $python circuit.py num_of_qubits num_of_gpus
if (len(sys.argv) != 2):
    print("Call using $python circuit.py n_qubits\n")
    exit()

n_qubits = int(sys.argv[1])

## Create simulator object
sim = dmsim.Simulation(n_qubits)

## Quantum ripple-carry adder from Cuccaro et al, quant-ph/0410184
## Define circuit module functions as below
def majority(sim, a, b, c):
    sim.CX(c, b)
    sim.CX(c, a)
    sim.CCX(a, b, c)

def unmaj(sim, a, b, c):
    sim.CCX(a, b, c)
    sim.CX(c, a)
    sim.CX(a, b)

## Add the gates to the circuit
sim.X(1)
sim.X(5)
sim.X(6)
sim.X(7)
sim.X(8)
## Add a new circuit for the current density matrix
majority(sim, 0, 5, 1)
majority(sim, 1, 6, 2)
majority(sim, 2, 7, 3)
majority(sim, 3, 8, 4)
sim.CX(4, 9)
unmaj(sim, 3, 8, 4)
unmaj(sim, 2, 7, 3)
unmaj(sim, 1, 6, 2)
unmaj(sim, 0, 5, 1)

#sim.set_noise_1q_sop(np.ones(16))
#sim.set_noise_1q_sop(np.ones(16))

res = sim.measure_all(10)
## Print measurement results
print ("\n===============  Measurement (tests=" + str(len(res)) + ") ================")
for i in range(len(res)):
    print ("Test-"+str(i)+": " + "{0:b}".format(res[i]).zfill(n_qubits))
