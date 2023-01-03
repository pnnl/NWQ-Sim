# ----------------------------------------------------------------------
# NWQBench: Northwest Quantum Proxy Application Suite 
# ----------------------------------------------------------------------
# Ang Li, Samuel Stein, James Ang.
# Pacific Northwest National Laboratory(PNNL), U.S.
# BSD Lincese.
# Created 05/21/2021.
# ----------------------------------------------------------------------


import numpy as np
from qiskit import QuantumCircuit
from qiskit import execute, Aer
from qiskit_nwqsim_provider import NWQSimProvider
import sys
import math
import random
from mpi4py import MPI

comm = MPI.COMM_WORLD

#n_qubits = int(sys.argv[1])
n_qubits = 10

def cx_chain(qc,n):
    for i in range(0,n-1):
        qc.cx(i,i+1)



qc = QuantumCircuit(n_qubits, n_qubits)
qc.h(0)
cx_chain(qc,n_qubits)
qc.measure_all()
#qasm_file = open("qft_n" + str(n_qubits) + ".qasm","w")
#qasm_file.write(qc.qasm())
#qasm_file.close()

#print (qc)

simulator = Aer.get_backend('statevector_simulator')
job1 = execute(qc,simulator,shots=1000)
result1 = job1.result()
counts1 = result1.get_counts(qc)
print (counts1)
#print (result.get_statevector(qc))

nwqsim = NWQSimProvider('NWSimSimulator')
svsim = nwqsim.backends['svsim_gpu']
job2 = execute(qc,svsim,shots=1000)
result2 = job2.result()
counts2 = result2.get_counts(qc)
print (counts2)
#print (result.get_statevector(qc))


