# ----------------------------------------------------------------------
# NWQBench: Northwest Quantum Proxy Application Suite 
# ----------------------------------------------------------------------
# Ang Li, Samuel Stein, James Ang.
# Pacific Northwest National Laboratory(PNNL), U.S.
# BSD Lincese.
# Created 05/21/2021.
# ----------------------------------------------------------------------


import numpy as np
from math import pi
from qiskit import QuantumCircuit
from qiskit import execute, Aer
from qiskit_nwqsim_provider import NWQSimProvider
import sys
import math
import random

#n_qubits = int(sys.argv[1])
n_qubits = 10

def cx_chain(qc,n):
    for i in range(0,n-1):
        qc.cx(i,i+1)
        #qc.rx(pi/float(i+1),i)



qc = QuantumCircuit(n_qubits)
qc.h(0)
cx_chain(qc,n_qubits)

qc.measure_all()
#qasm_file = open("qft_n" + str(n_qubits) + ".qasm","w")
#qasm_file.write(qc.qasm())
#qasm_file.close()

#print (qc)

#simulator = Aer.get_backend('statevector_simulator')
#job1 = execute(qc,simulator,shots=1000)
#result1 = job1.result()
#counts1 = result1.get_counts(qc)
#print (counts1)
#print (result.get_statevector(qc))

nwqsim = NWQSimProvider('NWSimSimulator')
svsim = nwqsim.backends['nwqsim']

#print (qc.qasm())

# default number of shots: 1024
print("Default number of shots: 1024")
job2 = svsim.run(qc)
result2 = job2.result()
counts2 = result2.get_counts(qc)
print (counts2)
print (result2.get_statevector(qc))

# custom number of shots
print("\nCustom number of shots: 5000")
job3 = svsim.run(qc, shots=5000)
result3 = job3.result()
counts3 = result3.get_counts(qc)
print (counts3)
print (result3.get_statevector(qc))
print ("Are two statevectors the same?", np.allclose(result2.get_statevector(qc), result3.get_statevector(qc)))

