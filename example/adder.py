import numpy as np
from qiskit import QuantumCircuit
#from qiskit import execute, Aer
from qiskit import QuantumCircuit, transpile
from qiskit_nwqsim_provider import NWQSimProvider
import sys

def majority(qc,a,b,c):
    qc.cx(c,b)
    qc.cx(c,a)
    qc.ccx(a,b,c)

def unmaj(qc,a,b,c):
    qc.ccx(a,b,c)
    qc.cx(c,a)
    qc.cx(a,b)

def adder4(qc,a0,a1,a2,a3,b0,b1,b2,b3,cin,cout):
    majority(qc,cin,b0,a0)
    majority(qc,a0,b1,a1)
    majority(qc,a1,b2,a2)
    majority(qc,a2,b3,a3)
    qc.cx(a3,cout)
    unmaj(qc,a2,b3,a3)
    unmaj(qc,a1,b2,a2)
    unmaj(qc,a0,b1,a1)
    unmaj(qc,cin,b0,a0)


n_bits = int(sys.argv[1])
if n_bits %  4 != 0 or n_bits <= 0:
    print ("Number of adder bits should be a multiply of 4.\n")
    exit(0)

n_qubits = n_bits * 2 + 2 + int(n_bits/4) - 1
qc = QuantumCircuit(n_qubits, n_qubits)


#===================== Initialization ====================
#we compute a=1110, b=0001, cin=1 => a=0000,b=0001,cout=1

#a[0:n_bits-1]=1110
for i in range(1,n_bits):
    qc.x(i)

#b=[n_bits:2*n_bits-1]=0001
qc.x(n_bits)

#cin[2*n_bits] = 1
qc.x(2*n_bits)

#===================== Adder ====================

for i in range(0,n_bits,4):
    adder4(qc,i,i+1,i+2,i+3,i+n_bits,i+n_bits+1,i+n_bits+2,i+n_bits+3,n_bits*2+int(i/4),n_bits*2+int(i/4)+1)

qc.measure_all()
#qasm_file = open("adder_n" + str(n_qubits) + ".qasm","w")
#qasm_file.write(qc.qasm())
#qasm_file.close()

#simulator = Aer.get_backend('qasm_simulator')
#job = execute(qc,simulator,shots=10)
#result = job.result()
#counts = result.get_counts(qc)
#print (counts)

nwqsim = NWQSimProvider('NWQSim')
backend = nwqsim.backends['svsim_cpu']
backend.set_n_cpus(8)
trans_qc = transpile(qc, backend)
job = backend.run(trans_qc)
print(job.get_counts())



