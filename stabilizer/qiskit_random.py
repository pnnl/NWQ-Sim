import qiskit
import time
import random
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator

def benchmark_qiskit(n_qubits):
    for each in n_qubits:
        circuit = QuantumCircuit(each, each)
        n_repeats = 1000

        for _ in range(n_repeats):
            for cntrl in range(each):
                gate = random.randint(0, 1)
                if gate:
                    circuit.h(cntrl)
                else:
                    circuit.s(cntrl)
            
                target = random.randint(0, each - 2)
                if target == cntrl:
                    target += 1

                circuit.cx(cntrl, target)

                
                circuit.measure(target, target) 
        
        simulator = AerSimulator(method='stabilizer')

        start = time.perf_counter()
        simulator.run(circuit).result()
        end = time.perf_counter()
        
        print("Time", end - start)
        
        filename = f"/Users/garn195/Project Repositories/NWQ-Sim/stabilizer/sim_bench/qiskit_{each}.txt"
        with open(filename, "w") as file:
            file.write("qiskit\n")
            file.write(f"{end - start}\n")
            file.write(f"{each}\n")

qubit_test = []
i = 2
while i < pow(2, 15):
    qubit_test.append(i)
    i *= 2
benchmark_qiskit(qubit_test)

