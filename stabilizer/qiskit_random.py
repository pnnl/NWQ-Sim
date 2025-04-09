import random
import time
import os

from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator

def benchmark_qiskit(qubit_test):
    for n_qubits in qubit_test:
        layers = n_qubits
        rng = random.Random()
        qc = QuantumCircuit(n_qubits, n_qubits)

        for _ in range(layers):
            for j in range(n_qubits):
                if rng.randint(0, 1):
                    qc.h(j)
                else:
                    qc.s(j)

            if n_qubits > 1:
                target = rng.randint(0, n_qubits - 2)
                qc.cx(target, target + 1)

            for j in range(n_qubits):
                if rng.randint(0, 19) == 0:
                    qc.measure(j, j)

        simulator = AerSimulator(method='stabilizer')

        start = time.perf_counter()
        result = simulator.run(qc).result()
        end = time.perf_counter()
        elapsed = end - start

        print(f"Time: {elapsed:.6f} seconds")

        # Write results to file
        output_dir = "/people/garn195/NWQ-Sim/stabilizer/sim_bench"
        os.makedirs(output_dir, exist_ok=True)
        filename = os.path.join(output_dir, f"qiskit_{n_qubits}.txt")
        with open(filename, "w") as f:
            f.write("qiskit\n")
            f.write(f"{elapsed}\n")
            f.write(f"{n_qubits}\n")

qubit_test = []
i = 2
while i < 2**20:
    qubit_test.append(i)
    i *= 2

benchmark_qiskit(qubit_test)
