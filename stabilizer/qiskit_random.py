import random
import time
import os
import math

from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator

def benchmark_qiskit(qubit_test):
    simulator = AerSimulator(method='stabilizer')  # Init once

    for n_qubits in qubit_test:
        layers = (int) (math.log2(n_qubits))
        qc = QuantumCircuit(n_qubits, n_qubits)

        for _ in range(layers):
            half = n_qubits // 2
            for j in range(half):
                # Fast random Clifford layer (only S or H)
                gate = random.getrandbits(1)
                if gate:
                    qc.h(j)
                else:
                    qc.s(j)

                # CX entanglement
                qc.cx(j, j + half)

                # 20% chance of measurement
                if random.random() < 0.2:
                    qc.measure(j+half, j+half)

        # Optional: measure remaining qubits once at end
        # qc.measure_all()

        start = time.perf_counter()
        simulator.run(qc).result()
        end = time.perf_counter()
        elapsed = end - start

        print(n_qubits)
        print(f"Time: {elapsed:.6f} seconds")

        output_dir = "/people/garn195/NWQ-Sim/stabilizer/sim_bench"
        os.makedirs(output_dir, exist_ok=True)
        filename = os.path.join(output_dir, f"qiskit_{n_qubits}.txt")
        with open(filename, "w") as f:
            f.write("qiskit\n")
            f.write(f"{elapsed}\n")
            f.write(f"{n_qubits}\n")

qubit_test = []
i = 2**10
while i < (1 << 21):  # same as pow(2, 15)
    qubit_test.append(i)
    i *= 2

benchmark_qiskit(qubit_test)
