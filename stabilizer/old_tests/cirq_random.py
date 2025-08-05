import cirq
import random
import time
import os
import math

def benchmark_cirq(qubit_test):
    for n_qubits in qubit_test:
        layers = int(math.log2(n_qubits))
        qubits = [cirq.NamedQubit(f'q{i}') for i in range(n_qubits)]
        circuit = cirq.Circuit()

        for _ in range(layers):
            half = n_qubits // 2
            for j in range(half):
                gate = random.getrandbits(1)
                if gate:
                    circuit.append(cirq.H(qubits[j]))
                else:
                    circuit.append(cirq.S(qubits[j]))

                circuit.append(cirq.CNOT(qubits[j], qubits[j + half]))

                if random.random() < 0.2:
                    circuit.append(cirq.measure(qubits[j + half], key=f'm{j + half}'))

        # Optional final measurement (commented like in original)
        # circuit.append(cirq.measure(*qubits, key='all'))

        simulator = cirq.CliffordSimulator()
        start = time.perf_counter()
        simulator.run(circuit)
        end = time.perf_counter()
        elapsed = end - start

        print(n_qubits)
        print(f"Time: {elapsed:.6f} seconds")

        output_dir = "/people/garn195/NWQ-Sim/stabilizer/sim_bench"
        os.makedirs(output_dir, exist_ok=True)
        filename = os.path.join(output_dir, f"cirq_{n_qubits}.txt")
        with open(filename, "w") as f:
            f.write("cirq\n")
            f.write(f"{elapsed}\n")
            f.write(f"{n_qubits}\n")

qubit_test = [i**2 for i in range(4, 200, 4)]
benchmark_cirq(qubit_test)
