import stim
import random
import time
import os
import math

def benchmark_stim(qubit_test):
    for n_qubits in qubit_test:
        layers = (int) (math.log2(n_qubits))
        circuit = stim.Circuit()

        half = n_qubits // 2
        for _ in range(layers):
            for j in range(half):
                # Random H or S
                if random.getrandbits(1):
                    circuit.append("H", [j])
                else:
                    circuit.append("S", [j])

                # Entangle with CX
                circuit.append("CX", [j, j + half])

                # 20% chance of measurement
                if random.random() < 0.2:
                    circuit.append("M", [j+half])

        # Optionally: final measurements
        # circuit.append("M", range(n_qubits))
        simulator = stim.TableauSimulator()
        start = time.perf_counter()
        simulator.do_circuit(circuit)
        end = time.perf_counter()
        elapsed = end - start

        print(n_qubits)
        print(f"Stim Time: {elapsed:.6f} seconds")

        # Output
        output_dir = "/people/garn195/NWQ-Sim/stabilizer/sim_bench"
        os.makedirs(output_dir, exist_ok=True)
        filename = os.path.join(output_dir, f"stim_{n_qubits}.txt")
        with open(filename, "w") as f:
            f.write("stim\n")
            f.write(f"{elapsed}\n")
            f.write(f"{n_qubits}\n")

qubit_test = []
i = 2**16
while i < (1 << 21):
    qubit_test.append(i)
    i *= 2

benchmark_stim(qubit_test)
