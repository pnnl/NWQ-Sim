import stim
import time
import random
import os

def benchmark_stim(qubit_test):
    for n_qubits in qubit_test:
        print(f"\nStarting program for {n_qubits} qubits")
        layers = n_qubits
        circuit = stim.Circuit()
        rng = random.Random()

        for k in range(layers):
            for j in range(n_qubits):
                if rng.randint(0, 1):
                    circuit.append("H", [j])
                else:
                    circuit.append("S", [j])

            if n_qubits > 1:
                target = rng.randint(0, n_qubits - 2)
                circuit.append("CX", [target, target + 1])

            for j in range(n_qubits):
                if rng.randint(0, 19) == 0:
                    circuit.append("M", [j])

        simulator = stim.TableauSimulator()

        start = time.perf_counter()
        simulator.do_circuit(circuit)
        end = time.perf_counter()
        elapsed = end - start

        print(f"Time: {elapsed:.6f} seconds")

        output_dir = "/people/garn195/NWQ-Sim/stabilizer/sim_bench"
        os.makedirs(output_dir, exist_ok=True)
        filename = os.path.join(output_dir, f"stim_{n_qubits}.txt")
        with open(filename, "w") as f:
            f.write("stim\n")
            f.write(f"{elapsed}\n")
            f.write(f"{n_qubits}\n")

qubit_test = []
i = 2
while i < 2**20:
    qubit_test.append(i)
    i *= 2

benchmark_stim(qubit_test)
