import stim
import time
import random

def benchmark_stim(n_qubits):
    for each in n_qubits:
        circuit = stim.Circuit()
        n_repeats = each

        for _ in range(n_repeats):
            for cntrl in range(n_repeats):
                gate = random.randint(0, 1)
                if gate:
                    circuit.append_operation("H", [cntrl])
                else:
                    circuit.append_operation("S", [cntrl])

                target = random.randint(0, each - 2)
                if target == cntrl:
                    target += 1  # ensure target != control

                circuit.append_operation("CX", [cntrl, target])
                circuit.append_operation("M", [target])

        simulator = stim.TableauSimulator()

        start = time.perf_counter()
        simulator.do_circuit(circuit)
        end = time.perf_counter()

        print("Time", end - start)

        filename = f"/Users/garn195/Project Repositories/NWQ-Sim/stabilizer/sim_bench/stim_{each}.txt"
        with open(filename, "w") as file:
            file.write("stim\n")
            file.write(f"{end - start}\n")
            file.write(f"{each}\n")

qubit_test = []
i = 2
while i < 16:
    qubit_test.append(i)
    i *= 2

benchmark_stim(qubit_test)
