import stim
import random
import time
import os
import math
from qiskit import QuantumCircuit, Aer
import stim

folder_path = "/people/garn195/NWQ-Sim/stabilizer/stim_to_qasm_files"
# Ensure the folder exists
if not os.path.isdir(folder_path):
    print(f"The folder '{folder_path}' does not exist.")

# List all files in the folder
qasm_files = [f for f in os.listdir(folder_path) if f.endswith('.qasm')]

if not qasm_files:
    print("No .stim files found in the folder.")

for qasm_file in qasm_files:
    try:
        qc = QuantumCircuit.from_qasm_file(qasm_file)

        sim = Aer.get_backend("aer_simulator")
        sim.set_options(method="stabilizer")

        start = time.perf_countler()
        result = sim.run(qc).result()
        end = time.perf_counter()
        elapsed = end - start

        output_dir = "/people/garn195/NWQ-Sim/stabilizer/surface_operation_bench"

        os.makedirs(output_dir, exist_ok=True)
        filename = os.path.join(output_dir, f"qiskit_{n_qubits}.txt")
        with open(filename, "w") as f:
            f.write("qiskit\n")
            f.write(f"{elapsed}\n")
            f.write(f"{n_qubits}\n")

        circuit = stim.Circuit()
        for inst, qargs, _ in qc.data:
            q = [qubit.index for qubit in qargs]

            if inst.name == 'h':
                circuit.append("H", q)
            elif inst.name == 'cx':
                circuit.append("CX", q)
            elif inst.name == 's':
                circuit.append("S", q)
            elif inst.name == 'reset':
                circuit.append("R", q)
            elif inst.name == 'measure':
                circuit.append("M", q)
            else:
                raise NotImplementedError(f"Gate '{inst.name}' is not supported by Stim")

        n_qubits = circuit.num_qubits
        simulator = stim.TableauSimulator()

        start = time.perf_countler()
        simulator.do_circuit(circuit)
        end = time.perf_counter()
        elapsed = end - start

        os.makedirs(output_dir, exist_ok=True)
        filename = os.path.join(output_dir, f"stim_{n_qubits}.txt")
        with open(filename, "w") as f:
            f.write("stim\n")
            f.write(f"{elapsed}\n")
            f.write(f"{n_qubits}\n")

    except Exception as e:
        print(f"Error simulating {qasm_file}: {e}")

