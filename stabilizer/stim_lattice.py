import stim
import random
import time
import os
import math

folder_path = "/people/garn195/NWQ-Sim/stabilizer/stim_to_stim_files"
# Ensure the folder exists
if not os.path.isdir(folder_path):
    print(f"The folder '{folder_path}' does not exist.")

# List all files in the folder
stim_files = [f for f in os.listdir(folder_path) if f.endswith('.stim')]

if not stim_files:
    print("No .stim files found in the folder.")

# Loop through each .stim file and simulate
for stim_file in stim_files:
    stim_file_path = os.path.join(folder_path, stim_file)
    print(f"Simulating: {stim_file_path}")

    try:
        # Read the stim file
        with open(stim_file_path, 'r') as file:
            circuit_text = file.read()

        circuit = stim.Circuit(circuit_text)
        n_qubits = circuit.num_qubits
        simulator = stim.TableauSimulator()

        start = time.perf_counter()
        simulator.do_circuit(circuit)
        end = time.perf_counter()
        elapsed = end - start

        output_dir = "/people/garn195/NWQ-Sim/stabilizer/surface_operation_bench"
        os.makedirs(output_dir, exist_ok=True)
        filename = os.path.join(output_dir, f"stim_{n_qubits}.txt")
        with open(filename, "w") as f:
            f.write("stim\n")
            f.write(f"{elapsed}\n")
            f.write(f"{n_qubits}\n")

    except Exception as e:
        print(f"Error simulating {stim_file}: {e}")

