import sys
import os
import time
import numpy as np
from qiskit_aer import AerSimulator
from qiskit.qasm2 import load
from qiskit.providers import JobTimeoutError

def benchmark_random_qiskit():
    """Benchmark random circuits using Qiskit with stabilizer method"""
    
    # Input and output directories
    input_dir = "/people/garn195/NWQ-Sim/stabilizer/random_bench_dense/"
    output_dir = "/people/garn195/NWQ-Sim/stabilizer/random_bench_dense/"
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize Qiskit stabilizer simulator
    simulator = AerSimulator(method='stabilizer')
    
    for qasm_file in os.listdir(input_dir):
        if qasm_file.endswith(".qasm") and qasm_file.startswith("random_"):
            qasm_filename = os.path.join(input_dir, qasm_file)
            try:
                print(f"Processing: {qasm_filename}")
                
                # Load circuit from QASM
                circuit = load(qasm_filename)
                
                # Benchmark Qiskit simulation
                start_time = time.time()
                
                # Run simulation with a 30-minute timeout
                job = simulator.run(circuit)
                job.result(timeout=1800)
                
                end_time = time.time()
                sim_time = end_time - start_time
                
                print(f"Sim time: {sim_time:.6f}s")
                
                # Extract details from filename for output
                parts = qasm_file.replace('.qasm', '').split('_')
                num_qubits = int(parts[1][1:])
                rounds = int(parts[2][1:])

    for qasm_file in qasm_files:
        # Extract number of qubits and skip if not > 22000
        # num_qubits = get_qubits_from_filename(qasm_file)
        # if num_qubits <= 22000:
        #     continue

        qasm_filename = os.path.join(input_dir, qasm_file)
        try:
            print(f"Processing: {qasm_filename}")
            
            # Load circuit from QASM
            circuit = load(qasm_filename)
            
            # Benchmark Qiskit simulation
            start_time = time.time()
            
            # Run simulation with a 30-minute timeout
            job = simulator.run(circuit)
            job.result(timeout=1200)
            
            end_time = time.time()
            sim_time = end_time - start_time
            
            print(f"Sim time: {sim_time:.6f}s")
            
            # Extract details from filename for output
            parts = qasm_file.replace('.qasm', '').split('_')
            # num_qubits is already extracted above
            rounds = int(parts[2][1:])

            # Write benchmark results
            output_filename = f"{output_dir}qiskit_random_q{circuit.num_qubits}_r{rounds}.txt"
            with open(output_filename, "w") as outfile:
                outfile.write("qiskit_random\n")
                outfile.write(f"{sim_time:.6f}\n")
                outfile.write(f"{circuit.num_qubits}\n")
                outfile.write(f"{len(circuit.data)}\n")  # Number of operations
                outfile.write(f"{rounds}\n")
        
        except JobTimeoutError:
            print(f"Timeout processing {qasm_filename} after 30 minutes.")
            continue
        except Exception as e:
            print(f"Error processing {qasm_filename}: {e}")
            continue

if __name__ == "__main__":
    benchmark_random_qiskit()
