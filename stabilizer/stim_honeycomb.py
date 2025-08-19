import sys
import os
import stim
import importlib
import time
import numpy as np

# Add the src directory to path
sys.path.insert(0, os.path.abspath('src'))

import honeycomb_code as hc
importlib.reload(hc)

def benchmark_honeycomb_stim():
    """Benchmark honeycomb code circuits using Stim"""
    
    # Create output directory if it doesn't exist

    output_dir = "/people/garn195/NWQ-Sim/stabilizer/honeycomb_bench/"

    os.makedirs(output_dir, exist_ok=True)
    
    # Test parameters - only varying distance
    distances = [3, 5, 7, 9, 11, 13]
    rounds = 1
    injection = '0'
    logical_basis = 'z'
    
    for distance in distances:
        try:
            print(f"Processing: d={distance}")
            
            # Generate honeycomb circuit
            code = hc.HoneycombQubits()
            code.select_active(distance=distance)
            
            circuit = code.inject_and_measure(
                rounds=rounds, 
                injection=injection, 
                logical_basis=logical_basis, 
                detector_reset=False
            )
            
            # Convert to QASM and save
            qasm_str = circuit.to_qasm(open_qasm_version=2, skip_dets_and_obs=True)
            qasm_filename = f"{output_dir}honeycomb_d{distance}.qasm"
            with open(qasm_filename, "w") as f:
                f.write(qasm_str)
            
            # Benchmark Stim simulation
            start_time = time.time()
            
            # Compile and run simulation
            sampler = circuit.compile_detector_sampler()
            shots = 1000
            dets, observables = sampler.sample(
                shots=shots,
                separate_observables=True,
                bit_packed=True,
            )
            
            end_time = time.time()
            sim_time = end_time - start_time
            
            print(f"Sim time: {sim_time:.6f}s")
            
            # Write benchmark results
            filename = f"{output_dir}stim_honeycomb_d{distance}.txt"
            with open(filename, "w") as outfile:
                outfile.write("stim\n")
                outfile.write(f"{sim_time:.6f}\n")
                outfile.write(f"{circuit.num_qubits}\n")
                outfile.write(f"{len(circuit)}\n")  # Number of operations
                outfile.write(f"{distance}\n")      # Distance parameter
            
        except Exception as e:
            print(f"Error processing d={distance}: {e}")
            continue

if __name__ == "__main__":
    benchmark_honeycomb_stim()