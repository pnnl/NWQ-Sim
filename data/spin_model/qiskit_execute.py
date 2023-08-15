import qiskit
from qiskit import QuantumCircuit, transpile, Aer
import sys
import time

def execute_on_gpu(qasm_filename, shots):
    # Load the QASM file into a Quantum Circuit
    qc = QuantumCircuit.from_qasm_file(qasm_filename)
    qc.measure_all() 
    simulator_gpu = Aer.get_backend('aer_simulator_statevector_gpu',  cuStateVec_enable=True)
    #simulator_gpu.set_options(device='GPU')
    
    # Record the start time
    start_time = time.time()

    # Execute the circuit on the simulator
    result = simulator_gpu.run(qc, shots=shots, fusion_enable=False, optimization_level=0).result()

    # Record the end time
    end_time = time.time()

    # Calculate the execution time
    execution_time = end_time - start_time
    print(result) 


    print(f"Execution Time: {execution_time:.4f} seconds, qiskit_reported: {result.to_dict()['time_taken']:.4f} seconds")

    return result

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Please provide the QASM filename as an argument.")
        sys.exit(1)
    
    shots = int(sys.argv[1])
    index = int(sys.argv[2])
    
    qasm_filename = f"spin_model_{index}.qasm"
    result = execute_on_gpu(qasm_filename, shots)
    


