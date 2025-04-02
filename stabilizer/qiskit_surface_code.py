from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
import time

def create_surface_code_circuit(distance, rounds):
    grid_size = distance + 1
    circuit = QuantumCircuit(grid_size ** 2, grid_size ** 2)
    
    def q(i, j):
        return i * grid_size + j
    
    for _ in range(rounds):
        for i in range(distance):
            for j in range(distance):
                if (i + j) % 2 == 1:
                    ancilla = q(i, j)
                    data_neighbors = [(i+1, j), (i, j+1), (i-1, j), (i, j-1)]
                    data_neighbors = [q(x, y) for x, y in data_neighbors if 0 <= x < grid_size and 0 <= y < grid_size]
                    
                    circuit.h(ancilla)
                    for data in data_neighbors:
                        circuit.cx(ancilla, data)
                    circuit.h(ancilla)
                    circuit.measure(ancilla, ancilla)
        
        for i in range(distance):
            for j in range(distance):
                if (i + j) % 2 == 0:
                    ancilla = q(i, j)
                    data_neighbors = [(i+1, j), (i, j+1), (i-1, j), (i, j-1)]
                    data_neighbors = [q(x, y) for x, y in data_neighbors if 0 <= x < grid_size and 0 <= y < grid_size]
                    
                    for data in data_neighbors:
                        circuit.cx(data, ancilla)
                    circuit.measure(ancilla, ancilla)
    
    return circuit

for d in range(1, 100, 2):
    distance = d
    rounds = 1
    circuit = create_surface_code_circuit(distance, rounds)
    simulator = AerSimulator(method="stabilizer")
    
    start = time.perf_counter()
    simulator.run(circuit).result()
    end = time.perf_counter()
    
    filename = f"/people/garn195/NWQ-Sim/stabilizer/surface_code_hpc/qiskit_{distance}.txt"
    with open(filename, "w") as file:
        file.write("qiskit\n")
        file.write(f"{end-start}\n")
        file.write(f"{distance}\n")
        file.write(f"{rounds}\n")
        file.write(f"{(2 * distance**2) + 1}\n")
