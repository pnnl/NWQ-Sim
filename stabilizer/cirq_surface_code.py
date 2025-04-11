import cirq
import time

def create_surface_code_circuit(distance, rounds):
    grid_size = (2 * distance) - 1
    qubits = [[cirq.GridQubit(i, j) for j in range(grid_size)] for i in range(grid_size)]
    circuit = cirq.Circuit()

    def q(i, j):
        return qubits[i][j]

    for _ in range(rounds):
        # X stabilizers (odd rows, even columns)
        for i in range(grid_size):
            for j in range(grid_size):
                if (i % 2 == 1) and (j % 2 == 0):
                    ancilla = q(i, j)
                    data_neighbors = [(i+1, j), (i-1, j), (i, j+1), (i, j-1)]
                    data_neighbors = [q(x, y) for x, y in data_neighbors if 0 <= x < grid_size and 0 <= y < grid_size]

                    circuit.append(cirq.H(ancilla))
                    for data in data_neighbors:
                        circuit.append(cirq.CNOT(ancilla, data))
                    circuit.append(cirq.H(ancilla))
                    circuit.append(cirq.measure(ancilla, key=f"m_{i}_{j}"))

        # Z stabilizers (even rows, odd columns)
        for i in range(grid_size):
            for j in range(grid_size):
                if (i % 2 == 0) and (j % 2 == 1):
                    ancilla = q(i, j)
                    data_neighbors = [(i, j+1), (i, j-1), (i+1, j), (i-1, j)]
                    data_neighbors = [q(x, y) for x, y in data_neighbors if 0 <= x < grid_size and 0 <= y < grid_size]

                    for data in data_neighbors:
                        circuit.append(cirq.CNOT(data, ancilla))
                    circuit.append(cirq.measure(ancilla, key=f"m_{i}_{j}"))

    return circuit, [q(i, j) for i in range(grid_size) for j in range(grid_size)]

for d in range(3, 152, 2):
    distance = d
    rounds = 1
    circuit, all_qubits = create_surface_code_circuit(distance, rounds)
    simulator = cirq.CliffordSimulator()

    start = time.perf_counter()
    simulator.run(circuit)
    end = time.perf_counter()

    filename = f"/people/garn195/NWQ-Sim/stabilizer/fowler_surface_code/cirq_{distance}.txt"
    with open(filename, "w") as file:
        file.write("cirq\n")
        file.write(f"{end - start}\n")
        file.write(f"{distance}\n")
        file.write(f"{rounds}\n")
        file.write(f"{(2 * distance - 1) ** 2}\n")
