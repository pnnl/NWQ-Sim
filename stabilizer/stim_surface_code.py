import stim
import time
import random
import numpy as np

for d in range(97, 121, 2):
    circuit = stim.Circuit()

    distance = d
    rounds = 1

    n = 2 * pow(distance, 2) + 1

    grid_size = distance + 1

    def q(i, j):
        return i * grid_size + j

    for i in range(distance):
        for j in range(distance):
            if (i + j) % 2 == 1:
                ancilla = q(i, j)
                data_neighbors = [q(i+1, j), q(i, j+1), q(i-1, j), q(i, j-1)]
                data_neighbors = [n for n in data_neighbors if 0 <= n < grid_size**2]

                circuit.append("H", [ancilla])
                for data in data_neighbors:
                    circuit.append("CX", [ancilla, data])
                circuit.append("H", [ancilla]) 
                circuit.append("M", [ancilla]) 

    for i in range(distance):
        for j in range(distance):
            if (i + j) % 2 == 0: 
                ancilla = q(i, j)
                data_neighbors = [q(i+1, j), q(i, j+1), q(i-1, j), q(i, j-1)]
                data_neighbors = [n for n in data_neighbors if 0 <= n < grid_size**2]

                for data in data_neighbors:
                    circuit.append("CX", [data, ancilla])
                circuit.append("M", [ancilla])

    #Inject an X error on a data qubit
    #circuit.append("X_ERROR(0.1)", [q(1, 1)]) 

    # samples = circuit.sample(num_samples=5)
    # print("\n Measurement Samples:")
    # print(len(circuit))

    simulator = stim.TableauSimulator()
    start = time.perf_counter()
    simulator.do(circuit)
    end = time.perf_counter()


    filename = "/people/garn195/NWQ-Sim/stabilizer/surface_code_hpc/stim_"+ (str)(distance) + ".txt"
    with open(filename, "w") as file:
        file.write("stim\n")
        file.write((str)(end-start)+"\n")
        file.write((str)(distance)+"\n")
        file.write((str)(rounds)+"\n")
        file.write((str)(n)+"\n")