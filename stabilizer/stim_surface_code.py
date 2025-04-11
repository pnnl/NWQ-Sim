import stim
import time
import random
import numpy as np

for d in range(63, 152, 2):
    circuit = stim.Circuit()

    distance = d
    rounds = 1
    grid_size = (2 * distance) - 1  
    n = pow(grid_size, 2)

    def q(i, j):
        return i * grid_size + j

    for i in range(grid_size):
        for j in range(grid_size):
            if ((i%2 == 1) and (j%2 == 0)):
                ancilla = q(i, j)
                data_neighbors = [(i+1, j), (i-1, j), (i, j+1), (i, j-1)]
                data_neighbors = [q(x, y) for x, y in data_neighbors if 0 <= x < grid_size and 0 <= y < grid_size]

                circuit.append("H", [ancilla])
                for data in data_neighbors:
                    circuit.append("CX", [ancilla, data])
                circuit.append("H", [ancilla]) 
                circuit.append("M", [ancilla]) 

    for i in range(grid_size):
        for j in range(grid_size):
            if ((i%2 == 0) and (j%2 == 1)):
                ancilla = q(i, j)
                data_neighbors = [(i+1, j), (i-1, j), (i, j+1), (i, j-1)]
                data_neighbors = [q(x, y) for x, y in data_neighbors if 0 <= x < grid_size and 0 <= y < grid_size]

                for data in data_neighbors:
                    circuit.append("CX", [data, ancilla])
                circuit.append("M", [ancilla])

    simulator = stim.TableauSimulator()
    start = time.perf_counter()
    simulator.do(circuit)
    end = time.perf_counter()


    filename = "/people/garn195/NWQ-Sim/stabilizer/fowler_surface_code/stim_"+ (str)(distance) + ".txt"
    with open(filename, "w") as file:
        file.write("stim\n")
        file.write((str)(end-start)+"\n")
        file.write((str)(distance)+"\n")
        file.write((str)(rounds)+"\n")
        file.write((str)(n)+"\n")