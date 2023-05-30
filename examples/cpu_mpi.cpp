#include <mpi.h>
#include <iostream>

#include "../include/svsim/sv_mpi.hpp" // Derived Class
#include "../include/circuit.hpp"      // Assuming this is the file that contains your Circuit class

#include "example_utils.hpp"

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);

    int num_qubits, num_operations;
    bool use_openmp = false;

    parse_arguments(argc, argv, use_openmp, num_qubits, num_operations);

    auto circuit = generate_random_circuit(num_qubits, num_operations);

    QuantumState *state = new SV_MPI(num_qubits);
    state->sim(circuit);

    delete state;
    delete circuit;

    MPI_Finalize();
    return 0;
}
