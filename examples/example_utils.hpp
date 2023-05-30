#pragma once

#include <bitset>
#include <climits>
#include <iostream>

#include <cstdlib>
#include <ctime>

#include <cmath> // For M_PI

#include <getopt.h>
#include <string>

#include "../include/state.hpp"        // Base Class
#include "../include/svsim/sv_cpu.hpp" // Derived Class
#include "../include/svsim/sv_omp.hpp" // Derived Class
#include "../include/circuit.hpp"      // Assuming this is the file that contains your Circuit class

using namespace NWQSim;

Circuit *get_adder_circuit()
{
    Circuit *circuit = new Circuit(4);
    circuit->X(0);
    circuit->X(1);
    circuit->H(3);

    circuit->CX(2, 3);
    circuit->T(0);
    circuit->T(1);
    circuit->T(2);
    circuit->TDG(3);

    circuit->CX(0, 1);
    circuit->CX(2, 3);
    circuit->CX(3, 0);
    circuit->CX(1, 2);
    circuit->CX(0, 1);
    circuit->CX(2, 3);

    circuit->TDG(0);
    circuit->TDG(1);
    circuit->TDG(2);

    circuit->T(3);
    circuit->CX(0, 1);
    circuit->CX(2, 3);
    circuit->S(3);
    circuit->CX(3, 0);
    circuit->H(3);
    return circuit;
}

Circuit *generate_random_circuit(int num_qubits, int num_operations)
{
    Circuit *circuit = new Circuit(num_qubits);

    srand(1024); // Seed random number generator with current time

    for (int i = 0; i < num_operations; i++)
    {
        int qubit = rand() % num_qubits; // Pick a random qubit index
        int gate_type = rand() % 4;      // Pick a random gate type

        switch (gate_type)
        {
        case 0:
        {
            double angle = ((double)rand() / (RAND_MAX)) * 2 * M_PI; // Generate random angle between 0 and 2π
            circuit->RZ(angle, qubit);
            break;
        }
        case 1:
            circuit->X(qubit);
            break;
        case 2:
            circuit->SX(qubit);
            break;
        case 3:
        {
            if (num_qubits > 1)
            {
                int qubit2;
                do
                {
                    qubit2 = rand() % num_qubits; // Pick a second, different random qubit index
                } while (qubit == qubit2);
                circuit->CX(qubit, qubit2);
            }
            break;
        }
        }
    }

    return circuit;
}

void parse_arguments(int argc, char **argv, bool &use_openmp, int &num_qubits, int &num_operations)
{
    int option_index = 0;
    static struct option long_options[] = {
        {"omp", no_argument, 0, 'o'},
        {0, 0, 0, 0}};

    int c;
    while ((c = getopt_long(argc, argv, "o", long_options, &option_index)) != -1)
    {
        switch (c)
        {
        case 'o':
            use_openmp = true;
            break;
        default:
            break;
        }
    }

    if (argc - optind != 2)
    {
        std::cerr << "Usage: " << argv[0] << " [--omp] num_qubits num_operations" << std::endl;
        exit(1);
    }

    num_qubits = std::stoi(argv[optind]);
    num_operations = std::stoi(argv[optind + 1]);
}