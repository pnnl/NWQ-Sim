#include <cassert>
#include <map>
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

#include "example_utils.hpp"

using namespace NWQSim;

int main(int argc, char **argv)
{
    int num_qubits, num_operations;
    bool use_openmp = false;

    parse_arguments(argc, argv, use_openmp, num_qubits, num_operations);

    auto circuit = generate_random_circuit(num_qubits, num_operations);

    QuantumState *state;

    if (use_openmp)
        state = new SV_OMP(num_qubits);
    else
        state = new SV_CPU(num_qubits);

    state->sim(circuit);

    // IdxType repetition = 1024;
    // auto res = state->measure_all(repetition);

    // std::map<IdxType, IdxType> result_dict;

    // for (IdxType i = 0; i < repetition; i++)
    // {
    //     if (result_dict.find(res[i]) != result_dict.end())
    //         result_dict[res[i]] += 1;
    //     else
    //         result_dict.insert({res[i], 1});
    // }

    // printf("\n===============  Measurement (tests=%lld) ================\n", repetition);

    // // Iterate over the map and print key-value pairs
    // for (const auto &pair : result_dict)
    // {
    //     // std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;

    //     printf("\"%s\" : %lld\n", pair.first, pair.second);
    // }

    // // Determine the minimum number of binary digits required
    // int minDigits = 1;
    // for (const auto &pair : result_dict)
    // {
    //     int numDigits = sizeof(int) * CHAR_BIT - __builtin_clz(pair.first);
    //     minDigits = std::max(minDigits, numDigits);
    // }

    // // Iterate over the map and print key-value pairs
    // for (const auto &pair : result_dict)
    // {
    //     std::cout << "Key: " << std::bitset<sizeof(int) * 1>(pair.first) << ", Value: " << pair.second << std::endl;
    // }

    delete state;
    delete circuit;

    return 0;
}