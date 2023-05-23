#include <cassert>
#include <map>
#include <bitset>
#include <climits>
#include <iostream>

#include <cmath> // For M_PI
// #include "../include/NWQSim.hpp"      // Base Class
#include "../include/svsim/svsim_cpu.hpp" // Derived Class
#include "../include/circuit.hpp"         // Assuming this is the file that contains your Circuit class

using namespace NWQSim;

int main()
{
    Circuit circuit;

    int n = 20;
    for (int i = 0; i < n; i++)
    {
        circuit.H(i);
    }
    for (int i = 0; i < n - 1; i++)
    {
        circuit.CX(0, i + 1);
    }
    for (int i = 0; i < n; i++)
    {
        circuit.H(i);
    }
    for (int i = 0; i < n - 1; i++)
    {
        circuit.CX(0, i + 1);
    }

    SVSIM_CPU sim(circuit.num_qubits());
    sim.sim(circuit);

    IdxType repetition = 1024;
    auto res = sim.measure_all(repetition);

    std::map<IdxType, IdxType> result_dict;

    for (IdxType i = 0; i < repetition; i++)
    {
        if (result_dict.find(res[i]) != result_dict.end())
            result_dict[res[i]] += 1;
        else
            result_dict.insert({res[i], 1});
    }

    // printf("\n===============  Measurement (tests=%lld) ================\n", repetition);

    // // Iterate over the map and print key-value pairs
    // for (const auto &pair : result_dict)
    // {
    //     // std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;

    //     printf("\"%s\" : %lld\n", pair.first, pair.second);
    // }

    // Determine the minimum number of binary digits required
    int minDigits = 1;
    for (const auto &pair : result_dict)
    {
        int numDigits = sizeof(int) * CHAR_BIT - __builtin_clz(pair.first);
        minDigits = std::max(minDigits, numDigits);
    }

    // Iterate over the map and print key-value pairs
    for (const auto &pair : result_dict)
    {
        std::cout << "Key: " << std::bitset<sizeof(int) * 1>(pair.first) << ", Value: " << pair.second << std::endl;
    }

    return 0;
}