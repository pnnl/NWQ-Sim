#include <cassert>
#include <cmath>                  // For M_PI
#include "../include/NWQSim.hpp"  // Assuming this is the file that contains your Circuit class
#include "../include/circuit.hpp" // Assuming this is the file that contains your Circuit class

using namespace NWQSim;

int main()
{
    Circuit circuit;
    circuit.X(0);
    circuit.H(1);

    Simulation sim(circuit.num_qubits());
    sim.sim(circuit);
    auto res = sim.measure_all(10);

    for (int i = 0; i < 10; i++)
    {
        std::cout << res[i] << std::endl;
    }
    return 0;
}