
#include <iostream>
#include "circuit.hpp"

int main(int argc, char const *argv[])
{
    /* code */

    auto c = NWQSim::Circuit(1);

    std::cout << c.circuitToString() << std::endl;
    return 0;
}
