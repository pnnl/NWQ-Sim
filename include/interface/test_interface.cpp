#include "interface_util.hpp"

int main()
{
    saveMapAsJson(NWQSim::getBsDMGateSP(), "customized_basis_gates.json");
    // NWQSim::getBsDMGateSP();
    return 0;
}