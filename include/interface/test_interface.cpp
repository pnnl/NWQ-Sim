#include "interface_util.hpp"

int main()
{
    saveMapAsJson(NWQSim::getBsDMGateSP(), "basis_gate_sp.json");
    // NWQSim::getBsDMGateSP();
    return 0;
}