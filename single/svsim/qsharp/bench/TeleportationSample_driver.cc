// ---------------------------------------------------------------------------
// NWQSim: Northwest Quantum Simulation Environment 
// ---------------------------------------------------------------------------
// Ang Li, Senior Computer Scientist
// Pacific Northwest National Laboratory(PNNL), U.S.
// Homepage: http://www.angliphd.com
// GitHub repo: http://www.github.com/pnnl/SV-Sim
// PNNL-IPID: 32166, ECCN: EAR99, IR: PNNL-SA-161181
// BSD Lincese.
// ---------------------------------------------------------------------------
#include <cassert> 
#include <iostream> 
#include <memory> 
#include "QirRuntimeApi_I.hpp" 
#include "QirContext.hpp"
//#include "config.h"

extern "C" void Microsoft__Quantum__Samples__Teleportation__RunProgram();
extern "C" Microsoft::Quantum::IRuntimeDriver* GetSVSim(); 

int main(int argc, char *argv[])
{
    Microsoft::Quantum::IRuntimeDriver* svsim = GetSVSim();
    Microsoft::Quantum::InitializeQirContext(svsim, false);
    Microsoft__Quantum__Samples__Teleportation__RunProgram();
    return 0;
}
