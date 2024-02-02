// ---------------------------------------------------------------------------
// NWQSim: Northwest Quantum Simulation Environment 
// ---------------------------------------------------------------------------
// Ang Li, Senior Computer Scientist
// Pacific Northwest National Laboratory(PNNL), U.S.
// Homepage: http://www.angliphd.com
// GitHub repo: http://www.github.com/pnnl/NWQ-Sim
// PNNL-IPID: 32166, ECCN: EAR99, IR: PNNL-SA-161181
// BSD Lincese.
// ---------------------------------------------------------------------------
#include <cassert> 
#include <iostream> 
#include <memory> 
#include "QirRuntimeApi_I.hpp" 
#include "QirContext.hpp"

class SVSimSimulator;

extern "C" void Qrng__SampleRandomNumber();
extern "C" Microsoft::Quantum::IRuntimeDriver* GetSVSim(); 

int main(int argc, char *argv[])
{
    Microsoft::Quantum::IRuntimeDriver* svsim = GetSVSim();
    Microsoft::Quantum::InitializeQirContext(svsim, false);
    Qrng__SampleRandomNumber();
    return 0;
}
