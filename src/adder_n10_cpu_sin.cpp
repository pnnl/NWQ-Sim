// ---------------------------------------------------------------------------
// NWQsim: Northwest Quantum Circuit Simulation Environment
// ---------------------------------------------------------------------------
// Ang Li, Senior Computer Scientist
// Pacific Northwest National Laboratory(PNNL), U.S.
// Homepage: http://www.angliphd.com
// GitHub repo: http://www.github.com/pnnl/DM-Sim
// PNNL-IPID: 31919-E, ECCN: EAR99, IR: PNNL-SA-143160
// BSD Lincese.
// ---------------------------------------------------------------------------
// File: adder_n10_cpu_omp.cpp
// A 10-qubit adder example based on OpenMP using CPU backend.
// ---------------------------------------------------------------------------
#include <stdio.h>
#include <time.h>
#include "util.h"
#include "svsim_cpu_sin.hpp"

//Use the SVSim namespace to enable C++/CUDA APIs
using namespace DMSim;

//You can define circuit module functions as below.
void majority(Simulation &sim, const IdxType a, const IdxType b, const IdxType c)
{
    //sim.CX(c,b);
    //sim.CX(c,a);
    //sim.append(Simulation::CCX(a, b, c));
}
void unmaj(Simulation &sim, const IdxType a, const IdxType b, const IdxType c)
{
    //sim.append(Simulation::CCX(a, b, c));
    //sim.append(Simulation::CX(c, a));
    //sim.append(Simulation::CX(a, b));
}

int main()
{
//=================================== Initialization =====================================
    srand(time(0));
    int n_qubits = 2;

    //Obtain a simulator object
    Simulation sim(n_qubits);
    
    //Simulation sim;
    //sim.AllocateQubit();

    //======== Method-1 ===========
    sim.X(0);
    sim.H(1);

    //sim.append(Simulation::H(1));
    //sim.append(Simulation::S(0));
    //sim.append(Simulation::T(1));
    //sim.append(Simulation::SDG(2));
    //sim.append(Simulation::T(3));

    //sim.append(Simulation::S(0));
    //sim.append(Simulation::T(1));
    //sim.append(Simulation::S(2));
    
    
    //sim.append(Simulation::H(1));
    //======== Method-2 ===========
    //sim.append(Simulation::C3X(0,2,3,1));

    //sim.append(Simulation::CCX(2,1,3));





    //Upload to GPU, ready for execution
    //sim.upload();

    //Run the simulation
    //sim.sim();
    
    //Measure
    //auto* res = sim.measure(5);
    //print_measurement(res, n_qubits, 5);


    auto* res = sim.measure_all(5);
    print_measurement(res, n_qubits, 5);
    //sim.print_res_sv();
    //delete res; 
    
    return 0;
}

