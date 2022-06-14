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
// File: test.cpp
// Unit test for single CPU.
// ---------------------------------------------------------------------------

#include <stdio.h>
#include <time.h>
#include "util.h"
#include "svsim_cpu_sin.hpp"

using namespace DMSim;
using namespace std;
#define TEST(X) pass = pass && X;

bool check_sv(Simulation& sim, ValType* sv_real_expected, 
        ValType* sv_imag_expected)
{
    bool pass = true;
    for (int i=0; i<sim.dim; i++)
    {
        ValType real_diff = abs(sv_real_expected[i] - sim.sv_real_cpu[i]);
        ValType imag_diff = abs(sv_imag_expected[i] - sim.sv_imag_cpu[i]);
        //printf("(%lf,%lf)",real_diff, imag_diff);
        if ( real_diff > ERROR_BAR || imag_diff > ERROR_BAR)
        {
            pass = false;
            break;
        }
    }
    return pass;
}


//============== X Gate ================
bool test_X()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.X(0);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0,0,0,1};
    ValType sv_imag_expected[dim] = {0,0,0,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "X gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

/*
//============== Y Gate ================
bool test_Y()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.Y(0);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0,0,0,1};
    ValType sv_imag_expected[dim] = {0,0,0,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "Y gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== Z Gate ================
bool test_Z()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.H(0);
    sim.Z(0);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,-0.5,-0.5,0.5};
    ValType sv_imag_expected[dim] = {0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "Z gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== H Gate ================
bool test_H()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.H(0);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.5,0.5,0.5};
    ValType sv_imag_expected[dim] = {0,0,0,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "H gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== S Gate ================
bool test_S()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.H(0);
    sim.S(0);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0,0,0.5};
    ValType sv_imag_expected[dim] = {0,0.5,-0.5,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "S gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}


//============== SDG Gate ================
bool test_SDG()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.H(0);
    sim.AdjointS(0);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0,0,0.5};
    ValType sv_imag_expected[dim] = {0,-0.5,0.5,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "SDG gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== T Gate ================
bool test_T()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.H(0);
    sim.T(0);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.353553,0.353553,0.5};
    ValType sv_imag_expected[dim] = {0,0.353553, -0.353553, 0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "T gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== TDG Gate ================
bool test_TDG()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.H(0);
    sim.AdjointT(0);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.353553,0.353553,0.5};
    ValType sv_imag_expected[dim] = {0,-0.353553,0.353553,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "TDG gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== RI Gate ================
bool test_RI()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.H(0);
    sim.T(0);
    sim.RI(PI/4.,0);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.353553,0.353553,0.5};
    ValType sv_imag_expected[dim] = {0,0.353553, -0.353553, 0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "RI gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== RX Gate ================
bool test_RX()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.RX(PI/3.,0);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.75,0,0,0.25};
    ValType sv_imag_expected[dim] = {0,-0.433013,0.433013,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "RX gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}


//============== RY Gate ================
bool test_RY()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.RY(PI/6.,0);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.933013, 0.25, 0.25, 0.066987};
    ValType sv_imag_expected[dim] = {0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "RY gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}



//============== RZ Gate ================
bool test_RZ()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.H(0);
    sim.RZ(PI/5.,0);
    sim.sim();

    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.404508,0.404508,0.5};
    ValType sv_imag_expected[dim] = {0,0.293893,-0.293893,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "RZ gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== EI Gate ================
bool test_EI()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.H(0);
    sim.EI(PI/5.,0);
    sim.sim();
    
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.5,0.5,0.5};
    ValType sv_imag_expected[dim] = {0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "EI gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== EX Gate ================
bool test_EX()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.H(0);
    sim.T(0);
    sim.EX(PI/3.,0);
    sim.sim();
    
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.193814,0.353553,0.353553,0.806186};
    ValType sv_imag_expected[dim] = {0,-0.176777,0.176777,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "EX gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== EY Gate ================
bool test_EY()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.H(0);
    sim.T(0);
    sim.EY(PI/3.,0);
    sim.sim();
    
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.806186, -0.176777, -0.176777, 0.193814};
    ValType sv_imag_expected[dim] = {0, 0.353553, -0.353553, 0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "EY gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}


//============== EZ Gate ================
bool test_EZ()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.H(0);
    sim.T(0);
    sim.EZ(PI/3.,0);
    sim.sim();
    
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.129410,0.129410,0.5};
    ValType sv_imag_expected[dim] = {0,-0.482963,0.482963,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "EZ gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}



//============== CX Gate ================
bool test_CX()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.AllocateQubit();
    sim.H(0);
    sim.ControlledX(1,((IdxType)1<<0));
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0,0,0.5, 0,0,0,0, 0,0,0,0, 0.5,0,0,0.5};
    ValType sv_imag_expected[dim] = {0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CX gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}


//============== CY Gate ================
bool test_CY()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.AllocateQubit();
    sim.H(0);
    sim.ControlledY(1,((IdxType)1<<0));
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0,0,0, 0,0,0,0, 0,0,0,0, 0.0,0,0,0.5};
    ValType sv_imag_expected[dim] = {0.0,0,0,0.5, 0,0,0,0, 0,0,0,0, -0.5,0,0,0.0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CY gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CZ Gate ================
bool test_CZ()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.AllocateQubit();
    sim.H(0);
    sim.ControlledZ(1,((IdxType)1<<0));
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.5,0,0, 0.5,0.5,0,0, 0,0,0,0, 0.0,0,0,0};
    ValType sv_imag_expected[dim] = {0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CZ gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CH Gate ================
bool test_CH()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.AllocateQubit();
    sim.H(0);
    sim.ControlledH(1,((IdxType)1<<0));
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.353553,0,0.353553, 0.353553,0.25,0,0.25, 0,0,0,0, 0.353553,0.25,0,0.25};
    ValType sv_imag_expected[dim] = {0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CH gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CS Gate ================
bool test_CS()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.AllocateQubit();
    sim.H(0);
    sim.ControlledS(1,((IdxType)1<<0));
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.5,0,0, 0.5,0.5,0,0, 0,0,0,0, 0,0,0,0};
    ValType sv_imag_expected[dim] = {0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CS gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CT Gate ================
bool test_CT()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.AllocateQubit();
    sim.H(0);
    sim.ControlledT(1,((IdxType)1<<0));
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.5,0,0, 0.5,0.5,0,0, 0,0,0,0, 0,0,0,0};
    ValType sv_imag_expected[dim] = {0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CT gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CAS Gate ================
bool test_CAS()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.AllocateQubit();
    sim.H(0);
    sim.ControlledAdjointS(1,((IdxType)1<<0));
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.5,0,0, 0.5,0.5,0,0, 0,0,0,0, 0,0,0,0};
    ValType sv_imag_expected[dim] = {0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CAS gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CAT Gate ================
bool test_CAT()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.AllocateQubit();
    sim.H(0);
    sim.ControlledAdjointT(1,((IdxType)1<<0));
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.5,0,0, 0.5,0.5,0,0, 0,0,0,0, 0,0,0,0};
    ValType sv_imag_expected[dim] = {0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CAT gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}


//============== CRI Gate ================
bool test_CRI()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.AllocateQubit();
    sim.H(0);
    sim.ControlledRI(PI/3.,1,((IdxType)1<<0));
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.433013,0,0, 0.433013,0.5,0,0, 0,0,0,0, 0,0,0,0};
    ValType sv_imag_expected[dim] = {0,-0.25,0,0, 0.25,0,0,0, 0,0,0,0, 0,0,0,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CRI gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CRX Gate ================
bool test_CRX()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.AllocateQubit();
    sim.H(0);
    sim.ControlledRX(PI/3.,1,((IdxType)1<<0));
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.433013,0,0, 0.433013,0.375,0,0, 0,0,0,0, 0,0,0,0.125};
    ValType sv_imag_expected[dim] = {0,0,0,-0.25, 0,0,0,-0.216506, 0,0,0,0, 0.25,0.216506,0,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CRX gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CRY Gate ================
bool test_CRY()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.AllocateQubit();
    sim.H(0);
    sim.ControlledRY(PI/3.,1,((IdxType)1<<0));
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.433013,0,0.25, 0.433013,0.375,0,0.216506, 0,0,0,0, 0.25,0.216506,0,0.125};
    ValType sv_imag_expected[dim] = {0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CRY gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CRZ Gate ================
bool test_CRZ()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.AllocateQubit();
    sim.H(0);
    sim.ControlledRZ(PI/3.,1,((IdxType)1<<0));
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.433013,0,0, 0.433013,0.5,0,0, 0,0,0,0, 0,0,0,0};
    ValType sv_imag_expected[dim] = {0,-0.25,0,0, 0.25,0,0,0, 0,0,0,0, 0,0,0,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CRZ gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CEI Gate ================
bool test_CEI()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.AllocateQubit();
    sim.H(0);
    sim.ControlledEI(PI/3.,1,((IdxType)1<<0));
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.25,0,0, 0.25,0.5,0,0, 0,0,0,0, 0,0,0,0};
    ValType sv_imag_expected[dim] = {0,0.433013,0,0, -0.433013,0,0,0, 0,0,0,0, 0,0,0,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CEI gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CEX Gate ================
bool test_CEX()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.AllocateQubit();
    sim.H(0);
    sim.ControlledEX(PI/3.,1,((IdxType)1<<0));
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.25,0,0, 0.25,0.125,0,0, 0,0,0,0, 0,0,0,0.375};
    ValType sv_imag_expected[dim] = {0,0,0,0.433013, 0,0,0,0.216506, 0,0,0,0, -0.433013,-0.216506,0,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CEX gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CEY Gate ================
bool test_CEY()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.AllocateQubit();
    sim.H(0);
    sim.ControlledEY(PI/3.,1,((IdxType)1<<0));
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.25,0,-0.433013, 0.25,0.125,0,-0.216506, 0,0,0,0, -0.433013,-0.216506,0,0.375};
    ValType sv_imag_expected[dim] = {0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CEY gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CEZ Gate ================
bool test_CEZ()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim;
    sim.AllocateQubit();
    sim.AllocateQubit();
    sim.H(0);
    sim.ControlledEZ(PI/3.,1,((IdxType)1<<0));
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.25,0,0, 0.25,0.50,0, 0,0,0,0, 0,0,0,0};
    ValType sv_imag_expected[dim] = {0,0.433013,0,0, -0.433013,0,0,0, 0,0,0,0, 0,0,0,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CEZ gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}
*/




int main()
{
//=================================== Initialization =====================================
    srand(RAND_SEED);
    bool pass = true;

    TEST(test_X());

    /*
    TEST(test_Y());
    TEST(test_Z());
    TEST(test_H());

    TEST(test_S());
    TEST(test_T());
    TEST(test_SDG());
    TEST(test_TDG());

    TEST(test_RI());
    TEST(test_RX());
    TEST(test_RY());
    TEST(test_RZ());
    
    TEST(test_EI());
    TEST(test_EX());
    TEST(test_EY());
    TEST(test_EZ());

    TEST(test_CX());
    TEST(test_CY());
    TEST(test_CZ());
    TEST(test_CH());
    TEST(test_CS());
    TEST(test_CT());

    TEST(test_CRI());
    TEST(test_CRX());
    TEST(test_CRY());
    TEST(test_CRZ());

    TEST(test_CEI());
    TEST(test_CEX());
    TEST(test_CEY());
    TEST(test_CEZ());

    TEST(test_CAS());
    TEST(test_CAT());
    */

    std::cout << "\nUnit Test for CPU SIN " 
        << (pass?"Success":"Failed") << " !!!" << std::endl;

    return 0;
}

