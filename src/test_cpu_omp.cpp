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
//#define USE_AVX512

#include <stdio.h>
#include <time.h>
#include "util.h"
#include "dmsim_cpu_omp.hpp"

using namespace DMSim;
using namespace std;
#define TEST(X) pass = pass && X;
const IdxType n_cpus = 2;

bool check_sv(Simulation& sim, ValType* sv_real_expected, 
        ValType* sv_imag_expected)
{
    bool pass = true;
    for (int i=0; i<sim.dim; i++)
    {
        ValType real_diff = abs(sv_real_expected[i] - sim.sv_real[i]);
        ValType imag_diff = abs(sv_imag_expected[i] - sim.sv_imag[i]);
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
    Simulation sim(n_qubits, n_cpus);
    sim.X(0);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0,0,0,1};
    ValType sv_imag_expected[dim] = {0,0,0,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "X gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== Y Gate ================
bool test_Y()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
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
    Simulation sim(n_qubits, n_cpus);
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
    Simulation sim(n_qubits, n_cpus);
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
    Simulation sim(n_qubits, n_cpus);
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
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.SDG(0);
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
    Simulation sim(n_qubits, n_cpus);
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
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.TDG(0);
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
    Simulation sim(n_qubits, n_cpus);
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
    Simulation sim(n_qubits, n_cpus);
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
    Simulation sim(n_qubits, n_cpus);
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
    Simulation sim(n_qubits, n_cpus);
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

//============== SX Gate ================
bool test_SX()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.SX(0);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0,0,0.5};
    ValType sv_imag_expected[dim] = {0,-0.5,0.5,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "SX gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== P Gate ================
bool test_P()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.P(PI/7.,0);
    sim.sim();

    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.45048,0.45048,0.5};
    ValType sv_imag_expected[dim] = {0,0.21694,-0.21694,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "P gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== U Gate ================
bool test_U()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.U(PI/3,PI/5,PI/9,0);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.75,0.35031,0.35031,0.25};
    ValType sv_imag_expected[dim] = {0,0.25452,-0.25452,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "U gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CX Gate ================
bool test_CX()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.CX(0,1);
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
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.CY(0,1);
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
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.CZ(0,1);
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
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.CH(0,1);
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
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.CS(0,1);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.5,0,0, 0.5,0.5,0,0, 0,0,0,0, 0,0,0,0};
    ValType sv_imag_expected[dim] = {0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CS gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CSDG Gate ================
bool test_CSDG()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.CSDG(0,1);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.5,0,0, 0.5,0.5,0,0, 0,0,0,0, 0,0,0,0};
    ValType sv_imag_expected[dim] = {0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CSDG gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}


//============== CT Gate ================
bool test_CT()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.CT(0,1);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.5,0,0, 0.5,0.5,0,0, 0,0,0,0, 0,0,0,0};
    ValType sv_imag_expected[dim] = {0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CT gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CTDG Gate ================
bool test_CTDG()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.CTDG(0,1);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.5,0,0, 0.5,0.5,0,0, 0,0,0,0, 0,0,0,0};
    ValType sv_imag_expected[dim] = {0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CTDG gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}


//============== CRI Gate ================
bool test_CRI()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.CRI(PI/3.,0,1);
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
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.CRX(PI/3.,0,1);
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
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.CRY(PI/3.,0,1);
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
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.CRZ(PI/3.,0,1);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.433013,0,0, 0.433013,0.5,0,0, 0,0,0,0, 0,0,0,0};
    ValType sv_imag_expected[dim] = {0,-0.25,0,0, 0.25,0,0,0, 0,0,0,0, 0,0,0,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CRZ gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CSX Gate ================
bool test_CSX()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.CSX(0,1);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.25,0,0.25, 0.25,0.25,0,0, 0,0,0,0, 0.25,0,0,0.25};
    ValType sv_imag_expected[dim] = {0,0.25,0,-0.25, -0.25,0,0,-0.25, 0,0,0,0, 0.25,0.25,0,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CSX gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CSX Gate ================
bool test_CP()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.H(1);
    sim.CP(PI/5,0,1);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.25,0.25,0.25,0.20225, 0.25,0.25,0.25,0.20225, 0.25,0.25,0.25,0.20225, 0.20225,0.20225,0.20225,0.25};
    ValType sv_imag_expected[dim] = {0,0,0,0.14695, 0,0,0,0.14695, 0,0,0,0.14695, -0.14695,-0.14695,-0.14695,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CP gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CU Gate ================
bool test_CU()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.CU(PI/5.,PI/6.,PI/7.,PI/8.,0,1);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.43933,0,0.09406,  0.43933,0.45225,0,0.12726, 0,0,0,0, 0.09406,0.12726,0,0.04775};
    ValType sv_imag_expected[dim] = {0,0.18198,0,0.12258,  -0.18198,0,0,0.07347, 0,0,0,0, -0.12258,-0.07347,0,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CU gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}
//============== ID Gate ================
bool test_ID()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.X(0);
    sim.ID(0);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0,0,0,1};
    ValType sv_imag_expected[dim] = {0,0,0,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "ID gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== SWAP Gate ================
bool test_SWAP()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.SWAP(0,1);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0,0.5,0, 0,0,0,0, 0.5,0,0.5,0, 0,0,0,0};
    ValType sv_imag_expected[dim] = {0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "SWAP gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}
//============== Reset Gate ================
bool test_RESET()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.RESET(0);
    sim.sim();
    sim.print_res_sv();
    ValType sv_real_expected[dim] = {1,0,0,0};
    ValType sv_imag_expected[dim] = {0,0,0,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "RESET gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CCX Gate ================
bool test_CCX()
{
    const int n_qubits = 3;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.H(1);
    sim.CCX(0,1,2);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.25,0.25,0.25,0,0,0,0,0.25, 0.25,0.25,0.25,0,0,0,0,0.25, 0.25,0.25,0.25,0,0,0,0,0.25, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0.25,0.25,0.25,0,0,0,0,0.25};
    ValType sv_imag_expected[dim] = {0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CCX gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== CSWAP Gate ================
bool test_CSWAP()
{
    const int n_qubits = 3;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.H(1);
    sim.CSWAP(0,1,2);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.25,0.25,0.25,0,0,0.25,0,0, 0.25,0.25,0.25,0,0,0.25,0,0, 0.25,0.25,0.25,0,0,0.25,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0, 0.25,0.25,0.25,0,0,0.25,0,0, 0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0};
    ValType sv_imag_expected[dim] = {0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "CSWAP gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== U1 Gate ================
bool test_U1()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.U1(PI/5.,0);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.5,0.40451,0.40451,0.5};
    ValType sv_imag_expected[dim] = {0,0.29389,-0.29389,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "U1 gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== U2 Gate ================
bool test_U2()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.U2(PI/5.,PI/6.,0);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.06699,-0.14695,-0.14695,0.93301};
    ValType sv_imag_expected[dim] = {0,0.20225,-0.20225,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "U2 gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== U3 Gate ================
bool test_U3()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.H(0);
    sim.U3(PI/5.,PI/6.,PI/7.,0);
    sim.sim();
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0.23521,0.20715,0.20715,0.76479};
    ValType sv_imag_expected[dim] = {0,0.3701,-0.3701,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "U3 gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== M ================
bool test_M()
{
    const int n_qubits = 1;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.X(0);
    IdxType res = sim.measure(0);
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0,0,0,1};
    ValType sv_imag_expected[dim] = {0,0,0,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "M gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

//============== MA ================
bool test_MA()
{
    const int n_qubits = 2;
    const int dim = ((IdxType)1<<(2*n_qubits));
    bool pass = true;
    Simulation sim(n_qubits, n_cpus);
    sim.X(0);
    sim.CX(0,1);
    IdxType* res = sim.measure_all(1);
    //sim.print_res_sv();
    ValType sv_real_expected[dim] = {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1};
    ValType sv_imag_expected[dim] = {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0};
    pass = check_sv(sim, sv_real_expected, sv_imag_expected);
    std::cout << "MA gate test: " << (pass?"Success":"Failed") << std::endl;
    return pass;
}

int main()
{
//=================================== Initialization =====================================
    srand(time(0));
    bool pass = true;

    TEST(test_X());
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
    TEST(test_SX());
    TEST(test_P());
    TEST(test_U());

    TEST(test_CX());
    TEST(test_CY());
    TEST(test_CZ());
    TEST(test_CH());
    TEST(test_CS());
    TEST(test_CSDG());
    TEST(test_CT());
    TEST(test_CTDG());
    TEST(test_CRI());
    TEST(test_CRX());
    TEST(test_CRY());
    TEST(test_CRZ());
    TEST(test_CSX());
    TEST(test_CP());
    TEST(test_CU());

    TEST(test_ID());
    TEST(test_SWAP());
    TEST(test_RESET());
    TEST(test_CCX());
    TEST(test_CSWAP());
    TEST(test_U1());
    TEST(test_U2());
    TEST(test_U3());

    TEST(test_M());
    TEST(test_MA());

    std::cout << "\nUnit Test for CPU SIN " 
        << (pass?"Success":"Failed") << " !!!" << std::endl;

    return 0;
}

