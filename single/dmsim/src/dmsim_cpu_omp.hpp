// ---------------------------------------------------------------------------
// NWQsim: Northwest Quantum Circuit Simulation Environment
// ---------------------------------------------------------------------------
// Ang Li, Senior Computer Scientist
// Pacific Northwest National Laboratory(PNNL), U.S.
// Homepage: http://www.angliphd.com
// GitHub repo: http://www.github.com/pnnl/SV-Sim
// PNNL-IPID: 32166, ECCN: EAR99, IR: PNNL-SA-161181
// BSD Lincese.
// ---------------------------------------------------------------------------
// File: dmsim_cpu_omp.hpp
// OpenMP CPU based density-matrix simulation 
// ---------------------------------------------------------------------------
#ifndef DMSIM_CPU_OMP_HPP_
#define DMSIM_CPU_OMP_HPP_
#include <assert.h>
#include <vector>
#include <math.h>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <cstring>
#include "config.h"
#include "util.h"
#include "gate.h"

using namespace NWQSim;

namespace DMSim
{
using namespace std;
class Gate;
class Simulation;
void Purity_Check(const Gate& g, const IdxType t, const Simulation* sim, ValType* sv_real, ValType* sv_imag);
//Unified gate functional pointer definition
using func_t = void (*)(const Gate*, const Simulation*, ValType*, ValType*);
//Declare major simulation kernel
void simulation_kernel(Simulation*);
/***********************************************
 * Declare gate function pointers
 ***********************************************/
extern func_t pX;
extern func_t pY;
extern func_t pZ;
extern func_t pH;
extern func_t pS;
extern func_t pSDG;
extern func_t pT;
extern func_t pTDG;
extern func_t pRI;
extern func_t pRX;
extern func_t pRY;
extern func_t pRZ;
extern func_t pSX;
extern func_t pP;
extern func_t pU;
extern func_t pCX;
extern func_t pCY;
extern func_t pCZ;
extern func_t pCH;
extern func_t pCS;
extern func_t pCSDG;
extern func_t pCT;
extern func_t pCTDG;
extern func_t pCRI;
extern func_t pCRX;
extern func_t pCRY;
extern func_t pCRZ;
extern func_t pCSX;
extern func_t pCP;
extern func_t pCU;
extern func_t pID;
extern func_t pSWAP;
extern func_t pM;
extern func_t pMA;
extern func_t pRESET;
/***********************************************
 * Gate Definition
 ***********************************************/
class Gate
{
public:
    Gate(enum OP _op_name, 
            func_t _op, 
            IdxType _qubit,
            IdxType _ctrl=0,
            ValType _theta=0, 
            ValType _phi=0, 
            ValType _lambda=0
        ) :
        op(_op), op_name(_op_name), qubit(_qubit), ctrl(_ctrl),
        theta(_theta), phi(_phi), lambda(_lambda) {}
    ~Gate() {}
    //applying the embedded gate operation 
    void exe_op(Simulation* sim, ValType* sv_real, ValType* sv_imag)
    {
        (*(this->op))(this, sim, sv_real, sv_imag);
    }
    //Gate operation
    func_t op;
    //Gate name
    enum OP op_name;
    //Qubit 
    IdxType qubit;
    //ctrl qubit
    IdxType ctrl;
    //Rotation parameter 0
    ValType theta;
    //Rotation parameter 1
    ValType phi;
    //Rotation parameter 2
    ValType lambda;
}; //end of Gate definition

/***********************************************
 * Circuit Definition
 ***********************************************/
class Circuit
{
public:
    Circuit(IdxType _n_qubits=0):
        n_qubits(_n_qubits), n_gates(0), circuit_cpu(NULL) {}
    ~Circuit() { clear(); }
    void append(Gate& g)
    {
#ifdef PRINT_GATE_TRACE
        printf("%s(theta:%lf,phi:%lf,lambda:%lf,ctrl:%llu,qubit:%llu)\n",
                OP_NAMES[g.op_name], g.theta, g.phi, g.lambda, g.ctrl, g.qubit);
#endif
        if (g.qubit >= n_qubits) 
        {
            printf("%s(theta:%lf,phi:%lf,lambda:%lf,ctrl:%llu,qubit:%llu)\n",
                    OP_NAMES[g.op_name], g.theta, g.phi, g.lambda, g.ctrl, g.qubit);
            throw std::logic_error("Qubit out of range!");
        }
        circuit.push_back(g);
        delete (&g);
        n_gates++;
    }
    void AllocateQubit() 
    { 
        ++n_qubits; 
#ifdef PRINT_QUBIT_TRACE
        printf("Allocate 1 qubit, in total: %llu\n",n_qubits);
#endif
    }
    void ReleaseQubit()
    {
        --n_qubits;
#ifdef PRINT_QUBIT_TRACE
        printf("Release 1 qubit, in total: %llu\n", n_qubits);
#endif
    }
    void clear()
    {
        circuit.clear();
        n_gates = 0;
        SAFE_FREE_HOST(circuit_cpu);
    }
    Gate* upload()
    {
        /* TODO: 
         * Here is where transpilation and opitmization is performed */
        SAFE_FREE_HOST(circuit_cpu);
        SAFE_ALOC_HOST(circuit_cpu, n_gates*sizeof(Gate));
        memcpy(circuit_cpu, circuit.data(), n_gates*sizeof(Gate));
        return circuit_cpu;
    }
public:
    // number of qubits
    IdxType n_qubits;
    // number of gates
    IdxType n_gates;
    // user input gate sequence
    vector<Gate> circuit;
    // gate sequence executed
    Gate* circuit_cpu;
};
/***********************************************
 * Simulation Definition
 ***********************************************/
class Simulation
{
public:
    Simulation(IdxType _n_qubits=0, IdxType _n_cpus = DEFAULT_OMP_PE) 
        : n_cpus(_n_cpus),
        n_gates(0), 
        sv_real(NULL),
        sv_imag(NULL),
        m_real(NULL),
        randoms(NULL),
        results(NULL),
        noise_1q_real(NULL),
        noise_1q_imag(NULL),
        noise_2q_real(NULL),
        noise_2q_imag(NULL)
    {
        if (_n_qubits == 0) update(N_QUBIT_SLOT/2, 0);
        else update(_n_qubits, 0);
        //initialization
        circuit_handle = new Circuit(_n_qubits);
        circuit_handle_cpu = NULL;
        SAFE_ALOC_HOST(sv_real, sv_size);
        SAFE_ALOC_HOST(sv_imag, sv_size);
        SAFE_ALOC_HOST(m_real, sv_size+1);
        memset(sv_real, 0, sv_size);
        memset(sv_imag, 0, sv_size);
        memset(m_real, 0, sv_size+1);
        //State-vector initial state [0..0] = 1
        sv_real[0] = 1.0;
        srand(time(0));
#ifdef PRINT_SIM_TRACE
        printf("DMSim_cpu is initialized!\n");
#endif
    }
    ~Simulation()
    {
        //Release circuit
        if (circuit_handle != NULL) delete circuit_handle;
        //Release circuit
        SAFE_FREE_HOST(sv_real);
        SAFE_FREE_HOST(sv_imag);
        SAFE_FREE_HOST(m_real);
        SAFE_FREE_HOST(randoms);
        SAFE_FREE_HOST(results);
        SAFE_FREE_HOST(noise_1q_real);
        SAFE_FREE_HOST(noise_1q_imag);
        SAFE_FREE_HOST(noise_2q_real);
        SAFE_FREE_HOST(noise_2q_imag);
#ifdef PRINT_SIM_TRACE
        printf("DMSim_cpu is finalized!\n\n");
#endif
    }
    void AllocateQubit()
    {
        circuit_handle->AllocateQubit();
    }
    void ReleaseQubit()
    {
        circuit_handle->ReleaseQubit();
    }
    // =============================== Gates ===================================
    // Basic gates
    void X(IdxType qubit)
    {
        Gate* G = new Gate(OP::X,pX,qubit);
        circuit_handle->append(*G);
    }
    void Y(IdxType qubit)
    {
        Gate* G = new Gate(OP::Y,pY,qubit);
        circuit_handle->append(*G);
    }
    void Z(IdxType qubit)
    {
        Gate* G = new Gate(OP::Z,pZ,qubit);
        circuit_handle->append(*G);
    }
    void H(IdxType qubit)
    {
        Gate* G = new Gate(OP::H,pH,qubit);
        circuit_handle->append(*G);
    }
    void S(IdxType qubit)
    {
        Gate* G = new Gate(OP::S,pS,qubit);
        circuit_handle->append(*G);
    }
    void SDG(IdxType qubit)
    {
        Gate* G = new Gate(OP::SDG,pSDG,qubit);
        circuit_handle->append(*G);
    }
    void T(IdxType qubit)
    {
        Gate* G = new Gate(OP::T,pT,qubit);
        circuit_handle->append(*G);
    }
    void TDG(IdxType qubit)
    {
        Gate* G = new Gate(OP::TDG,pTDG,qubit);
        circuit_handle->append(*G);
    }
    void RI(ValType theta, IdxType qubit)
    {
        Gate* G = new Gate(OP::RI,pRI,qubit,0,theta);
        circuit_handle->append(*G);
    }
    void RX(ValType theta, IdxType qubit)
    {
        Gate* G = new Gate(OP::RX,pRX,qubit,0,theta);
        circuit_handle->append(*G);
    }
    void RY(ValType theta, IdxType qubit)
    {
        Gate* G = new Gate(OP::RY,pRY,qubit,0,theta);
        circuit_handle->append(*G);
    }
    void RZ(ValType theta, IdxType qubit)
    {
        Gate* G = new Gate(OP::RZ,pRZ,qubit,0,theta);
        circuit_handle->append(*G);
    }
    void SX(IdxType qubit)
    {
        Gate* G = new Gate(OP::SX,pSX,qubit);
        circuit_handle->append(*G);
    }
    void P(ValType theta, IdxType qubit)
    {
        Gate* G = new Gate(OP::P,pP,qubit,0,theta);
        circuit_handle->append(*G);
    }
    void U(ValType theta, ValType phi, ValType lambda, IdxType qubit)
    {
        Gate* G = new Gate(OP::U,pU,qubit,0,theta,phi,lambda);
        circuit_handle->append(*G);
    }
    // ctrlled gates
    void CX(IdxType ctrl, IdxType qubit)
    {
        Gate* G = new Gate(OP::CX,pCX,qubit,ctrl);
        circuit_handle->append(*G);
    }
    void CY(IdxType ctrl, IdxType qubit)
    {
        Gate* G = new Gate(OP::CY,pCY,qubit,ctrl);
        circuit_handle->append(*G);
    }
    void CZ(IdxType ctrl, IdxType qubit)
    {
        Gate* G = new Gate(OP::CZ,pCZ,qubit,ctrl);
        circuit_handle->append(*G);
    }
    void CH(IdxType ctrl, IdxType qubit)
    {
        Gate* G = new Gate(OP::CH,pCH,qubit,ctrl);
        circuit_handle->append(*G);
    }
    void CS(IdxType ctrl, IdxType qubit)
    {
        Gate* G = new Gate(OP::CS,pCS,qubit,ctrl);
        circuit_handle->append(*G);
    }
    void CSDG(IdxType ctrl, IdxType qubit)
    {
        Gate* G = new Gate(OP::CSDG,pCSDG,qubit,ctrl);
        circuit_handle->append(*G);
    }
    void CT(IdxType ctrl, IdxType qubit)
    {
        Gate* G = new Gate(OP::CT,pCT,qubit,ctrl);
        circuit_handle->append(*G);
    }
    void CTDG(IdxType ctrl, IdxType qubit)
    {
        Gate* G = new Gate(OP::CTDG,pCTDG,qubit,ctrl);
        circuit_handle->append(*G);
    }
    void CRI(ValType theta, IdxType ctrl, IdxType qubit)
    {
        Gate* G = new Gate(OP::CRI,pCRI,qubit,ctrl,theta);
        circuit_handle->append(*G);
    }
    void CRX(ValType theta, IdxType ctrl, IdxType qubit)
    {
        Gate* G = new Gate(OP::CRX,pCRX,qubit,ctrl,theta);
        circuit_handle->append(*G);
    }
    void CRY(ValType theta, IdxType ctrl, IdxType qubit)
    {
        Gate* G = new Gate(OP::CRY,pCRY,qubit,ctrl,theta);
        circuit_handle->append(*G);
    }
    void CRZ(ValType theta, IdxType ctrl, IdxType qubit)
    {
        Gate* G = new Gate(OP::CRZ,pCRZ,qubit,ctrl,theta);
        circuit_handle->append(*G);
    }
    void CSX(IdxType ctrl, IdxType qubit)
    {
        Gate* G = new Gate(OP::CSX,pCSX,qubit,ctrl);
        circuit_handle->append(*G);
    }
    void CP(ValType theta, IdxType ctrl, IdxType qubit)
    {
        Gate* G = new Gate(OP::CP,pCP,qubit,ctrl,theta);
        circuit_handle->append(*G);
    }
    void CU(ValType theta, ValType phi, ValType lambda, ValType gamma, IdxType ctrl, IdxType qubit)
    {
        P(gamma,ctrl);
        Gate* G = new Gate(OP::CU,pCU,qubit,ctrl,theta,phi,lambda);
        circuit_handle->append(*G);
    }
    //Other
    void ID(IdxType qubit)
    {
        Gate* G = new Gate(OP::ID,pID,qubit);
        circuit_handle->append(*G);
    }
    void SWAP(IdxType qubit0, IdxType qubit1)
    {
        Gate* G = new Gate(OP::SWAP,pSWAP,qubit0,qubit1);
        circuit_handle->append(*G);
    }
    void M(IdxType qubit, IdxType pauli) //default is pauli-Z
    {
        SAFE_FREE_HOST(results);
        SAFE_ALOC_HOST(results, sizeof(IdxType));
        memset(results, 0, sizeof(IdxType));
        ValType rand = randomval();
        Gate* G = new Gate(OP::M,pM,qubit,pauli,rand);
        circuit_handle->append(*G);
    }
    void MA(IdxType repetition) //default is pauli-Z
    {
        SAFE_FREE_HOST(results);
        SAFE_ALOC_HOST(results, sizeof(IdxType)*repetition);
        memset(results, 0, sizeof(IdxType)*repetition);
        SAFE_FREE_HOST(randoms);
        SAFE_ALOC_HOST(randoms, sizeof(ValType)*repetition);
        for (IdxType i=0; i<repetition; i++) 
            randoms[i] = (ValType)rand()/(ValType)RAND_MAX;
        Gate* G = new Gate(OP::MA,pMA,0,repetition);
        circuit_handle->append(*G);
    }
    void RESET(IdxType qubit)
    {
        Gate* G = new Gate(OP::RESET,pRESET,qubit);
        circuit_handle->append(*G);
    }
    //Composition Gates
    void MZ(IdxType qubit)
    {
        M(qubit,0);
    }
    void CCX(IdxType ctrl0, IdxType ctrl1, IdxType qubit)
    {
        H(qubit);
        CX(ctrl1,qubit);
        TDG(qubit);
        CX(ctrl0,qubit);
        T(qubit);
        CX(ctrl1,qubit);
        TDG(qubit);
        CX(ctrl0,qubit);
        T(ctrl1);
        T(qubit);
        H(qubit);
        CX(ctrl0,ctrl1);
        T(ctrl0);
        TDG(ctrl1);
        CX(ctrl0,ctrl1);
    }
    void CSWAP(IdxType ctrl, IdxType qubit0, IdxType qubit1)
    {
        CX(qubit1, qubit0);
        CCX(ctrl, qubit0, qubit1);
        CX(qubit1, qubit0);
    }
    void U1(ValType lambda, IdxType qubit)
    {
        U(0,0,lambda,qubit);
    }
    void U2(ValType phi, ValType lambda, IdxType qubit) 
    {
        RI(-HALF*(phi+lambda),qubit);
        U(PI/2,phi,lambda,qubit);
    }
    void U3(ValType theta, ValType phi, ValType lambda, IdxType qubit)
    {
        RI(-HALF*(phi+lambda),qubit);
        U(theta,phi,lambda,qubit);
    }
    // =============================== Other functions ===================================
    void reset_sim()
    {
        memset(sv_real, 0, sv_size);
        memset(sv_imag, 0, sv_size);
        memset(m_real, 0, sv_size+1);
        SAFE_FREE_HOST(results);
        SAFE_FREE_HOST(randoms);
        //State Vector initial state [0..0] = 1
        sv_real[0] = 1.0;
        reset_circuit();        
    }
    void reset_circuit()
    {
        circuit_handle->clear();
#ifdef PRINT_CIRCUIT_TRACE
        printf("Circuit is reset!\n");
#endif
    }
    IdxType get_n_qubits()
    {
        return circuit_handle->n_qubits;
    }
    IdxType get_n_gates()
    {
        return circuit_handle->n_gates;
    }
    void set_seed(IdxType seed)
    {
        srand(seed);
    }
    void update(const IdxType _n_qubits, const IdxType _n_gates)
    {
        if (_n_qubits > (N_QUBIT_SLOT/2)) 
            throw std::logic_error("Request more qubits than slots!");
#ifdef USE_AVX512
        if (_n_qubits < 3) 
            throw std::logic_error("Number of qubits is too small for AVX512!");
#endif
        this->n_qubits = _n_qubits;
        this->n_gates = _n_gates;
        this->dim = ((IdxType)1<<(2*_n_qubits));
        this->half_dim = (IdxType)1<<(2*_n_qubits-1);
        this->sv_size = dim*(IdxType)sizeof(ValType);
    }
    void sim()
    {
        update(circuit_handle->n_qubits, circuit_handle->n_gates);
        circuit_handle_cpu = circuit_handle->upload();
#ifdef PRINT_SIM_TRACE
        printf("DMSim_cpu is running!\n");
#endif
#ifdef PRINT_KERNEL_TRACE
        double* sim_times;
        SAFE_ALOC_HOST(sim_times, sizeof(double)*n_cpus);
#endif

#pragma omp parallel num_threads(n_cpus) 
        {
            int d = omp_get_thread_num();
#ifdef PRINT_KERNEL_TRACE
            cpu_timer sim_timer;
            #pragma omp barrier
            sim_timer.start_timer();
#endif
            //================================
            simulation_kernel(this);
            //================================
#ifdef PRINT_KERNEL_TRACE
            sim_timer.stop_timer();
            sim_times[d] = sim_timer.measure();
#endif
        }
#ifdef PRINT_KERNEL_TRACE
        double avg_sim_time = 0;
        for (unsigned d=0; d<n_cpus; d++)
        {
            avg_sim_time += sim_times[d];
        }
        avg_sim_time /= (double)n_cpus;
        printf("============== DMSim_cpu_omp ===============\n");
#ifdef USE_AVX512
        printf("nqubits:%llu, ngates:%llu, ncores(avx512):%llu, sim:%.3lf ms.\n",
                n_qubits, n_gates, n_cpus, avg_sim_time);
#else
        printf("nqubits:%llu, ngates:%llu, ncores:%llu, sim:%.3lf ms.\n",
                n_qubits, n_gates, n_cpus, avg_sim_time);
#endif
        printf("=====================================\n");
#endif
        reset_circuit();
    }
    //measure one qubit, default in Pauli_Z
    //return 0 or 1 for that qubit 
    IdxType measure(IdxType qubit, IdxType pauli=0) 
    {
        this->M(qubit, pauli);
        this->sim();
        return this->results[0];
    }
    IdxType measureZ(IdxType qubit)
    {
        this->MZ(qubit);
        this->sim();
        return this->results[0];
    }
    IdxType* measure_all(IdxType repetition=DEFAULT_REPETITIONS)
    {
        this->MA(repetition);
        this->sim();
        return this->results;
    }
    void print_res_sv()
    {
        printf("----- Real SV ------\n");
        for (IdxType i=0; i<dim; i++) 
            printf("%lf ", sv_real[i]);
        printf("\n");
        printf("----- Imag SV ------\n");
        for (IdxType i=0; i<dim; i++) 
            printf("%lf ", sv_imag[i]);
        printf("\n");
    }
public:
    IdxType n_cpus;
    // n_qubits is the number of qubits
    IdxType n_qubits;
    IdxType dim;
    IdxType half_dim;
    IdxType sv_size;
    IdxType n_gates;
    //cpu arrays, they are used as alias of particular pointers
    ValType* sv_real;
    ValType* sv_imag;
    //For measurement interal usage
    ValType* m_real;
    //For measurement randoms
    ValType* randoms;
    //For measurement result
    IdxType* results;
    //circuit
    Circuit* circuit_handle;
    //circuit cpu
    Gate* circuit_handle_cpu;
    //noise 1q
    ValType* noise_1q_real;
    ValType* noise_1q_imag;
    //noise 2q
    ValType* noise_2q_real;
    ValType* noise_2q_imag;
};

//============== Purity Check  ================
void Purity_Check(const Gate& g, const IdxType t, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    if (omp_get_thread_num() == 0)
    {
        ValType purity = 0; 
        for (IdxType i=0; i<(((IdxType)1<<(sim->n_qubits))); i++)
            purity += abs(sv_real[(i<<(sim->n_qubits))+i]); //diagonal element
        if ( abs(purity - 1.0) > ERROR_BAR )
            printf("Purity Check fails after Gate-%llu=>%s(theta:%lf,phi:%lf,lambda:%lf,ctrl:%llu,qubit:%llu) with %lf\n", 
                    t, OP_NAMES[g.op_name], g.theta, g.phi, g.lambda, g.ctrl, g.qubit, purity);
    }
}
/***********************************************
 * Simulation Main Kernel
 ***********************************************/
void simulation_kernel(Simulation* sim)
{
    for (IdxType t=0; t<(sim->n_gates); t++)
    {
        Gate& g = (sim->circuit_handle_cpu)[t];
        g.exe_op(sim, sim->sv_real, sim->sv_imag);
#ifdef PURITY_CHECK 
        Purity_Check(g, t, sim, sim->sv_real, sim->sv_imag);
#endif
    }
}

#ifdef USE_AVX512
#include "nwqsim_cpu_avx512.hpp"
//#include "nwqsim_cpu_avx2.hpp"

#else

/***********************************************
 * Key Macros
 ***********************************************/
//Common
#define PUT(arr,i,val) (arr[(i)]=(val))
#define GET(arr,i) (arr[(i)])
//Load 1st qubit coefficient
#define LOAD_Q0 \
    const ValType el0_real=GET(sv_real,pos0);\
    const ValType el0_imag=GET(sv_imag,pos0);\
    const ValType el1_real=GET(sv_real,pos1);\
    const ValType el1_imag=GET(sv_imag,pos1);
//Save 1st qubit coefficient
#define STORE_Q0 \
    PUT(sv_real,pos0,sv_real_pos0);\
    PUT(sv_imag,pos0,sv_imag_pos0);\
    PUT(sv_real,pos1,sv_real_pos1);\
    PUT(sv_imag,pos1,sv_imag_pos1);
//Load 2nd qubit coefficient
#define LOAD_Q1 \
    const ValType el2_real=GET(sv_real,pos2);\
    const ValType el2_imag=GET(sv_imag,pos2);\
    const ValType el3_real=GET(sv_real,pos3);\
    const ValType el3_imag=GET(sv_imag,pos3);
//Save 2nd qubit coefficient
#define STORE_Q1 \
    PUT(sv_real,pos2,sv_real_pos2);\
    PUT(sv_imag,pos2,sv_imag_pos2);\
    PUT(sv_real,pos3,sv_real_pos3);\
    PUT(sv_imag,pos3,sv_imag_pos3);
//Define MG-BSP machine operation header (Optimized version)
#define OP_HEAD \
    _Pragma("omp for schedule(auto)") \
    for (IdxType i=0; i<(sim->half_dim);i++){ \
        IdxType outer = (i >> qubit); \
        IdxType inner = (i & (((IdxType)1<<qubit)-1)); \
        IdxType offset = (outer << (qubit+1)); \
        IdxType pos0 = offset + inner; \
        IdxType pos1 = pos0 + ((IdxType)1<<qubit); 
//Define operation header for 2-qubit
#define OP_HEAD_2Q \
    const IdxType q0dim = ((IdxType)1<<max(ctrl,qubit));\
    const IdxType q1dim = ((IdxType)1<<min(ctrl,qubit));\
    assert (ctrl != qubit);\
    const IdxType outer_factor=((sim->dim)+q0dim+q0dim-1)>>(max(ctrl,qubit)+1);\
    const IdxType mider_factor=(q0dim+q1dim+q1dim-1)>>(min(ctrl,qubit)+1);\
    const IdxType inner_factor = q1dim;\
    const IdxType qubit1_dim = ((IdxType)1<<ctrl);\
    const IdxType qubit2_dim = ((IdxType)1<<qubit);\
    _Pragma("omp for schedule(auto)") \
    for (IdxType i=0; i<outer_factor*mider_factor*inner_factor; i++){ \
        IdxType outer = ((i/inner_factor)/(mider_factor))*(q0dim+q0dim);\
        IdxType mider = ((i/inner_factor)%(mider_factor))*(q1dim+q1dim);\
        IdxType inner = i%inner_factor;\
        IdxType pos0 = outer+mider+inner;\
        IdxType pos1 = outer+mider+inner+qubit2_dim;\
        IdxType pos2 = outer+mider+inner+qubit1_dim;\
        IdxType pos3 = outer+mider+inner+q0dim+q1dim;

//Define operation tail
#define OP_TAIL } _Pragma("omp barrier")  

//For C2 and C4 gates
#define DIV2E(x,y) ((x)>>(y))
#define MOD2E(x,y) ((x)&(((IdxType)1<<(y))-(IdxType)1)) 
#define EXP2E(x) ((IdxType)1<<(x))
#define SV16IDX(x) ( ((x>>3)&1)*EXP2E(qubit0) + ((x>>2)&1)*EXP2E(qubit1) + ((x>>1)&1)*EXP2E(qubit2) + ((x&1)*EXP2E(qubit3)) )

/***********************************************
 * Gate Implementation
 ***********************************************/
//============== X Gate ================
inline void X_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    OP_HEAD;
    LOAD_Q0;
    ValType sv_real_pos0 = el1_real;
    ValType sv_imag_pos0 = el1_imag;
    ValType sv_real_pos1 = el0_real;
    ValType sv_imag_pos1 = el0_imag;
    STORE_Q0;
    OP_TAIL;
}
//============== Y Gate ================
inline void Y_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    OP_HEAD;
    LOAD_Q0;
    ValType sv_real_pos0 =  el1_imag; 
    ValType sv_imag_pos0 = -el1_real;
    ValType sv_real_pos1 = -el0_imag;
    ValType sv_imag_pos1 =  el0_real;
    STORE_Q0;
    OP_TAIL;
}
//============== ConjugateY Gate ================
inline void ConjugateY_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    OP_HEAD;
    LOAD_Q0;
    ValType sv_real_pos0 = -el1_imag; 
    ValType sv_imag_pos0 =  el1_real;
    ValType sv_real_pos1 =  el0_imag;
    ValType sv_imag_pos1 = -el0_real;
    STORE_Q0;
    OP_TAIL;
}
//============== Z Gate ================
inline void Z_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    OP_HEAD;
    const ValType el1_real = GET(sv_real, pos1);
    const ValType el1_imag = GET(sv_imag, pos1);
    ValType sv_real_pos1 = -el1_real;
    ValType sv_imag_pos1 = -el1_imag;
    PUT(sv_real, pos1, sv_real_pos1);
    PUT(sv_imag, pos1, sv_imag_pos1);
    OP_TAIL;
}
//============== H Gate ================
inline void H_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    OP_HEAD;
    LOAD_Q0;
    ValType sv_real_pos0 = S2I*(el0_real + el1_real); 
    ValType sv_imag_pos0 = S2I*(el0_imag + el1_imag);
    ValType sv_real_pos1 = S2I*(el0_real - el1_real);
    ValType sv_imag_pos1 = S2I*(el0_imag - el1_imag);
    STORE_Q0;
    OP_TAIL;
}
//============== S Gate ================
inline void S_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    OP_HEAD;
    const ValType el1_real = GET(sv_real, pos1);
    const ValType el1_imag = GET(sv_imag, pos1);
    ValType sv_real_pos1 = -el1_imag;
    ValType sv_imag_pos1 =  el1_real;
    PUT(sv_real, pos1, sv_real_pos1);
    PUT(sv_imag, pos1, sv_imag_pos1);
    OP_TAIL;
}
//============== SDG Gate ================
inline void SDG_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    OP_HEAD;
    const ValType el1_real = GET(sv_real, pos1);
    const ValType el1_imag = GET(sv_imag, pos1);
    ValType sv_real_pos1 =  el1_imag;
    ValType sv_imag_pos1 = -el1_real;
    PUT(sv_real, pos1, sv_real_pos1);
    PUT(sv_imag, pos1, sv_imag_pos1);
    OP_TAIL;
}
//============== T Gate ================
inline void T_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    OP_HEAD;
    const ValType el1_real = GET(sv_real, pos1);
    const ValType el1_imag = GET(sv_imag, pos1);
    ValType sv_real_pos1 = S2I*(el1_real-el1_imag);
    ValType sv_imag_pos1 = S2I*(el1_real+el1_imag);
    PUT(sv_real, pos1, sv_real_pos1);
    PUT(sv_imag, pos1, sv_imag_pos1);
    OP_TAIL;
}
//============== TDG Gate ================
inline void TDG_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    OP_HEAD;
    const ValType el1_real = GET(sv_real, pos1);
    const ValType el1_imag = GET(sv_imag, pos1);
    ValType sv_real_pos1 = S2I*( el1_real+el1_imag);
    ValType sv_imag_pos1 = S2I*(-el1_real+el1_imag);
    PUT(sv_real, pos1, sv_real_pos1);
    PUT(sv_imag, pos1, sv_imag_pos1);
    OP_TAIL;
}
//============== RI Gate ================
inline void RI_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const IdxType qubit)
{
    ValType ri_real =  cos(HALF*theta);
    ValType ri_imag = -sin(HALF*theta);
    OP_HEAD;
    LOAD_Q0;
    ValType sv_real_pos0 = (el0_real * ri_real) - (el0_imag * ri_imag);
    ValType sv_imag_pos0 = (el0_real * ri_imag) + (el0_imag * ri_real);
    ValType sv_real_pos1 = (el1_real * ri_real) - (el1_imag * ri_imag);
    ValType sv_imag_pos1 = (el1_real * ri_imag) + (el1_imag * ri_real);
    STORE_Q0;
    OP_TAIL;
}
//============== ConjugateRI Gate ================
inline void ConjugateRI_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const IdxType qubit)
{
    ValType ri_real = cos(HALF*theta);
    ValType ri_imag = sin(HALF*theta);
    OP_HEAD;
    LOAD_Q0;
    ValType sv_real_pos0 = (el0_real * ri_real) - (el0_imag * ri_imag);
    ValType sv_imag_pos0 = (el0_real * ri_imag) + (el0_imag * ri_real);
    ValType sv_real_pos1 = (el1_real * ri_real) - (el1_imag * ri_imag);
    ValType sv_imag_pos1 = (el1_real * ri_imag) + (el1_imag * ri_real);
    STORE_Q0;
    OP_TAIL;
}
//============== RX Gate ================
inline void RX_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const IdxType qubit)
{
    ValType rx_real =  cos(HALF*theta);
    ValType rx_imag = -sin(HALF*theta);
    OP_HEAD;
    LOAD_Q0;
    ValType sv_real_pos0 =  (rx_real * el0_real) - (rx_imag * el1_imag);
    ValType sv_imag_pos0 =  (rx_real * el0_imag) + (rx_imag * el1_real);
    ValType sv_real_pos1 = -(rx_imag * el0_imag) + (rx_real * el1_real);
    ValType sv_imag_pos1 =  (rx_imag * el0_real) + (rx_real * el1_imag);
    STORE_Q0;
    OP_TAIL;
}
//============== ConjugateRX Gate ================
inline void ConjugateRX_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const IdxType qubit)
{
    ValType rx_real = cos(HALF*theta);
    ValType rx_imag = sin(HALF*theta);
    OP_HEAD;
    LOAD_Q0;
    ValType sv_real_pos0 =  (rx_real * el0_real) - (rx_imag * el1_imag);
    ValType sv_imag_pos0 =  (rx_real * el0_imag) + (rx_imag * el1_real);
    ValType sv_real_pos1 = -(rx_imag * el0_imag) + (rx_real * el1_real);
    ValType sv_imag_pos1 =  (rx_imag * el0_real) + (rx_real * el1_imag);
    STORE_Q0;
    OP_TAIL;
}
//============== RY Gate ================
inline void RY_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const IdxType qubit)
{
    ValType e0_real =  cos(HALF*theta);
    ValType e1_real = -sin(HALF*theta);
    ValType e2_real =  sin(HALF*theta);
    ValType e3_real =  cos(HALF*theta);
    OP_HEAD;
    LOAD_Q0;
    ValType sv_real_pos0 = (e0_real * el0_real) + (e1_real * el1_real);
    ValType sv_imag_pos0 = (e0_real * el0_imag) + (e1_real * el1_imag);
    ValType sv_real_pos1 = (e2_real * el0_real) + (e3_real * el1_real);
    ValType sv_imag_pos1 = (e2_real * el0_imag) + (e3_real * el1_imag);
    STORE_Q0;
    OP_TAIL;
}
//============== RZ Gate ================
inline void RZ_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const IdxType qubit)
{
    ValType e0_real =  cos(HALF*theta);
    ValType e0_imag = -sin(HALF*theta);
    ValType e3_real =  cos(HALF*theta);
    ValType e3_imag =  sin(HALF*theta);
    OP_HEAD;
    LOAD_Q0;
    ValType sv_real_pos0 = (el0_real * e0_real) - (el0_imag * e0_imag);
    ValType sv_imag_pos0 = (el0_real * e0_imag) + (el0_imag * e0_real);
    ValType sv_real_pos1 = (el1_real * e3_real) - (el1_imag * e3_imag);
    ValType sv_imag_pos1 = (el1_real * e3_imag) + (el1_imag * e3_real);
    STORE_Q0;
    OP_TAIL;
}
//============== ConjugateRZ Gate ================
inline void ConjugateRZ_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const IdxType qubit)
{
    ValType e0_real =  cos(HALF*theta);
    ValType e0_imag =  sin(HALF*theta);
    ValType e3_real =  cos(HALF*theta);
    ValType e3_imag = -sin(HALF*theta);
    OP_HEAD;
    LOAD_Q0;
    ValType sv_real_pos0 = (el0_real * e0_real) - (el0_imag * e0_imag);
    ValType sv_imag_pos0 = (el0_real * e0_imag) + (el0_imag * e0_real);
    ValType sv_real_pos1 = (el1_real * e3_real) - (el1_imag * e3_imag);
    ValType sv_imag_pos1 = (el1_real * e3_imag) + (el1_imag * e3_real);
    STORE_Q0;
    OP_TAIL;
}
//============== SX Gate ================
inline void SX_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    OP_HEAD;
    LOAD_Q0;
    ValType sv_real_pos0 = HALF * (el0_real - el0_imag + el1_real + el1_imag);
    ValType sv_imag_pos0 = HALF * (el0_imag + el0_real + el1_imag - el1_real);
    ValType sv_real_pos1 = HALF * (el0_real + el0_imag + el1_real - el1_imag);
    ValType sv_imag_pos1 = HALF * (el0_imag - el0_real + el1_imag + el1_real);
    STORE_Q0;
    OP_TAIL;
}
//============== ConjugateSX Gate ================
inline void ConjugateSX_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    OP_HEAD;
    LOAD_Q0;
    ValType sv_real_pos0 = HALF * (el0_real + el0_imag + el1_real - el1_imag);
    ValType sv_imag_pos0 = HALF * (el0_imag - el0_real + el1_imag + el1_real);
    ValType sv_real_pos1 = HALF * (el0_real - el0_imag + el1_real + el1_imag);
    ValType sv_imag_pos1 = HALF * (el0_imag + el0_real + el1_imag - el1_real);
    STORE_Q0;
    OP_TAIL;
}

//============== P Gate ================
inline void P_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const IdxType qubit)
{
    ValType e3_real = cos(theta);
    ValType e3_imag = sin(theta);
    OP_HEAD;
    const ValType el1_real = GET(sv_real, pos1);
    const ValType el1_imag = GET(sv_imag, pos1);
    ValType sv_real_pos1 = (e3_real * el1_real) - (e3_imag * el1_imag);
    ValType sv_imag_pos1 = (e3_real * el1_imag) + (e3_imag * el1_real);
    PUT(sv_real, pos1, sv_real_pos1);
    PUT(sv_imag, pos1, sv_imag_pos1);
    OP_TAIL;
}

//============== ConjugateP Gate ================
inline void ConjugateP_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag,
        const ValType theta, const IdxType qubit)
{
    ValType e3_real = cos(theta);
    ValType e3_imag = -sin(theta);
    OP_HEAD;
    const ValType el1_real = GET(sv_real, pos1);
    const ValType el1_imag = GET(sv_imag, pos1);
    ValType sv_real_pos1 = (e3_real * el1_real) - (e3_imag * el1_imag);
    ValType sv_imag_pos1 = (e3_real * el1_imag) + (e3_imag * el1_real);
    PUT(sv_real, pos1, sv_real_pos1);
    PUT(sv_imag, pos1, sv_imag_pos1);
    OP_TAIL;
}


//============== U Gate ================
inline void U_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const ValType phi, const ValType lambda, const IdxType qubit)
{
    ValType e0_real = cos(HALF*theta);
    ValType e0_imag = 0;
    ValType e1_real = -cos(lambda)*sin(HALF*theta);
    ValType e1_imag = -sin(lambda)*sin(HALF*theta);
    ValType e2_real = cos(phi)*sin(HALF*theta);
    ValType e2_imag = sin(phi)*sin(HALF*theta);
    ValType e3_real = cos(phi+lambda)*cos(HALF*theta);
    ValType e3_imag = sin(phi+lambda)*cos(HALF*theta);
    OP_HEAD;
    LOAD_Q0;
    ValType sv_real_pos0 =  (e0_real * el0_real) - (e0_imag * el0_imag)
                           +(e1_real * el1_real) - (e1_imag * el1_imag);
    ValType sv_imag_pos0 =  (e0_real * el0_imag) + (e0_imag * el0_real)
                           +(e1_real * el1_imag) + (e1_imag * el1_real);
    ValType sv_real_pos1 =  (e2_real * el0_real) - (e2_imag * el0_imag)
                           +(e3_real * el1_real) - (e3_imag * el1_imag);
    ValType sv_imag_pos1 =  (e2_real * el0_imag) + (e2_imag * el0_real)
                           +(e3_real * el1_imag) + (e3_imag * el1_real);
    STORE_Q0;
    OP_TAIL;
}
//============== ConjugateU Gate ================
inline void ConjugateU_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta, const ValType phi, const ValType lambda, const IdxType qubit)
{
    ValType e0_real = cos(HALF*theta);
    ValType e0_imag = 0;
    ValType e1_real = -cos(lambda)*sin(HALF*theta);
    ValType e1_imag = sin(lambda)*sin(HALF*theta);
    ValType e2_real = cos(phi)*sin(HALF*theta);
    ValType e2_imag = -sin(phi)*sin(HALF*theta);
    ValType e3_real = cos(phi+lambda)*cos(HALF*theta);
    ValType e3_imag = -sin(phi+lambda)*cos(HALF*theta);
    OP_HEAD;
    LOAD_Q0;
    ValType sv_real_pos0 =  (e0_real * el0_real) - (e0_imag * el0_imag)
                           +(e1_real * el1_real) - (e1_imag * el1_imag);
    ValType sv_imag_pos0 =  (e0_real * el0_imag) + (e0_imag * el0_real)
                           +(e1_real * el1_imag) + (e1_imag * el1_real);
    ValType sv_real_pos1 =  (e2_real * el0_real) - (e2_imag * el0_imag)
                           +(e3_real * el1_real) - (e3_imag * el1_imag);
    ValType sv_imag_pos1 =  (e2_real * el0_imag) + (e2_imag * el0_real)
                           +(e3_real * el1_imag) + (e3_imag * el1_real);
    STORE_Q0;
    OP_TAIL;
}
//============== CX Gate ================
inline void CX_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    LOAD_Q1;
    ValType sv_real_pos2 = el3_real;
    ValType sv_imag_pos2 = el3_imag;
    ValType sv_real_pos3 = el2_real;
    ValType sv_imag_pos3 = el2_imag;
    STORE_Q1;
    OP_TAIL;
}
//============== CY Gate ================
inline void CY_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    LOAD_Q1;
    ValType sv_real_pos2 =  el3_imag; 
    ValType sv_imag_pos2 = -el3_real;
    ValType sv_real_pos3 = -el2_imag;
    ValType sv_imag_pos3 =  el2_real;
    STORE_Q1;
    OP_TAIL;
}
//============== Conjugate CY Gate ================
inline void ConjugateCY_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    LOAD_Q1;
    ValType sv_real_pos2 = -el3_imag; 
    ValType sv_imag_pos2 =  el3_real;
    ValType sv_real_pos3 =  el2_imag;
    ValType sv_imag_pos3 = -el2_real;
    STORE_Q1;
    OP_TAIL;
}
//============== CZ Gate ================
inline void CZ_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    const ValType el3_real = GET(sv_real, pos3);
    const ValType el3_imag = GET(sv_imag, pos3);
    ValType sv_real_pos3 = -el3_real;
    ValType sv_imag_pos3 = -el3_imag;
    PUT(sv_real, pos3, sv_real_pos3);
    PUT(sv_imag, pos3, sv_imag_pos3);
    OP_TAIL;
}
//============== CH Gate ================
inline void CH_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    LOAD_Q1;
    ValType sv_real_pos2 = S2I*(el2_real + el3_real); 
    ValType sv_imag_pos2 = S2I*(el2_imag + el3_imag);
    ValType sv_real_pos3 = S2I*(el2_real - el3_real);
    ValType sv_imag_pos3 = S2I*(el2_imag - el3_imag);
    STORE_Q1;
    OP_TAIL;
}
//============== CS Gate ================
inline void CS_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    const ValType el3_real = GET(sv_real, pos3);
    const ValType el3_imag = GET(sv_imag, pos3);
    ValType sv_real_pos3 = -el3_imag;
    ValType sv_imag_pos3 =  el3_real;
    PUT(sv_real, pos3, sv_real_pos3);
    PUT(sv_imag, pos3, sv_imag_pos3);
    OP_TAIL;
}
//============== CSDG Gate ================
inline void CSDG_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    const ValType el3_real = GET(sv_real, pos3);
    const ValType el3_imag = GET(sv_imag, pos3);
    ValType sv_real_pos3 =  el3_imag;
    ValType sv_imag_pos3 = -el3_real;
    PUT(sv_real, pos3, sv_real_pos3);
    PUT(sv_imag, pos3, sv_imag_pos3);
    OP_TAIL;
}
//============== CT Gate ================
inline void CT_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    const ValType el3_real = GET(sv_real, pos3);
    const ValType el3_imag = GET(sv_imag, pos3);
    ValType sv_real_pos3 = S2I*(el3_real-el3_imag);
    ValType sv_imag_pos3 = S2I*(el3_real+el3_imag);
    PUT(sv_real, pos3, sv_real_pos3);
    PUT(sv_imag, pos3, sv_imag_pos3);
    OP_TAIL;
}
//============== CTDG Gate ================
inline void CTDG_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    const ValType el3_real = GET(sv_real, pos3);
    const ValType el3_imag = GET(sv_imag, pos3);
    ValType sv_real_pos3 = S2I*( el3_real+el3_imag);
    ValType sv_imag_pos3 = S2I*(-el3_real+el3_imag);
    PUT(sv_real, pos3, sv_real_pos3);
    PUT(sv_imag, pos3, sv_imag_pos3);
    OP_TAIL;
}
//============== CRI Gate ================
inline void CRI_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const IdxType ctrl, const IdxType qubit)
{
    ValType ri_real =  cos(HALF*theta);
    ValType ri_imag = -sin(HALF*theta);
    OP_HEAD_2Q;
    LOAD_Q1;
    ValType sv_real_pos2 = (el2_real * ri_real) - (el2_imag * ri_imag);
    ValType sv_imag_pos2 = (el2_real * ri_imag) + (el2_imag * ri_real);
    ValType sv_real_pos3 = (el3_real * ri_real) - (el3_imag * ri_imag);
    ValType sv_imag_pos3 = (el3_real * ri_imag) + (el3_imag * ri_real);
    STORE_Q1;
    OP_TAIL;
}
//============== ConjugateCRI Gate ================
inline void ConjugateCRI_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const IdxType ctrl, const IdxType qubit)
{
    ValType ri_real = cos(HALF*theta);
    ValType ri_imag = sin(HALF*theta);
    OP_HEAD_2Q;
    LOAD_Q1;
    ValType sv_real_pos2 = (el2_real * ri_real) - (el2_imag * ri_imag);
    ValType sv_imag_pos2 = (el2_real * ri_imag) + (el2_imag * ri_real);
    ValType sv_real_pos3 = (el3_real * ri_real) - (el3_imag * ri_imag);
    ValType sv_imag_pos3 = (el3_real * ri_imag) + (el3_imag * ri_real);
    STORE_Q1;
    OP_TAIL;
}
//============== CRX Gate ================
inline void CRX_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const IdxType ctrl, const IdxType qubit)
{
    ValType rx_real =  cos(HALF*theta);
    ValType rx_imag = -sin(HALF*theta);
    OP_HEAD_2Q;
    LOAD_Q1;
    ValType sv_real_pos2 =  (rx_real * el2_real) - (rx_imag * el3_imag);
    ValType sv_imag_pos2 =  (rx_real * el2_imag) + (rx_imag * el3_real);
    ValType sv_real_pos3 = -(rx_imag * el2_imag) + (rx_real * el3_real);
    ValType sv_imag_pos3 =  (rx_imag * el2_real) + (rx_real * el3_imag);
    STORE_Q1;
    OP_TAIL;
}
//============== ConjugateCRX Gate ================
inline void ConjugateCRX_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const IdxType ctrl, const IdxType qubit)
{
    ValType rx_real = cos(HALF*theta);
    ValType rx_imag = sin(HALF*theta);
    OP_HEAD_2Q;
    LOAD_Q1;
    ValType sv_real_pos2 =  (rx_real * el2_real) - (rx_imag * el3_imag);
    ValType sv_imag_pos2 =  (rx_real * el2_imag) + (rx_imag * el3_real);
    ValType sv_real_pos3 = -(rx_imag * el2_imag) + (rx_real * el3_real);
    ValType sv_imag_pos3 =  (rx_imag * el2_real) + (rx_real * el3_imag);
    STORE_Q1;
    OP_TAIL;
}
//============== CRY Gate ================
inline void CRY_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const IdxType ctrl, const IdxType qubit)
{
    ValType e0_real =  cos(HALF*theta);
    ValType e1_real = -sin(HALF*theta);
    ValType e2_real =  sin(HALF*theta);
    ValType e3_real =  cos(HALF*theta);
    OP_HEAD_2Q;
    LOAD_Q1;
    ValType sv_real_pos2 = (e0_real * el2_real) + (e1_real * el3_real);
    ValType sv_imag_pos2 = (e0_real * el2_imag) + (e1_real * el3_imag);
    ValType sv_real_pos3 = (e2_real * el2_real) + (e3_real * el3_real);
    ValType sv_imag_pos3 = (e2_real * el2_imag) + (e3_real * el3_imag);
    STORE_Q1;
    OP_TAIL;
}
//============== CRZ Gate ================
inline void CRZ_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const IdxType ctrl, const IdxType qubit)
{
    ValType e0_real =  cos(HALF*theta);
    ValType e0_imag = -sin(HALF*theta);
    ValType e3_real =  cos(HALF*theta);
    ValType e3_imag =  sin(HALF*theta);
    OP_HEAD_2Q;
    LOAD_Q1;
    ValType sv_real_pos2 = (e0_real * el2_real) - (e0_imag * el2_imag);
    ValType sv_imag_pos2 = (e0_imag * el2_real) + (e0_real * el2_imag);
    ValType sv_real_pos3 = (e3_real * el3_real) - (e3_imag * el3_imag);
    ValType sv_imag_pos3 = (e3_imag * el3_real) + (e3_real * el3_imag);
    STORE_Q1;
    OP_TAIL;
}
//============== ConjugateCRZ Gate ================
inline void ConjugateCRZ_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const IdxType ctrl, const IdxType qubit)
{
    ValType e0_real =  cos(HALF*theta);
    ValType e0_imag =  sin(HALF*theta);
    ValType e3_real =  cos(HALF*theta);
    ValType e3_imag = -sin(HALF*theta);
    OP_HEAD_2Q;
    LOAD_Q1;
    ValType sv_real_pos2 = (e0_real * el2_real) - (e0_imag * el2_imag);
    ValType sv_imag_pos2 = (e0_imag * el2_real) + (e0_real * el2_imag);
    ValType sv_real_pos3 = (e3_real * el3_real) - (e3_imag * el3_imag);
    ValType sv_imag_pos3 = (e3_imag * el3_real) + (e3_real * el3_imag);
    STORE_Q1;
    OP_TAIL;
}
//============== CSX Gate ================
inline void CSX_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    LOAD_Q1;
    ValType sv_real_pos2 = HALF * (el2_real - el2_imag + el3_real + el3_imag);
    ValType sv_imag_pos2 = HALF * (el2_imag + el2_real + el3_imag - el3_real);
    ValType sv_real_pos3 = HALF * (el2_real + el2_imag + el3_real - el3_imag);
    ValType sv_imag_pos3 = HALF * (el2_imag - el2_real + el3_imag + el3_real);
    STORE_Q1;
    OP_TAIL;
}
//============== ConjugateCSX Gate ================
inline void ConjugateCSX_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    LOAD_Q1;
    ValType sv_real_pos2 = HALF * (el2_real + el2_imag + el3_real - el3_imag);
    ValType sv_imag_pos2 = HALF * (el2_imag - el2_real + el3_imag + el3_real);
    ValType sv_real_pos3 = HALF * (el2_real - el2_imag + el3_real + el3_imag);
    ValType sv_imag_pos3 = HALF * (el2_imag + el2_real + el3_imag - el3_real);
    STORE_Q1;
    OP_TAIL;
}
//============== CP Gate ================
inline void CP_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType theta,
        const IdxType ctrl, const IdxType qubit)
{
    ValType e3_real = cos(theta);
    ValType e3_imag = sin(theta);
    OP_HEAD_2Q;
    const ValType el3_real = GET(sv_real, pos1);
    const ValType el3_imag = GET(sv_imag, pos1);
    ValType sv_real_pos3 = (e3_real * el3_real) - (e3_imag * el3_imag);
    ValType sv_imag_pos3 = (e3_real * el3_imag) + (e3_imag * el3_real);
    PUT(sv_real, pos3, sv_real_pos3);
    PUT(sv_imag, pos3, sv_imag_pos3);
    OP_TAIL;
}
//============== ConjugateCP Gate ================
inline void ConjugateCP_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const IdxType ctrl, const IdxType qubit)
{
    ValType e3_real = cos(theta);
    ValType e3_imag = -sin(theta);
    OP_HEAD_2Q;
    const ValType el3_real = GET(sv_real, pos1);
    const ValType el3_imag = GET(sv_imag, pos1);
    ValType sv_real_pos3 = (e3_real * el3_real) - (e3_imag * el3_imag);
    ValType sv_imag_pos3 = (e3_real * el3_imag) + (e3_imag * el3_real);
    PUT(sv_real, pos3, sv_real_pos3);
    PUT(sv_imag, pos3, sv_imag_pos3);
    OP_TAIL;
}

//============== CU Gate ================
inline void CU_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const ValType phi, const ValType lambda, 
        const IdxType ctrl, const IdxType qubit)
{
    ValType e0_real = cos(HALF*theta);
    ValType e0_imag = 0;
    ValType e1_real = -cos(lambda)*sin(HALF*theta);
    ValType e1_imag = -sin(lambda)*sin(HALF*theta);
    ValType e2_real = cos(phi)*sin(HALF*theta);
    ValType e2_imag = sin(phi)*sin(HALF*theta);
    ValType e3_real = cos(phi+lambda)*cos(HALF*theta);
    ValType e3_imag = sin(phi+lambda)*cos(HALF*theta);
    OP_HEAD_2Q;
    LOAD_Q1;
    ValType sv_real_pos2 =  (e0_real * el2_real) - (e0_imag * el2_imag)
                           +(e1_real * el3_real) - (e1_imag * el3_imag);
    ValType sv_imag_pos2 =  (e0_real * el2_imag) + (e0_imag * el2_real)
                           +(e1_real * el3_imag) + (e1_imag * el3_real);
    ValType sv_real_pos3 =  (e2_real * el2_real) - (e2_imag * el2_imag)
                           +(e3_real * el3_real) - (e3_imag * el3_imag);
    ValType sv_imag_pos3 =  (e2_real * el2_imag) + (e2_imag * el2_real)
                           +(e3_real * el3_imag) + (e3_imag * el3_real);
    STORE_Q1;
    OP_TAIL;
}
//============== ConjugateCU Gate ================
inline void ConjugateCU_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType theta, const ValType phi, const ValType lambda, 
        const IdxType ctrl, const IdxType qubit)
{
    ValType e0_real = cos(HALF*theta);
    ValType e0_imag = 0;
    ValType e1_real = -cos(lambda)*sin(HALF*theta);
    ValType e1_imag = sin(lambda)*sin(HALF*theta);
    ValType e2_real = cos(phi)*sin(HALF*theta);
    ValType e2_imag = -sin(phi)*sin(HALF*theta);
    ValType e3_real = cos(phi+lambda)*cos(HALF*theta);
    ValType e3_imag = -sin(phi+lambda)*cos(HALF*theta);
    OP_HEAD_2Q;
    LOAD_Q1;
    ValType sv_real_pos2 =  (e0_real * el2_real) - (e0_imag * el2_imag)
                           +(e1_real * el3_real) - (e1_imag * el3_imag);
    ValType sv_imag_pos2 =  (e0_real * el2_imag) + (e0_imag * el2_real)
                           +(e1_real * el3_imag) + (e1_imag * el3_real);
    ValType sv_real_pos3 =  (e2_real * el2_real) - (e2_imag * el2_imag)
                           +(e3_real * el3_real) - (e3_imag * el3_imag);
    ValType sv_imag_pos3 =  (e2_real * el2_imag) + (e2_imag * el2_real)
                           +(e3_real * el3_imag) + (e3_imag * el3_real);
    STORE_Q1;
    OP_TAIL;
}
//============== ID Gate ================
inline void ID_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
}
//============== SWAP Gate ================
inline void SWAP_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    const ValType el1_real = GET(sv_real, pos1);
    const ValType el1_imag = GET(sv_imag, pos1);
    const ValType el2_real = GET(sv_real, pos2);
    const ValType el2_imag = GET(sv_imag, pos2);
    ValType sv_real_pos1 = el2_real; 
    ValType sv_imag_pos1 = el2_imag; 
    ValType sv_real_pos2 = el1_real;
    ValType sv_imag_pos2 = el1_imag;
    PUT(sv_real, pos1, sv_real_pos1);
    PUT(sv_imag, pos1, sv_imag_pos1);
    PUT(sv_real, pos2, sv_real_pos2);
    PUT(sv_imag, pos2, sv_imag_pos2);
    OP_TAIL;
}
//============== C2 Gate ================
//Arbitrary 2-qubit gate
inline void C2_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag,
        const IdxType ctrl, const IdxType qubit)
{
    OP_HEAD_2Q;
    const ValType el0_real = GET(sv_real, pos0);
    const ValType el0_imag = GET(sv_imag, pos0);
    const ValType el1_real = GET(sv_real, pos1);
    const ValType el1_imag = GET(sv_imag, pos1);
    const ValType el2_real = GET(sv_real, pos2);
    const ValType el2_imag = GET(sv_imag, pos2);
    const ValType el3_real = GET(sv_real, pos3);
    const ValType el3_imag = GET(sv_imag, pos3);
    //Real part
    ValType sv_real_pos0 = (gm_real[ 0] * el0_real) - (gm_imag[ 0] * el0_imag)
                          +(gm_real[ 1] * el1_real) - (gm_imag[ 1] * el1_imag)
                          +(gm_real[ 2] * el2_real) - (gm_imag[ 2] * el2_imag)
                          +(gm_real[ 3] * el3_real) - (gm_imag[ 3] * el3_imag);
    ValType sv_real_pos1 = (gm_real[ 4] * el0_real) - (gm_imag[ 4] * el0_imag)
                          +(gm_real[ 5] * el1_real) - (gm_imag[ 5] * el1_imag)
                          +(gm_real[ 6] * el2_real) - (gm_imag[ 6] * el2_imag)
                          +(gm_real[ 7] * el3_real) - (gm_imag[ 7] * el3_imag);
    ValType sv_real_pos2 = (gm_real[ 8] * el0_real) - (gm_imag[ 8] * el0_imag)
                          +(gm_real[ 9] * el1_real) - (gm_imag[ 9] * el1_imag)
                          +(gm_real[10] * el2_real) - (gm_imag[10] * el2_imag)
                          +(gm_real[11] * el3_real) - (gm_imag[11] * el3_imag);
    ValType sv_real_pos3 = (gm_real[12] * el0_real) - (gm_imag[12] * el0_imag)
                          +(gm_real[13] * el1_real) - (gm_imag[13] * el1_imag)
                          +(gm_real[14] * el2_real) - (gm_imag[14] * el2_imag)
                          +(gm_real[15] * el3_real) - (gm_imag[15] * el3_imag);
    //Imag part
    ValType sv_imag_pos0 = (gm_real[ 0] * el0_imag) + (gm_imag[ 0] * el0_real)
                          +(gm_real[ 1] * el1_imag) + (gm_imag[ 1] * el1_real)
                          +(gm_real[ 2] * el2_imag) + (gm_imag[ 2] * el2_real)
                          +(gm_real[ 3] * el3_imag) + (gm_imag[ 3] * el3_real);
    ValType sv_imag_pos1 = (gm_real[ 4] * el0_imag) + (gm_imag[ 4] * el0_real)
                          +(gm_real[ 5] * el1_imag) + (gm_imag[ 5] * el1_real)
                          +(gm_real[ 6] * el2_imag) + (gm_imag[ 6] * el2_real)
                          +(gm_real[ 7] * el3_imag) + (gm_imag[ 7] * el3_real);
    ValType sv_imag_pos2 = (gm_real[ 8] * el0_imag) + (gm_imag[ 8] * el0_real)
                          +(gm_real[ 9] * el1_imag) + (gm_imag[ 9] * el1_real)
                          +(gm_real[10] * el2_imag) + (gm_imag[10] * el2_real)
                          +(gm_real[11] * el3_imag) + (gm_imag[11] * el3_real);
    ValType sv_imag_pos3 = (gm_real[12] * el0_imag) + (gm_imag[12] * el0_real)
                          +(gm_real[13] * el1_imag) + (gm_imag[13] * el1_real)
                          +(gm_real[14] * el2_imag) + (gm_imag[14] * el2_real)
                          +(gm_real[15] * el3_imag) + (gm_imag[15] * el3_real);
    PUT(sv_real, pos0, sv_real_pos0);
    PUT(sv_imag, pos0, sv_imag_pos0);
    PUT(sv_real, pos1, sv_real_pos1);
    PUT(sv_imag, pos1, sv_imag_pos1);
    PUT(sv_real, pos2, sv_real_pos2);
    PUT(sv_imag, pos2, sv_imag_pos2);
    PUT(sv_real, pos3, sv_real_pos3);
    PUT(sv_imag, pos3, sv_imag_pos3);
    OP_TAIL;
}
//============== C4 Gate ================
//Arbitrary 4-qubit gate
inline void C4_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag,
        const IdxType qubit0, const IdxType qubit1,
        const IdxType qubit2, const IdxType qubit3)
{
    assert (qubit0 != qubit1); //Non-cloning
    assert (qubit0 != qubit2); //Non-cloning
    assert (qubit0 != qubit3); //Non-cloning
    assert (qubit1 != qubit2); //Non-cloning
    assert (qubit1 != qubit3); //Non-cloning
    assert (qubit2 != qubit3); //Non-cloning
    //need to sort qubits: min->max: p, q, r, s
    const IdxType v0 = min(qubit0, qubit1);
    const IdxType v1 = min(qubit2, qubit3);
    const IdxType v2 = max(qubit0, qubit1);
    const IdxType v3 = max(qubit2, qubit3);
    const IdxType p = min(v0,v1); 
    const IdxType q = min(min(v2,v3),max(v0,v1)); 
    const IdxType r = max(min(v2,v3),max(v0,v1)); 
    const IdxType s = max(v2,v3);
    _Pragma("omp for schedule(auto)") 
    for (IdxType i=0; i<((sim->dim)>>4); i++)
    {
        const IdxType term0 = MOD2E(i,p);
        const IdxType term1 = MOD2E(DIV2E(i,p),q-p-1)*EXP2E(p+1);
        const IdxType term2 = MOD2E(DIV2E(DIV2E(i,p),q-p-1),r-q-1)*EXP2E(q+1);
        const IdxType term3 = MOD2E(DIV2E(DIV2E(DIV2E(i,p),q-p-1),r-q-1),s-r-1)*EXP2E(r+1);
        const IdxType term4 = DIV2E(DIV2E(DIV2E(DIV2E(i,p),q-p-1),r-q-1),s-r-1)*EXP2E(s+1);
        const IdxType term = term4 + term3 + term2 + term1 + term0;
        const ValType el_real[16] = { 
            GET(sv_real,term+SV16IDX(0)),  GET(sv_real,term+SV16IDX(1)),
            GET(sv_real,term+SV16IDX(2)),  GET(sv_real,term+SV16IDX(3)),
            GET(sv_real,term+SV16IDX(4)),  GET(sv_real,term+SV16IDX(5)),
            GET(sv_real,term+SV16IDX(6)),  GET(sv_real,term+SV16IDX(7)),
            GET(sv_real,term+SV16IDX(8)),  GET(sv_real,term+SV16IDX(9)),
            GET(sv_real,term+SV16IDX(10)), GET(sv_real,term+SV16IDX(11)),
            GET(sv_real,term+SV16IDX(12)), GET(sv_real,term+SV16IDX(13)),
            GET(sv_real,term+SV16IDX(14)), GET(sv_real,term+SV16IDX(15))
        };
        const ValType el_imag[16] = { 
            GET(sv_imag,term+SV16IDX(0)),  GET(sv_imag,term+SV16IDX(1)),
            GET(sv_imag,term+SV16IDX(2)),  GET(sv_imag,term+SV16IDX(3)),
            GET(sv_imag,term+SV16IDX(4)),  GET(sv_imag,term+SV16IDX(5)),
            GET(sv_imag,term+SV16IDX(6)),  GET(sv_imag,term+SV16IDX(7)),
            GET(sv_imag,term+SV16IDX(8)),  GET(sv_imag,term+SV16IDX(9)),
            GET(sv_imag,term+SV16IDX(10)), GET(sv_imag,term+SV16IDX(11)),
            GET(sv_imag,term+SV16IDX(12)), GET(sv_imag,term+SV16IDX(13)),
            GET(sv_imag,term+SV16IDX(14)), GET(sv_imag,term+SV16IDX(15))
        };
        #pragma unroll
        for (unsigned j=0; j<16; j++)
        {
            ValType res_real = 0;
            ValType res_imag = 0;
            #pragma unroll
            for (unsigned k=0; k<16; k++)
            {
                res_real += (el_real[k] * gm_real[j*16+k]) - (el_imag[k] * gm_imag[j*16+k]);
                res_imag += (el_real[k] * gm_imag[j*16+k]) + (el_imag[k] * gm_real[j*16+k]);
            }
            PUT(sv_real, term+SV16IDX(j), res_real);
            PUT(sv_imag, term+SV16IDX(j), res_imag);
        }
    }
    _Pragma("omp barrier")  
}
#endif

//============== M Gate (Measure 1 qubit in Pauli-Z) ================
inline void M_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const ValType rand, const IdxType qubit)
{
    ValType * m_real = sim->m_real;
    IdxType mask = ((IdxType)1<<qubit);
    #pragma omp for schedule(auto) 
    for (IdxType i=0; i<((IdxType)1<<(sim->n_qubits)); i++)
    {
        if ((i&mask)==0) m_real[i] = 0.;
        else m_real[i] = abs(sv_real[(i<<(sim->n_qubits))+i]); //diagonal element
    }
    #pragma omp barrier
    for (IdxType k=(sim->half_dim); k>0; k>>=1)
    {
        #pragma omp for schedule(auto) 
        for (IdxType i=0; i<k; i++)
        {
            m_real[i] += m_real[i+k];
        }
        #pragma omp barrier
    }
    ValType prob_of_one = m_real[0];
    #pragma omp barrier
    if (rand < prob_of_one)
    {
        ValType normalize_factor = prob_of_one;
        #pragma omp for schedule(auto) 
        for (IdxType i=0; i<sim->dim; i++)
        {
            if ((i&mask) == 0)
            {
                sv_real[i] = 0.;
                sv_imag[i] = 0.;
            }
            else
            {
                sv_real[i] /= normalize_factor;
                sv_imag[i] /= normalize_factor;
            }
        }
    }
    else
    {
        #pragma omp for schedule(auto) 
        for (IdxType i=0; i<sim->dim; i++)
        {
            ValType normalize_factor = 1.0-prob_of_one;
            if ((i&mask) == 0)
            {
                sv_real[i] /= normalize_factor;
                sv_imag[i] /= normalize_factor;
            }
            else
            {
                sv_real[i] = 0;
                sv_imag[i] = 0;
            }
        }
    }
    #pragma omp barrier
    if (omp_get_thread_num()==0)
        sim->results[0] = (rand<prob_of_one?1:0);
}
//============== MA Gate (Measure all qubits in Pauli-Z) ================
inline void MA_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType repetition)
{
    ValType * m_real = sim->m_real;
    IdxType n_size = (IdxType)1<<(sim->n_qubits);
    #pragma omp for schedule(auto)
    for (IdxType i=0; i<n_size; i++)
        m_real[i] = abs(sv_real[(i<<(sim->n_qubits))+i]); //diagonal element
    //Parallel prefix sum
    for (IdxType d=0; d<(sim->n_qubits); d++)
    {
        #pragma omp for schedule(auto)
        for (IdxType k=0; k<n_size; k+=((IdxType)1<<(d+1)))
            m_real[k+(1<<(d+1))-1] = m_real[k+(1<<d)-1] + m_real[k+(1<<(d+1))-1];
        #pragma omp barrier
    }
    if (omp_get_thread_num() == 0)
    {
        m_real[n_size] = m_real[n_size-1];
        m_real[n_size-1] = 0;
        ValType purity = fabs(m_real[n_size]);
        if ( abs(purity - 1.0) > ERROR_BAR )
            printf("MA: Purity Check fails with %lf\n", purity);
    }
    #pragma omp barrier
    for (long long d=(sim->n_qubits)-1; d>=0; d--)
    {
        #pragma omp for schedule(auto)
        for (IdxType k=0; k<n_size-1; k+=((IdxType)1<<(d+1)))
        {
            ValType tmp = m_real[k+((IdxType)1<<d)-1];
            m_real[k+((IdxType)1<<d)-1] = m_real[k+((IdxType)1<<(d+1))-1];
            m_real[k+((IdxType)1<<(d+1))-1] = tmp + m_real[k+((IdxType)1<<(d+1))-1];
        }
        #pragma omp barrier
    }
    #pragma omp for schedule(auto)
    for (IdxType i=0; i<repetition; i++)
    {
        ValType r = sim->randoms[i];
        for (IdxType j=0; j<((IdxType)1<<(sim->n_qubits)); j++)
            if (m_real[j]<=r && r<m_real[j+1])
                sim->results[i] = j;
    }
}
 //============== Reset ================
inline void RESET_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, const IdxType qubit)
{
    ValType * m_real = sim->m_real;
    IdxType mask = ((IdxType)1<<qubit);
    #pragma omp for schedule(auto) 
    for (IdxType i=0; i<((IdxType)1<<(sim->n_qubits)); i++)
    {
        if ((i&mask)==0) m_real[i] = 0.;
        else m_real[i] = abs(sv_real[(i<<(sim->n_qubits))+i]); //diagonal element
    }
    #pragma omp barrier
    for (IdxType k=(sim->half_dim); k>0; k>>=1)
    {
        #pragma omp for schedule(auto) 
        for (IdxType i=0; i<k; i++)
        {
            m_real[i] += m_real[i+k];
        }
        #pragma omp barrier
    }
    ValType prob = m_real[0];
    #pragma omp barrier
    ValType normalize_factor = 1.0-prob;
    #pragma omp for schedule(auto) 
    for (IdxType i=0; i<sim->dim; i++)
    {
        if ((i&mask) == 0)
        {
            sv_real[i] /= normalize_factor;
            sv_imag[i] /= normalize_factor;
        }
        else
        {
            sv_real[i] = 0;
            sv_imag[i] = 0;
        }
    }
}
/***********************************************
 * Operation Implementation
 ***********************************************/
inline void X_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    X_GATE(sim, sv_real, sv_imag, g->qubit); 
    X_GATE(sim, sv_real, sv_imag, (g->qubit)+(sim->n_qubits));
    if (sim->noise_1q_real != NULL) C2_GATE(sim, sv_real, sv_imag, sim->noise_1q_real, sim->noise_1q_imag, g->qubit, (g->qubit)+(sim->n_qubits));
}
inline void Y_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    Y_GATE(sim, sv_real, sv_imag, g->qubit); 
    ConjugateY_GATE(sim, sv_real, sv_imag, (g->qubit)+(sim->n_qubits));
    if (sim->noise_1q_real != NULL) C2_GATE(sim, sv_real, sv_imag, sim->noise_1q_real, sim->noise_1q_imag, g->qubit, (g->qubit)+(sim->n_qubits));
}
inline void Z_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    Z_GATE(sim, sv_real, sv_imag, g->qubit); 
    Z_GATE(sim, sv_real, sv_imag, (g->qubit)+(sim->n_qubits));
    if (sim->noise_1q_real != NULL) C2_GATE(sim, sv_real, sv_imag, sim->noise_1q_real, sim->noise_1q_imag, g->qubit, (g->qubit)+(sim->n_qubits));
}
inline void H_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    H_GATE(sim, sv_real, sv_imag, g->qubit); 
    H_GATE(sim, sv_real, sv_imag, (g->qubit)+(sim->n_qubits));
    if (sim->noise_1q_real != NULL) C2_GATE(sim, sv_real, sv_imag, sim->noise_1q_real, sim->noise_1q_imag, g->qubit, (g->qubit)+(sim->n_qubits));
}
inline void S_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    S_GATE(sim, sv_real, sv_imag, g->qubit); 
    SDG_GATE(sim, sv_real, sv_imag, (g->qubit)+(sim->n_qubits));
    if (sim->noise_1q_real != NULL) C2_GATE(sim, sv_real, sv_imag, sim->noise_1q_real, sim->noise_1q_imag, g->qubit, (g->qubit)+(sim->n_qubits));
}
inline void SDG_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    SDG_GATE(sim, sv_real, sv_imag, g->qubit); 
    S_GATE(sim, sv_real, sv_imag, (g->qubit)+(sim->n_qubits));
    if (sim->noise_1q_real != NULL) C2_GATE(sim, sv_real, sv_imag, sim->noise_1q_real, sim->noise_1q_imag, g->qubit, (g->qubit)+(sim->n_qubits));
}
inline void T_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    T_GATE(sim, sv_real, sv_imag, g->qubit); 
    TDG_GATE(sim, sv_real, sv_imag, (g->qubit)+(sim->n_qubits));
    if (sim->noise_1q_real != NULL) C2_GATE(sim, sv_real, sv_imag, sim->noise_1q_real, sim->noise_1q_imag, g->qubit, (g->qubit)+(sim->n_qubits));
}
inline void TDG_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    TDG_GATE(sim, sv_real, sv_imag, g->qubit); 
    T_GATE(sim, sv_real, sv_imag, (g->qubit)+(sim->n_qubits));
    if (sim->noise_1q_real != NULL) C2_GATE(sim, sv_real, sv_imag, sim->noise_1q_real, sim->noise_1q_imag, g->qubit, (g->qubit)+(sim->n_qubits));
}
inline void RI_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    RI_GATE(sim, sv_real, sv_imag, g->theta, g->qubit); 
    ConjugateRI_GATE(sim, sv_real, sv_imag, g->theta, (g->qubit)+(sim->n_qubits));
    if (sim->noise_1q_real != NULL) C2_GATE(sim, sv_real, sv_imag, sim->noise_1q_real, sim->noise_1q_imag, g->qubit, (g->qubit)+(sim->n_qubits));
}
inline void RX_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    RX_GATE(sim, sv_real, sv_imag, g->theta, g->qubit); 
    ConjugateRX_GATE(sim, sv_real, sv_imag, g->theta, (g->qubit)+(sim->n_qubits));
    if (sim->noise_1q_real != NULL) C2_GATE(sim, sv_real, sv_imag, sim->noise_1q_real, sim->noise_1q_imag, g->qubit, (g->qubit)+(sim->n_qubits));
}
inline void RY_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    RY_GATE(sim, sv_real, sv_imag, g->theta, g->qubit); 
    RY_GATE(sim, sv_real, sv_imag, g->theta, (g->qubit)+(sim->n_qubits));
    if (sim->noise_1q_real != NULL) C2_GATE(sim, sv_real, sv_imag, sim->noise_1q_real, sim->noise_1q_imag, g->qubit, (g->qubit)+(sim->n_qubits));
}
inline void RZ_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    RZ_GATE(sim, sv_real, sv_imag, g->theta, g->qubit); 
    ConjugateRZ_GATE(sim, sv_real, sv_imag, g->theta, (g->qubit)+(sim->n_qubits));
    if (sim->noise_1q_real != NULL) C2_GATE(sim, sv_real, sv_imag, sim->noise_1q_real, sim->noise_1q_imag, g->qubit, (g->qubit)+(sim->n_qubits));
}
inline void SX_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    SX_GATE(sim, sv_real, sv_imag, g->qubit); 
    ConjugateSX_GATE(sim, sv_real, sv_imag, (g->qubit)+(sim->n_qubits));
    if (sim->noise_1q_real != NULL) C2_GATE(sim, sv_real, sv_imag, sim->noise_1q_real, sim->noise_1q_imag, g->qubit, (g->qubit)+(sim->n_qubits));
}
inline void P_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    P_GATE(sim, sv_real, sv_imag, g->theta, g->qubit); 
    ConjugateP_GATE(sim, sv_real, sv_imag, g->theta, (g->qubit)+(sim->n_qubits));
    if (sim->noise_1q_real != NULL) C2_GATE(sim, sv_real, sv_imag, sim->noise_1q_real, sim->noise_1q_imag, g->qubit, (g->qubit)+(sim->n_qubits));
}
inline void U_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    U_GATE(sim, sv_real, sv_imag, g->theta, g->phi, g->lambda, g->qubit); 
    ConjugateU_GATE(sim, sv_real, sv_imag, g->theta, g->phi, g->lambda, (g->qubit)+(sim->n_qubits));
    if (sim->noise_1q_real != NULL) C2_GATE(sim, sv_real, sv_imag, sim->noise_1q_real, sim->noise_1q_imag, g->qubit, (g->qubit)+(sim->n_qubits));
}
inline void CX_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    CX_GATE(sim, sv_real, sv_imag, g->ctrl, g->qubit); 
    CX_GATE(sim, sv_real, sv_imag, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
    if (sim->noise_2q_real != NULL) C4_GATE(sim, sv_real, sv_imag, sim->noise_2q_real, sim->noise_2q_imag, g->ctrl, g->qubit, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
}
inline void CY_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    CY_GATE(sim, sv_real, sv_imag, g->ctrl, g->qubit); 
    ConjugateCY_GATE(sim, sv_real, sv_imag, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
    if (sim->noise_2q_real != NULL) C4_GATE(sim, sv_real, sv_imag, sim->noise_2q_real, sim->noise_2q_imag, g->ctrl, g->qubit, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
}
inline void CZ_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    CZ_GATE(sim, sv_real, sv_imag, g->ctrl, g->qubit); 
    CZ_GATE(sim, sv_real, sv_imag, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
    if (sim->noise_2q_real != NULL) C4_GATE(sim, sv_real, sv_imag, sim->noise_2q_real, sim->noise_2q_imag, g->ctrl, g->qubit, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
}
inline void CH_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    CH_GATE(sim, sv_real, sv_imag, g->ctrl, g->qubit); 
    CH_GATE(sim, sv_real, sv_imag, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
    if (sim->noise_2q_real != NULL) C4_GATE(sim, sv_real, sv_imag, sim->noise_2q_real, sim->noise_2q_imag, g->ctrl, g->qubit, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
} 
inline void CS_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    CS_GATE(sim, sv_real, sv_imag, g->ctrl, g->qubit); 
    CSDG_GATE(sim, sv_real, sv_imag, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
    if (sim->noise_2q_real != NULL) C4_GATE(sim, sv_real, sv_imag, sim->noise_2q_real, sim->noise_2q_imag, g->ctrl, g->qubit, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
} 
inline void CSDG_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    CSDG_GATE(sim, sv_real, sv_imag, g->ctrl, g->qubit); 
    CS_GATE(sim, sv_real, sv_imag, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
    if (sim->noise_2q_real != NULL) C4_GATE(sim, sv_real, sv_imag, sim->noise_2q_real, sim->noise_2q_imag, g->ctrl, g->qubit, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
} 
inline void CT_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    CT_GATE(sim, sv_real, sv_imag, g->ctrl, g->qubit); 
    CTDG_GATE(sim, sv_real, sv_imag, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
    if (sim->noise_2q_real != NULL) C4_GATE(sim, sv_real, sv_imag, sim->noise_2q_real, sim->noise_2q_imag, g->ctrl, g->qubit, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
} 
inline void CTDG_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    CTDG_GATE(sim, sv_real, sv_imag, g->ctrl, g->qubit); 
    CT_GATE(sim, sv_real, sv_imag, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
    if (sim->noise_2q_real != NULL) C4_GATE(sim, sv_real, sv_imag, sim->noise_2q_real, sim->noise_2q_imag, g->ctrl, g->qubit, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
} 
inline void CRI_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    CRI_GATE(sim, sv_real, sv_imag, g->theta, g->ctrl, g->qubit);
    ConjugateCRI_GATE(sim, sv_real, sv_imag, g->theta, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
    if (sim->noise_2q_real != NULL) C4_GATE(sim, sv_real, sv_imag, sim->noise_2q_real, sim->noise_2q_imag, g->ctrl, g->qubit, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
} 
inline void CRX_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    CRX_GATE(sim, sv_real, sv_imag, g->theta, g->ctrl, g->qubit);
    ConjugateCRX_GATE(sim, sv_real, sv_imag, g->theta, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
    if (sim->noise_2q_real != NULL) C4_GATE(sim, sv_real, sv_imag, sim->noise_2q_real, sim->noise_2q_imag, g->ctrl, g->qubit, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
} 
inline void CRY_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    CRY_GATE(sim, sv_real, sv_imag, g->theta, g->ctrl, g->qubit);
    CRY_GATE(sim, sv_real, sv_imag, g->theta, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
    if (sim->noise_2q_real != NULL) C4_GATE(sim, sv_real, sv_imag, sim->noise_2q_real, sim->noise_2q_imag, g->ctrl, g->qubit, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
} 
inline void CRZ_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    CRZ_GATE(sim, sv_real, sv_imag, g->theta, g->ctrl, g->qubit);
    ConjugateCRZ_GATE(sim, sv_real, sv_imag, g->theta, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
    if (sim->noise_2q_real != NULL) C4_GATE(sim, sv_real, sv_imag, sim->noise_2q_real, sim->noise_2q_imag, g->ctrl, g->qubit, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
} 
inline void CSX_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    CSX_GATE(sim, sv_real, sv_imag, g->ctrl, g->qubit); 
    ConjugateCSX_GATE(sim, sv_real, sv_imag, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
    if (sim->noise_2q_real != NULL) C4_GATE(sim, sv_real, sv_imag, sim->noise_2q_real, sim->noise_2q_imag, g->ctrl, g->qubit, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
}
inline void CP_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    CP_GATE(sim, sv_real, sv_imag, g->theta, g->ctrl, g->qubit); 
    ConjugateCP_GATE(sim, sv_real, sv_imag, g->theta, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
    if (sim->noise_2q_real != NULL) C4_GATE(sim, sv_real, sv_imag, sim->noise_2q_real, sim->noise_2q_imag, g->ctrl, g->qubit, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
}
inline void CU_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    CU_GATE(sim, sv_real, sv_imag, g->theta, g->phi, g->lambda, g->ctrl, g->qubit); 
    ConjugateCU_GATE(sim, sv_real, sv_imag, g->theta, g->phi, g->lambda, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
    if (sim->noise_2q_real != NULL) C4_GATE(sim, sv_real, sv_imag, sim->noise_2q_real, sim->noise_2q_imag, g->ctrl, g->qubit, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
} 
inline void ID_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    ID_GATE(sim, sv_real, sv_imag, g->qubit); 
    ID_GATE(sim, sv_real, sv_imag, (g->qubit)+(sim->n_qubits));
}
inline void SWAP_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    SWAP_GATE(sim, sv_real, sv_imag, g->ctrl, g->qubit);
    SWAP_GATE(sim, sv_real, sv_imag, (g->ctrl)+(sim->n_qubits),  (g->qubit)+(sim->n_qubits));
    if (sim->noise_2q_real != NULL) C4_GATE(sim, sv_real, sv_imag, sim->noise_2q_real, sim->noise_2q_imag, g->ctrl, g->qubit, (g->ctrl)+(sim->n_qubits), (g->qubit)+(sim->n_qubits));
} 
inline void M_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    IdxType qubit = g->qubit; 
    IdxType pauli = g->ctrl;
    ValType rand = g->theta;
    if (pauli == 1)
    {
        H_GATE(sim, sv_real, sv_imag, qubit);
        H_GATE(sim, sv_real, sv_imag, (sim->n_qubits)+qubit);
    }
    if (pauli == 2)
    {
        SDG_GATE(sim, sv_real, sv_imag, qubit);
        S_GATE(sim, sv_real, sv_imag, (sim->n_qubits)+qubit);
        H_GATE(sim, sv_real, sv_imag, qubit);
        H_GATE(sim, sv_real, sv_imag, (sim->n_qubits)+qubit);
    }
    M_GATE(sim, sv_real, sv_imag, rand, g->qubit); 
    if (pauli == 1)
    {
        H_GATE(sim, sv_real, sv_imag, qubit);
        H_GATE(sim, sv_real, sv_imag, (sim->n_qubits)+qubit);
    }
    if (pauli == 2)
    {
        H_GATE(sim, sv_real, sv_imag, qubit);
        H_GATE(sim, sv_real, sv_imag, (sim->n_qubits)+qubit);
        S_GATE(sim, sv_real, sv_imag, qubit);
        SDG_GATE(sim, sv_real, sv_imag, (sim->n_qubits)+qubit);
    }
} 
inline void MA_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    IdxType repetition = g->ctrl;
    MA_GATE(sim, sv_real, sv_imag, repetition);
} 
inline void RESET_OP(const Gate* g, const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    RESET_GATE(sim, sv_real, sv_imag, g->qubit);
    RESET_GATE(sim, sv_real, sv_imag, (g->qubit)+(sim->n_qubits));
} 

// ============================ Function Pointers ================================
func_t pX = X_OP;
func_t pY = Y_OP; 
func_t pZ = Z_OP; 
func_t pH = H_OP; 
func_t pS = S_OP;
func_t pSDG = SDG_OP;
func_t pT = T_OP;
func_t pTDG = TDG_OP;
func_t pRI = RI_OP;
func_t pRX = RX_OP;
func_t pRY = RY_OP;
func_t pRZ = RZ_OP;
func_t pSX = SX_OP;
func_t pP = P_OP;
func_t pU = U_OP;
func_t pCX = CX_OP;
func_t pCY = CY_OP;
func_t pCZ = CZ_OP;
func_t pCH = CH_OP;
func_t pCS = CS_OP;
func_t pCSDG = CSDG_OP;
func_t pCT = CT_OP;
func_t pCTDG = CTDG_OP;
func_t pCRI = CRI_OP;
func_t pCRX = CRX_OP;
func_t pCRY = CRY_OP;
func_t pCRZ = CRZ_OP;
func_t pCSX = CSX_OP;
func_t pCP = CP_OP;
func_t pCU = CU_OP;
func_t pID = ID_OP;
func_t pSWAP = SWAP_OP;
func_t pM = M_OP;
func_t pMA = MA_OP;
func_t pRESET = RESET_OP;

}; //namespace DMSim
#endif
