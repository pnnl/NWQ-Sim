// ---------------------------------------------------------------------------
// NWQsim: Northwest Quantum Circuit Simulation Environment
// ---------------------------------------------------------------------------
// Ang Li, Senior Computer Scientist
// Pacific Northwest National Laboratory(PNNL), U.S.
// Homepage: http://www.angliphd.com
// GitHub repo: http://www.github.com/pnnl/NWQ-Sim
// PNNL-IPID: 32166, ECCN: EAR99, IR: PNNL-SA-161181
// BSD Lincese.
// ---------------------------------------------------------------------------
// File: svsim_cpu_mpi.hpp
// MPI based implementation of the scale-out DV-Sim gates and 
// simulation runtime using CPU backend.
// ---------------------------------------------------------------------------
#ifndef DMSIM_CPU_SIN_HPP
#define DMSIM_CPU_SIN_HPP
#include <assert.h>
#include <random>
#include <complex.h>
#include <vector>
#include <sstream>
#include <string>
#include <iostream>
#include <cuda.h>
#include <mma.h>
#include "config.h"
#include "gate.h"
#include "device_noise.hpp"

#ifndef DISABLE_GATE_FUSION
#include "fusion.h"
#endif


namespace NWQSim
{
using namespace std;
    
class Gate;
class Simulation;
void Purity_Check(const Simulation* sim, const IdxType t, ValType* sv_real, ValType* sv_imag);

//Declare major simulation kernel
void simulation_kernel(Simulation*);

/***********************************************
 * Circuit Definition
 ***********************************************/
class Circuit
{
public:
    Circuit(IdxType _n_qubits=0):
        n_qubits(_n_qubits), n_gates(0), circuit_cpu(NULL)
    {}
    ~Circuit() { clear(); }
    void append(Gate& g)
    {
#ifdef PRINT_GATE_TRACE
        std::cout << g.gateToString() << std::flush;
#endif
        if (g.qubit >= n_qubits || ((g.op_name != MA) && (g.ctrl >= n_qubits))) 
        {
            std::cerr << g.gateToString() << std::flush;
            string msg = "Gate uses qubit out of range of " + to_string(n_qubits) + " !\n"; 
            throw std::logic_error(msg.c_str());
        }
        circuit.push_back(g);
        delete (&g);
        n_gates++;
        
    }
    void AllocateQubit() 
    { 
        n_qubits++; 
#ifdef PRINT_QUBIT_ALC_AND_RLS
        std::cout << "allocate 1 qubit, now in total:" << n_qubits << std::endl;
#endif
    }
    void ReleaseQubit()
    {
        n_qubits--;
#ifdef PRINT_QUBIT_ALC_AND_RLS
        std::cout << "release 1 qubit, now in total:" << n_qubits << std::endl;
#endif
    }
    void clear()
    {
        circuit.clear();
        n_gates = 0;
        SAFE_FREE_HOST(circuit_cpu);
    }
    void reset()
    {
        clear();
    }

    Gate* upload()
    {
#ifdef DISABLE_GATE_FUSION
        //====================== No Fuse =====================
        SAFE_FREE_HOST(circuit_cpu);
        SAFE_ALOC_HOST(circuit_cpu, n_gates*sizeof(Gate));
        memcpy(circuit_cpu, circuit.data(), n_gates*sizeof(Gate));
        //====================================================
#else
        //====================== Fuse ========================
        vector<Gate> tmp_circuit;
        tmp_circuit.clear();
        gate_fusion_1q(circuit, tmp_circuit, n_qubits);
        gate_fusion_2q(tmp_circuit, fused_circuit, n_qubits);

        this->n_gates = fused_circuit.size();
        SAFE_FREE_HOST(circuit_cpu);
        SAFE_ALOC_HOST(circuit_cpu, n_gates*sizeof(Gate));
        memcpy(circuit_cpu, fused_circuit.data(), n_gates*sizeof(Gate));
        //====================================================
#endif
        return circuit_cpu;
    }
    std::string circuitToString()
    {
        stringstream ss;
        for (IdxType t=0; t<n_gates; t++)
            ss << circuit[t].gateToString();
        return ss.str();
    }
public:
    // number of logical qubits in the circuit
    IdxType n_qubits;
    // number of physical qubits in the device
    IdxType n_phy_qubits;
    // number of gates
    IdxType n_gates;
    // user input gate sequence
    vector<Gate> circuit;
    // fused gate sequence
    vector<Gate> fused_circuit;
    // gate sequence executed
    Gate* circuit_cpu;
};
/***********************************************
 * Simulation Definition
 ***********************************************/
class Simulation
{
public:
    Simulation(IdxType _n_qubits=N_QUBIT_SLOT) : 
        n_qubits(_n_qubits),
        dim((IdxType)1<<(2*n_qubits)), 
        half_dim((IdxType)1<<(2*n_qubits-1)),
        sv_size(dim*(IdxType)sizeof(ValType)),
        diag_num(((IdxType)1<<n_qubits)+1),
        diag_size(diag_num*sizeof(ValType)),
        n_gates(0), 
        cpu_mem(0),
        sv_real(NULL),
        sv_imag(NULL),
        m_real(NULL),
        randoms(NULL),
        results(NULL)
    {
        //initialization
        device_name = DEVICE_CONFIG_NAME;
        device_config = readConfigFile(DEVICE_CONFIG_PATH, device_name);
        IdxType device_qubits = device_config["num_qubits"];
        if (device_qubits < n_qubits)
        {
            string msg = "Error: Circuit uses " + to_string(n_qubits) + " qubits, more than " 
                + to_string(device_qubits) + "qubits in the device!!\n"; 
            throw std::logic_error(msg.c_str());
        }
        if (n_qubits > N_QUBIT_SLOT) 
            throw std::invalid_argument("Requesting more qubits than threshold!");

        //CPU side initialization
        SAFE_ALOC_HOST(sv_real, sv_size);
        SAFE_ALOC_HOST(sv_imag, sv_size);
        SAFE_ALOC_HOST(m_real, diag_size);
        memset(sv_real, 0, sv_size);
        memset(sv_imag, 0, sv_size);
        memset(m_real, 0, diag_size);
        cpu_mem += sv_size*2+diag_size;
        //State-vector initial state [0..0] = 1
        sv_real[0] = 1.;
        circuit_handle = new Circuit(n_qubits);
        circuit_handle->n_phy_qubits = n_qubits;
        circuit_handle_cpu = NULL;
        //Use time as random seed
        rng.seed(time(0));
#ifdef PRINT_SIM_TRACE
        printf("DMSim_cpu is initialized!\n");
#endif
    }
    ~Simulation()
    {
        //Release circuit
        if (circuit_handle != NULL)
            delete circuit_handle;
        //Release for CPU side
        SAFE_FREE_HOST(sv_real);
        SAFE_FREE_HOST(sv_imag);
        SAFE_FREE_HOST(m_real);

        SAFE_FREE_HOST(randoms);
        SAFE_FREE_HOST(results);
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
    // =============================== Input Gates ===================================
    void X(IdxType qubit)
    {
        Gate* G = new Gate(OP::C2, qubit, -1, 0);
        //if (i_cpu == 0) printf("x q[%lld];\n",qubit);

#ifdef DISABLE_NOISE
        std::complex<double> perfect_sp1q[4][4] = {};
        backendGateNoiseSp("x", perfect_sp1q[0], qubit, -1);
        G->set_gm(perfect_sp1q[0],4);
#else
        std::complex<double> noise_sp1q[4][4] = {};
        backendGateNoiseSp(device_config, "x", noise_sp1q[0], qubit, -1);
        //if (i_cpu == 0) 
        //{
        //printf("\n======= X =======\n");
        //printGate(noise_sp1q[0], 4);
        //printf("==================\n");
        //}
        G->set_gm(noise_sp1q[0],4);
#endif
        circuit_handle->append(*G);
    }
    void ID(IdxType qubit)
    {
        Gate* G = new Gate(OP::C2, qubit, -1, 0);
#ifdef DISABLE_NOISE
        std::complex<double> perfect_sp1q[4][4] = {};
        backendGateNoiseSp("id", perfect_sp1q[0], qubit, -1);
        G->set_gm(perfect_sp1q[0],4);
#else
        std::complex<double> noise_sp1q[4][4] = {};
        backendGateNoiseSp(device_config, "id", noise_sp1q[0], qubit, -1);
        G->set_gm(noise_sp1q[0],4);
#endif
        circuit_handle->append(*G);
    }
    void RZ(ValType theta, IdxType qubit)
    {
        Gate* G = new Gate(OP::C2, qubit, -1, theta);
        //if (i_cpu == 0) printf("rz(%lf) q[%lld];\n",qubit, theta);
#ifdef DISABLE_NOISE
        std::complex<double> perfect_sp1q[4][4] = {};
        backendGateNoiseSp("rz", perfect_sp1q[0], qubit, -1, theta);
        G->set_gm(perfect_sp1q[0],4);
#else
        std::complex<double> noise_sp1q[4][4] = {};
        backendGateNoiseSp(device_config, "rz", noise_sp1q[0], qubit, -1, theta);
        G->set_gm(noise_sp1q[0],4);
#endif
        circuit_handle->append(*G);
    }
    void SX(IdxType qubit)
    {
        Gate* G = new Gate(OP::C2, qubit, -1, 0);
        //if (i_cpu == 0) printf("sx q[%lld];\n",qubit);
#ifdef DISABLE_NOISE
        std::complex<double> perfect_sp1q[4][4] = {};
        backendGateNoiseSp("sx", perfect_sp1q[0], qubit, -1);
        G->set_gm(perfect_sp1q[0],4);
#else
        std::complex<double> noise_sp1q[4][4] = {};
        backendGateNoiseSp(device_config, "sx", noise_sp1q[0], qubit, -1);
        G->set_gm(noise_sp1q[0],4);
#endif
        circuit_handle->append(*G);
    }
    void CX(IdxType ctrl, IdxType qubit)
    {
        Gate* G = new Gate(OP::C4, qubit, ctrl, 0);
        //if (i_cpu == 0) printf("cx q[%lld], q[%lld];\n", ctrl, qubit);
#ifdef DISABLE_NOISE
        std::complex<double> perfect_sp2q[16][16] = {};
        backendGateNoiseSp("cx", perfect_sp2q[0], ctrl, qubit);
        G->set_gm(perfect_sp2q[0],16);
#else
        std::complex<double> noise_sp2q[16][16] = {};
        backendGateNoiseSp(device_config, "cx", noise_sp2q[0], ctrl, qubit);
        G->set_gm(noise_sp2q[0],16);
#endif
        circuit_handle->append(*G);

    }
    void M(IdxType qubit) //default is pauli-Z
    {
        SAFE_FREE_HOST(results);
        SAFE_ALOC_HOST(results, sizeof(IdxType));
        memset(results, 0, sizeof(IdxType));
        ValType rand = uni_dist(rng);
        Gate* G = new Gate(OP::M,qubit,-1,rand);
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
            randoms[i] = uni_dist(rng);
        Gate* G = new Gate(OP::MA,0,repetition,0);
        circuit_handle->append(*G);
    }
    void RESET(IdxType qubit)
    {
        Gate* G = new Gate(OP::RESET,qubit,-1,0);
        circuit_handle->append(*G);
    }
    void reset_sim()
    {
        //Reset CPU input & output
        memset(sv_real, 0, sv_size);
        memset(sv_imag, 0, sv_size);
        memset(m_real, 0, diag_size);
        //State Vector initial state [0..0] = 1
        sv_real[0] = 1.;
        reset_circuit();
    }
    void reset_circuit()
    {
        circuit_handle->reset();
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
    void clear_circuit()
    {
        circuit_handle->clear();
    }
    void update(const IdxType _n_qubits, const IdxType _n_gates)
    {
        this->n_gates = _n_gates;
    }
    std::string circuitToString()
    {
        return circuit_handle->circuitToString();
    }
    void sim()
    {
        IdxType input_gates = circuit_handle->n_gates;
        circuit_handle_cpu = circuit_handle->upload();
        //update should be put after upload where gate fusion is applied
        //which may change the number of gates in the circuit
        update(circuit_handle->n_qubits, circuit_handle->n_gates);
#ifdef PRINT_SIM_TRACE
        printf("DMSim_CPU is running! Requesting %lld qubits.\n", circuit_handle->n_qubits);
#endif
#ifdef PRINT_KERNEL_TRACE
        double sim_time;
        cpu_timer sim_timer;
        sim_timer.start_timer();
#endif
        //=========================================
        simulation_kernel(this);
        //=========================================

#ifdef PRINT_KERNEL_TRACE
        sim_timer.stop_timer();
        sim_time = sim_timer.measure();
        printf("\n============== DM-Sim ===============\n");
#ifdef USE_AVX512
        printf("nqubits:%lld, ngates:%lld, sim_gates:%lld, ncpus(avx512):%lld, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms, mem:%.3lf MB, mem_per_cpu:%.3lf MB\n",
                n_qubits, input_gates, n_gates, 1, sim_time, 0., 
                sim_time, cpu_mem/1024/1024, cpu_mem/1024/1024);
#else
        printf("nqubits:%lld, ngates:%lld, sim_gates:%lld, ncpus:%lld, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms, mem:%.3lf MB, mem_per_cpu:%.3lf MB\n",
                n_qubits, input_gates, n_gates, 1, sim_time, 0., 
                sim_time, cpu_mem/1024/1024, cpu_mem/1024/1024);
#endif
        printf("=====================================\n");
        fflush(stdout);
#endif
        reset_circuit();
    }
    IdxType measure(IdxType qubit) 
    {
        this->M(qubit);
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
        IdxType num = ((IdxType)1<<n_qubits);
        printf("----- Real DM diag ------\n");
        for (IdxType i=0; i<num; i++) 
        {
            printf("%lf ", sv_real[i*num+i]);
            if ((i+1)%8==0) printf("\n");
        }
        printf("\n");
    }
public:
    // n_qubits is the number of qubits
    IdxType n_qubits;
    // cpu_scale is 2^x of the number of CPUs, e.g., with 8 CPUs the cpu_scale is 3 (2^3=8)
    IdxType dim;
    IdxType half_dim;
    IdxType sv_size;
    IdxType n_gates;
    // for measurement along the dm diagnol
    IdxType diag_num;
    IdxType diag_size;
    //CPU arrays
    ValType* sv_real;
    ValType* sv_imag;
    //For joint measurement
    ValType* m_real;
    //For measurement randoms
    ValType* randoms;
    //For measurement result
    IdxType* results;
    //Device name
    std::string device_name;
    //Device config dict
    json device_config;
    //Random
    std::mt19937 rng;
    std::uniform_real_distribution<ValType> uni_dist;
    //CPU memory usage
    ValType cpu_mem;
    //cricuit
    Circuit* circuit_handle;
    //circuit cpu
    Gate* circuit_handle_cpu;
};

/***********************************************
 * Simulation Main Kernel
 ***********************************************/
void simulation_kernel(Simulation* sim)
{
    for (IdxType t=0; t<(sim->n_gates); t++)
    {
        //IdxType t0;
        //if (sim->i_cpu == 0 && blockIdx.x == 0 && threadIdx.x == 0) t0 = clock64();

        ((sim->circuit_handle_cpu)[t]).exe_op(sim, sim->sv_real, sim->sv_imag);

        /*
        if (sim->i_cpu == 0 && blockIdx.x == 0 && threadIdx.x == 0)
        {
            IdxType time = clock64() - t0;
            printf("%s(%lld) ticks: %lld\n",OP_NAMES_NVGPU[sim->circuit_handle_cpu[t].op_name],sim->circuit_handle_cpu[t].qubit, time);
        }
         */
#ifdef PURITY_CHECK
        Purity_Check(sim, t, sim->sv_real, sim->sv_imag);
#endif
    }
}

//================================= Gate Definition ========================================

#ifdef USE_AVX512
#include "svsim_cpu_avx512.hpp"
#else


/***********************************************
 * Key Macros
 ***********************************************/
#define LOCAL_G(arr,i) (arr[(i)])
#define LOCAL_P(arr,i,val) (arr[(i)]=(val))
#define BARR while(0){};

//============== Local 2-qubit Gate  ================
void C2_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag, const IdxType qubit0, const IdxType qubit1)
{
    const IdxType per_pe_work = ((sim->dim)>>2);
    assert (qubit0 != qubit1); //Non-cloning
    const IdxType q0dim = ((IdxType)1 << max(qubit0, qubit1) );
    const IdxType q1dim = ((IdxType)1 << min(qubit0, qubit1) );
    const IdxType outer_factor = ((sim->dim) + q0dim + q0dim - 1) >> (max(qubit0,qubit1)+1);
    const IdxType mider_factor = (q0dim + q1dim + q1dim - 1) >> (min(qubit0,qubit1)+1);
    const IdxType inner_factor = q1dim;
    const IdxType qubit0_dim = ((IdxType)1 << qubit0);
    const IdxType qubit1_dim = ((IdxType)1 << qubit1);

    for (IdxType i=0; i<per_pe_work; i++)
    {
        IdxType outer = ((i/inner_factor) / (mider_factor)) * (q0dim+q0dim);
        IdxType mider = ((i/inner_factor) % (mider_factor)) * (q1dim+q1dim);
        IdxType inner = i % inner_factor;
        IdxType pos0 = outer + mider + inner;
        IdxType pos1 = outer + mider + inner + qubit1_dim;
        IdxType pos2 = outer + mider + inner + qubit0_dim;
        IdxType pos3 = outer + mider + inner + q0dim + q1dim;

        const ValType el0_real = LOCAL_G(sv_real, pos0);
        const ValType el0_imag = LOCAL_G(sv_imag, pos0);
        const ValType el1_real = LOCAL_G(sv_real, pos1);
        const ValType el1_imag = LOCAL_G(sv_imag, pos1);
        const ValType el2_real = LOCAL_G(sv_real, pos2);
        const ValType el2_imag = LOCAL_G(sv_imag, pos2);
        const ValType el3_real = LOCAL_G(sv_real, pos3);
        const ValType el3_imag = LOCAL_G(sv_imag, pos3);

        //Real part
        ValType sv_real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag)
                              +(gm_real[1] * el1_real) - (gm_imag[1] * el1_imag)
                              +(gm_real[2] * el2_real) - (gm_imag[2] * el2_imag)
                              +(gm_real[3] * el3_real) - (gm_imag[3] * el3_imag);
        ValType sv_real_pos1 = (gm_real[4] * el0_real) - (gm_imag[4] * el0_imag)
                              +(gm_real[5] * el1_real) - (gm_imag[5] * el1_imag)
                              +(gm_real[6] * el2_real) - (gm_imag[6] * el2_imag)
                              +(gm_real[7] * el3_real) - (gm_imag[7] * el3_imag);
        ValType sv_real_pos2 = (gm_real[8] * el0_real) - (gm_imag[8] * el0_imag)
                              +(gm_real[9] * el1_real) - (gm_imag[9] * el1_imag)
                              +(gm_real[10] * el2_real) - (gm_imag[10] * el2_imag)
                              +(gm_real[11] * el3_real) - (gm_imag[11] * el3_imag);
        ValType sv_real_pos3 = (gm_real[12] * el0_real) - (gm_imag[12] * el0_imag)
                              +(gm_real[13] * el1_real) - (gm_imag[13] * el1_imag)
                              +(gm_real[14] * el2_real) - (gm_imag[14] * el2_imag)
                              +(gm_real[15] * el3_real) - (gm_imag[15] * el3_imag);

        //Imag part
        ValType sv_imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real)
                              +(gm_real[1] * el1_imag) + (gm_imag[1] * el1_real)
                              +(gm_real[2] * el2_imag) + (gm_imag[2] * el2_real)
                              +(gm_real[3] * el3_imag) + (gm_imag[3] * el3_real);
        ValType sv_imag_pos1 = (gm_real[4] * el0_imag) + (gm_imag[4] * el0_real)
                              +(gm_real[5] * el1_imag) + (gm_imag[5] * el1_real)
                              +(gm_real[6] * el2_imag) + (gm_imag[6] * el2_real)
                              +(gm_real[7] * el3_imag) + (gm_imag[7] * el3_real);
        ValType sv_imag_pos2 = (gm_real[8] * el0_imag) + (gm_imag[8] * el0_real)
                              +(gm_real[9] * el1_imag) + (gm_imag[9] * el1_real)
                              +(gm_real[10] * el2_imag) + (gm_imag[10] * el2_real)
                              +(gm_real[11] * el3_imag) + (gm_imag[11] * el3_real);
        ValType sv_imag_pos3 = (gm_real[12] * el0_imag) + (gm_imag[12] * el0_real)
                              +(gm_real[13] * el1_imag) + (gm_imag[13] * el1_real)
                              +(gm_real[14] * el2_imag) + (gm_imag[14] * el2_real)
                              +(gm_real[15] * el3_imag) + (gm_imag[15] * el3_real);

        LOCAL_P(sv_real, pos0, sv_real_pos0); 
        LOCAL_P(sv_real, pos1, sv_real_pos1); 
        LOCAL_P(sv_real, pos2, sv_real_pos2); 
        LOCAL_P(sv_real, pos3, sv_real_pos3); 

        LOCAL_P(sv_imag, pos0, sv_imag_pos0); 
        LOCAL_P(sv_imag, pos1, sv_imag_pos1); 
        LOCAL_P(sv_imag, pos2, sv_imag_pos2); 
        LOCAL_P(sv_imag, pos3, sv_imag_pos3); 
    }
    //BARR;
}


#define SV4IDX(x) (((x>>1)&1)*EXP2E(qubit0) + ((x&1)*EXP2E(qubit1)) )

#define DIV2E(x,y) ((x)>>(y))
#define MOD2E(x,y) ((x)&(((IdxType)1<<(y))-(IdxType)1)) 
#define EXP2E(x) ((IdxType)1<<(x))
#define SV16IDX(x) ( ((x>>3)&1)*EXP2E(qubit0) + ((x>>2)&1)*EXP2E(qubit1) + ((x>>1)&1)*EXP2E(qubit2) + ((x&1)*EXP2E(qubit3)) )

//============== Local 4-qubit Gate ================
void C4_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag, const IdxType qubit0, const IdxType qubit1,
        const IdxType qubit2, const IdxType qubit3)
{
    const IdxType per_pe_work = ((sim->dim)>>4);
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

    for (IdxType i=0; i<per_pe_work; i++)
    {
        const IdxType term0 = MOD2E(i,p);
        const IdxType term1 = MOD2E(DIV2E(i,p),q-p-1)*EXP2E(p+1);
        const IdxType term2 = MOD2E(DIV2E(DIV2E(i,p),q-p-1),r-q-1)*EXP2E(q+1);
        const IdxType term3 = MOD2E(DIV2E(DIV2E(DIV2E(i,p),q-p-1),r-q-1),s-r-1)*EXP2E(r+1);
        const IdxType term4 = DIV2E(DIV2E(DIV2E(DIV2E(i,p),q-p-1),r-q-1),s-r-1)*EXP2E(s+1);
        const IdxType term = term4 + term3 + term2 + term1 + term0;
        const ValType el_real[16] = { 
            LOCAL_G(sv_real,term+SV16IDX(0)),  LOCAL_G(sv_real,term+SV16IDX(1)),
            LOCAL_G(sv_real,term+SV16IDX(2)),  LOCAL_G(sv_real,term+SV16IDX(3)),
            LOCAL_G(sv_real,term+SV16IDX(4)),  LOCAL_G(sv_real,term+SV16IDX(5)),
            LOCAL_G(sv_real,term+SV16IDX(6)),  LOCAL_G(sv_real,term+SV16IDX(7)),
            LOCAL_G(sv_real,term+SV16IDX(8)),  LOCAL_G(sv_real,term+SV16IDX(9)),
            LOCAL_G(sv_real,term+SV16IDX(10)), LOCAL_G(sv_real,term+SV16IDX(11)),
            LOCAL_G(sv_real,term+SV16IDX(12)), LOCAL_G(sv_real,term+SV16IDX(13)),
            LOCAL_G(sv_real,term+SV16IDX(14)), LOCAL_G(sv_real,term+SV16IDX(15))
        };
        const ValType el_imag[16] = { 
            LOCAL_G(sv_imag,term+SV16IDX(0)),  LOCAL_G(sv_imag,term+SV16IDX(1)),
            LOCAL_G(sv_imag,term+SV16IDX(2)),  LOCAL_G(sv_imag,term+SV16IDX(3)),
            LOCAL_G(sv_imag,term+SV16IDX(4)),  LOCAL_G(sv_imag,term+SV16IDX(5)),
            LOCAL_G(sv_imag,term+SV16IDX(6)),  LOCAL_G(sv_imag,term+SV16IDX(7)),
            LOCAL_G(sv_imag,term+SV16IDX(8)),  LOCAL_G(sv_imag,term+SV16IDX(9)),
            LOCAL_G(sv_imag,term+SV16IDX(10)), LOCAL_G(sv_imag,term+SV16IDX(11)),
            LOCAL_G(sv_imag,term+SV16IDX(12)), LOCAL_G(sv_imag,term+SV16IDX(13)),
            LOCAL_G(sv_imag,term+SV16IDX(14)), LOCAL_G(sv_imag,term+SV16IDX(15))
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
            LOCAL_P(sv_real, term+SV16IDX(j), res_real);
            LOCAL_P(sv_imag, term+SV16IDX(j), res_imag);
        }
    }
    //BARR;
}
#endif

void M_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, ValType* gm_real, 
        ValType* gm_imag, const IdxType qubit, const ValType rand)
{
    IdxType mask = ((IdxType)1<<qubit);
    ValType prob_of_one = 0;
    for (IdxType i=0; i<((IdxType)1<<(sim->n_qubits)); i++)
    {
        if ((i&mask)!=0) prob_of_one += fabs(sv_real[(i<<(sim->n_qubits))+i]);
    }
    if (rand <= prob_of_one) gm_real[15] = 1.0/prob_of_one;
    else gm_real[0] = 1.0/(1.0-prob_of_one);
    BARR;
    C2_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, qubit, qubit+sim->n_qubits);
    BARR;
    sim->results[0] = (rand<=prob_of_one?1:0);
}

void MA_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag,
        const IdxType repetition)
{
    ValType* m_real = sim->m_real;
    const IdxType n_size = (IdxType)1<<(sim->n_qubits);
    m_real[0] = 0;
    for (IdxType i=1; i<=(n_size); i++)
        m_real[i] = m_real[i-1] + fabs(sv_real[((i-1)<<(sim->n_qubits))+i-1]);
    ValType purity = fabs(m_real[n_size]);

    if ( fabs(purity - 1.0) > ERROR_BAR )
        printf("MA: Purity Check fails with %lf\n", purity);
    if (repetition < n_size)
    {
        for (IdxType j=0; j<n_size; j++)
        {
            ValType lower = m_real[j];
            ValType upper = (j+1==n_size)? 1:m_real[j+1];
            for (IdxType i=0; i<repetition; i++)
            {
                ValType r = sim->randoms[i];
                if (lower<=r && r<upper) sim->results[i] = j; 
            }
        }
    }
    else
    {
        for (IdxType i=0; i<repetition; i++)
        {
            IdxType lo = 0;
            IdxType hi = ((IdxType)1<<(sim->n_qubits));
            IdxType mid;
            ValType r = sim->randoms[i];
            while (hi-lo>1)
            {
                mid = lo + (hi-lo)/2;
                if (r >= m_real[mid])
                    lo = mid;
                else
                    hi = mid;
            }
            sim->results[i] = lo;
        }
    }
}

void RESET_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag,
        const IdxType qubit)
{
    IdxType mask = ((IdxType)1<<qubit);
    ValType prob_of_one = 0;
    for (IdxType i=0; i<((IdxType)1<<(sim->n_qubits)); i++)
    {
        if ((i&mask)!=0) prob_of_one += fabs(sv_real[(i<<(sim->n_qubits))+i]);
    }
    assert(prob_of_one <= 1.0);
    
    if (prob_of_one < 1.0) //still possible to normalize
    {
        ValType factor = 1.0/(1.0-prob_of_one);
        for (IdxType i=0; i<sim->dim; i++)
        {
            if ((i&mask) == 0)
            {
                sv_real[i] *= factor;
                sv_imag[i] *= factor;
            }
            else
            {
                sv_real[i] = 0;
                sv_imag[i] = 0;
            }
        }
    }
    else
    {
        for (IdxType i=0; i<sim->dim; i++)
        {
            if ( (i & mask) == 0)
            {
                IdxType dual_i = i^mask;
                sv_real[i] = sv_real[dual_i];
                sv_imag[i] = sv_imag[dual_i];
                sv_real[dual_i] = 0;
                sv_imag[dual_i] = 0;
            }
        }
    }
}

//============== Purity Check  ================
void Purity_Check(const Simulation* sim, const IdxType t, ValType* sv_real, ValType* sv_imag)
{
    ValType purity = 0; 
    for (IdxType i=0; i<(((IdxType)1<<(sim->n_qubits))); i++)
        purity += fabs(sv_real[(i<<(sim->n_qubits))+i]);
    if ( fabs(purity - 1.0) > ERROR_BAR )
    {
        Gate* g = &sim->circuit_handle_cpu[t];
        printf("Purity Check fails after Gate-%lld=>%s(ctrl:%lld,qubit:%lld,theta:%lf) with %lf\n",t,OP_NAMES[g->op_name],g->ctrl,g->qubit,g->theta,purity);
    }
}




//=====================================================================================
//Per-gate execution function
void Gate::exe_op(Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    if (op_name == RESET)
    {
        RESET_GATE(sim, sv_real, sv_imag, qubit);
    }
    else if (op_name == M)
    {
        M_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, qubit, theta);
    }
    else if (op_name == MA)
    {
        MA_GATE(sim, sv_real, sv_imag, ctrl); 
    }
    else if (op_name == C2) 
    {
        C2_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, qubit, qubit+(sim->n_qubits));
    }
    else if (op_name == C4) 
    {
        C4_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, ctrl, qubit, ctrl+(sim->n_qubits), qubit+(sim->n_qubits));
    }
}

}; //namespace NWQSim

#endif
