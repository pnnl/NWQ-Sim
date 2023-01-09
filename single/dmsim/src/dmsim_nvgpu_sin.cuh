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
// File: svsim_nvgpu_sin.cuh
// MPI and NVSHMEM based implementation of the scale-out SV-Sim gates and 
// simulation runtime using NVIDIA GPU backend.
// ---------------------------------------------------------------------------
#ifndef DMSIM_NVGPU_SIN_CUH
#define DMSIM_NVGPU_SIN_CUH
#include <assert.h>
#include <random>
#include <complex.h>
#include <cooperative_groups.h>
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
using namespace cooperative_groups;
using namespace std;
using namespace nvcuda;
    
class Gate;
class Simulation;
__device__ void Purity_Check(const Simulation* sim, const IdxType t, ValType* sv_real, ValType* sv_imag);


//Declare major simulation kernel
__global__ void simulation_kernel(Simulation*);

/***********************************************
 * Circuit Definition
 ***********************************************/
class Circuit
{
public:
    Circuit(IdxType _n_qubits=0):
        n_qubits(_n_qubits), n_gates(0), circuit_gpu(NULL)
    {}
    ~Circuit() { clear(); }
    void append(Gate& g)
    {
#ifdef PRINT_GATE_TRACE
        if (nvshmem_my_pe() == 0) std::cout << g.gateToString() << std::flush;
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
        SAFE_FREE_GPU(circuit_gpu);
    }
    void reset()
    {
        clear();
    }
    Gate* upload()
    {
#ifdef DISABLE_GATE_FUSION
        //====================== No Fuse =====================
        SAFE_FREE_GPU(circuit_gpu);
        SAFE_ALOC_GPU(circuit_gpu, n_gates*sizeof(Gate));
        cudaSafeCall(cudaMemcpy(circuit_gpu, circuit.data(), n_gates*sizeof(Gate), cudaMemcpyHostToDevice));
        //====================================================
#else
        //====================== Fuse ========================
        vector<Gate> tmp_circuit;
        tmp_circuit.clear();
        gate_fusion_1q(circuit, tmp_circuit, n_qubits);
        gate_fusion_2q(tmp_circuit, fused_circuit, n_qubits);

        this->n_gates = fused_circuit.size();
        SAFE_FREE_GPU(circuit_gpu);
        SAFE_ALOC_GPU(circuit_gpu, n_gates*sizeof(Gate));
        cudaSafeCall(cudaMemcpy(circuit_gpu, fused_circuit.data(), n_gates*sizeof(Gate), cudaMemcpyHostToDevice));
        //====================================================
#endif
        return circuit_gpu;
    }
    std::string circuitToString()
    {
        stringstream ss;
        for (IdxType t=0; t<n_gates; t++)
            ss << circuit[t].gateToString();
        return ss.str();
    }
public:
    // number of qubits
    IdxType n_qubits;
    // number of gates
    IdxType n_gates;
    // user input gate sequence
    vector<Gate> circuit;
    // fused gate sequence
    vector<Gate> fused_circuit;
    // gate sequence executed
    Gate* circuit_gpu;
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
        n_gates(0), 
        gpu_mem(0),
        sim_gpu(NULL),
        sv_real(NULL),
        sv_imag(NULL),
        m_real(NULL),
        m_imag(NULL),
        randoms(NULL),
        randoms_gpu(NULL),
        results(NULL),
        results_gpu(NULL)
    {
        i_gpu = DEFAULT_SIN_GPU;
        //load device configuration file
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
        
        //always be 0 since 1-MPI maps to 1-GPU
        cudaSafeCall(cudaSetDevice(i_gpu));
        //CPU side initialization
        SAFE_ALOC_HOST(sv_real_cpu, sv_size);
        SAFE_ALOC_HOST(sv_imag_cpu, sv_size);
        memset(sv_real_cpu, 0, sv_size);
        memset(sv_imag_cpu, 0, sv_size);
        //State-vector initial state [0..0] = 1
        if (i_gpu == 0) sv_real_cpu[0] = 1.;
        //NVSHMEM GPU memory allocation
        SAFE_ALOC_GPU(sv_real, sv_size);
        SAFE_ALOC_GPU(sv_imag, sv_size);
        SAFE_ALOC_GPU(m_real, sv_size);
        SAFE_ALOC_GPU(m_imag, sv_size);

        cudaCheckError(); 
        gpu_mem += sv_size*4;
        //Initialize Circuit 
        circuit_handle = new Circuit(n_qubits);
        circuit_handle_gpu = NULL;
        //GPU memory initilization
        cudaSafeCall(cudaMemcpy(sv_real, sv_real_cpu, sv_size, 
                    cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(sv_imag, sv_imag_cpu, sv_size, 
                    cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemset(m_real, 0, sv_size));
        cudaSafeCall(cudaMemset(m_imag, 0, sv_size));
        SAFE_ALOC_GPU(sim_gpu, sizeof(Simulation));
        rng.seed(time(0));
    }
    ~Simulation()
    {
        //Release circuit
        if (circuit_handle != NULL)
            delete circuit_handle;
        //Release for CPU side
        SAFE_FREE_HOST(sv_real_cpu);
        SAFE_FREE_HOST(sv_imag_cpu);
        SAFE_FREE_HOST(randoms);
        SAFE_FREE_HOST(results);
        //Release for GPU side
        SAFE_FREE_GPU(sv_real);
        SAFE_FREE_GPU(sv_imag);
        SAFE_FREE_GPU(m_real);
        SAFE_FREE_GPU(m_imag);
        SAFE_FREE_GPU(sim_gpu);
        SAFE_FREE_GPU(randoms_gpu);
        SAFE_FREE_GPU(results_gpu);
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
        //if (i_gpu == 0) printf("x q[%lld];\n",qubit);

#ifdef DISABLE_NOISE
        std::complex<double> perfect_sp1q[4][4] = {};
        backendGateNoiseSp("x", perfect_sp1q[0], qubit, -1);
        G->set_gm(perfect_sp1q[0],4);
#else
        std::complex<double> noise_sp1q[4][4] = {};
        backendGateNoiseSp(device_config, "x", noise_sp1q[0], qubit, -1);
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
        //if (i_gpu == 0) printf("rz(%lf) q[%lld];\n",qubit, theta);
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
        //if (i_gpu == 0) printf("sx q[%lld];\n",qubit);
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
        //if (i_gpu == 0) printf("cx q[%lld], q[%lld];\n", ctrl, qubit);
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
        SAFE_FREE_GPU(results_gpu);
        SAFE_ALOC_GPU(results_gpu, sizeof(IdxType));
        cudaSafeCall(cudaMemset(results_gpu, 0, sizeof(IdxType)));
        ValType rand = uni_dist(rng);
        Gate* G = new Gate(OP::M,qubit,-1,rand);
        circuit_handle->append(*G);
    }
    void MA(IdxType repetition) //default is pauli-Z
    {
        SAFE_FREE_HOST(results);
        SAFE_ALOC_HOST(results, sizeof(IdxType)*repetition);
        memset(results, 0, sizeof(IdxType)*repetition);
        SAFE_FREE_GPU(results_gpu);
        SAFE_ALOC_GPU(results_gpu, sizeof(IdxType)*repetition);
        cudaSafeCall(cudaMemset(results_gpu, 0, sizeof(IdxType)*repetition));
        SAFE_FREE_HOST(randoms);
        SAFE_ALOC_HOST(randoms, sizeof(ValType)*repetition);
        for (IdxType i=0; i<repetition; i++) 
            randoms[i] = uni_dist(rng);
        SAFE_FREE_GPU(randoms_gpu);
        SAFE_ALOC_GPU(randoms_gpu, sizeof(ValType)*repetition);
        cudaSafeCall(cudaMemcpy(randoms_gpu, randoms, 
                    sizeof(ValType)*repetition, cudaMemcpyHostToDevice));
        Gate* G = new Gate(OP::MA,0,repetition,0);
        circuit_handle->append(*G);
    }
    void RESET(IdxType qubit)
    {
        Gate* G = new Gate(OP::RESET,qubit);
        circuit_handle->append(*G);
    }
    void reset_sim()
    {
        //Reset CPU input & output
        memset(sv_real_cpu, 0, sv_size);
        memset(sv_imag_cpu, 0, sv_size);
        //State Vector initial state [0..0] = 1
        if (i_gpu == 0) sv_real_cpu[0] = 1.;
        //GPU side initialization
        cudaSafeCall(cudaMemcpy(sv_real, sv_real_cpu, 
                    sv_size, cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(sv_imag, sv_imag_cpu, 
                    sv_size, cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemset(m_real, 0, sv_size));
        cudaSafeCall(cudaMemset(m_imag, 0, sv_size));
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
        
        //assert(_n_qubits <= (N_QUBIT_SLOT/2));
        //this->n_qubits = _n_qubits;
        //this->dim = ((IdxType)1<<(2*n_qubits));
        //this->half_dim = (IdxType)1<<(2*n_qubits-1);
        //this->sv_size = dim*(IdxType)sizeof(ValType);
        //this->lg2_m_gpu = 2*n_qubits-gpu_scale;
        //this->m_gpu = ((IdxType)1<<(lg2_m_gpu));
        //this->sv_size = sv_size/n_gpus;
    }
    std::string circuitToString()
    {
        return circuit_handle->circuitToString();
    }
    void sim()
    {
        cudaSafeCall(cudaSetDevice(i_gpu));
        IdxType input_gates = circuit_handle->n_gates;
        circuit_handle_gpu = circuit_handle->upload();
        //update should be put after upload where gate fusion is applied
        //which may change the number of gates in the circuit
        update(circuit_handle->n_qubits, circuit_handle->n_gates);
        cudaSafeCall(cudaMemcpy(sim_gpu, this, 
                    sizeof(Simulation), cudaMemcpyHostToDevice));

#ifdef PRINT_SIM_TRACE
        printf("DMSim_gpu is running! Requesting %lld qubits.\n", circuit_handle->n_qubits);
#endif

#ifdef PRINT_KERNEL_TRACE
        double sim_time;
        gpu_timer sim_timer;
#endif
        dim3 gridDim(1,1,1);
        cudaDeviceProp deviceProp;
        cudaSafeCall(cudaGetDeviceProperties(&deviceProp, 0));
        //8*16 is per warp shared-memory usage, with real and imag
        unsigned smem_size = THREADS_CTA_NVGPU/32*8*16*2*sizeof(ValType);
        int numBlocksPerSm;
        cudaSafeCall(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm, 
                    simulation_kernel, THREADS_CTA_NVGPU, smem_size));
        gridDim.x = numBlocksPerSm * deviceProp.multiProcessorCount;
        void* args[] = {&sim_gpu};
        cudaSafeCall(cudaDeviceSynchronize());
#ifdef PRINT_KERNEL_TRACE
        sim_timer.start_timer();
#endif
        cudaLaunchCooperativeKernel((void*)simulation_kernel,gridDim,
                THREADS_CTA_NVGPU,args,smem_size);
        cudaSafeCall(cudaDeviceSynchronize());
#ifdef PRINT_KERNEL_TRACE
        sim_timer.stop_timer();
        sim_time = sim_timer.measure();
#endif
        cudaCheckError();
#ifdef PRINT_KERNEL_TRACE
        printf("\n============== DM-Sim ===============\n");
        printf("nqubits:%lld, ngates:%lld, sim_gates:%lld, ngpus:%lld, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms, mem:%.3lf MB, mem_per_gpu:%.3lf MB\n",
                n_qubits, input_gates, n_gates, 1, sim_time, 0., 
                sim_time, gpu_mem/1024/1024, gpu_mem/1024/1024);
        printf("=====================================\n");
#endif
        clear_circuit();
    }
    IdxType measure(IdxType qubit) 
    {
        this->M(qubit);
        this->sim();
        cudaSafeCall(cudaMemcpy(results, results_gpu, sizeof(IdxType), cudaMemcpyDeviceToHost));

        return this->results[0];
    }
    IdxType* measure_all(IdxType repetition=DEFAULT_REPETITIONS)
    {
        this->MA(repetition);
        this->sim();
        cudaSafeCall(cudaMemcpy(results, results_gpu, sizeof(IdxType)*repetition, cudaMemcpyDeviceToHost));
        return this->results;
    }
    void print_res_sv()
    {
        cudaCheckError();
        cudaSafeCall(cudaMemcpy(sv_real_cpu, sv_real, sv_size, cudaMemcpyDeviceToHost));
        cudaSafeCall(cudaMemcpy(sv_imag_cpu, sv_imag, sv_size, cudaMemcpyDeviceToHost));

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
    // which gpu
    IdxType i_gpu;
    // gpu_scale is 2^x of the number of GPUs, e.g., with 8 GPUs the gpu_scale is 3 (2^3=8)
    IdxType dim;
    IdxType half_dim;
    IdxType sv_size;
    IdxType n_gates;
    //CPU arrays
    ValType* sv_real_cpu;
    ValType* sv_imag_cpu;
    //GPU arrays
    ValType* sv_real;
    ValType* sv_imag;
    //For joint measurement
    ValType* m_real;
    ValType* m_imag;
    //For measurement randoms
    ValType* randoms;
    ValType* randoms_gpu;
    //For measurement result
    IdxType* results;
    IdxType* results_gpu;
    //Device name
    std::string device_name;
    //Device config dict
    json device_config;
    //Random
    std::mt19937 rng;
    std::uniform_real_distribution<ValType> uni_dist;
    //GPU memory usage
    ValType gpu_mem;
    //cricuit
    Circuit* circuit_handle;
    //circuit gpu
    Gate* circuit_handle_gpu;
    //hold the GPU-side simulator instances
    Simulation* sim_gpu;
};

#define LOCAL_G(arr,i) arr[(i)]
#define LOCAL_P(arr,i,val) arr[(i)] = val;
#define BARR grid.sync();

//============== Check Trace (debug purpose) ================
__device__ __inline__ void CHECK_TRACE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, IdxType t)
{
    grid_group grid = this_grid();
    if (blockIdx.x ==0 && threadIdx.x == 0)
    {
        ValType trace = 0;
        for (IdxType i=0; i<((IdxType)1<<(sim->n_qubits)); i++)
        {
            const ValType val = LOCAL_G(sv_real, (i<<(sim->n_qubits))+i);
            trace += abs(val);
        }
        printf("%s: Trace is: %lf\n", OP_NAMES_NVGPU[sim->circuit_handle_gpu[t].op_name], trace);
    }
    BARR;
}

__global__ void simulation_kernel(Simulation* sim)
{
    for (IdxType t=0; t<(sim->n_gates); t++)
    {
        //CHECK_TRACE(sim, sim->sv_real, sim->sv_imag, t);
        //if (blockIdx.x == 0 && threadIdx.x == 0) 
        //printf("==%lld== %s(qubit:%lld, ctrl:%lld, theta:%lf)\n",sim->circuit_handle_gpu[t].op_name, OP_NAMES_NVGPU[sim->circuit_handle_gpu[t].op_name], sim->circuit_handle_gpu[t].qubit, sim->circuit_handle_gpu[t].ctrl, sim->circuit_handle_gpu[t].theta);

        //IdxType t0;
        //if (sim->i_gpu == 0 && blockIdx.x == 0 && threadIdx.x == 0) t0 = clock64();

        ((sim->circuit_handle_gpu)[t]).exe_op(sim, sim->sv_real, sim->sv_imag);

        /*
        if (blockIdx.x == 0 && threadIdx.x == 0)
        {
            IdxType time = clock64() - t0;
            printf("%s(%lld) ticks: %lld\n",OP_NAMES_NVGPU[sim->circuit_handle_gpu[t].op_name],sim->circuit_handle_gpu[t].qubit, time);
        }
         */
        //CHECK_TRACE(sim, sim->sv_real, sim->sv_imag, t);
#ifdef PURITY_CHECK
        Purity_Check(sim, t, sim->sv_real, sim->sv_imag);
#endif
    }
}

//================================= Gate Definition ========================================
//============== Unified 2-qubit Gate without comm optimization ================
__device__ __inline__ void C2_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag, const IdxType qubit0, const IdxType qubit1)
{
    grid_group grid = this_grid(); 
    const int tid = blockDim.x * blockIdx.x + threadIdx.x; 
    const IdxType per_pe_work = ((sim->dim)>>2);
    assert (qubit0 != qubit1); //Non-cloning
    
    const IdxType q0dim = ((IdxType)1 << max(qubit0, qubit1) );
    const IdxType q1dim = ((IdxType)1 << min(qubit0, qubit1) );
    const IdxType outer_factor = ((sim->dim) + q0dim + q0dim - 1) >> (max(qubit0,qubit1)+1);
    const IdxType mider_factor = (q0dim + q1dim + q1dim - 1) >> (min(qubit0,qubit1)+1);
    const IdxType inner_factor = q1dim;
    const IdxType qubit0_dim = ((IdxType)1 << qubit0);
    const IdxType qubit1_dim = ((IdxType)1 << qubit1);

    for (IdxType i=tid; i<per_pe_work; i+=blockDim.x*gridDim.x) 
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
}


#define SV4IDX(x) (((x>>1)&1)*EXP2E(qubit0) + ((x&1)*EXP2E(qubit1)) )
#define DIV2E(x,y) ((x)>>(y))
#define MOD2E(x,y) ((x)&(((IdxType)1<<(y))-(IdxType)1)) 
#define EXP2E(x) ((IdxType)1<<(x))
#define SV16IDX(x) ( ((x>>3)&1)*EXP2E(qubit0) + ((x>>2)&1)*EXP2E(qubit1) + ((x>>1)&1)*EXP2E(qubit2) + ((x&1)*EXP2E(qubit3)) )

//============== Unified 4-qubit Gate ================
__device__ __inline__ void C4_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag, const IdxType qubit0, const IdxType qubit1,
        const IdxType qubit2, const IdxType qubit3)
{
    grid_group grid = this_grid(); 
    const int tid = blockDim.x * blockIdx.x + threadIdx.x; 
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

    for (IdxType i=(sim->i_gpu)*per_pe_work+tid; i<(sim->i_gpu+1)*per_pe_work;
            i+=blockDim.x*gridDim.x) 
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

#ifdef ENABLE_TENSOR_CORES

__device__ __inline__ IdxType get_term(IdxType idx, IdxType p, IdxType q, IdxType r, IdxType s)
{
    const IdxType term0 = MOD2E(idx,p);
    const IdxType term1 = MOD2E(DIV2E(idx,p),q-p-1)*EXP2E(p+1);
    const IdxType term2 = MOD2E(DIV2E(DIV2E(idx,p),q-p-1),r-q-1)*EXP2E(q+1);
    const IdxType term3 = MOD2E(DIV2E(DIV2E(DIV2E(idx,p),q-p-1),r-q-1),s-r-1)*EXP2E(r+1);
    const IdxType term4 = DIV2E(DIV2E(DIV2E(DIV2E(idx,p),q-p-1),r-q-1),s-r-1)*EXP2E(s+1);
    const IdxType term = term4 + term3 + term2 + term1 + term0;
    return term;
}

//============== Unified 4-qubit Gate with TC V1 ================
//This is with tensorcore optimization
__device__ __inline__ void C4TCV1_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag, const IdxType qubit0, const IdxType qubit1,
        const IdxType qubit2, const IdxType qubit3)
{
    grid_group grid = this_grid(); 
    const int tid = blockDim.x * blockIdx.x + threadIdx.x; 
    const int laneid = (threadIdx.x & 31);
    const int warpid = (tid >> 5);
    const int wid = (threadIdx.x >> 5);
    const int nwarps = ((blockDim.x * gridDim.x)>>5);
    const IdxType per_pe_work = ((sim->dim)>>4);
    //Our design require at least (gpu_scale + 4 + log(8)) qubits, because
    //tensorcore demands [8,4]*[4,8], we need at least 8 on B's row dim.
    //This implies for 2 GPUs, we need at least 8 qubit for SV, thus 4 qubit for DM.
    //The smaller condition is difficult to handle due to the complex get_term function,
    //We need to figure out, after the get_term, which row/col it maps to and filter out.
    //But the row and col will be quite difficult to obtain.

    extern __shared__ ValType els[];
    ValType* el_real_s = &els[wid*8*16*2];//per warp 8*16 for real 
    ValType* el_imag_s = &els[wid*8*16*2+8*16]; //per warp 8*16 for imag

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

    wmma::fragment<wmma::matrix_a, 8, 8, 4, ValType, wmma::row_major> a_frag_real_up;
    wmma::fragment<wmma::matrix_a, 8, 8, 4, ValType, wmma::row_major> a_frag_imag_up;
    wmma::fragment<wmma::matrix_a, 8, 8, 4, ValType, wmma::row_major> a_frag_real_dn;
    wmma::fragment<wmma::matrix_a, 8, 8, 4, ValType, wmma::row_major> a_frag_imag_dn;
    wmma::fragment<wmma::matrix_b, 8, 8, 4, ValType, wmma::col_major> b_frag_real;
    wmma::fragment<wmma::matrix_b, 8, 8, 4, ValType, wmma::col_major> b_frag_imag;
    wmma::fragment<wmma::accumulator, 8, 8, 4, ValType> c_frag_real_up;
    wmma::fragment<wmma::accumulator, 8, 8, 4, ValType> c_frag_imag_up;;
    wmma::fragment<wmma::accumulator, 8, 8, 4, ValType> c_frag_real_dn;
    wmma::fragment<wmma::accumulator, 8, 8, 4, ValType> c_frag_imag_dn;;

    //if (threadIdx.x == 0 && blockIdx.x == 0) 
    //printf("GPU-%lld: dim:%lld, per_pe_work:%lld, nwarps:%d, start:%lld, end:%lld, step:%lld\n",sim->i_gpu, sim->dim,per_pe_work, nwarps,
    //(sim->i_gpu)*per_pe_work+(warpid<<3), (sim->i_gpu+1)*per_pe_work, (nwarps<<3));

    for (IdxType i=(warpid<<3); i<per_pe_work; i+=(nwarps<<3)) 
    {
        //if (laneid == 0) printf("=== warpid: %d ==== \n",warpid);
#pragma unroll
        for (int j=0; j<8; j++)
        {
            const IdxType term = get_term(i+j,p,q,r,s);
            if (laneid < 16) 
            {
                IdxType addr = term + SV16IDX(laneid);
                el_real_s[j*16+laneid] = LOCAL_G(sv_real, addr);
            }
            else 
            {
                IdxType addr = term + SV16IDX(laneid-16);
                el_imag_s[j*16+laneid-16] = LOCAL_G(sv_imag, addr);
            }
        }
        __syncwarp();
        wmma::fill_fragment(c_frag_real_up, 0.0);
        wmma::fill_fragment(c_frag_imag_up, 0.0);
        wmma::fill_fragment(c_frag_real_dn, 0.0);
        wmma::fill_fragment(c_frag_imag_dn, 0.0);
       #pragma unroll
        for (unsigned c = 0; c < 4; c++)
        {
            //load A from const mem
            load_matrix_sync(a_frag_real_up, &gm_real[c*4], 16);
            load_matrix_sync(a_frag_imag_up, &gm_imag[c*4], 16);
            load_matrix_sync(a_frag_real_dn, &gm_real[16*8+c*4], 16);
            load_matrix_sync(a_frag_imag_dn, &gm_imag[16*8+c*4], 16);

            //load B from shared mem
            load_matrix_sync(b_frag_real, &el_real_s[c*4], 16);
            load_matrix_sync(b_frag_imag, &el_imag_s[c*4], 16);
            //complex multiplication
            mma_sync(c_frag_imag_up,  a_frag_real_up, b_frag_imag, c_frag_imag_up);
            mma_sync(c_frag_imag_up,  a_frag_imag_up, b_frag_real, c_frag_imag_up);
            mma_sync(c_frag_real_up,  a_frag_real_up, b_frag_real, c_frag_real_up);
            a_frag_imag_up.x[0] = -a_frag_imag_up.x[0];
            mma_sync(c_frag_real_up, a_frag_imag_up, b_frag_imag, c_frag_real_up);

            mma_sync(c_frag_imag_dn,  a_frag_real_dn, b_frag_imag, c_frag_imag_dn);
            mma_sync(c_frag_imag_dn,  a_frag_imag_dn, b_frag_real, c_frag_imag_dn);
            mma_sync(c_frag_real_dn,  a_frag_real_dn, b_frag_real, c_frag_real_dn);
            a_frag_imag_dn.x[0] = -a_frag_imag_dn.x[0];
            mma_sync(c_frag_real_dn, a_frag_imag_dn, b_frag_imag, c_frag_real_dn);
        }
        //Store first result per segment-C
        IdxType j0 = ((2*laneid)&7);
        IdxType k0 = ((2*laneid)>>3);
        const IdxType term0 = get_term(i+j0,p,q,r,s) + SV16IDX(k0);
        LOCAL_P(sv_real, term0, c_frag_real_up.x[0]);
        LOCAL_P(sv_imag, term0, c_frag_imag_up.x[0]);

        //Store second result per segment-C
        IdxType j1 = ((2*laneid+1)&7);
        IdxType k1 = ((2*laneid+1)>>3);
        const IdxType term1 = get_term(i+j1,p,q,r,s) + SV16IDX(k1);
        LOCAL_P(sv_real, term1, c_frag_real_up.x[1]);
        LOCAL_P(sv_imag, term1, c_frag_imag_up.x[1]);

        //Store first result per segment-C
        IdxType j2 = ((2*laneid)&7);
        IdxType k2 = ((2*laneid)>>3)+8;
        const IdxType term2 = get_term(i+j2,p,q,r,s) + SV16IDX(k2);
        LOCAL_P(sv_real, term2, c_frag_real_dn.x[0]);
        LOCAL_P(sv_imag, term2, c_frag_imag_dn.x[0]);

        //Store second result per segment-C
        IdxType j3 = ((2*laneid+1)&7);
        IdxType k3 = ((2*laneid+1)>>3)+8;
        const IdxType term3 = get_term(i+j3,p,q,r,s) + SV16IDX(k3);
        LOCAL_P(sv_real, term3, c_frag_real_dn.x[1]);
        LOCAL_P(sv_imag, term3, c_frag_imag_dn.x[1]);
    }
    //BARR;
}

//============== Unified 4-qubit Gate with TC V2 ================
//This is with further tensorcore optimization
__device__ __inline__ void C4TCV2_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag, const IdxType qubit0, const IdxType qubit1,
        const IdxType qubit2, const IdxType qubit3)
{
    grid_group grid = this_grid(); 
    const int tid = blockDim.x * blockIdx.x + threadIdx.x; 
    const int laneid = (threadIdx.x & 31);
    const int warpid = (tid >> 5);
    const int wid = (threadIdx.x >> 5);
    const int nwarps = ((blockDim.x * gridDim.x)>>5);
    const IdxType per_pe_work = ((sim->dim)>>4);
    //Our design require at least (gpu_scale + 4 + log(8)) qubits, because
    //tensorcore demands [8,4]*[4,8], we need at least 8 on B's row dim.
    //This implies for 2 GPUs, we need at least 8 qubit for SV, thus 4 qubit for DM.
    //The smaller condition is difficult to handle due to the complex get_term function,
    //We need to figure out, after the get_term, which row/col it maps to and filter out.
    //But the row and col will be quite difficult to obtain.

    extern __shared__ ValType els[];
    ValType* el_real_s = &els[wid*8*16*2];//per warp 8*16 for real 
    ValType* el_imag_s = &els[wid*8*16*2+8*16]; //per warp 8*16 for imag

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

    wmma::fragment<wmma::matrix_a, 8, 8, 4, ValType, wmma::row_major> a_frag_real;
    wmma::fragment<wmma::matrix_a, 8, 8, 4, ValType, wmma::row_major> a_frag_imag;
    wmma::fragment<wmma::matrix_b, 8, 8, 4, ValType, wmma::col_major> b_frag_real;
    wmma::fragment<wmma::matrix_b, 8, 8, 4, ValType, wmma::col_major> b_frag_imag;
    wmma::fragment<wmma::accumulator, 8, 8, 4, ValType> c_frag_real;
    wmma::fragment<wmma::accumulator, 8, 8, 4, ValType> c_frag_imag;

    for (IdxType i=(warpid<<3); i<per_pe_work; i+=(nwarps<<3)) 
    {
        //load from matrix-B
        for (int j=0; j<8; j++)
        {
            const IdxType term = get_term(i+j,p,q,r,s);
            if (laneid < 16) el_real_s[j*16+laneid] = LOCAL_G(sv_real, term+SV16IDX(laneid));
            else el_imag_s[j*16+laneid-16] = LOCAL_G(sv_imag, term+SV16IDX(laneid-16));
        }
        __syncwarp();

        //============= For matrix-A upper-part ============
        wmma::fill_fragment(c_frag_real, 0.0);
        wmma::fill_fragment(c_frag_imag, 0.0);
        for (unsigned c = 0; c < 4; c++)
        {
            //load A from const mem
            load_matrix_sync(a_frag_real, &gm_real[c*4], 16);
            load_matrix_sync(a_frag_imag, &gm_imag[c*4], 16);
            //load B from shared mem
            load_matrix_sync(b_frag_real, &el_real_s[c*4], 16);
            load_matrix_sync(b_frag_imag, &el_imag_s[c*4], 16);
            //complex multiplication
            mma_sync(c_frag_imag,  a_frag_real, b_frag_imag, c_frag_imag);
            mma_sync(c_frag_imag,  a_frag_imag, b_frag_real, c_frag_imag);
            mma_sync(c_frag_real,  a_frag_real, b_frag_real, c_frag_real);
            a_frag_imag.x[0] = -a_frag_imag.x[0];
            mma_sync(c_frag_real, a_frag_imag, b_frag_imag, c_frag_real);
        }
        //Store first result per segment-C
        IdxType j0 = ((2*laneid)&7);
        IdxType k0 = ((2*laneid)>>3);
        const IdxType term0 = get_term(i+j0,p,q,r,s) + SV16IDX(k0);
        LOCAL_P(sv_real, term0, c_frag_real.x[0]);
        LOCAL_P(sv_imag, term0, c_frag_imag.x[0]);
        //Store second result per segment-C
        IdxType j1 = ((2*laneid+1)&7);
        IdxType k1 = ((2*laneid+1)>>3);
        const IdxType term1 = get_term(i+j1,p,q,r,s) + SV16IDX(k1);
        LOCAL_P(sv_real, term1, c_frag_real.x[1]);
        LOCAL_P(sv_imag, term1, c_frag_imag.x[1]);

        //============= For matrix-A lower-part ============
        wmma::fill_fragment(c_frag_real, 0.0);
        wmma::fill_fragment(c_frag_imag, 0.0);
        for (unsigned c = 0; c < 4; c++)
        {
            //load A from const mem
            load_matrix_sync(a_frag_real, &gm_real[16*8+c*4], 16);
            load_matrix_sync(a_frag_imag, &gm_imag[16*8+c*4], 16);
            //load B from shared mem
            load_matrix_sync(b_frag_real, &el_real_s[c*4], 16);
            load_matrix_sync(b_frag_imag, &el_imag_s[c*4], 16);
            //complex multiplication
            mma_sync(c_frag_imag,  a_frag_real, b_frag_imag, c_frag_imag);
            mma_sync(c_frag_imag,  a_frag_imag, b_frag_real, c_frag_imag);
            mma_sync(c_frag_real,  a_frag_real, b_frag_real, c_frag_real);
            a_frag_imag.x[0] = -a_frag_imag.x[0];
            mma_sync(c_frag_real, a_frag_imag, b_frag_imag, c_frag_real);
        }
        //Store first result per segment-C
        IdxType j2 = ((2*laneid)&7);
        IdxType k2 = ((2*laneid)>>3)+8;
        const IdxType term2 = get_term(i+j2,p,q,r,s) + SV16IDX(k2);
        LOCAL_P(sv_real, term2, c_frag_real.x[0]);
        LOCAL_P(sv_imag, term2, c_frag_imag.x[0]);

        //Store second result per segment-C
        IdxType j3 = ((2*laneid+1)&7);
        IdxType k3 = ((2*laneid+1)>>3)+8;
        const IdxType term3 = get_term(i+j3,p,q,r,s) + SV16IDX(k3);
        LOCAL_P(sv_real, term3, c_frag_real.x[1]);
        LOCAL_P(sv_imag, term3, c_frag_imag.x[1]);
    }
    //BARR;
}

#endif

__device__ __inline__ void M_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, ValType* gm_real, ValType* gm_imag,
        const IdxType qubit, const ValType rand)
{
    grid_group grid = this_grid();
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    ValType * m_real = sim->m_real;
    IdxType mask = ((IdxType)1<<qubit);

    for (IdxType i=tid; i<((IdxType)1<<(sim->n_qubits)); i+=blockDim.x*gridDim.x)
    {
        const ValType val = LOCAL_G(sv_real, (i<<(sim->n_qubits))+i);
        if ( (i & mask) == 0) m_real[i] = 0;
        else m_real[i] = abs(val);
    }
    grid.sync();
    for (IdxType k=((IdxType)1<<(sim->n_qubits-1)); k>0; k>>=1)
    {
        for (IdxType i=tid; i<k; i+=blockDim.x*gridDim.x) m_real[i] += m_real[i+k];
        grid.sync();
    }
    ValType prob_of_one = m_real[0];
    grid.sync();

    if (tid==0)
    {
        if (rand <= prob_of_one)
        {
            gm_real[15] = 1.0/prob_of_one;
        }
        else
        {
            gm_real[0] = 1.0/(1.0-prob_of_one);
        }
    }
    BARR;
    C2_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, qubit, qubit+sim->n_qubits);
    BARR;
    if (tid==0) sim->results_gpu[0] = (rand<=prob_of_one?1:0);
}

__device__ __inline__ void MA_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag,
        const IdxType repetition)
{
    grid_group grid = this_grid();
    const IdxType n_size = (IdxType)1<<(sim->n_qubits);
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    ValType * m_real = sim->m_real;

    BARR;

    for (IdxType i=tid; i<n_size; i+=blockDim.x*gridDim.x)
        m_real[i] = abs(LOCAL_G(sv_real, (i<<(sim->n_qubits))+i));

    BARR;

    //Parallel prefix sum
    for (IdxType d=0; d<(sim->n_qubits); d++)
    {
        IdxType step = (IdxType)1<<(d+1);
        for (IdxType k=tid*step; k<n_size; k+=step*blockDim.x*gridDim.x)
            m_real[k+(1<<(d+1))-1] = m_real[k+(1<<d)-1] + m_real[k+(1<<(d+1))-1];
        grid.sync();
    }


    if (tid == 0)
    {
        ValType purity = fabs(m_real[n_size-1]);
        m_real[n_size-1] = 0;
        if ( abs(purity - 1.0) > ERROR_BAR )
            printf("MA: Purity Check fails with %lf\n", purity);
    }

    grid.sync();
    for (IdxType d=(sim->n_qubits)-1; d>=0; d--)
    {
        IdxType step = (IdxType)1<<(d+1);
        for (IdxType k=tid*step; k<n_size-1; k+=step*blockDim.x*gridDim.x)
        {
            ValType tmp = m_real[k+((IdxType)1<<d)-1];
            m_real[k+((IdxType)1<<d)-1] = m_real[k+((IdxType)1<<(d+1))-1];
            m_real[k+((IdxType)1<<(d+1))-1] = tmp + m_real[k+((IdxType)1<<(d+1))-1];
        }
        grid.sync();
    }

    for (IdxType j=tid; j<n_size; j+=blockDim.x*gridDim.x)
    {
        ValType lower = m_real[j];
        ValType upper = (j+1==n_size)? 1:m_real[j+1];
        for (IdxType i=0; i<repetition; i++)
        {
            ValType r = sim->randoms_gpu[i];
            if (lower<=r && r<upper) sim->results_gpu[i] = j;
        }
    }
    BARR;
}


__device__ __inline__ void RESET_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag,
        const IdxType qubit)
{
    grid_group grid = this_grid();
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    const IdxType per_pe_work_dm = (sim->dim);
    ValType * m_real = sim->m_real;
    IdxType mask = ((IdxType)1<<qubit);

    for (IdxType i=tid; i<((IdxType)1<<(sim->n_qubits)); i+=blockDim.x*gridDim.x)
    {
        if ( (i & mask) == 0) 
            m_real[i] = 0;
        else
            m_real[i] = abs(LOCAL_G(sv_real, (i<<(sim->n_qubits))+i));
    }
    grid.sync();
    for (IdxType k=((IdxType)1<<(sim->n_qubits-1)); k>0; k>>=1)
    {
        for (IdxType i=tid; i<k; i+=blockDim.x*gridDim.x) m_real[i] += m_real[i+k];
        grid.sync();
    }


    BARR;
    ValType prob_of_one = m_real[0];
    grid.sync();

    if (prob_of_one < 1.0) //still possible to normalize
    {
        ValType factor = 1.0/(1.0-prob_of_one);
        for (IdxType i=tid; i<per_pe_work_dm; i+=blockDim.x*gridDim.x)
        {
            if ( (i & mask) == 0)
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
        for (IdxType i=tid; i<per_pe_work_dm; i+=blockDim.x*gridDim.x)
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
    BARR;
}


__device__ __inline__ void Purity_Check(const Simulation* sim, const IdxType t, ValType* sv_real, ValType* sv_imag)
{
    grid_group grid = this_grid();
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    ValType * m_real = sim->m_real;

    BARR;

    for (IdxType i=tid; i<((IdxType)1<<(sim->n_qubits)); i+=blockDim.x*gridDim.x)
    {
        const ValType val = LOCAL_G(sv_real, (i<<(sim->n_qubits))+i);
        m_real[i] = abs(val);
    }
    grid.sync();
    for (IdxType k=((IdxType)1<<(sim->n_qubits-1)); k>0; k>>=1)
    {
        for (IdxType i=tid; i<k; i+=blockDim.x*gridDim.x) m_real[i] += m_real[i+k];
        grid.sync();
    }
    if (threadIdx.x==0 && blockIdx.x==0)
    {
        ValType purity = m_real[0];
        if (abs(purity-1.0) > ERROR_BAR)
        {
            Gate* g = &sim->circuit_handle_gpu[t];
            printf("Purity Check fails after Gate-%lld=>%s(ctrl:%lld,qubit:%lld,theta:%lf) with %lf\n",t,OP_NAMES_NVGPU[g->op_name],g->ctrl,g->qubit,g->theta,purity);
        }
    }
    BARR;
}




//=====================================================================================
//Per-gate execution function
__device__ void Gate::exe_op(Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    grid_group grid = this_grid(); 
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
#ifdef ENABLE_TENSOR_CORES
        C4TCV1_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, ctrl, qubit, ctrl+(sim->n_qubits), qubit+(sim->n_qubits));
#else
        C4_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, ctrl, qubit, ctrl+(sim->n_qubits), qubit+(sim->n_qubits));
#endif

    }
    grid.sync();
}

}; //namespace NWQSim

#endif
