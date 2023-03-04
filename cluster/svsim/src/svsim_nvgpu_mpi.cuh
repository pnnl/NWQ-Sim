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
// File: svsim_nvgpu_mpi.cuh
// MPI and NVSHMEM based implementation of the scale-out SV-Sim gates and 
// simulation runtime using NVIDIA GPU backend.
// ---------------------------------------------------------------------------

#ifndef SVSIM_NVGPU_MPI_CUH
#define SVSIM_NVGPU_MPI_CUH
#include <assert.h>
#include <random>
#include <complex.h>
#include <cooperative_groups.h>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include <string>
#include <iostream>
#include <cuda.h>
#include <nvshmem.h>
#include <nvshmemx.h>
#include <nvshmemx_error.h>
#include <mma.h>
#include "config.h"
#include "gate.h"
#include "metric.hpp"
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

//Simulation runtime
__global__ void simulation_kernel(Simulation*);


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
#ifdef PRINT_CIRCUIT_METRICS
        circuit_metrics(circuit, n_qubits);
#endif
#ifdef DISABLE_GATE_FUSION
        //====================== No Fuse =====================
        SAFE_FREE_GPU(circuit_gpu);
        SAFE_ALOC_GPU(circuit_gpu, n_gates*sizeof(Gate));
        cudaSafeCall(cudaMemcpy(circuit_gpu, circuit.data(), n_gates*sizeof(Gate), cudaMemcpyHostToDevice));
        //====================================================
#else
        //====================== Fuse ========================
        vector<Gate> tmp1_circuit;
        vector<Gate> tmp2_circuit;
        vector<Gate> tmp3_circuit;
        tmp1_circuit.clear();
        tmp2_circuit.clear();

        gate_fusion_1q(circuit, tmp1_circuit, n_qubits);
        gate_fusion_2q_absorb_1q_forward(tmp1_circuit, tmp2_circuit, n_qubits);
        gate_fusion_2q_absorb_1q_backward(tmp2_circuit, tmp3_circuit, n_qubits);
        gate_fusion_2q(tmp3_circuit, fused_circuit, n_qubits);

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
    // number of gpus
    IdxType n_gpus;
    // number of gates
    IdxType n_gates;
    // user input gate sequence
    vector<Gate> circuit;
    // fused gate sequence
    vector<Gate> fused_circuit;
    Gate* circuit_gpu;
};


class Simulation
{
public:
    Simulation(IdxType _n_qubits=N_QUBIT_SLOT) : 
        comm_global(MPI_COMM_WORLD),
        n_qubits(_n_qubits),
        dim((IdxType)1<<(n_qubits)), 
        half_dim((IdxType)1<<(n_qubits-1)),
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
        //set GPUs and communication
        nvshmemx_init_attr_t attr;
        MPI_Comm comm = comm_global;
        attr.mpi_comm = &comm;
        nvshmemx_init_attr(NVSHMEMX_INIT_WITH_MPI_COMM, &attr);
        n_gpus = nvshmem_n_pes();
        i_gpu = nvshmem_my_pe();
        //always be 0 since 1-MPI maps to 1-GPU
        cudaSafeCall(cudaSetDevice(0));
        gpu_scale = floor(log((double)n_gpus+0.5)/log(2.0));
        lg2_m_gpu = n_qubits-gpu_scale;
        m_gpu = ((IdxType)1<<(lg2_m_gpu));
        sv_size_per_gpu = sv_size/n_gpus;
        //CPU side initialization
        if (!is_power_of_2(n_gpus))
            throw std::logic_error("Error: Number of GPUs should be power of 2.");
        if (dim % n_gpus != 0)
            throw std::logic_error("Error: Number of GPUs is too large or too small.");
        if (lg2_m_gpu < 5)
            throw std::logic_error("Error: Each GPU should have at least 5 qubits for multi-node version. Please increase qubits or reduce GPUs");
        //CPU side initialization
        SAFE_ALOC_HOST(sv_real_cpu, sv_size_per_gpu);
        SAFE_ALOC_HOST(sv_imag_cpu, sv_size_per_gpu);
        memset(sv_real_cpu, 0, sv_size_per_gpu);
        memset(sv_imag_cpu, 0, sv_size_per_gpu);
        //State-vector initial state [0..0] = 1
        if (i_gpu == 0) sv_real_cpu[0] = 1.;
        //NVSHMEM GPU memory allocation
        sv_real = (ValType*)nvshmem_malloc(sv_size_per_gpu);
        sv_imag = (ValType*)nvshmem_malloc(sv_size_per_gpu);
        m_real = (ValType*)nvshmem_malloc(sv_size_per_gpu);
        m_imag = (ValType*)nvshmem_malloc(sv_size_per_gpu);
        cudaCheckError(); 
        gpu_mem += sv_size_per_gpu*4;
        //Initialize Circuit 
        circuit_handle = new Circuit(n_qubits);
        circuit_handle->n_gpus = n_gpus;
        circuit_handle_gpu = NULL;
        //GPU memory initilization
        cudaSafeCall(cudaMemcpy(sv_real, sv_real_cpu, sv_size_per_gpu, 
                    cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(sv_imag, sv_imag_cpu, sv_size_per_gpu, 
                    cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemset(m_real, 0, sv_size_per_gpu));
        cudaSafeCall(cudaMemset(m_imag, 0, sv_size_per_gpu));
        SAFE_ALOC_GPU(sim_gpu, sizeof(Simulation));
        rng.seed(time(0));
    }

    ~Simulation()
    {
        //Release circuit
        if (circuit_handle != NULL)
            delete circuit_handle;
        //Release for GPU side
        nvshmem_free(sv_real);
        nvshmem_free(sv_imag);
        nvshmem_free(m_real);
        nvshmem_free(m_imag);
        //Release for CPU side
        SAFE_FREE_HOST(sv_real_cpu);
        SAFE_FREE_HOST(sv_imag_cpu);
        SAFE_FREE_HOST(randoms);
        SAFE_FREE_HOST(results);
        SAFE_FREE_GPU(sim_gpu);
        SAFE_FREE_GPU(randoms_gpu);
        
        //SAFE_FREE_GPU(results_gpu);
        nvshmem_free(results_gpu);
        
        nvshmem_finalize();

    }
    void AllocateQubit()
    {
        circuit_handle->AllocateQubit();
    }
    void ReleaseQubit()
    {
        circuit_handle->ReleaseQubit();
    }
    // =============================== Standard Gates ===================================
    void X(IdxType qubit)
    {
        //Pauli X-gate: bit flip
        /** X = [0 1]
                [1 0]
         */
        Gate* G = new Gate(OP::C1,qubit);
        ValType gm_real[4] = {0,1,1,0};
        ValType gm_imag[4] = {0,0,0,0};
        G->set_gm(gm_real, gm_imag, 2);
        circuit_handle->append(*G);
    }
    void Y(IdxType qubit)
    {
        //Pauli-Y gate: bit and phase flip
        /** Y = [0 -i]
                [i  0]
         */
        Gate* G = new Gate(OP::C1,qubit);
        ValType gm_real[4] = {0,0,0,0};
        ValType gm_imag[4] = {0,-1,1,0};
        G->set_gm(gm_real, gm_imag, 2);
        circuit_handle->append(*G);
    }
    void Z(IdxType qubit)
    {
        //Pauli-Z gate: phase flip
        /** Z = [1  0]
                [0 -1]
         */
        Gate* G = new Gate(OP::C1,qubit);
        ValType gm_real[4] = {1,0,0,-1};
        ValType gm_imag[4] = {0,0,0,0};
        G->set_gm(gm_real, gm_imag, 2);
        circuit_handle->append(*G);
    }
    void H(IdxType qubit)
    {
        //Clifford gate: Hadamard
        /** H = 1/sqrt(2) * [1  1]
                            [1 -1]
         */
        Gate* G = new Gate(OP::C1,qubit);
        ValType gm_real[4] = {S2I,S2I,S2I,-S2I};
        ValType gm_imag[4] = {0,0,0,0};
        G->set_gm(gm_real, gm_imag, 2);
        circuit_handle->append(*G);
    }
    void S(IdxType qubit)
    {
        //Clifford gate: sqrt(Z) phase gate
        /** S = [1 0]
                [0 i]
        */
        Gate* G = new Gate(OP::C1,qubit);
        ValType gm_real[4] = {1,0,0,0};
        ValType gm_imag[4] = {0,0,0,1};
        G->set_gm(gm_real, gm_imag, 2);
        circuit_handle->append(*G);
    }
    void SDG(IdxType qubit)
    {
        //Clifford gate: conjugate of sqrt(Z) phase gate
        /** SDG = [1  0]
                  [0 -i]
        */
        Gate* G = new Gate(OP::C1,qubit);
        ValType gm_real[4] = {1,0,0,0};
        ValType gm_imag[4] = {0,0,0,-1};
        G->set_gm(gm_real, gm_imag, 2);
        circuit_handle->append(*G);
    }
    void T(IdxType qubit)
    {
        //C3 gate: sqrt(S) phase gate
        /** T = [1 0]
                [0 s2i+s2i*i]
        */
        Gate* G = new Gate(OP::C1,qubit);
        ValType gm_real[4] = {1,0,0,S2I};
        ValType gm_imag[4] = {0,0,0,S2I};
        G->set_gm(gm_real, gm_imag, 2);
        circuit_handle->append(*G);
    }
    void TDG(IdxType qubit)
    {
        //C3 gate: conjugate of sqrt(S) phase gate
        /** TDG = [1 0]
                  [0 s2i-s2i*i]
        */
        Gate* G = new Gate(OP::C1,qubit);
        ValType gm_real[4] = {1,0,0,S2I};
        ValType gm_imag[4] = {0,0,0,-S2I};
        G->set_gm(gm_real, gm_imag, 2);
        circuit_handle->append(*G);
    }
    void RI(ValType theta, IdxType qubit)
    {
        //Global phase gate
        /** RI = [e^(ia) 0] = [cos(a)+i*sin(a) 0]
                 [0 e^(ia)]   [0 cos(a)+i*sin(a)]
        */
        Gate* G = new Gate(OP::C1,qubit,-1,theta);
        ValType gm_real[4] = {cos(theta),0,0,cos(theta)};
        ValType gm_imag[4] = {sin(theta),0,0,sin(theta)};
        G->set_gm(gm_real, gm_imag, 2);
        circuit_handle->append(*G);
    }
    void RX(ValType theta, IdxType qubit)
    {
        //Rotation around X axis
        /** RX = [cos(a/2) -i*sin(a/2)]
                 [-i*sin(a/2) cos(a/2)]
        */
        Gate* G = new Gate(OP::C1,qubit,-1,theta);
        ValType gm_real[4] = {cos(HALF*theta),0,0,cos(HALF*theta)};
        ValType gm_imag[4] = {0,-sin(HALF*theta),-sin(HALF*theta),0};
        G->set_gm(gm_real, gm_imag, 2);
        circuit_handle->append(*G);
    }
    void RY(ValType theta, IdxType qubit)
    {
        //Rotation around Y axis
        /** RY = [cos(a/2) -sin(a/2)]
                 [sin(a/2)  cos(a/2)]
        */
        Gate* G = new Gate(OP::C1,qubit,-1,theta);
        ValType gm_real[4] = {cos(HALF*theta),-sin(HALF*theta),sin(HALF*theta),cos(HALF*theta)};
        ValType gm_imag[4] = {0};
        G->set_gm(gm_real, gm_imag, 2);
        circuit_handle->append(*G);
    }
    void RZ(ValType theta, IdxType qubit)
    {
        //Rotation around Z axis
        /** RZ = [cos(a/2)-i*sin(a/2)  0]
                 [0  cos(a/2)+i*sin(a/2)]
        */
        Gate* G = new Gate(OP::C1,qubit,-1,theta);
        ValType gm_real[4] = {cos(HALF*theta),0,0,cos(HALF*theta)};
        ValType gm_imag[4] = {-sin(HALF*theta),0,0,sin(HALF*theta)};
        G->set_gm(gm_real, gm_imag, 2);
        circuit_handle->append(*G);
    }
    void SX(IdxType qubit)
    {
        //sqrt(X) gate, basis gate for IBMQ
        /** SX = 1/2 [1+i 1-i]
                     [1-i 1+i]
        */
        Gate* G = new Gate(OP::C1,qubit);
        ValType gm_real[4] = {HALF,HALF,HALF,HALF};
        ValType gm_imag[4] = {HALF,-HALF,-HALF,HALF};
        G->set_gm(gm_real, gm_imag, 2);
        circuit_handle->append(*G);
    }
    void P(ValType theta, IdxType qubit)
    {
        //Phase gate defined by Qiskit
        /** P = [1, 0     ]  = [1,0]
                [0, e^(ia)]    [0,cos(a)+i*sin(a)]
        */
        Gate* G = new Gate(OP::C1,qubit,-1,theta);
        ValType gm_real[4] = {1,0,0,cos(theta)};
        ValType gm_imag[4] = {0,0,0,sin(theta)};
        G->set_gm(gm_real, gm_imag, 2);
        circuit_handle->append(*G);
    }
    void U(ValType theta, ValType phi, ValType lam, IdxType qubit)
    {
        //Generic single-qubit rotation gate with 3 Euler angles
        /** U = [cos(theta/2), -e^(i*lam)sin(theta/2)]
                [e^(i*phi)sin(theta/2), e^(i*(phi+lam))cos(theta/2)]
        */
        Gate* G = new Gate(OP::C1,qubit,-1,theta);
        ValType gm_real[4] = { cos(HALF*theta), 
                              -cos(lam)*sin(HALF*theta),
                               cos(phi)*sin(HALF*theta),
                               cos(phi+lam)*cos(HALF*theta)};
        ValType gm_imag[4] = { 0,
                              -sin(lam)*sin(HALF*theta),
                               sin(phi)*sin(HALF*theta),
                               sin(lam+phi)*cos(HALF*theta)};
        G->set_gm(gm_real, gm_imag, 2);
        circuit_handle->append(*G);
    }
    void CX(IdxType ctrl, IdxType qubit)
    {
        //Controlled-NOT or CNOT
        /**  CX   = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 0 1]
                    [0 0 1 0]
        */
        Gate* G = new Gate(OP::C2,qubit,ctrl);
        ValType gm_real[16] = {1,0,0,0,
                               0,1,0,0,
                               0,0,0,1,
                               0,0,1,0};
        ValType gm_imag[16] = {0};
        G->set_gm(gm_real, gm_imag, 4);
        circuit_handle->append(*G);
    }
    void CY(IdxType ctrl, IdxType qubit)
    {
        //Controlled-Y
        /**  CY   = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 0 -i]
                    [0 0 i 0]
        */
        Gate* G = new Gate(OP::C2,qubit,ctrl);
        ValType gm_real[16] = {1,0,0,0,
                               0,1,0,0,
                               0,0,0,0,
                               0,0,0,0};
        ValType gm_imag[16] = {0,0,0,0,
                               0,0,0,0,
                               0,0,0,-1,
                               0,0,1,0};
        G->set_gm(gm_real, gm_imag, 4);
        circuit_handle->append(*G);
    }
    void CZ(IdxType ctrl, IdxType qubit)
    {
        //Controlled-Z
        /**  CZ   = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 1 0]
                    [0 0 0 -1]
        */
        Gate* G = new Gate(OP::C2,qubit,ctrl);
        ValType gm_real[16] = {1,0,0,0,
                               0,1,0,0,
                               0,0,1,0,
                               0,0,0,-1};
        ValType gm_imag[16] = {0};
        G->set_gm(gm_real, gm_imag, 4);
        circuit_handle->append(*G);
    }
    void CH(IdxType ctrl, IdxType qubit)
    {
        //Controlled-H
        /**  CH   = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 s2i s2i]
                    [0 0 s2i -s2i]
        */
        Gate* G = new Gate(OP::C2,qubit,ctrl);
        ValType gm_real[16] = {1,0,0,0,
                               0,1,0,0,
                               0,0,S2I,S2I,
                               0,0,S2I,-S2I};
        ValType gm_imag[16] = {0};
        G->set_gm(gm_real, gm_imag, 4);
        circuit_handle->append(*G);
    }
    void CS(IdxType ctrl, IdxType qubit)
    {
        //Controlled-S
        /**  CS   = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 1 0]
                    [0 0 0 i]
        */
        Gate* G = new Gate(OP::C2,qubit,ctrl);
        ValType gm_real[16] = {1,0,0,0,
                               0,1,0,0,
                               0,0,1,0,
                               0,0,0,0};
        ValType gm_imag[16] = {0,0,0,0,
                               0,0,0,0,
                               0,0,0,0,
                               0,0,0,1};
        G->set_gm(gm_real, gm_imag, 4);
        circuit_handle->append(*G);
    }
    void CSDG(IdxType ctrl, IdxType qubit)
    {
        //Controlled-SDG
        /**  CSDG = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 1 0]
                    [0 0 0 -i]
        */
        Gate* G = new Gate(OP::C2,qubit,ctrl);
        ValType gm_real[16] = {1,0,0,0,
                               0,1,0,0,
                               0,0,1,0,
                               0,0,0,0};
        ValType gm_imag[16] = {0,0,0,0,
                               0,0,0,0,
                               0,0,0,0,
                               0,0,0,-1};
        G->set_gm(gm_real, gm_imag, 4);
        circuit_handle->append(*G);
    }
    void CT(IdxType ctrl, IdxType qubit)
    {
        //Controlled-T
        /**  CT   = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 1 0]
                    [0 0 0 s2i+si2*i]
        */
        Gate* G = new Gate(OP::C2,qubit,ctrl);
        ValType gm_real[16] = {1,0,0,0,
                               0,1,0,0,
                               0,0,1,0,
                               0,0,0,S2I};
        ValType gm_imag[16] = {0,0,0,0,
                               0,0,0,0,
                               0,0,0,0,
                               0,0,0,S2I};
        G->set_gm(gm_real, gm_imag, 4);
        circuit_handle->append(*G);
    }
    void CTDG(IdxType ctrl, IdxType qubit)
    {
        //Controlled-TDG
        /**  CTDG = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 1 0]
                    [0 0 0 s2i-si2*i]
        */
        Gate* G = new Gate(OP::C2,qubit,ctrl);
        ValType gm_real[16] = {1,0,0,0,
                               0,1,0,0,
                               0,0,1,0,
                               0,0,0,S2I};
        ValType gm_imag[16] = {0,0,0,0,
                               0,0,0,0,
                               0,0,0,0,
                               0,0,0,-S2I};
        G->set_gm(gm_real, gm_imag, 4);
        circuit_handle->append(*G);
    }
    void CRX(ValType theta, IdxType ctrl, IdxType qubit)
    {
        //Controlled-RX
        /**  CRX  = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 cos(a/2) -i*sin(a/2)]
                    [0 0 -i*sin(a/2) cos(a/2)]
        */
        Gate* G = new Gate(OP::C2,qubit,ctrl,theta);
        ValType gm_real[16] = {1,0,0,0,
                               0,1,0,0,
                               0,0,cos(HALF*theta),0,
                               0,0,0,cos(HALF*theta)};
        ValType gm_imag[16] = {0,0,0,0,
                               0,0,0,0,
                               0,0,0,-sin(HALF*theta),
                               0,0,-sin(HALF*theta),0};
        G->set_gm(gm_real, gm_imag, 4);
        circuit_handle->append(*G);
    }
    void CRY(ValType theta, IdxType ctrl, IdxType qubit)
    {
        //Controlled-RY
        /**  CRY  = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 cos(a/2) -sin(a/2)]
                    [0 0 sin(a/2)  cos(a/2)]
        */
        Gate* G = new Gate(OP::C2,qubit,ctrl,theta);
        ValType gm_real[16] = {1,0,0,0,
                               0,1,0,0,
                               0,0,cos(HALF*theta),-sin(HALF*theta),
                               0,0,sin(HALF*theta),cos(HALF*theta)};
        ValType gm_imag[16] = {0};
        G->set_gm(gm_real, gm_imag, 4);
        circuit_handle->append(*G);
    }
    void CRZ(ValType theta, IdxType ctrl, IdxType qubit)
    {
        //Controlled-RZ
        /**  CRZ  = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 cos(a/2)-i*sin(a/2)  0]
                    [0 0 0  cos(a/2)+i*sin(a/2)]
        */
        Gate* G = new Gate(OP::C2,qubit,ctrl,theta);
        ValType gm_real[16] = {1,0,0,0,
                               0,1,0,0,
                               0,0,cos(HALF*theta),0,
                               0,0,0,cos(HALF*theta)};
        ValType gm_imag[16] = {0,0,0,0,
                               0,0,0,0,
                               0,0,-sin(HALF*theta),0,
                               0,0,0,sin(HALF*theta)};
        G->set_gm(gm_real, gm_imag, 4);
        circuit_handle->append(*G);
    }
    void CSX(IdxType ctrl, IdxType qubit)
    {
        //Controlled-SX
        /**  CSX  = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 (1+i)/2 (1-i)/2]
                    [0 0 (1-i)/2 (1+i)/2]
        */
        Gate* G = new Gate(OP::C2,qubit,ctrl);
        ValType gm_real[16] = {1,0,0,0,
                               0,1,0,0,
                               0,0,HALF,HALF,
                               0,0,HALF,HALF};
        ValType gm_imag[16] = {0,0,0,0,
                               0,0,0,0,
                               0,0,HALF,-HALF,
                               0,0,-HALF,HALF};
        G->set_gm(gm_real, gm_imag, 4);
        circuit_handle->append(*G);
    }
    void CP(ValType theta, IdxType ctrl, IdxType qubit)
    {
        //Controlled-P
        /**  CP   = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 1 0]
                    [0 0 0 cos(a)+i*sin(a)]
        */
        Gate* G = new Gate(OP::C2,qubit,ctrl,theta);
        ValType gm_real[16] = {1,0,0,0,
                               0,1,0,0,
                               0,0,1,0,
                               0,0,0,cos(theta)};
        ValType gm_imag[16] = {0,0,0,0,
                               0,0,0,0,
                               0,0,0,0,
                               0,0,0,sin(theta)};
        G->set_gm(gm_real, gm_imag, 4);
        circuit_handle->append(*G);
    }
    void CU(ValType theta, ValType phi, ValType lam, ValType gamma,
            IdxType ctrl, IdxType qubit)
    {
        //Controlled-U, w.s.p. to Qiksik: https://qiskit.org/documentation/stubs/qiskit.circuit.library.CUGate.html
        /**  CU   = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 e^(i*gamma)cos(theta/2), -e^(i*(gamma+lam))sin(theta/2)]
                    [0 0 e^(i*(gamma+phi))sin(theta/2), e^(i*(gamma+phi+lam))cos(theta/2)]
        */
        Gate* G = new Gate(OP::C2,qubit,ctrl,theta);
        ValType gm_real[16] = {1,0,0,0,
                               0,1,0,0,
                               0,0,cos(gamma)*cos(HALF*theta), -cos(gamma+lam)*sin(HALF*theta),
                               0,0,cos(gamma+phi)*sin(HALF*theta),cos(gamma+phi+lam)*cos(HALF*theta)};
        ValType gm_imag[16] = {0,0,0,0,
                               0,0,0,0,
                               0,0,sin(gamma)*cos(HALF*theta), -sin(gamma+lam)*sin(HALF*theta),
                               0,0,sin(gamma+phi)*sin(HALF*theta),sin(gamma+phi+lam)*cos(HALF*theta)};
        G->set_gm(gm_real, gm_imag, 4);
        circuit_handle->append(*G);
    }
    void ID(IdxType qubit)
    {
        //Identity gate
        /** ID  = [1 0]
                  [0 1]
        */
        Gate* G = new Gate(OP::C1,qubit);
        ValType gm_real[4] = {1,0,0,1};
        ValType gm_imag[4] = {0,0,0,0};
        G->set_gm(gm_real, gm_imag, 2);
        circuit_handle->append(*G);
    }
    void SWAP(IdxType ctrl, IdxType qubit)
    {
        //SWAP gate
        /**  SWAP = [1 0 0 0]
                    [0 0 1 0]
                    [0 1 0 0]
                    [0 0 0 1]
        */
        Gate* G = new Gate(OP::C2,ctrl,qubit);
        ValType gm_real[16] = {1,0,0,0,
                               0,0,1,0,
                               0,1,0,0,
                               0,0,0,1};
        ValType gm_imag[16] = {0,0,0,0,
                               0,0,0,0,
                               0,0,0,0,
                               0,0,0,0};
        G->set_gm(gm_real, gm_imag, 4);
        circuit_handle->append(*G);
    }
    void M(IdxType qubit) //default is pauli-Z
    {
        SAFE_FREE_HOST(results);
        SAFE_ALOC_HOST(results, sizeof(IdxType));
        memset(results, 0, sizeof(IdxType));
        //SAFE_FREE_GPU(results_gpu);
        //SAFE_ALOC_GPU(results_gpu, sizeof(IdxType));

        nvshmem_free(results_gpu);
        results_gpu = (IdxType*)nvshmem_malloc(sizeof(IdxType));


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
        //SAFE_FREE_GPU(results_gpu);
        //SAFE_ALOC_GPU(results_gpu, sizeof(IdxType)*repetition);
        
        nvshmem_free(results_gpu);
        results_gpu = (IdxType*)nvshmem_malloc(sizeof(IdxType)*repetition);

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

    // ============================== Other Gate Definition ================================
    void U3(ValType theta, ValType phi, ValType lam, IdxType qubit)
    {
        U(theta, phi, lam, qubit);
    }
    void U2(ValType phi, ValType lam, IdxType qubit)
    {
        U(PI/2, phi, lam, qubit);
    }
    void U1(ValType lam, IdxType qubit)
    {
        U(0, 0, lam, qubit);
    }
    void CCX(IdxType qubit0, IdxType qubit1, IdxType qubit2)
    {
        H(qubit2);
        CX(qubit1, qubit2);
        TDG(qubit2);
        CX(qubit0, qubit2);
        T(qubit2);
        CX(qubit1, qubit2);
        T(qubit1);
        TDG(qubit2);
        CX(qubit0, qubit2);
        T(qubit2);
        CX(qubit0, qubit1);
        T(qubit0);
        TDG(qubit1);
        H(qubit2);
        CX(qubit0, qubit1);
    }
    void CSWAP(IdxType qubit0, IdxType qubit1, IdxType qubit2)
    {
        CX(qubit2, qubit1);
        CCX(qubit0, qubit1, qubit2);
        CX(qubit2, qubit1);
    }
    void RXX(ValType theta, IdxType qubit0, IdxType qubit1)
    {
        H(qubit0);
        H(qubit1);
        CX(qubit0, qubit1);
        RZ(theta, qubit1);
        CX(qubit0, qubit1);
        H(qubit0);
        H(qubit1);
    }
    void RYY(ValType theta, IdxType qubit0, IdxType qubit1)
    {
        RX(PI/2, qubit0);
        RX(PI/2, qubit1);
        CX(qubit0, qubit1);
        RZ(theta, qubit1);
        CX(qubit0, qubit1);
        RX(PI/2, qubit0);
        RX(PI/2, qubit1);
    }
    void RZZ(ValType theta, IdxType qubit0, IdxType qubit1)
    {
        CX(qubit0, qubit1);
        RZ(theta, qubit1);
        CX(qubit0, qubit1);
    }



    // =============================== End of Gate Define ===================================
    void reset_sim()
    {
        //for (IdxType i=0; i<circuit_handle->n_qubits; i++)
        //circuit_handle->ReleaseQubit();
        
        //Reset CPU input & output
        memset(sv_real_cpu, 0, sv_size_per_gpu);
        memset(sv_imag_cpu, 0, sv_size_per_gpu);
        //State Vector initial state [0..0] = 1
        if (i_gpu == 0) sv_real_cpu[0] = 1.;
        //GPU side initialization
        cudaSafeCall(cudaMemcpy(sv_real, sv_real_cpu, 
                    sv_size_per_gpu, cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(sv_imag, sv_imag_cpu, 
                    sv_size_per_gpu, cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemset(m_real, 0, sv_size_per_gpu));
        cudaSafeCall(cudaMemset(m_imag, 0, sv_size_per_gpu));
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
        rng.seed(seed);
    }
    void clear_circuit()
    {
        circuit_handle->clear();
    }
    void update(const IdxType _n_qubits, const IdxType _n_gates)
    {
        this->n_gates = _n_gates;
        
        //assert(_n_qubits <= (N_QUBIT_SLOT));
        //this->n_qubits = _n_qubits;
        //this->dim = ((IdxType)1<<(n_qubits));
        //this->half_dim = (IdxType)1<<(n_qubits-1);
        //this->sv_size = dim*(IdxType)sizeof(ValType);
        //this->lg2_m_gpu = n_qubits-gpu_scale;
        //this->m_gpu = ((IdxType)1<<(lg2_m_gpu));
        //this->sv_size_per_gpu = sv_size/n_gpus;
    }
    std::string circuitToString()
    {
        return circuit_handle->circuitToString();
    }
    void sim()
    {
        cudaSafeCall(cudaSetDevice(0));
#ifdef PRINT_KERNEL_TRACE
        double* sim_times;
        double sim_time;
        gpu_timer sim_timer;
        if (i_gpu == 0)
        {
            SAFE_ALOC_HOST(sim_times, sizeof(double)*n_gpus);
            memset(sim_times, 0, sizeof(double)*n_gpus);
        }
        IdxType input_gates = circuit_handle->n_gates;
#endif
        circuit_handle_gpu = circuit_handle->upload();
        //update should be put after upload where gate fusion is applied
        //which may change the number of gates in the circuit
        update(circuit_handle->n_qubits, circuit_handle->n_gates);
        cudaSafeCall(cudaMemcpy(sim_gpu, this, 
                    sizeof(Simulation), cudaMemcpyHostToDevice));
#ifdef PRINT_SIM_TRACE
        printf("SVSim_gpu is running! Requesting %lld qubits.\n", circuit_handle->n_qubits);
#endif
        dim3 gridDim(1,1,1);
        cudaDeviceProp deviceProp;
        cudaSafeCall(cudaGetDeviceProperties(&deviceProp, 0));
        //8*16 is per warp shared-memory usage for C4 TC, with real and imag
        //unsigned smem_size = THREADS_CTA_NVGPU/32*8*16*2*sizeof(ValType);
        unsigned smem_size = 0*sizeof(ValType);
        int numBlocksPerSm;
        cudaSafeCall(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm, 
                    simulation_kernel, THREADS_CTA_NVGPU, smem_size));
        gridDim.x = numBlocksPerSm * deviceProp.multiProcessorCount;
        void* args[] = {&sim_gpu};
        cudaSafeCall(cudaDeviceSynchronize());
#ifdef PRINT_KERNEL_TRACE
        MPI_Barrier(MPI_COMM_WORLD);
        sim_timer.start_timer();
#endif
        NVSHMEM_CHECK(nvshmemx_collective_launch((const void*)simulation_kernel,gridDim,
                    THREADS_CTA_NVGPU,args,smem_size,0));
        cudaSafeCall(cudaDeviceSynchronize());
#ifdef PRINT_KERNEL_TRACE
        sim_timer.stop_timer();
        sim_time = sim_timer.measure();
#endif
        cudaCheckError();
#ifdef PRINT_KERNEL_TRACE
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(&sim_time, 1, MPI_DOUBLE,
                &sim_times[i_gpu], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (i_gpu ==0)
        {
            double avg_sim_time = 0;
            for (unsigned d=0; d<n_gpus; d++)
            {
                avg_sim_time += sim_times[d];
            }
            avg_sim_time /= (double)n_gpus;
            printf("\n============== SV-Sim ===============\n");
            printf("nqubits:%lld, ngates:%lld, sim_gates:%lld, ngpus:%lld, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms, mem:%.3lf MB, mem_per_gpu:%.3lf MB\n",
                    n_qubits, input_gates, n_gates, n_gpus, avg_sim_time, 0., 
                    avg_sim_time, gpu_mem/1024/1024*n_gpus, gpu_mem/1024/1024);
            printf("=====================================\n");
            SAFE_FREE_HOST(sim_times);
        }
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
        cudaSafeCall(cudaMemcpy(sv_real_cpu, sv_real, sv_size_per_gpu, cudaMemcpyDeviceToHost));
        cudaSafeCall(cudaMemcpy(sv_imag_cpu, sv_imag, sv_size_per_gpu, cudaMemcpyDeviceToHost));

        ValType* sv_diag_real = NULL;
        ValType* sv_diag_imag = NULL;
        if (i_gpu == 0) SAFE_ALOC_HOST(sv_diag_real, dim*sizeof(ValType));
        if (i_gpu == 0) SAFE_ALOC_HOST(sv_diag_imag, dim*sizeof(ValType));

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gather(sv_real_cpu, m_gpu, MPI_DOUBLE,
                &sv_diag_real[i_gpu*m_gpu], m_gpu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather(sv_imag_cpu, m_gpu, MPI_DOUBLE,
                &sv_diag_imag[i_gpu*m_gpu], m_gpu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        if (i_gpu == 0) 
        {
            IdxType num = ((IdxType)1<<n_qubits);
            printf("----- SVSim ------\n");
            for (IdxType i=0; i<num; i++) 
            {
                printf("(%.3lf,%.3lfj) ", sv_diag_real[i], sv_diag_imag[i]);
                if ((i+1)%8==0) printf("\n");
            }
            printf("\n");
            SAFE_FREE_HOST(sv_diag_real);
            SAFE_FREE_HOST(sv_diag_imag);
        }
    }
public:
    // n_qubits is the number of qubits
    IdxType n_qubits;
    // which gpu
    IdxType i_gpu;
    IdxType gpu_scale;
    IdxType n_gpus;
    IdxType lg2_m_gpu;
    IdxType m_gpu;
    IdxType sv_size_per_gpu;

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
    MPI_Comm comm_global;
};

#define PGAS_P(arr,i,val) nvshmem_double_p(&(arr)[(i)&((sim->m_gpu)-1)], (val), ((i)>>(sim->lg2_m_gpu)) )
#define PGAS_G(arr,i) nvshmem_double_g(&(arr)[(i)&((sim->m_gpu)-1)], ((i)>>(sim->lg2_m_gpu)) )
#define BARR if(threadIdx.x==0 && blockIdx.x==0) nvshmem_barrier_all(); grid.sync();

//============== Check Trace (debug purpose) ================
__device__ __inline__ void CHECK_TRACE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, IdxType t)
{
    grid_group grid = this_grid();
    if (sim->i_gpu == 0 && blockIdx.x ==0 && threadIdx.x == 0)
    {
        ValType trace = 0;
        for (IdxType i=0; i<((IdxType)1<<(sim->n_qubits)); i++)
        {
            const ValType val = PGAS_G(sv_real, i);
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
        //if (sim->i_gpu == 0 && blockIdx.x == 0 && threadIdx.x == 0) 
        //printf("==%lld== %s(qubit:%lld, ctrl:%lld, theta:%lf)\n",sim->circuit_handle_gpu[t].op_name, OP_NAMES_NVGPU[sim->circuit_handle_gpu[t].op_name], sim->circuit_handle_gpu[t].qubit, sim->circuit_handle_gpu[t].ctrl, sim->circuit_handle_gpu[t].theta);

        //IdxType t0;
        //if (sim->i_gpu == 0 && blockIdx.x == 0 && threadIdx.x == 0) t0 = clock64();

        ((sim->circuit_handle_gpu)[t]).exe_op(sim, sim->sv_real, sim->sv_imag);

        /*
        if (sim->i_gpu == 0 && blockIdx.x == 0 && threadIdx.x == 0)
        {
            IdxType time = clock64() - t0;
            printf("G%lld:%s(ctrl:%lld,qubit:%lld) ticks: %lld\n",t,
                    OP_NAMES_NVGPU[sim->circuit_handle_gpu[t].op_name],
                    sim->circuit_handle_gpu[t].ctrl, 
                    sim->circuit_handle_gpu[t].qubit, 
                    time);
        }
         */
        //CHECK_TRACE(sim, sim->sv_real, sim->sv_imag, t);
#ifdef PURITY_CHECK
        Purity_Check(sim, t, sim->sv_real, sim->sv_imag);
#endif
    }
}

//================================= Gate Definition ========================================

#define LOCAL_G(arr,i) arr[(i)&(sim->m_gpu-1)]
#define LOCAL_P(arr,i,val) arr[(i)&(sim->m_gpu-1)] = val;

//============== Local Unified 1-qubit Gate ================
__device__ __inline__ void C1_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag, const IdxType qubit)
{
    grid_group grid = this_grid(); 
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x; 
    const IdxType per_pe_work = ((sim->half_dim)>>(sim->gpu_scale));

    for (IdxType i=(sim->i_gpu)*per_pe_work+tid; i<(sim->i_gpu+1)*per_pe_work; i+=blockDim.x*gridDim.x)
    { 
        IdxType outer = (i >> qubit); 
        IdxType inner =  (i & (((IdxType)1<<qubit)-1)); 
        IdxType offset = (outer << (qubit+1)); 
        IdxType pos0 = ((offset + inner) & (sim->m_gpu-1)); 
        IdxType pos1 = ((offset + inner + ((IdxType)1<<qubit)) & (sim->m_gpu-1)); 
        const ValType el0_real = LOCAL_G(sv_real, pos0);
        const ValType el0_imag = LOCAL_G(sv_imag, pos0);
        const ValType el1_real = LOCAL_G(sv_real, pos1);
        const ValType el1_imag = LOCAL_G(sv_imag, pos1);
        ValType sv_real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag)
                              +(gm_real[1] * el1_real) - (gm_imag[1] * el1_imag);
        ValType sv_imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real)
                              +(gm_real[1] * el1_imag) + (gm_imag[1] * el1_real);
        ValType sv_real_pos1 = (gm_real[2] * el0_real) - (gm_imag[2] * el0_imag)
                              +(gm_real[3] * el1_real) - (gm_imag[3] * el1_imag);
        ValType sv_imag_pos1 = (gm_real[2] * el0_imag) + (gm_imag[2] * el0_real)
                              +(gm_real[3] * el1_imag) + (gm_imag[3] * el1_real);
        LOCAL_P(sv_real, pos0, sv_real_pos0);
        LOCAL_P(sv_imag, pos0, sv_imag_pos0);
        LOCAL_P(sv_real, pos1, sv_real_pos1);
        LOCAL_P(sv_imag, pos1, sv_imag_pos1);
    } 
    grid.sync(); 
}

/** The C1V1 version relies on bidirectional communication between each GPU pair. However, on Perlmutter ss11 with the recent slingshot-V2 upgrade, bidirectional communication
  encounters a bug. If that is the case, we use C1V2.
*/
__device__ __inline__ void C1V1_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag, const IdxType qubit)
{
    grid_group grid = this_grid(); 
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x; 
    const IdxType per_pe_work = ((sim->half_dim)>>(sim->gpu_scale));

    if (qubit < sim->lg2_m_gpu) 
    {
        C1_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, qubit);
    }
    else
    {
        IdxType pair_gpu = (sim->i_gpu)^((IdxType)1<<(qubit-(sim->lg2_m_gpu)));
        IdxType flag = (sim->i_gpu<pair_gpu) ? 0:1; 
        ValType* sv_real_pos0 = NULL;
        ValType* sv_imag_pos0 = NULL;
        ValType* sv_real_pos1 = NULL;
        ValType* sv_imag_pos1 = NULL;
        if (flag == 0)//copy right to local
        {
            sv_real_pos0 = sv_real;
            sv_imag_pos0 = sv_imag;
            sv_real_pos1 = sim->m_real;
            sv_imag_pos1 = sim->m_imag;
            if (tid == 0) nvshmem_double_get(sv_real_pos1, sv_real, per_pe_work, pair_gpu);
            if (tid == 0) nvshmem_double_get(sv_imag_pos1, sv_imag, per_pe_work, pair_gpu);
        }
        else //copy left to local
        {
            sv_real_pos0 = sim->m_real;
            sv_imag_pos0 = sim->m_imag;
            sv_real_pos1 = &sv_real[per_pe_work];
            sv_imag_pos1 = &sv_imag[per_pe_work];
            if (tid == 0) nvshmem_double_get(sv_real_pos0, &sv_real[per_pe_work], per_pe_work, pair_gpu );
            if (tid == 0) nvshmem_double_get(sv_imag_pos0, &sv_imag[per_pe_work], per_pe_work, pair_gpu );
        }
        grid.sync();
        for (IdxType i=tid; i<per_pe_work; i+=blockDim.x*gridDim.x)
        {
            const ValType el0_real = LOCAL_G(sv_real_pos0, i);
            const ValType el0_imag = LOCAL_G(sv_imag_pos0, i);
            const ValType el1_real = LOCAL_G(sv_real_pos1, i);
            const ValType el1_imag = LOCAL_G(sv_imag_pos1, i);
            ValType real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag)
                               +(gm_real[1] * el1_real) - (gm_imag[1] * el1_imag);
            ValType imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real)
                               +(gm_real[1] * el1_imag) + (gm_imag[1] * el1_real);
            ValType real_pos1 = (gm_real[2] * el0_real) - (gm_imag[2] * el0_imag)
                               +(gm_real[3] * el1_real) - (gm_imag[3] * el1_imag);
            ValType imag_pos1 = (gm_real[2] * el0_imag) + (gm_imag[2] * el0_real)
                               +(gm_real[3] * el1_imag) + (gm_imag[3] * el1_real);
            LOCAL_P(sv_real_pos0, i, real_pos0);
            LOCAL_P(sv_imag_pos0, i, imag_pos0);
            LOCAL_P(sv_real_pos1, i, real_pos1);
            LOCAL_P(sv_imag_pos1, i, imag_pos1);
        }
        grid.sync();
        if (flag == 0)//copy local to right
        {
            if (tid == 0) nvshmem_double_put(sv_real, sv_real_pos1, per_pe_work, pair_gpu);
            if (tid == 0) nvshmem_double_put(sv_imag, sv_imag_pos1, per_pe_work, pair_gpu);
        }
        else //copy local to left
        {
            if (tid == 0) nvshmem_double_put(&sv_real[per_pe_work], sv_real_pos0, per_pe_work, pair_gpu);
            if (tid == 0) nvshmem_double_put(&sv_imag[per_pe_work], sv_imag_pos0, per_pe_work, pair_gpu);
        }
    }
}


/** This is a less optimized version. Half of the GPUs stay idle for better communication efficiency. Use this version only when bidirectional NVSHMEM communicaiton is not well-supported.
*/
__device__ __inline__ void C1V2_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag, const IdxType qubit)
{
    grid_group grid = this_grid(); 
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x; 
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale));

    if (qubit < sim->lg2_m_gpu) 
    {
        C1_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, qubit);
    }
    else
    {
        IdxType pair_gpu = (sim->i_gpu)^((IdxType)1<<(qubit-(sim->lg2_m_gpu)));
        if (sim->i_gpu>pair_gpu) return;
        ValType* sv_real_remote = sim->m_real;
        ValType* sv_imag_remote = sim->m_imag;
        if (tid == 0) nvshmem_double_get(sv_real_remote, sv_real, per_pe_work, pair_gpu );
        if (tid == 0) nvshmem_double_get(sv_imag_remote, sv_imag, per_pe_work, pair_gpu );
        grid.sync();

        for (IdxType i=tid; i<per_pe_work; i+=blockDim.x*gridDim.x)
        {
            const ValType el0_real = LOCAL_G(sv_real, i);
            const ValType el0_imag = LOCAL_G(sv_imag, i);
            const ValType el1_real = LOCAL_G(sv_real_remote, i);
            const ValType el1_imag = LOCAL_G(sv_imag_remote, i);

            ValType real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag)
                               +(gm_real[1] * el1_real) - (gm_imag[1] * el1_imag);
            ValType imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real)
                               +(gm_real[1] * el1_imag) + (gm_imag[1] * el1_real);
            ValType real_pos1 = (gm_real[2] * el0_real) - (gm_imag[2] * el0_imag)
                               +(gm_real[3] * el1_real) - (gm_imag[3] * el1_imag);
            ValType imag_pos1 = (gm_real[2] * el0_imag) + (gm_imag[2] * el0_real)
                               +(gm_real[3] * el1_imag) + (gm_imag[3] * el1_real);
            LOCAL_P(sv_real, i, real_pos0);
            LOCAL_P(sv_imag, i, imag_pos0);
            LOCAL_P(sv_real_remote, i, real_pos1);
            LOCAL_P(sv_imag_remote, i, imag_pos1);
        }
        grid.sync();
        if (tid == 0) nvshmem_double_put(sv_real, sv_real_remote, per_pe_work, pair_gpu);
        if (tid == 0) nvshmem_double_put(sv_imag, sv_imag_remote, per_pe_work, pair_gpu);
    }
}

//============== Local 2-qubit Gate  ================
__device__ __inline__ void C2_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag, const IdxType qubit0, const IdxType qubit1)
{
    grid_group grid = this_grid(); 
    const int tid = blockDim.x * blockIdx.x + threadIdx.x; 
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale+2));
    assert (qubit0 != qubit1); //Non-cloning
    
    const IdxType q0dim = ((IdxType)1 << max(qubit0, qubit1) );
    const IdxType q1dim = ((IdxType)1 << min(qubit0, qubit1) );
    const IdxType outer_factor = ((sim->dim) + q0dim + q0dim - 1) >> (max(qubit0,qubit1)+1);
    const IdxType mider_factor = (q0dim + q1dim + q1dim - 1) >> (min(qubit0,qubit1)+1);
    const IdxType inner_factor = q1dim;
    const IdxType qubit0_dim = ((IdxType)1 << qubit0);
    const IdxType qubit1_dim = ((IdxType)1 << qubit1);

    for (IdxType i=(sim->i_gpu)*per_pe_work+tid; i<(sim->i_gpu+1)*per_pe_work;
            i+=blockDim.x*gridDim.x) 
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

//============== Unified 2-qubit Gate ================
//Perform communication optimization here
__device__ __inline__ void C2V1_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag, const IdxType qubit0, const IdxType qubit1)
{
    assert (qubit0 != qubit1); //Non-cloning

    if (qubit0 < sim->lg2_m_gpu && qubit1 < sim->lg2_m_gpu) 
    {
        C2_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, qubit0, qubit1);
    }
    else
    {
        grid_group grid = this_grid(); 
        const int tid = blockDim.x * blockIdx.x + threadIdx.x; 
        //For processing a remote qubit, half of the GPUs will be idle for
        //better communication efficiency. Depending on qubit position, 
        //for a GPU pair, we only use GPUs with smaller ids. 
        //Consequently, each GPU should take double of the workload than before
        //Therefore, here it is gpu_scale+1 not gpu_scale+2
        const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale+1)); 
        const IdxType per_pe_num = ((sim->dim)>>(sim->gpu_scale));
        const IdxType p = min(qubit0, qubit1);
        const IdxType q = max(qubit0, qubit1);

        //load data from pair GPU
        IdxType pair_gpu = (sim->i_gpu)^((IdxType)1<<(q-(sim->lg2_m_gpu)));
        if (sim->i_gpu>pair_gpu) return;

        ValType* sv_real_remote = sim->m_real;
        ValType* sv_imag_remote = sim->m_imag;

        if (tid == 0) nvshmem_double_get(sv_real_remote, sv_real, per_pe_num, pair_gpu );
        if (tid == 0) nvshmem_double_get(sv_imag_remote, sv_imag, per_pe_num, pair_gpu );
        grid.sync();

        for (IdxType i=(sim->i_gpu)*per_pe_work+tid; i<(sim->i_gpu+1)*per_pe_work;
                i+=blockDim.x*gridDim.x) 
        {
            ValType el_real[4];
            ValType el_imag[4];
            ValType res_real[4] = {0};
            ValType res_imag[4] = {0};
            const IdxType term0 = MOD2E(i,p);
            const IdxType term1 = MOD2E(DIV2E(i,p),q-p-1)*EXP2E(p+1);
            const IdxType term2 = DIV2E(DIV2E(i,p),q-p-1)*EXP2E(q+1);
            const IdxType term = term2 + term1 + term0;
            el_real[0] = LOCAL_G(sv_real, term+SV4IDX(0));
            el_imag[0] = LOCAL_G(sv_imag, term+SV4IDX(0));
            el_real[3] = LOCAL_G(sv_real_remote, term+SV4IDX(3));
            el_imag[3] = LOCAL_G(sv_imag_remote, term+SV4IDX(3));

            if (qubit0 == q) //qubit0 is the remote qubit
            {
                el_real[1] = LOCAL_G(sv_real, term+SV4IDX(1));
                el_imag[1] = LOCAL_G(sv_imag, term+SV4IDX(1));
                el_real[2] = LOCAL_G(sv_real_remote, term+SV4IDX(2));
                el_imag[2] = LOCAL_G(sv_imag_remote, term+SV4IDX(2));
            }
            else //qubit1 is the remote qubit
            {
                el_real[1] = LOCAL_G(sv_real_remote, term+SV4IDX(1));
                el_imag[1] = LOCAL_G(sv_imag_remote, term+SV4IDX(1));
                el_real[2] = LOCAL_G(sv_real, term+SV4IDX(2));
                el_imag[2] = LOCAL_G(sv_imag, term+SV4IDX(2));
            }
            #pragma unroll
            for (unsigned j=0; j<4; j++)
            {
                #pragma unroll
                for (unsigned k=0; k<4; k++)
                {
                    res_real[j] += (el_real[k] * gm_real[j*4+k]) - (el_imag[k] * gm_imag[j*4+k]);
                    res_imag[j] += (el_real[k] * gm_imag[j*4+k]) + (el_imag[k] * gm_real[j*4+k]);
                }
            }
            LOCAL_P(sv_real, term+SV4IDX(0), res_real[0]);
            LOCAL_P(sv_imag, term+SV4IDX(0), res_imag[0]);
            LOCAL_P(sv_real_remote, term+SV4IDX(3), res_real[3]);
            LOCAL_P(sv_imag_remote, term+SV4IDX(3), res_imag[3]);

            if (qubit0 == q) //qubit0 is the remote qubit
            {
                LOCAL_P(sv_real, term+SV4IDX(1), res_real[1]);
                LOCAL_P(sv_imag, term+SV4IDX(1), res_imag[1]);
                LOCAL_P(sv_real_remote, term+SV4IDX(2), res_real[2]);
                LOCAL_P(sv_imag_remote, term+SV4IDX(2), res_imag[2]);
            }
            else //qubit1 is the remote qubit
            {
                LOCAL_P(sv_real_remote, term+SV4IDX(1), res_real[1]);
                LOCAL_P(sv_imag_remote, term+SV4IDX(1), res_imag[1]);
                LOCAL_P(sv_real, term+SV4IDX(2), res_real[2]);
                LOCAL_P(sv_imag, term+SV4IDX(2), res_imag[2]);
            }
        }
        grid.sync();
        if (tid == 0) nvshmem_double_put(sv_real, sv_real_remote, per_pe_num, pair_gpu);
        if (tid == 0) nvshmem_double_put(sv_imag, sv_imag_remote, per_pe_num, pair_gpu);
        //BARR;
    }
}

__device__ __inline__ void M_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        ValType* gm_real, ValType* gm_imag, const IdxType qubit, const ValType rand)
{
    grid_group grid = this_grid();
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale));
    ValType * m_real = sim->m_real;
    IdxType mask = ((IdxType)1<<qubit);

    for (IdxType i=tid; i<sim->m_gpu; i+=blockDim.x*gridDim.x)
    {
        IdxType idx = (sim->i_gpu)*per_pe_work+i;
        if ( (idx & mask) == 0) 
            m_real[i] = 0;
        else 
            m_real[i] = sv_real[i]*sv_real[i] + sv_imag[i]*sv_imag[i];
    }
    BARR;

    //Parallel reduction
    for (IdxType k=((IdxType)1<<(sim->n_qubits-1)); k>0; k>>=1)
    {
        for (IdxType i=tid; i<k; i+=blockDim.x*gridDim.x) 
        {
            IdxType local_gpu = i>>(sim->lg2_m_gpu);
            if (local_gpu == sim->i_gpu)
            {
                IdxType local_idx = i&(sim->m_gpu-1);
                m_real[local_idx] += PGAS_G(m_real, i+k);
            }
        }
        BARR;
    }

    if (tid==0 && sim->i_gpu!=0) m_real[0] = PGAS_G(m_real,0);
    grid.sync();
    ValType prob_of_one = m_real[0];
    grid.sync();
 
    if (rand < prob_of_one)
    {
        ValType factor = 1./sqrt(prob_of_one);
        for (IdxType i=tid; i<sim->m_gpu; i+=blockDim.x*gridDim.x)
        {
            IdxType idx = (sim->i_gpu)*per_pe_work+i;
            if ( (idx & mask) == 0)
            {
                sv_real[i] = 0;
                sv_imag[i] = 0;
            }
            else
            {
                sv_real[i] *= factor;
                sv_imag[i] *= factor;
            }
        }
    }
    else
    {
        ValType factor = 1./sqrt(1.-prob_of_one); 
        for (IdxType i=tid; i<sim->m_gpu; i+=blockDim.x*gridDim.x)
        {
            IdxType idx = (sim->i_gpu)*per_pe_work+i;
            if ( (idx & mask) == 0)
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
    if (tid==0) sim->results_gpu[0] = (rand<=prob_of_one?1:0);
    BARR;
}

__device__ __inline__ void MA_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag,
        const IdxType repetition)
{
    grid_group grid = this_grid();
    const IdxType n_size = (IdxType)1<<(sim->n_qubits);
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    ValType * m_real = sim->m_real;
    for (IdxType i=tid; i<sim->m_gpu; i+=blockDim.x*gridDim.x)
    {
        m_real[i] = sv_real[i]*sv_real[i] + sv_imag[i]*sv_imag[i];
    }
    BARR;

    //Parallel prefix sum
    for (IdxType d=0; d<(sim->n_qubits); d++)
    {
        IdxType step = (IdxType)1<<(d+1);
        for (IdxType k=tid*step; k<n_size; k+=step*blockDim.x*gridDim.x)
        {
            IdxType idx = k+((IdxType)1<<(d+1))-1;
            IdxType local_gpu = idx >> (sim->lg2_m_gpu);
            if (local_gpu == sim->i_gpu)
            {
                IdxType local_idx = idx & ((sim->m_gpu)-1);
                m_real[local_idx] += PGAS_G(m_real, k+((IdxType)1<<d)-1);
            }
        }
        BARR;
    }

    if (sim->i_gpu == (sim->n_gpus-1) && tid == 0) //last GPU
    {
        ValType val = LOCAL_G(m_real, n_size-1);
        LOCAL_P(m_real, n_size-1, 0);
        ValType purity = fabs(val);
        if ( fabs(purity - 1.0) > ERROR_BAR )
            printf("MA: Purity Check fails with %lf\n", purity);
    }
    BARR;

    for (IdxType d=(sim->n_qubits)-1; d>=0; d--)
    {
        IdxType step = (IdxType)1<<(d+1);
        for (IdxType k=tid*step; k<n_size-1; k+=step*blockDim.x*gridDim.x)
        {
            IdxType idx = k+((IdxType)1<<(d+1))-1;
            IdxType local_gpu = idx >> (sim->lg2_m_gpu);
            if (local_gpu == sim->i_gpu)
            {
                IdxType local_idx = idx & ((sim->m_gpu)-1);
                IdxType remote_idx = k+((IdxType)1<<d)-1;
                ValType tmp = PGAS_G(m_real, remote_idx);
                ValType tmp2 = m_real[local_idx];
                PGAS_P(m_real, remote_idx, tmp2);
                m_real[local_idx] = tmp + tmp2;
            }
        }
        BARR;
    }

    for (IdxType j=tid; j<n_size; j+=blockDim.x*gridDim.x)
    {
        IdxType local_gpu = j>>(sim->lg2_m_gpu);
        if (local_gpu == sim->i_gpu)
        {
            ValType lower = LOCAL_G(m_real,j);
            ValType upper = (j+1==n_size)? 1:PGAS_G(m_real,j+1);
            for (IdxType i=0; i<repetition; i++)
            {
                ValType r = sim->randoms_gpu[i];
                if (lower<=r && r<upper) 
                    nvshmem_longlong_p(&sim->results_gpu[i], j, 0);
            }
        }
    }

    BARR;
    
    if (sim->i_gpu != 0 && tid == 0) 
        nvshmem_longlong_get(sim->results_gpu, sim->results_gpu, repetition, 0);
    BARR;
}

__device__ __inline__ void MAV1_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag,
        const IdxType repetition)
{
    grid_group grid = this_grid();
    const IdxType n_size = (IdxType)1<<(sim->n_qubits);
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    ValType * m_real = sim->m_real;
    ValType * m_imag = sim->m_imag;
    for (IdxType i=tid; i<sim->m_gpu; i+=blockDim.x*gridDim.x)
    {
        m_real[i] = sv_real[i]*sv_real[i] + sv_imag[i]*sv_imag[i];
    }
    grid.sync();

    //local parallel prefix sum
    for (IdxType d=0; d<(sim->lg2_m_gpu); d++)
    {
        IdxType step = (IdxType)1<<(d+1);
        for (IdxType k=tid*step; k<sim->m_gpu; k+=step*blockDim.x*gridDim.x)
        {
            m_real[(k+((IdxType)1<<(d+1))-1)] += m_real[k+((IdxType)1<<d)-1];
        }
        grid.sync();
    }
    if (tid == 0)
    {
        m_imag[0] = m_real[sim->m_gpu-1]; //local sum
        m_real[sim->m_gpu-1] = 0;
    }
    grid.sync();
    for (IdxType d=(sim->lg2_m_gpu)-1; d>=0; d--)
    {
        IdxType step = (IdxType)1<<(d+1);
        for (IdxType k=tid*step; k<sim->m_gpu-1; k+=step*blockDim.x*gridDim.x)
        {
            ValType tmp = m_real[k+((IdxType)1<<d)-1];
            ValType tmp2 = m_real[(k+((IdxType)1<<(d+1))-1)];
            m_real[k+((IdxType)1<<d)-1] = tmp2;
            m_real[(k+((IdxType)1<<(d+1))-1)] = tmp + tmp2;
        }
        grid.sync();
    }
    BARR;

    if (sim->i_gpu == 0 && tid == 0) //first GPU
    {
        ValType partial = 0;
        for (IdxType g=0; g<sim->n_gpus; g++)
        {
            nvshmem_double_p(&m_imag[1], partial, g);
            ValType inc = nvshmem_double_g(&m_imag[0], g);
            partial += inc;
        }
        ValType purity = fabs(partial);
        if ( fabs(purity - 1.0) > ERROR_BAR )
            printf("MA: Purity Check fails with %lf\n", purity);
    }

    BARR;
    for (IdxType i=tid; i<sim->m_gpu; i+=blockDim.x*gridDim.x)
    {
        m_real[i] += m_imag[1];
    }

    BARR;
    for (IdxType j=tid; j<n_size; j+=blockDim.x*gridDim.x)
    {
        IdxType local_gpu = j>>(sim->lg2_m_gpu);
        if (local_gpu == sim->i_gpu)
        {
            ValType lower = LOCAL_G(m_real,j);
            ValType upper = (j+1==n_size)? 1:PGAS_G(m_real,j+1);
            for (IdxType i=0; i<repetition; i++)
            {
                ValType r = sim->randoms_gpu[i];
                if (lower<=r && r<upper) 
                    nvshmem_longlong_p(&sim->results_gpu[i], j, 0);
            }
        }
    }
    BARR;
    if (sim->i_gpu != 0 && tid == 0) 
        nvshmem_longlong_get(sim->results_gpu, sim->results_gpu, repetition, 0);
    BARR;
}




__device__ __inline__ void RESET_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType qubit)
{
    grid_group grid = this_grid();
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale));
    ValType * m_real = sim->m_real;
    IdxType mask = ((IdxType)1<<qubit);

    for (IdxType i=tid; i<sim->m_gpu; i+=blockDim.x*gridDim.x)
    {
        IdxType idx = (sim->i_gpu)*per_pe_work+i;
        if ( (idx & mask) == 0) 
            m_real[i] = 0;
        else 
            m_real[i] = sv_real[i]*sv_real[i] + sv_imag[i]*sv_imag[i];
    }
    BARR;

    //Parallel reduction
    for (IdxType k=((IdxType)1<<(sim->n_qubits-1)); k>0; k>>=1)
    {
        for (IdxType i=tid; i<k; i+=blockDim.x*gridDim.x) 
        {
            IdxType local_gpu = i>>(sim->lg2_m_gpu);
            if (local_gpu == sim->i_gpu)
            {
                IdxType local_idx = i&(sim->m_gpu-1);
                m_real[local_idx] += PGAS_G(m_real, i+k);
            }
        }
        BARR;
    }

    if (tid==0 && sim->i_gpu!=0) m_real[0] = PGAS_G(m_real,0);
    grid.sync();
    ValType prob_of_one = m_real[0];
    grid.sync();

    if (prob_of_one < 1.0) //still possible to normalize
    {
        ValType factor = 1.0/sqrt(1.0-prob_of_one);
        for (IdxType i=tid; i<sim->m_gpu; i+=blockDim.x*gridDim.x)
        {
            IdxType idx = (sim->i_gpu)*per_pe_work+i;
            if ( (idx & mask) == 0)
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
        if ((qubit+sim->n_qubits)>=sim->lg2_m_gpu) //remote qubit, need switch
        {
            IdxType pair_gpu = (sim->i_gpu)^((IdxType)1<<(qubit-(sim->lg2_m_gpu)));
            assert(pair_gpu != sim->i_gpu);
            ValType* sv_real_remote = sim->m_real;
            ValType* sv_imag_remote = sim->m_imag;
            if (tid == 0) nvshmem_double_get(sv_real_remote, sv_real, per_pe_work, pair_gpu );
            if (tid == 0) nvshmem_double_get(sv_imag_remote, sv_imag, per_pe_work, pair_gpu );
            BARR;
            for (IdxType i=tid; i<sim->m_gpu; i+=blockDim.x*gridDim.x)
            {
                sv_real[i] = sv_real_remote[i];
                sv_imag[i] = sv_imag_remote[i];
            }
        }
        else
        {
            for (IdxType i=tid; i<sim->m_gpu; i+=blockDim.x*gridDim.x)
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
    BARR;
}

//============== SWAP Gate ================
//This gate is for internal usage. It is used
//when ctrl and target qubits are remote qubis, we then
//swap one of them to a local qubit,
//perform the C2 gate, and then swap back
//It is assumed qubit0 is local, qubit1 is remote
__device__ __inline__ void SWAP_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType qubit0, const IdxType qubit1)
{
    assert (qubit0 != qubit1); //Non-cloning
    grid_group grid = this_grid(); 
    const int tid = blockDim.x * blockIdx.x + threadIdx.x; 

    //For processing a remote qubit, half of the GPUs will be idle for
    //better communication efficiency. Depending on qubit position, 
    //for a GPU pair, we only use GPUs with smaller ids. 
    //Consequently, each GPU should take double of the workload than before
    //Therefore, here it is gpu_scale+1 not gpu_scale+2
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale+1)); 
    const IdxType per_pe_num = ((sim->dim)>>(sim->gpu_scale));
    const IdxType p = min(qubit0, qubit1);
    const IdxType q = max(qubit0, qubit1);

    //load data from pair GPU
    IdxType pair_gpu = (sim->i_gpu)^((IdxType)1<<(q-(sim->lg2_m_gpu)));
    if (sim->i_gpu>pair_gpu) return;

    ValType* sv_real_remote = sim->m_real;
    ValType* sv_imag_remote = sim->m_imag;

    if (tid == 0) nvshmem_double_get(sv_real_remote, sv_real, per_pe_num, pair_gpu );
    if (tid == 0) nvshmem_double_get(sv_imag_remote, sv_imag, per_pe_num, pair_gpu );
    grid.sync();

    for (IdxType i=(sim->i_gpu)*per_pe_work+tid; i<(sim->i_gpu+1)*per_pe_work;
            i+=blockDim.x*gridDim.x) 
    {
        ValType el_real[4];
        ValType el_imag[4];
        ValType res_real[4] = {0};
        ValType res_imag[4] = {0};
        const IdxType term0 = MOD2E(i,p);
        const IdxType term1 = MOD2E(DIV2E(i,p),q-p-1)*EXP2E(p+1);
        const IdxType term2 = DIV2E(DIV2E(i,p),q-p-1)*EXP2E(q+1);
        const IdxType term = term2 + term1 + term0;

        el_real[1] = LOCAL_G(sv_real_remote, term+SV4IDX(1));
        el_imag[1] = LOCAL_G(sv_imag_remote, term+SV4IDX(1));
        el_real[2] = LOCAL_G(sv_real, term+SV4IDX(2));
        el_imag[2] = LOCAL_G(sv_imag, term+SV4IDX(2));

        res_real[1] = el_real[2];
        res_imag[1] = el_imag[2];
        res_real[2] = el_real[1];
        res_imag[2] = el_imag[1];

        LOCAL_P(sv_real_remote, term+SV4IDX(1), res_real[1]);
        LOCAL_P(sv_imag_remote, term+SV4IDX(1), res_imag[1]);
        LOCAL_P(sv_real, term+SV4IDX(2), res_real[2]);
        LOCAL_P(sv_imag, term+SV4IDX(2), res_imag[2]);
    }
    grid.sync();
    if (tid == 0) nvshmem_double_put(sv_real, sv_real_remote, per_pe_num, pair_gpu);
    if (tid == 0) nvshmem_double_put(sv_imag, sv_imag_remote, per_pe_num, pair_gpu);
    //BARR;
}


__device__ __inline__ void Purity_Check(const Simulation* sim, const IdxType t, ValType* sv_real, ValType* sv_imag)
{
    grid_group grid = this_grid();
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale));
    ValType * m_real = sim->m_real;

    for (IdxType i=(sim->i_gpu)*per_pe_work+tid; i<(sim->i_gpu+1)*per_pe_work; i+=blockDim.x*gridDim.x)
    {
        ValType val_real = PGAS_G(sv_real,i);
        ValType val_imag = PGAS_G(sv_imag,i);
        PGAS_P(m_real,i, (val_real*val_real)+(val_imag*val_imag));
    }
    BARR;

    if (sim->i_gpu == 0)
    {
        for (IdxType k=((IdxType)1<<(sim->n_qubits-1)); k>0; k>>=1)
        {
            for (IdxType i=tid; i<k; i+=blockDim.x*gridDim.x) 
            {
                ValType a = PGAS_G(m_real, i);
                ValType b = PGAS_G(m_real, i+k);
                PGAS_P(m_real,i,a+b);
            }
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
    }
    BARR;
}


//=====================================================================================
//Per-gate execution function
__device__ void Gate::exe_op(Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    grid_group grid = this_grid(); 
    //only need sync when operating on remote qubits
    if ((ctrl>=sim->lg2_m_gpu) || (qubit>=sim->lg2_m_gpu))
    {
        grid.sync();
        if( threadIdx.x==0 && blockIdx.x==0 ) nvshmem_barrier_all(); 
    }
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
        //MAV1_GATE(sim, sv_real, sv_imag, ctrl); 
    }
    else if (op_name == C1) 
    {
        C1V2_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, qubit);
    }
    else if (op_name == C2) 
    {
        if ((ctrl>=sim->lg2_m_gpu) && (qubit>=sim->lg2_m_gpu))
        {
            SWAP_GATE(sim, sv_real, sv_imag, 0, ctrl);
            BARR;
            C2V1_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, 0, qubit);
            BARR;
            SWAP_GATE(sim, sv_real, sv_imag, 0, ctrl);
        }
        else
        {
            C2V1_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, ctrl, qubit);
        }
    }
    //only need sync when operating on remote qubits
    if ((ctrl>=sim->lg2_m_gpu) || (qubit>=sim->lg2_m_gpu))
        if( threadIdx.x==0 && blockIdx.x==0 ) nvshmem_barrier_all(); 
    grid.sync();
    //BARR;
}

}; //namespace NWQSim

#endif

