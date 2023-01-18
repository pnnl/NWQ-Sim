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
// Single GPU SV-Sim simulation runtime using NVIDIA GPU backend.
// ---------------------------------------------------------------------------

#ifndef SVSIM_NVGPU_SIN_CUH
#define SVSIM_NVGPU_SIN_CUH
#include <assert.h>
#include <random>
#include <complex.h>
#include <cooperative_groups.h>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <cuda.h>
#include <mma.h>
#include "config.h"
#include "gate.h"

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
        n_qubits(_n_qubits),
        dim((IdxType)1<<(n_qubits)), 
        half_dim((IdxType)1<<(n_qubits-1)),
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
        sv_size = dim*(IdxType)sizeof(ValType);
        i_gpu = DEFAULT_SIN_GPU;
        //always be 0 since 1-MPI maps to 1-GPU
        cudaSafeCall(cudaSetDevice(i_gpu));
        //CPU side initialization
        SAFE_ALOC_HOST(sv_real_cpu, sv_size);
        SAFE_ALOC_HOST(sv_imag_cpu, sv_size);
        memset(sv_real_cpu, 0, sv_size);
        memset(sv_imag_cpu, 0, sv_size);
        //State-vector initial state [0..0] = 1
        sv_real_cpu[0] = 1.;
        //NVSHMEM GPU memory allocation
        SAFE_ALOC_GPU(sv_real, sv_size);
        SAFE_ALOC_GPU(sv_imag, sv_size);
        SAFE_ALOC_GPU(m_real, sv_size+sizeof(ValType));
        SAFE_ALOC_GPU(m_imag, sv_size+sizeof(ValType));
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
        cudaSafeCall(cudaMemset(m_real, 0, sv_size+sizeof(ValType)));
        cudaSafeCall(cudaMemset(m_imag, 0, sv_size+sizeof(ValType)));
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
        //Reset CPU input & output
        memset(sv_real_cpu, 0, sv_size);
        memset(sv_imag_cpu, 0, sv_size);
        //State Vector initial state [0..0] = 1
        sv_real_cpu[0] = 1.;
        //GPU side initialization
        cudaSafeCall(cudaMemcpy(sv_real, sv_real_cpu, 
                    sv_size, cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(sv_imag, sv_imag_cpu, 
                    sv_size, cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemset(m_real, 0, sv_size+sizeof(ValType)));
        cudaSafeCall(cudaMemset(m_imag, 0, sv_size+sizeof(ValType)));
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
        printf("SVSim_gpu is running! Requesting %lld qubits.\n", circuit_handle->n_qubits);
#endif
#ifdef PRINT_KERNEL_TRACE
        double sim_time;
        gpu_timer sim_timer;
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
        printf("\n============== SV-Sim ===============\n");
        printf("nqubits:%lld, ngates:%lld, sim_gates:%lld, ngpus:%d, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms, mem:%.3lf MB, mem_per_gpu:%.3lf MB\n",
                n_qubits, input_gates, n_gates, 1, sim_time, 0., 
                sim_time, gpu_mem/1024/1024, gpu_mem/1024/1024);
        printf("=====================================\n");
        fflush(stdout);
#endif
        clear_circuit();
    }
    IdxType measure(IdxType qubit) 
    {
        this->M(qubit);
        this->sim();
        cudaSafeCall(cudaMemcpy(results, results_gpu, sizeof(IdxType), cudaMemcpyDeviceToHost));
        //cout << "Measure result: " << this->results[0] << endl;

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
        cudaSafeCall(cudaMemcpy(sv_real_cpu, sv_real, sv_size, cudaMemcpyDeviceToHost));
        cudaSafeCall(cudaMemcpy(sv_imag_cpu, sv_imag, sv_size, cudaMemcpyDeviceToHost));

        IdxType num = ((IdxType)1<<n_qubits);
        printf("----- SVSim ------\n");
        for (IdxType i=0; i<num; i++) 
        {
            printf("(%.3lf,%.3lfj) ", sv_real_cpu[i], sv_imag_cpu[i]);
            if ((i+1)%8==0) printf("\n");
        }
        printf("\n");
    }
public:
    // n_qubits is the number of qubits
    IdxType n_qubits;
    // which gpu
    IdxType i_gpu;
    IdxType sv_size;

    IdxType dim;
    IdxType half_dim;
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
            const ValType val = LOCAL_G(sv_real, i);
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
        //if (blockIdx.x == 0 && threadIdx.x == 0) t0 = clock64();

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


//============== Local Unified 1-qubit Gate ================
__device__ __inline__ void C1_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag, const IdxType qubit)
{
    grid_group grid = this_grid(); 
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x; 
    for (IdxType i=tid; i<sim->half_dim; i+=blockDim.x*gridDim.x)
    { 
        IdxType outer = (i >> qubit); 
        IdxType inner =  (i & (((IdxType)1<<qubit)-1)); 
        IdxType offset = (outer << (qubit+1)); 
        IdxType pos0 = (offset + inner);
        IdxType pos1 = (offset + inner + ((IdxType)1<<qubit));
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

//============== Local 2-qubit Gate  ================
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
    //BARR;
}

__device__ __inline__ void M_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        ValType* gm_real, ValType* gm_imag, const IdxType qubit, const ValType rand)
{
    grid_group grid = this_grid();
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    ValType * m_real = sim->m_real;
    IdxType mask = ((IdxType)1<<qubit);

    for (IdxType i=tid; i<sim->dim; i+=blockDim.x*gridDim.x)
    {
        if ( (i & mask) == 0) m_real[i] = 0;
        else m_real[i] = sv_real[i]*sv_real[i] + sv_imag[i]*sv_imag[i];
    }
    BARR;

    //Parallel reduction
    for (IdxType k=((IdxType)1<<(sim->n_qubits-1)); k>0; k>>=1)
    {
        for (IdxType i=tid; i<k; i+=blockDim.x*gridDim.x) 
        {
            m_real[i] += m_real[i+k];
        }
        BARR;
    }
    ValType prob_of_one = m_real[0];
    grid.sync();
 
    if (rand < prob_of_one)
    {
        ValType factor = 1./sqrt(prob_of_one);
        for (IdxType i=tid; i<sim->dim; i+=blockDim.x*gridDim.x)
        {
            if ( (i & mask) == 0)
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
        for (IdxType i=tid; i<sim->dim; i+=blockDim.x*gridDim.x)
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
    
    for (IdxType i=tid; i<sim->dim; i+=blockDim.x*gridDim.x)
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
            m_real[(k+((IdxType)1<<(d+1))-1)] += LOCAL_G(m_real, k+((IdxType)1<<d)-1);
        }
        BARR;
    }

    if (tid == 0)
    {
        ValType val = LOCAL_G(m_real, n_size-1);
        m_real[n_size] = val;
        LOCAL_P(m_real, n_size-1, 0);
        ValType purity = fabs(val);
        if ( abs(purity - 1.0) > ERROR_BAR )
            printf("MA: Purity Check fails with %lf\n", purity);

    }

    BARR;

    for (IdxType d=(sim->n_qubits)-1; d>=0; d--)
    {
        IdxType step = (IdxType)1<<(d+1);
        for (IdxType k=tid*step; k<n_size-1; k+=step*blockDim.x*gridDim.x)
        {
            ValType tmp = LOCAL_G(m_real, k+((IdxType)1<<d)-1);
            ValType tmp2 = LOCAL_G(m_real, (k+((IdxType)1<<(d+1))-1));
            LOCAL_P(m_real, k+((IdxType)1<<d)-1, tmp2);
            m_real[(k+((IdxType)1<<(d+1))-1)] = tmp + tmp2;
        }
        BARR;
    }

    for (IdxType j=tid; j<n_size; j+=blockDim.x*gridDim.x)
    {
        ValType lower = LOCAL_G(m_real,j);
        ValType upper = (j+1==n_size)? 1:LOCAL_G(m_real,j+1);
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
    ValType * m_real = sim->m_real;
    IdxType mask = ((IdxType)1<<qubit);

    for (IdxType i=tid; i<sim->dim; i+=blockDim.x*gridDim.x)
    {
        if ( (i & mask) == 0) m_real[i] = 0;
        else m_real[i] = sv_real[i]*sv_real[i] + sv_imag[i]*sv_imag[i];
    }
    BARR;
    for (IdxType k=((IdxType)1<<(sim->n_qubits-1)); k>0; k>>=1)
    {
        for (IdxType i=tid; i<k; i+=blockDim.x*gridDim.x) 
        {
            m_real[i] += m_real[i+k];
        }
        BARR;
    }
    
    BARR;
    ValType prob_of_one = m_real[0];
    grid.sync();

    if (prob_of_one < 1.0) //still possible to normalize
    {
        ValType factor = 1.0/sqrt(1.0-prob_of_one);
        for (IdxType i=tid; i<sim->dim; i+=blockDim.x*gridDim.x)
        {
            if ( (i & mask) == 0)
            {
                sv_real[i] *= factor;
                sv_imag[i] *= factor;
                m_real[i] = sv_real[i]*sv_real[i] + sv_imag[i]*sv_imag[i];
            }
            else
            {
                sv_real[i] = 0;
                sv_imag[i] = 0;
                m_real[i] = 0;
            }
        }
        BARR;
        for (IdxType k=((IdxType)1<<(sim->n_qubits-1)); k>0; k>>=1)
        {
            for (IdxType i=tid; i<k; i+=blockDim.x*gridDim.x) 
            {
                m_real[i] += m_real[i+k];
            }
            BARR;
        }
        ValType norm = m_real[0];
        grid.sync();
        if ( abs(norm-1.0) > 0)
        {
            ValType factor = 1.0/sqrt(norm);
            for (IdxType i=tid; i<sim->dim; i+=blockDim.x*gridDim.x)
            {
                sv_real[i] *= factor;
                sv_imag[i] *= factor;
            }
        }
    }
    //becuase qubit=0 probability is 0, we can't simply normalize
    else 
    {
        for (IdxType i=tid; i<sim->dim; i+=blockDim.x*gridDim.x)
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
    const IdxType per_pe_work = (sim->dim);
    ValType * m_real = sim->m_real;

    for (IdxType i=tid; i<sim->dim; i+=blockDim.x*gridDim.x)
    {
        m_real[i] = sv_real[i]*sv_real[i] + sv_imag[i]*sv_imag[i];
    }
    BARR;

    for (IdxType k=((IdxType)1<<(sim->n_qubits-1)); k>0; k>>=1)
    {
        for (IdxType i=tid; i<k; i+=blockDim.x*gridDim.x) 
        {
            m_real[i] += m_real[i+k];
        }
        grid.sync();
    }
    BARR;
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
    else if (op_name == C1) 
    {
        C1_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, qubit);
    }
    else if (op_name == C2) 
    {
        C2_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, ctrl, qubit);
    }
    grid.sync();
}

}; //namespace NWQSim

#endif

