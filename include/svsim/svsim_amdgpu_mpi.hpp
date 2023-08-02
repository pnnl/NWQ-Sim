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
// File: svsim_amdgpu_mpi.cuh
// MPI based implementation of the scale-out SV-Sim gates and 
// simulation runtime using AMD GPU backend.
// ---------------------------------------------------------------------------

#ifndef SVSIM_AMDGPU_MPI_HPP
#define SVSIM_AMDGPU_MPI_HPP
#include <hip/hip_runtime.h>
#include <assert.h>
#include <random>
#include <complex.h>
#include <hip/hip_cooperative_groups.h>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include <string>
#include <iostream>
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

class Gate;
class Simulation;
void Purity_Check(const Simulation* sim, const IdxType t);

//Simulation runtime
void simulation_kernel(Simulation*);

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
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) std::cout << g.gateToString() << std::flush;
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
        hipSafeCall(hipMemcpy(circuit_gpu, circuit.data(), n_gates*sizeof(Gate), hipMemcpyHostToDevice));
        fused_circuit.assign(circuit.begin(), circuit.end());
        //====================================================
#else
        //====================== Fuse ========================
        vector<Gate> tmp1_circuit;
        vector<Gate> tmp2_circuit;
        vector<Gate> tmp3_circuit;
        tmp1_circuit.clear();
        tmp2_circuit.clear();
        fused_circuit.clear();

        gate_fusion_1q(circuit, tmp1_circuit, n_qubits);
        gate_fusion_2q_absorb_1q_forward(tmp1_circuit, tmp2_circuit, n_qubits);
        gate_fusion_2q_absorb_1q_backward(tmp2_circuit, tmp3_circuit, n_qubits);
        gate_fusion_2q(tmp3_circuit, fused_circuit, n_qubits);

        this->n_gates = fused_circuit.size();
        SAFE_FREE_GPU(circuit_gpu);
        SAFE_ALOC_GPU(circuit_gpu, n_gates*sizeof(Gate));
        hipSafeCall(hipMemcpy(circuit_gpu, fused_circuit.data(), n_gates*sizeof(Gate), hipMemcpyHostToDevice));
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
        results(NULL),
        results_local(NULL)
    {
        int mpi_size;
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
        n_gpus = (IdxType)mpi_size;
        int mpi_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        i_gpu = (IdxType)mpi_rank;

        //always be 0 since 1-MPI maps to 1-GPU
        hipSafeCall(hipSetDevice(0));
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
        //GPU memory allocation
        SAFE_ALOC_GPU(sv_real, sv_size_per_gpu);
        SAFE_ALOC_GPU(sv_imag, sv_size_per_gpu);
        SAFE_ALOC_GPU(m_real, sv_size_per_gpu);
        SAFE_ALOC_GPU(m_imag, sv_size_per_gpu);
        hipCheckError();
        gpu_mem += sv_size_per_gpu*4;
        //Initialize Circuit 
        circuit_handle = new Circuit(n_qubits);
        circuit_handle->n_gpus = n_gpus;
        circuit_handle_gpu = NULL;
        //GPU memory initilization
        hipSafeCall(hipMemcpy(sv_real, sv_real_cpu, sv_size_per_gpu, 
                    hipMemcpyHostToDevice));
        hipSafeCall(hipMemcpy(sv_imag, sv_imag_cpu, sv_size_per_gpu, 
                    hipMemcpyHostToDevice));
        hipSafeCall(hipMemset(m_real, 0, sv_size_per_gpu));
        hipSafeCall(hipMemset(m_imag, 0, sv_size_per_gpu));
        SAFE_ALOC_GPU(sim_gpu, sizeof(Simulation));
        rng.seed(time(0));
    }

    ~Simulation()
    {
        //Release circuit
        if (circuit_handle != NULL)
            delete circuit_handle;
        //Release for GPU side
        SAFE_FREE_GPU(sv_real);
        SAFE_FREE_GPU(sv_imag);
        SAFE_FREE_GPU(m_real);
        SAFE_FREE_GPU(m_imag);
        SAFE_FREE_GPU(sim_gpu);

        //Release for CPU side
        SAFE_FREE_HOST(sv_real_cpu);
        SAFE_FREE_HOST(sv_imag_cpu);
        SAFE_FREE_HOST(randoms);
        SAFE_FREE_HOST(results);
        SAFE_FREE_HOST(results_local);
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
        ValType rand = uni_dist(rng);
        Gate* G = new Gate(OP::M,qubit,-1,rand);
        circuit_handle->append(*G);
    }
    void MA(IdxType repetition) //default is pauli-Z
    {
        SAFE_FREE_HOST(results);
        SAFE_ALOC_HOST(results, sizeof(IdxType)*repetition);
        memset(results, 0, sizeof(IdxType)*repetition);
        SAFE_FREE_HOST(results_local);
        SAFE_ALOC_HOST(results_local, sizeof(IdxType)*repetition);
        memset(results_local, 0, sizeof(IdxType)*repetition);
        SAFE_FREE_HOST(randoms);
        SAFE_ALOC_HOST(randoms, sizeof(ValType)*repetition);
        for (IdxType i=0; i<repetition; i++) 
            randoms[i] = uni_dist(rng);
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
        memset(sv_real_cpu, 0, sv_size_per_gpu);
        memset(sv_imag_cpu, 0, sv_size_per_gpu);
        //State Vector initial state [0..0] = 1
        if (i_gpu == 0) sv_real_cpu[0] = 1.;
        //GPU side initialization
        hipSafeCall(hipMemcpy(sv_real, sv_real_cpu, 
                    sv_size_per_gpu, hipMemcpyHostToDevice));
        hipSafeCall(hipMemcpy(sv_imag, sv_imag_cpu, 
                    sv_size_per_gpu, hipMemcpyHostToDevice));
        hipSafeCall(hipMemset(m_real, 0, sv_size_per_gpu));
        hipSafeCall(hipMemset(m_imag, 0, sv_size_per_gpu));
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
    }
    std::string circuitToString()
    {
        return circuit_handle->circuitToString();
    }
    void sim()
    {
        hipSafeCall(hipSetDevice(0));
#ifdef PRINT_KERNEL_TRACE
        double* sim_times;
        double sim_time;
        cpu_timer sim_timer;
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
        hipSafeCall(hipMemcpy(sim_gpu, this, 
                    sizeof(Simulation), hipMemcpyHostToDevice));
#ifdef PRINT_SIM_TRACE
        printf("SVSim_GPU is running! Requesting %lld qubits.\n", circuit_handle->n_qubits);
#endif
#ifdef PRINT_KERNEL_TRACE
        MPI_Barrier(MPI_COMM_WORLD);
        sim_timer.start_timer();
#endif
        //=========================================
        simulation_kernel(this);
        //=========================================

#ifdef PRINT_KERNEL_TRACE
        sim_timer.stop_timer();
        sim_time = sim_timer.measure();
        
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
        hipSafeCall(hipMemcpy(sv_real_cpu, sv_real, sv_size_per_gpu, hipMemcpyDeviceToHost));
        hipSafeCall(hipMemcpy(sv_imag_cpu, sv_imag, sv_size_per_gpu, hipMemcpyDeviceToHost));

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
    //GPU Buffer
    ValType* m_real;
    ValType* m_imag;
    //For measurement randoms
    ValType* randoms;
    //For measurement result
    IdxType* results;
    //Local measurement result for MA
    IdxType* results_local;
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

#define LAUNCH_KERNEL(X) dim3 gridDim(1,1,1); \
    hipDeviceProp_t deviceProp; \
    unsigned smem_size = 0; \
    hipSafeCall(hipGetDeviceProperties(&deviceProp, 0)); \
    int numBlocksPerSm; \
    hipSafeCall(hipOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm, \
        (X), THREADS_CTA_AMDGPU, smem_size)); \
    gridDim.x = numBlocksPerSm * deviceProp.multiProcessorCount; \
    hipSafeCall(hipLaunchCooperativeKernel((void*)(X),gridDim, \
                        THREADS_CTA_AMDGPU,args,smem_size,0)); \
    hipCheckError(); \
    hipSafeCall(hipDeviceSynchronize()); 


//================================= Gate Definition ========================================

#define LOCAL_G(arr,i) arr[(i)&(sim->m_gpu-1)]
#define LOCAL_P(arr,i,val) arr[(i)&(sim->m_gpu-1)] = val;

//============== Local Unified 1-qubit Gate ================
__global__ void C1LC_GATE(Simulation* sim, const IdxType t)
{
    ValType* sv_real = sim->sv_real;
    ValType* sv_imag = sim->sv_imag;
    const ValType* gm_real = sim->circuit_handle_gpu[t].gm_real;
    const ValType* gm_imag = sim->circuit_handle_gpu[t].gm_imag;
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x; 
    const IdxType per_pe_work = ((sim->half_dim)>>(sim->gpu_scale));
    IdxType qubit = sim->circuit_handle_gpu[t].qubit;
    
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
}

__global__ void C1RM_GATE(Simulation* sim, const IdxType t)
{
    ValType* sv_real = sim->sv_real;
    ValType* sv_imag = sim->sv_imag;
    ValType* sv_real_remote = sim->m_real;
    ValType* sv_imag_remote = sim->m_imag;
    const ValType* gm_real = sim->circuit_handle_gpu[t].gm_real;
    const ValType* gm_imag = sim->circuit_handle_gpu[t].gm_imag;
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x; 
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale));
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
}
 
void C1_GATE(const Simulation* sim, const IdxType t)
{
    Gate& g = (sim->circuit_handle->fused_circuit[t]);
    ValType* sv_real = sim->sv_real;
    ValType* sv_imag = sim->sv_imag;
    ValType* m_real= sim->m_real;
    ValType* m_imag= sim->m_imag;

    if (g.qubit < sim->lg2_m_gpu) 
    {
        void* args[] = {(void*)(&sim->sim_gpu), (void*)&t};
        LAUNCH_KERNEL(C1LC_GATE);
    }
    else
    {
        const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale));
        IdxType pair_gpu = (sim->i_gpu)^((IdxType)1<<(g.qubit-(sim->lg2_m_gpu)));
        MPI_Barrier(MPI_COMM_WORLD);

        if (sim->i_gpu>pair_gpu) 
        {
            MPI_Send(sv_real, per_pe_work, MPI_DOUBLE, pair_gpu, 0, MPI_COMM_WORLD);
            MPI_Send(sv_imag, per_pe_work, MPI_DOUBLE, pair_gpu, 1, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Recv(m_real, per_pe_work, MPI_DOUBLE, pair_gpu, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
            MPI_Recv(m_imag, per_pe_work, MPI_DOUBLE, pair_gpu, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
        }
        MPI_Barrier(MPI_COMM_WORLD);

        if (sim->i_gpu<pair_gpu) 
        {
            void* args[] = {(void*)(&sim->sim_gpu), (void*)&t};
            LAUNCH_KERNEL(C1RM_GATE);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        
        if (sim->i_gpu>pair_gpu) 
        {
            MPI_Recv(sv_real, per_pe_work, MPI_DOUBLE, pair_gpu, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(sv_imag, per_pe_work, MPI_DOUBLE, pair_gpu, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        else
        {
            MPI_Send(sim->m_real, per_pe_work, MPI_DOUBLE, pair_gpu, 2, MPI_COMM_WORLD);
            MPI_Send(sim->m_imag, per_pe_work, MPI_DOUBLE, pair_gpu, 3, MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

//============== Local 2-qubit Gate  ================
__global__ void C2LC_GATE(const Simulation* sim, const IdxType qubit0, 
        const IdxType qubit1, const IdxType t) 
{
    ValType* sv_real = sim->sv_real;
    ValType* sv_imag = sim->sv_imag;
    const ValType* gm_real = sim->circuit_handle_gpu[t].gm_real;
    const ValType* gm_imag = sim->circuit_handle_gpu[t].gm_imag;
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
}


#define SV4IDX(x) (((x>>1)&1)*EXP2E(qubit0) + ((x&1)*EXP2E(qubit1)) )

#define DIV2E(x,y) ((x)>>(y))
#define MOD2E(x,y) ((x)&(((IdxType)1<<(y))-(IdxType)1)) 
#define EXP2E(x) ((IdxType)1<<(x))
#define SV16IDX(x) ( ((x>>3)&1)*EXP2E(qubit0) + ((x>>2)&1)*EXP2E(qubit1) + ((x>>1)&1)*EXP2E(qubit2) + ((x&1)*EXP2E(qubit3)) )

//============== Unified 2-qubit Gate ================
//Perform communication optimization here
__global__ void C2RM_GATE(const Simulation* sim, const IdxType qubit0, 
        const IdxType qubit1, const IdxType t) 
{
    ValType* sv_real = sim->sv_real;
    ValType* sv_imag = sim->sv_imag;
    ValType* sv_real_remote = sim->m_real;
    ValType* sv_imag_remote = sim->m_imag;
    const ValType* gm_real = sim->circuit_handle_gpu[t].gm_real;
    const ValType* gm_imag = sim->circuit_handle_gpu[t].gm_imag;
    const int tid = blockDim.x * blockIdx.x + threadIdx.x; 
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale+1)); 
    const IdxType per_pe_num = ((sim->dim)>>(sim->gpu_scale));
    const IdxType p = min(qubit0, qubit1);
    const IdxType q = max(qubit0, qubit1);

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
}

void C2RMM_GATE(const Simulation* sim, const IdxType qubit0, const IdxType qubit1, const IdxType t)
{
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale)); 
    const IdxType p = min(qubit0, qubit1);
    const IdxType q = max(qubit0, qubit1);
    //load data from pair GPU
    IdxType pair_gpu = (sim->i_gpu)^((IdxType)1<<(q-(sim->lg2_m_gpu)));
    MPI_Barrier(MPI_COMM_WORLD);

    if (sim->i_gpu>pair_gpu) 
    {
        MPI_Send(sim->sv_real, per_pe_work, MPI_DOUBLE, pair_gpu, 0, MPI_COMM_WORLD);
        MPI_Send(sim->sv_imag, per_pe_work, MPI_DOUBLE, pair_gpu, 1, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Recv(sim->m_real, per_pe_work, MPI_DOUBLE, pair_gpu, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
        MPI_Recv(sim->m_imag, per_pe_work, MPI_DOUBLE, pair_gpu, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (sim->i_gpu<pair_gpu) 
    {
        void* args[] = {(void*)(&sim->sim_gpu), (void*)&qubit0, (void*)&qubit1, (void*)&t};
        LAUNCH_KERNEL(C2RM_GATE);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (sim->i_gpu>pair_gpu) 
    {
        MPI_Recv(sim->sv_real, per_pe_work, MPI_DOUBLE, pair_gpu, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(sim->sv_imag, per_pe_work, MPI_DOUBLE, pair_gpu, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
        MPI_Send(sim->m_real, per_pe_work, MPI_DOUBLE, pair_gpu, 2, MPI_COMM_WORLD);
        MPI_Send(sim->m_imag, per_pe_work, MPI_DOUBLE, pair_gpu, 3, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

//============== SWAP Gate ================
//This gate is for internal usage. It is used
//when ctrl and target qubits are remote qubis, we then
//swap one of them to a local qubit,
//perform the C2 gate, and then swap back
//It is assumed qubit0 is local, qubit1 is remote

__global__ void SWAPLC_GATE(const Simulation* sim, const IdxType qubit0, const IdxType qubit1) 
{
    assert (qubit0 != qubit1); //Non-cloning
    const int tid = blockDim.x * blockIdx.x + threadIdx.x; 
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale+1)); 
    ValType* sv_real = sim->sv_real;
    ValType* sv_imag = sim->sv_imag;
    ValType* sv_real_remote = sim->m_real;
    ValType* sv_imag_remote = sim->m_imag;
    const IdxType p = min(qubit0, qubit1);
    const IdxType q = max(qubit0, qubit1);

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

}




void SWAP_GATE(const Simulation* sim, const IdxType qubit0, const IdxType qubit1)
{
    assert (qubit0 != qubit1); //Non-cloning
    const IdxType per_pe_num = ((sim->dim)>>(sim->gpu_scale));
    const IdxType p = min(qubit0, qubit1);
    const IdxType q = max(qubit0, qubit1);

    //load data from pair GPU
    IdxType pair_gpu = (sim->i_gpu)^((IdxType)1<<(q-(sim->lg2_m_gpu)));
    MPI_Barrier(MPI_COMM_WORLD);

    if (sim->i_gpu>pair_gpu) 
    {
        MPI_Send(sim->sv_real, per_pe_num, MPI_DOUBLE, pair_gpu, 0, MPI_COMM_WORLD);
        MPI_Send(sim->sv_imag, per_pe_num, MPI_DOUBLE, pair_gpu, 1, MPI_COMM_WORLD);
    }
    else
    {
        MPI_Recv(sim->m_real, per_pe_num, MPI_DOUBLE, pair_gpu, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
        MPI_Recv(sim->m_imag, per_pe_num, MPI_DOUBLE, pair_gpu, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (sim->i_gpu<pair_gpu) 
    {
        void* args[] = {(void*)(&sim->sim_gpu), (void*)&qubit0, (void*)&qubit1};
        LAUNCH_KERNEL(SWAPLC_GATE);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (sim->i_gpu>pair_gpu) 
    {
        MPI_Recv(sim->sv_real, per_pe_num, MPI_DOUBLE, pair_gpu, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(sim->sv_imag, per_pe_num, MPI_DOUBLE, pair_gpu, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
        MPI_Send(sim->m_real, per_pe_num, MPI_DOUBLE, pair_gpu, 2, MPI_COMM_WORLD);
        MPI_Send(sim->m_imag, per_pe_num, MPI_DOUBLE, pair_gpu, 3, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}




void C2_GATE(const Simulation* sim, const IdxType t)
{
    IdxType ctrl = sim->circuit_handle_gpu[t].ctrl;
    IdxType qubit = sim->circuit_handle_gpu[t].qubit;
    assert (ctrl != qubit); //Non-cloning
    //Case-1: both ctrl and qubit are local
    
    if (ctrl < sim->lg2_m_gpu && qubit < sim->lg2_m_gpu) 
    {
        void* args[] = {(void*)(&sim->sim_gpu), &ctrl, &qubit, (void*)&t};
        LAUNCH_KERNEL(C2LC_GATE);
    }
    //Case-2: both ctrl and qubit are global
    else if (ctrl >= sim->lg2_m_gpu && qubit >= sim->lg2_m_gpu) 
    {
        SWAP_GATE(sim, 0, ctrl);
        MPI_Barrier(MPI_COMM_WORLD);
        C2RMM_GATE(sim, 0, qubit, t);
        MPI_Barrier(MPI_COMM_WORLD);
        SWAP_GATE(sim, 0, ctrl);
    }
    //Case-3: one of ctrl and qubit is global 
    else
    {
        C2RMM_GATE(sim, ctrl, qubit, t);
    }
}

__global__ void MLC1_GATE(const Simulation* sim, const IdxType t)
{
    grid_group grid = this_grid();
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    const IdxType qubit = sim->circuit_handle_gpu[t].qubit;
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale));
    IdxType mask = ((IdxType)1<<qubit);
    ValType* sv_real = sim->sv_real;
    ValType* sv_imag = sim->sv_imag;
    ValType * m_real = sim->m_real;
    for (IdxType i=tid; i<sim->m_gpu; i+=blockDim.x*gridDim.x)
    {
        IdxType idx = (sim->i_gpu)*per_pe_work+i;
        if ( (idx & mask) == 0) 
            m_real[i] = 0;
        else 
            m_real[i] = sv_real[i]*sv_real[i] + sv_imag[i]*sv_imag[i];
    }
    grid.sync();
    //Parallel reduction
    for (IdxType k=(sim->m_gpu>>(IdxType)1); k>0; k>>=1)
    {
        for (IdxType i=tid; i<k; i+=blockDim.x*gridDim.x) 
        {
            m_real[i] += m_real[i+k];
        }
        grid.sync();
    }
}

__global__ void MLC2_GATE(const Simulation* sim, const ValType prob_of_one, const IdxType t)
{
    grid_group grid = this_grid();
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    const IdxType qubit = sim->circuit_handle_gpu[t].qubit;
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale));
    ValType* sv_real = sim->sv_real;
    ValType* sv_imag = sim->sv_imag;
    ValType * m_real = sim->m_real;
    IdxType mask = ((IdxType)1<<qubit);
    const ValType rand = sim->circuit_handle_gpu[t].theta;
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
}

void M_GATE(const Simulation* sim, const IdxType t)
{
    {
        void* args[] = {(void*)(&sim->sim_gpu), (void*)&t};
        LAUNCH_KERNEL(MLC1_GATE);
    }
    ValType sum = 0;
    ValType prob_of_one = 0;
    hipSafeCall(hipMemcpy(&sum, &(sim->m_real[0]), sizeof(ValType), hipMemcpyHostToDevice));
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&sum, &prob_of_one, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    {
        void* args[] = {(void*)(&sim->sim_gpu), &prob_of_one, (void*)&t};
        LAUNCH_KERNEL(MLC2_GATE);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

void MA_GATE(const Simulation* sim, const IdxType t) 
{
    const IdxType sv_size_per_gpu = sim->sv_size_per_gpu;
    const IdxType repetition = sim->circuit_handle->fused_circuit[t].ctrl;
    ValType* sv_real_cpu = sim->sv_real_cpu;
    ValType* sv_imag_cpu = sim->sv_imag_cpu;
    hipSafeCall(hipMemcpy(sv_real_cpu, sim->sv_real, sv_size_per_gpu, hipMemcpyDeviceToHost));
    hipSafeCall(hipMemcpy(sv_imag_cpu, sim->sv_imag, sv_size_per_gpu, hipMemcpyDeviceToHost));

    ValType* m_buff = NULL;
    SAFE_ALOC_HOST(m_buff, sv_size_per_gpu);
    memset(m_buff, 0, sv_size_per_gpu);
    const IdxType n_size = (IdxType)1<<(sim->n_qubits);

    ValType reduce = 0;
    ValType partial = 0;
    for (IdxType i=0; i<sim->m_gpu; i++)
        reduce += sv_real_cpu[i]*sv_real_cpu[i] + sv_imag_cpu[i]*sv_imag_cpu[i];
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Scan(&reduce, &partial, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if (sim->i_gpu == (sim->n_gpus-1)) //last node
    {
        ValType purity = fabs(partial);
        if ( fabs(purity - 1.0) > ERROR_BAR )
            printf("MA: Purity Check fails with %lf\n", purity);
    }

    m_buff[0] = (partial - reduce); //per-node incremental val
    for (IdxType i=1; i<sim->m_gpu; i++)
        m_buff[i] = m_buff[i-1] + ((sv_real_cpu[i-1]*sv_real_cpu[i-1]) + (sv_imag_cpu[i-1]*sv_imag_cpu[i-1])); 

    for (IdxType j=0; j<n_size; j++)
    {
        IdxType local_cpu = j>>(sim->lg2_m_gpu);
        IdxType local_j = j & (sim->m_gpu-1);
        if (local_cpu == sim->i_gpu)
        {
            ValType lower = LOCAL_G(m_buff,j);
            ValType upper = 0;
            if (j+1 == n_size) upper = 1.0; //last element
            else upper = (local_j+1==sim->m_gpu)?partial:m_buff[local_j+1]; //last element per node
            
            for (IdxType i=0; i<repetition; i++)
            {
                ValType r = sim->randoms[i];
                //all nodes store partial results locally. since other entires are all
                //zeros, we can do all-reduce to merge
                if (lower<=r && r<upper) sim->results_local[i] = j;
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(sim->results_local, sim->results, repetition, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    SAFE_FREE_HOST(m_buff);
}

__global__ void RESETLC1_GATE(const Simulation* sim, const IdxType t)
{
    grid_group grid = this_grid();
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    const IdxType qubit = sim->circuit_handle_gpu[t].qubit;
    IdxType mask = ((IdxType)1<<qubit);
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale));
    ValType* sv_real = sim->sv_real;
    ValType* sv_imag = sim->sv_imag;
    ValType * m_real = sim->m_real;
    for (IdxType i=tid; i<sim->m_gpu; i+=blockDim.x*gridDim.x)
    {
        IdxType idx = (sim->i_gpu)*per_pe_work+i;
        if ( (idx & mask) == 0) 
            m_real[i] = 0;
        else 
            m_real[i] = sv_real[i]*sv_real[i] + sv_imag[i]*sv_imag[i];
    }
    grid.sync();
    //Parallel reduction
    for (IdxType k=(sim->m_gpu>>(IdxType)1); k>0; k>>=1)
    {
        for (IdxType i=tid; i<k; i+=blockDim.x*gridDim.x) 
        {
            m_real[i] += m_real[i+k];
        }
        grid.sync();
    }
}

__global__ void RESETLC2_GATE(const Simulation* sim, const ValType prob_of_one, const IdxType t)
{
    const IdxType qubit = sim->circuit_handle_gpu[t].qubit;
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale));
    ValType* sv_real = sim->sv_real;
    ValType* sv_imag = sim->sv_imag;
    ValType factor = 1.0/sqrt(1.0-prob_of_one);
    IdxType mask = ((IdxType)1<<qubit);
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

__global__ void RESETLC3_GATE(const Simulation* sim, const IdxType t)
{
    const IdxType qubit = sim->circuit_handle_gpu[t].qubit;
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale));
    ValType* sv_real = sim->sv_real;
    ValType* sv_imag = sim->sv_imag;
    IdxType mask = ((IdxType)1<<qubit);
    
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

void RESET_GATE(const Simulation* sim, const IdxType t) 
{
    const IdxType qubit = sim->circuit_handle->fused_circuit[t].qubit;
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale));
    ValType* sv_real = sim->sv_real;
    ValType* sv_imag = sim->sv_imag;

    void* args[] = {(void*)(&sim->sim_gpu),(void*)&t};
    LAUNCH_KERNEL(RESETLC1_GATE);
    
    ValType sum = 0;
    ValType prob_of_one = 0;
    hipSafeCall(hipMemcpy(&sum, &(sim->m_real[0]), sizeof(ValType), hipMemcpyHostToDevice));
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allreduce(&sum, &prob_of_one, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if (prob_of_one < 1.0) //still possible to normalize
    {
        void* args[] = {(void*)(&sim->sim_gpu), (void*)&prob_of_one, (void*)&t};
        LAUNCH_KERNEL(RESETLC2_GATE);
    }
    else
    {
        if ((qubit+sim->n_qubits)>=sim->lg2_m_gpu) //remote qubit, need switch
        {
            IdxType pair_gpu = (sim->i_gpu)^((IdxType)1<<(qubit-(sim->lg2_m_gpu)));
            assert(pair_gpu != sim->i_gpu);

            ValType* sv_real_remote = sim->m_real;
            ValType* sv_imag_remote = sim->m_imag;

            if (sim->i_gpu>pair_gpu) 
            {
                MPI_Send(sv_real, per_pe_work, MPI_DOUBLE, pair_gpu, 0, MPI_COMM_WORLD);
                MPI_Send(sv_imag, per_pe_work, MPI_DOUBLE, pair_gpu, 1, MPI_COMM_WORLD);
                MPI_Recv(sv_real_remote, per_pe_work, MPI_DOUBLE, pair_gpu, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(sv_imag_remote, per_pe_work, MPI_DOUBLE, pair_gpu, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else
            {
                MPI_Recv(sv_real_remote, per_pe_work, MPI_DOUBLE, pair_gpu, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
                MPI_Recv(sv_imag_remote, per_pe_work, MPI_DOUBLE, pair_gpu, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 

                MPI_Send(sv_real, per_pe_work, MPI_DOUBLE, pair_gpu, 2, MPI_COMM_WORLD);
                MPI_Send(sv_imag, per_pe_work, MPI_DOUBLE, pair_gpu, 3, MPI_COMM_WORLD);
            }
            MPI_Barrier(MPI_COMM_WORLD);

            hipSafeCall(hipMemcpy(sv_real, sv_real_remote, per_pe_work*sizeof(ValType), hipMemcpyDeviceToDevice));
            hipSafeCall(hipMemcpy(sv_imag, sv_imag_remote, per_pe_work*sizeof(ValType), hipMemcpyDeviceToDevice));

        }
        else
        {
            void* args[] = {(void*)(&sim->sim_gpu),(void*)&t};
            LAUNCH_KERNEL(RESETLC3_GATE);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
}


__global__ void PURITYLC_GATE(Simulation* sim, const IdxType t)
{
    grid_group grid = this_grid();
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale));
    ValType* sv_real = sim->sv_real;
    ValType* sv_imag = sim->sv_imag;
    ValType* m_real = sim->m_real;
    for (IdxType i=tid; i<sim->m_gpu; i+=blockDim.x*gridDim.x)
    {
        m_real[i] = sv_real[i]*sv_real[i] + sv_imag[i]*sv_imag[i];
    }
    grid.sync();
    //Parallel reduction
    for (IdxType k=(sim->m_gpu>>(IdxType)1); k>0; k>>=1)
    {
        for (IdxType i=tid; i<k; i+=blockDim.x*gridDim.x) 
        {
            m_real[i] += m_real[i+k];
        }
        grid.sync();
    }
}


void Purity_Check(const Simulation* sim, const IdxType t)
{
    MPI_Barrier(MPI_COMM_WORLD);
    void* args[] = {(void*)(&sim->sim_gpu), (void*)&t};
    LAUNCH_KERNEL(PURITYLC_GATE);
    ValType sum = 0;
    ValType purity = 0;
    hipSafeCall(hipMemcpy(&sum, &(sim->m_real[0]), sizeof(ValType), hipMemcpyDeviceToHost));
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&sum, &purity, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if (sim->i_gpu==0 && abs(purity-1.0) > ERROR_BAR)
    {
        Gate* g = &sim->circuit_handle->fused_circuit[t];
        printf("Purity Check fails after Gate-%lld=>%s(ctrl:%lld,qubit:%lld,theta:%lf) with %lf\n",t,OP_NAMES[g->op_name],g->ctrl,g->qubit,g->theta,purity);
    }
}

//=====================================================================================
//Per-gate execution function
void Gate::exe_op(Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
}

void simulation_kernel(Simulation* sim)
{
    for (IdxType t=0; t<(sim->n_gates); t++)
    {
        enum OP op_name = sim->circuit_handle->fused_circuit[t].op_name;
        if (op_name == RESET)
        {
            //RESET_GATE(sim, t);
        }
        else if (op_name == M)
        {
            //M_GATE(sim, t);
        }
        else if (op_name == MA)
        {
            MA_GATE(sim, t);
        }
        else if (op_name == C1) 
        {
            C1_GATE(sim, t);
        }
        else if (op_name == C2) 
        {
            C2_GATE(sim, t);
        }
#ifdef PURITY_CHECK
        Purity_Check(sim, t);
#endif
    }
}


}; //namespace NWQSim

#endif

