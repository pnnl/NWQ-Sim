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
#ifndef DMSIM_NVGPU_MPI_CUH
#define DMSIM_NVGPU_MPI_CUH
#include <assert.h>
#include <random>
#include <complex.h>
#include <cooperative_groups.h>
#include <vector>
#include <mpi.h>
#include <sstream>
#include <string>
#include <iostream>
#include <cuda.h>
#include <nvshmem.h>
#include <nvshmemx.h>
#include <nvshmemx_error.h>
#include <mma.h>
#include "config.h"
#include "gate.h"
#include "device_noise.hpp"

namespace NWQSim
{
using namespace cooperative_groups;
using namespace std;
using namespace nvcuda;
    
class Gate;
class Simulation;
__device__ void Purity_Check(const Simulation* sim, const IdxType t, ValType* sv_real, ValType* sv_imag);
__device__ void Normalization(const Simulation* sim, ValType* sv_real, ValType* sv_imag);


//Declare major simulation kernel
__global__ void simulation_kernel(Simulation*);

/***********************************************
 * Gate Definition
 ***********************************************/
class Gate
{
public:
    Gate(enum OP _op_name, IdxType _qubit, IdxType _ctrl, ValType _theta) :
        op_name(_op_name), qubit(_qubit), ctrl(_ctrl), theta(_theta) {
            //if (_op_name == C4)
            //printf("%s(theta:%lf,ctrl:%lld,qubit:%lld)\n",OP_NAMES[op_name], theta, ctrl, qubit);
            memset(gm_real, 0, sizeof(ValType)*256);
            memset(gm_imag, 0, sizeof(ValType)*256);
        }
    ~Gate() {}

    Gate(const Gate& g):op_name(g.op_name), qubit(g.qubit), ctrl(g.ctrl), theta(g.theta)
    {
        memcpy(gm_real, g.gm_real, 256*sizeof(ValType));
        memcpy(gm_imag, g.gm_imag, 256*sizeof(ValType));
    }

    //applying the embedded gate operation on GPU side
    __device__ void exe_op(Simulation* sim, ValType* sv_real, ValType* sv_imag);

    //set gm
    void set_gm(std::complex<ValType>* gm, IdxType dim)
    {
        if (dim == 4 || dim==16)
        {
            for (IdxType i=0; i<dim; i++)
                for (IdxType j=0; j<dim; j++)
                {
                    gm_real[i*dim+j] = gm[i*dim+j].real();
                    gm_imag[i*dim+j] = gm[i*dim+j].imag();
                }
        }
        else
        {
            throw std::logic_error("Dim should be 4 (1-qubit gate) or 16 (2-qubit gate)!");
        }
    }
    //Gate name
    enum OP op_name;
    //Qubits 
    IdxType qubit;
    IdxType ctrl;
    ValType theta;
    //4-qubit gate parameters (after fusion)
    ValType gm_real[256];
    ValType gm_imag[256];
}; //end of Gate definition

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
        int i_gpu = nvshmem_my_pe();
        //if (i_gpu == 0)
        {
            printf("%s(theta:%lf,ctrl:%lld,qubit:%lld)\n",OP_NAMES[g.op_name], g.theta, g.ctrl, g.qubit);
            fflush(stdout);
        }
#endif
        if (g.qubit >= n_qubits || ((g.op_name != MA) && (g.ctrl >= n_qubits))) 
        {
            printf("%s(theta:%lf,ctrl:%lld,qubit:%lld)\n",OP_NAMES[g.op_name], g.theta, g.ctrl, g.qubit);
            throw std::logic_error("Gate uses qubit out of range!");
        }

        circuit.push_back(g);
        delete (&g);
        n_gates++;
        
    }
    void AllocateQubit() 
    { 
        n_qubits++; 
#ifdef PRINT_QUBIT_ALC_AND_RLS
        printf("allocate 1 qubit, now in total: %lld\n",n_qubits);
#endif
    }
    void ReleaseQubit()
    {
        n_qubits--;
#ifdef PRINT_QUBIT_ALC_AND_RLS
        printf("release 1 qubit, now in total: %lld\n", n_qubits);
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
    //This is to fuse two gates (both 1-qubit gates or 2-qubit gates)
    //This is essentially matrix multiplication
    void fuse_gate(const Gate& g0, const Gate& g1, Gate& g, const IdxType dim)
    {
        for (IdxType m=0; m<dim; m++)
        {
            for (IdxType n=0; n<dim; n++)
            {
                g.gm_real[m*dim+n] = 0;
                g.gm_imag[m*dim+n] = 0;
                for (IdxType k=0; k<dim; k++)
                {
                    g.gm_real[m*dim+n] += g0.gm_real[m*dim+k] * g1.gm_real[k*dim+n] - 
                        g0.gm_imag[m*dim+k] * g1.gm_imag[k*dim+n];
                    g.gm_imag[m*dim+n] += g0.gm_real[m*dim+k] * g1.gm_imag[k*dim+n] +
                        g0.gm_imag[m*dim+k] * g1.gm_real[k*dim+n];
                }
            }
        }
    }
    //This is kronecker product
    void kron(const Gate& g0, const Gate& g1, Gate& g, const IdxType dim)
    {
        for (IdxType r = 0; r < dim; r++) 
        {
            for (IdxType s = 0; s < dim; s++) 
            {
                for (IdxType v = 0; v < dim; v++) 
                {
                    for (IdxType w = 0; w < dim; w++) 
                    {
                        g.gm_real[(dim*r+v)*dim*dim+(dim*s+w)] = g0.gm_real[r*dim+s] * g1.gm_real[v*dim+w] - 
                            g0.gm_imag[r*dim+s] * g1.gm_imag[v*dim+w];
                        g.gm_imag[(dim*r+v)*dim*dim+(dim*s+w)] = g0.gm_real[r*dim+s] * g1.gm_imag[v*dim+w] + 
                            g0.gm_imag[r*dim+s] * g1.gm_real[v*dim+w];
                    }
                }
            }
        }
    }
    //This function is used for switch ctrl and qubit in C4 gate matrix
    void reverse_ctrl_qubit(const Gate& g0, Gate& g)
    {
        Gate SWAP(OP::C2, -1, -1, 0);
        SWAP.gm_real[0] = SWAP.gm_real[6] = SWAP.gm_real[9] = SWAP.gm_real[15] = 1;
        Gate Super_SWAP(OP::C4, -1, -1, 0);
        kron(SWAP, SWAP, Super_SWAP, 4);

        Gate tmp_g(OP::C4, -1, -1, 0);
        fuse_gate(g0, Super_SWAP, tmp_g, 16);
        fuse_gate(Super_SWAP, tmp_g, g, 16);
    }
    //This function is to fuse C4 gates (allows switching ctrl/qubit, e.g., in a SWAP gate)
    void gate_fusion_2q(const vector<Gate>& circuit_in, vector<Gate>& circuit_out)
    {
        //prepare
        IdxType* table = new IdxType[n_qubits*n_qubits];
        bool* canfuse = new bool[n_qubits*n_qubits];
        for (IdxType i=0; i<n_qubits*n_qubits; i++) table[i] = -1;
        for (IdxType i=0; i<n_qubits*n_qubits; i++) canfuse[i] = false;
        circuit_out.clear();

        for (IdxType i=0; i<circuit_in.size(); i++)
        {
            if (circuit_in[i].op_name == C2) //1-qubit gate
            {
                IdxType qubit = circuit_in[i].qubit;
                Gate g(circuit_in[i]);
                for (IdxType j=0; j<n_qubits; j++)
                {
                    canfuse[j*n_qubits+qubit] = false;
                    canfuse[qubit*n_qubits+j] = false;
                }
                circuit_out.push_back(g);
            }
            if (circuit_in[i].op_name == C4) //2-qubit gate
            {
                IdxType qubit = circuit_in[i].qubit;
                IdxType ctrl = circuit_in[i].ctrl;
                Gate g0(circuit_in[i]);
                
                //we reverse ctrl-qubit to harvest more fusion opportunities
                // M = SWAP N SWAP
                /*
                if (ctrl > qubit) 
                {
                    reverse_ctrl_qubit(circuit_in[i], g0);
                    g0.ctrl = circuit_in[i].qubit;
                    ctrl = circuit_in[i].qubit;
                    g0.qubit = circuit_in[i].ctrl;
                    qubit = circuit_in[i].ctrl;
                }
                */

                if (canfuse[ctrl*n_qubits+qubit] == false) //cannot fuse
                {
                    Gate g(g0);
                    circuit_out.push_back(g);
                    table[ctrl*n_qubits+qubit] = circuit_out.size()-1; //point for later fusion
                    for (IdxType j=0; j<n_qubits; j++)
                    {
                        canfuse[qubit*n_qubits+j] = false;
                        canfuse[ctrl*n_qubits+j] = false;
                        canfuse[j*n_qubits+qubit] = false;
                        canfuse[j*n_qubits+ctrl] = false;
                    }
                    canfuse[ctrl*n_qubits+qubit] = true;
                }
                else //able to fuse
                {
                    Gate& g1 = circuit_out[table[ctrl*n_qubits+qubit]];
                    Gate final_g(OP::C4, -1, -1, 0);
                    fuse_gate(g0, g1, final_g, 16);
                    memcpy(g1.gm_real, final_g.gm_real, 256*sizeof(ValType));
                    memcpy(g1.gm_imag, final_g.gm_imag, 256*sizeof(ValType));
                }

            }
            if (circuit_in[i].op_name == M) //1-qubit measure gate
            {
                IdxType qubit = circuit_in[i].qubit;
                Gate g(circuit_in[i]);
                for (IdxType j=0; j<n_qubits; j++)
                {
                    canfuse[j*n_qubits+qubit] = false;
                    canfuse[qubit*n_qubits+j] = false;
                }
                circuit_out.push_back(g);
            }
            if (circuit_in[i].op_name == RESET) //1-qubit reset gate
            {
                IdxType qubit = circuit_in[i].qubit;
                Gate g(circuit_in[i]);
                for (IdxType j=0; j<n_qubits; j++)
                {
                    canfuse[j*n_qubits+qubit] = false;
                    canfuse[qubit*n_qubits+j] = false;
                }
                circuit_out.push_back(g);
            }
            if (circuit_in[i].op_name == MA) //all-qubit measure gate
            {
                Gate g(circuit_in[i]);
                for (IdxType j=0; j<n_qubits; j++)
                    for (IdxType k=0; k<n_qubits; k++)
                        canfuse[j*n_qubits+k] = false;
                circuit_out.push_back(g);
            }
            //*/
        }
        //clean
        delete[] table;
        delete[] canfuse;
    }
    //This function is to fuse C2 gates in a circuit
    void gate_fusion_1q(const vector<Gate>& circuit_in, vector<Gate>& circuit_out)
    {
        //prepare
        IdxType* table = new IdxType[n_qubits];
        bool* canfuse = new bool[n_qubits];
        for (IdxType i=0; i<n_qubits; i++) table[i] = -1;
        for (IdxType i=0; i<n_qubits; i++) canfuse[i] = false;
        circuit_out.clear();
        //parse
        for (IdxType i=0; i<circuit_in.size(); i++)
        {
            if (circuit_in[i].op_name == C2) //1-qubit gate
            {
                IdxType qubit = circuit_in[i].qubit;
                if (canfuse[qubit] == false) //cannot fuse
                {
                    Gate g(circuit_in[i]);
                    circuit_out.push_back(g);
                    canfuse[qubit] = true;
                    table[qubit] = circuit_out.size()-1; //point to this gate for later fusion
                }
                else //able to fuse
                {
                    //start to fuse circuit_in[i] and circuit_out[table[qubit]]
                    const Gate& g0 = circuit_in[i];
                    Gate& g1 = circuit_out[table[qubit]];

                    ValType res_real[16] = {0};
                    ValType res_imag[16] = {0};
                    for (int m=0; m<4; m++)
                        for (int n=0; n<4; n++)
                            for (int k=0; k<4; k++)
                            {
                                res_real[m*4+n] += g0.gm_real[m*4+k] * g1.gm_real[k*4+n] - 
                                                   g0.gm_imag[m*4+k] * g1.gm_imag[k*4+n];
                                res_imag[m*4+n] += g0.gm_real[m*4+k] * g1.gm_imag[k*4+n] +
                                                   g0.gm_imag[m*4+k] * g1.gm_real[k*4+n];
                            }
                    memcpy(g1.gm_real, res_real, 16*sizeof(ValType));
                    memcpy(g1.gm_imag, res_imag, 16*sizeof(ValType));
                }
            }
            if (circuit_in[i].op_name == C4) //2-qubit gate
            {
                canfuse[circuit_in[i].qubit] = false;
                canfuse[circuit_in[i].ctrl] = false;
                Gate g(circuit_in[i]);
                circuit_out.push_back(g);
            }
            if (circuit_in[i].op_name == M) //1-qubit measure gate
            {
                canfuse[circuit_in[i].qubit] = false;
                Gate g(circuit_in[i]);
                circuit_out.push_back(g);
            }
            if (circuit_in[i].op_name == RESET) //1-qubit reset gate
            {
                canfuse[circuit_in[i].qubit] = false;
                Gate g(circuit_in[i]);
                circuit_out.push_back(g);
            }
            if (circuit_in[i].op_name == MA) //all-qubit measure gate
            {
                Gate g(circuit_in[i]);
                circuit_out.push_back(g);
                for (IdxType q=0; q<n_qubits; q++) canfuse[q] = false;
            }
        }
        //clean
        delete[] table;
        delete[] canfuse;
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
        gate_fusion_1q(circuit, tmp_circuit);
        gate_fusion_2q(tmp_circuit, fused_circuit);

        //gate_fusion_2q(circuit, fused_circuit);

        this->n_gates = fused_circuit.size();
        SAFE_FREE_GPU(circuit_gpu);
        SAFE_ALOC_GPU(circuit_gpu, n_gates*sizeof(Gate));
        cudaSafeCall(cudaMemcpy(circuit_gpu, fused_circuit.data(), n_gates*sizeof(Gate), cudaMemcpyHostToDevice));
        //tmp_circuit.clear();
        //====================================================
#endif

        return circuit_gpu;
    }
public:
    // number of logical qubits in the circuit
    IdxType n_qubits;
    // number of physical qubits in the device
    IdxType n_phy_qubits;
    // number of gpus
    IdxType n_gpus;
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
    Simulation() : 
        comm_global(MPI_COMM_WORLD),
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
        //initialization
        device_name = DEVICE_CONFIG_NAME;
        device_config = readConfigFile(DEVICE_CONFIG_PATH, device_name);
        n_qubits = device_config["num_qubits"];
        if (n_qubits > N_QUBIT_SLOT) 
            throw std::invalid_argument("Requesting more qubits than available slots!");
        dim = ((IdxType)1<<(2*n_qubits)); 
        half_dim = ((IdxType)1<<(2*n_qubits-1));
        sv_size = (dim*(IdxType)sizeof(ValType));
        //always be 0 since 1-MPI maps to 1-GPU
        cudaSafeCall(cudaSetDevice(0));
        gpu_scale = floor(log((double)n_gpus+0.5)/log(2.0));
        lg2_m_gpu = 2*n_qubits-gpu_scale;
        m_gpu = ((IdxType)1<<(lg2_m_gpu));
        sv_size_per_gpu = sv_size/n_gpus;
        //CPU side initialization
        assert(is_power_of_2(n_gpus));
        assert(dim % n_gpus == 0);
        if (!is_power_of_2(n_gpus))
            throw std::logic_error("Error: Number of GPUs should be power of 2.");
        if (dim % n_gpus != 0)
            throw std::logic_error("Error: Number of GPUs is too large or too small.");
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
        m_real = (ValType*)nvshmem_malloc(sv_size_per_gpu+1);
        m_imag = (ValType*)nvshmem_malloc(sv_size_per_gpu+1);
        cudaCheckError(); 
        gpu_mem += sv_size_per_gpu*4;
        //Initialize Circuit 
        circuit_handle = new Circuit(n_qubits);
        circuit_handle->n_phy_qubits = n_qubits;
        circuit_handle->n_gpus = n_gpus;
        circuit_handle_gpu = NULL;
        //GPU memory initilization
        cudaSafeCall(cudaMemcpy(sv_real, sv_real_cpu, sv_size_per_gpu, 
                    cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(sv_imag, sv_imag_cpu, sv_size_per_gpu, 
                    cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemset(m_real, 0, sv_size_per_gpu+1));
        cudaSafeCall(cudaMemset(m_imag, 0, sv_size_per_gpu+1));
        SAFE_ALOC_GPU(sim_gpu, sizeof(Simulation));
        //Use time as random seed

        //IdxType seed = time(0); 
        //rng.seed(1657327996);
        //rng.seed(16573279);
        rng.seed(time(0));
        //printf("Seed is %lld\n",seed);
        //rng.seed(1657421151);
    }
    ~Simulation()
    {
        //Release circuit
        if (circuit_handle != NULL)
            delete circuit_handle;
        //Release for NVSHMEM
        nvshmem_free(sv_real);
        nvshmem_free(sv_imag);
        nvshmem_free(m_real);
        nvshmem_free(m_imag);
        //Release for CPU side
        SAFE_FREE_HOST(sv_real_cpu);
        SAFE_FREE_HOST(sv_imag_cpu);
        SAFE_FREE_HOST(randoms);
        SAFE_FREE_HOST(results);
        //Release for GPU side
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
        //if (i_gpu == 0) 
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
        Gate* G = new Gate(OP::RESET,qubit,-1,0);
        circuit_handle->append(*G);
    }
    void reset_sim()
    {
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
        cudaSafeCall(cudaMemset(m_real, 0, sv_size_per_gpu+1));
        cudaSafeCall(cudaMemset(m_imag, 0, sv_size_per_gpu+1));
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
    void update(const IdxType _n_qubits, const IdxType _n_gates)
    {
        //assert(_n_qubits <= (N_QUBIT_SLOT/2));
        //this->n_qubits = _n_qubits;
        //this->dim = ((IdxType)1<<(2*n_qubits));
        //this->half_dim = (IdxType)1<<(2*n_qubits-1);
        //this->sv_size = dim*(IdxType)sizeof(ValType);
        //this->lg2_m_gpu = 2*n_qubits-gpu_scale;
        //this->m_gpu = ((IdxType)1<<(lg2_m_gpu));
        //this->sv_size_per_gpu = sv_size/n_gpus;
        this->n_gates = _n_gates;
    }
    void sim()
    {
        cudaSafeCall(cudaSetDevice(0));
        double* sim_times;
#ifdef PRINT_KERNEL_TRACE
        if (i_gpu == 0)
        {
            SAFE_ALOC_HOST(sim_times, sizeof(double)*n_gpus);
            memset(sim_times, 0, sizeof(double)*n_gpus);
        }
#endif
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


        //if (i_gpu ==0)
        //{
            //printf("\n============== DM-Sim ===============\n");
            //printf("nqubits:%lld, ngates:%lld, sim_gates:%lld, ngpus:%lld\n",
                    //n_qubits, input_gates, n_gates, n_gpus);
            //printf("=====================================\n");
        //}
 

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
            printf("\n============== DM-Sim ===============\n");
            printf("nqubits:%lld, ngates:%lld, sim_gates:%lld, ngpus:%lld, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms, mem:%.3lf MB, mem_per_gpu:%.3lf MB\n",
                    n_qubits, input_gates, n_gates, n_gpus, avg_sim_time, 0., 
                    avg_sim_time, gpu_mem/1024/1024, gpu_mem/1024/1024);
            printf("=====================================\n");

            SAFE_FREE_HOST(sim_times);
        }
#endif
        reset_circuit();
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
            printf("----- Real DM diag ------\n");
            for (IdxType i=0; i<num; i++) 
            {
                printf("%lf ", sv_diag_real[i*num+i]);
                if ((i+1)%8==0) printf("\n");
            }
            printf("\n");
            //printf("----- Imag DM diag ------\n");
            //for (IdxType i=0; i<num; i++) 
            //printf("%lf ", sv_diag_imag[i*num+i]);
            //printf("\n");
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
    MPI_Comm comm_global;
};

#define PGAS_P(arr,i,val) nvshmem_double_p(&(arr)[(i)&((sim->m_gpu)-1)], (val), ((i)>>(sim->lg2_m_gpu)) )
#define PGAS_G(arr,i) nvshmem_double_g(&(arr)[(i)&((sim->m_gpu)-1)], ((i)>>(sim->lg2_m_gpu)) )
#define BARRIER if(threadIdx.x==0 && blockIdx.x==0) nvshmem_barrier_all(); grid.sync();

//============== Check Trace (debug purpose) ================
__device__ __inline__ void CHECK_TRACE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, IdxType t)
{
    grid_group grid = this_grid();
    if (sim->i_gpu == 0 && blockIdx.x ==0 && threadIdx.x == 0)
    {
        ValType trace = 0;
        for (IdxType i=0; i<((IdxType)1<<(sim->n_qubits)); i++)
        {
            const ValType val = PGAS_G(sv_real, (i<<(sim->n_qubits))+i);
            trace += abs(val);
        }
        printf("%s: Trace is: %lf\n", OP_NAMES_NVGPU[sim->circuit_handle_gpu[t].op_name], trace);
    }
    BARRIER;
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
            printf("%s(%lld) ticks: %lld\n",OP_NAMES_NVGPU[sim->circuit_handle_gpu[t].op_name],sim->circuit_handle_gpu[t].qubit, time);
        }
         */
        //CHECK_TRACE(sim, sim->sv_real, sim->sv_imag, t);

#ifdef PURITY_CHECK
        Purity_Check(sim, t, sim->sv_real, sim->sv_imag);
#endif
    }
    //Normalization(sim, sim->sv_real, sim->sv_imag);
}

//================================= Gate Definition ========================================
//============== Unified 2-qubit Gate without comm optimization ================
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

        const ValType el0_real = PGAS_G(sv_real, pos0);
        const ValType el0_imag = PGAS_G(sv_imag, pos0);
        const ValType el1_real = PGAS_G(sv_real, pos1);
        const ValType el1_imag = PGAS_G(sv_imag, pos1);
        const ValType el2_real = PGAS_G(sv_real, pos2);
        const ValType el2_imag = PGAS_G(sv_imag, pos2);
        const ValType el3_real = PGAS_G(sv_real, pos3);
        const ValType el3_imag = PGAS_G(sv_imag, pos3);

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

        PGAS_P(sv_real, pos0, sv_real_pos0); 
        PGAS_P(sv_real, pos1, sv_real_pos1); 
        PGAS_P(sv_real, pos2, sv_real_pos2); 
        PGAS_P(sv_real, pos3, sv_real_pos3); 

        PGAS_P(sv_imag, pos0, sv_imag_pos0); 
        PGAS_P(sv_imag, pos1, sv_imag_pos1); 
        PGAS_P(sv_imag, pos2, sv_imag_pos2); 
        PGAS_P(sv_imag, pos3, sv_imag_pos3); 
    }
    //BARRIER;
}


#define LOCAL_G(arr,i) arr[(i)&(sim->m_gpu-1)]
#define LOCAL_P(arr,i,val) arr[(i)&(sim->m_gpu-1)] = val;
#define SV4IDX(x) (((x>>1)&1)*EXP2E(qubit0) + ((x&1)*EXP2E(qubit1)) )

#define DIV2E(x,y) ((x)>>(y))
#define MOD2E(x,y) ((x)&(((IdxType)1<<(y))-(IdxType)1)) 
#define EXP2E(x) ((IdxType)1<<(x))
#define SV16IDX(x) ( ((x>>3)&1)*EXP2E(qubit0) + ((x>>2)&1)*EXP2E(qubit1) + ((x>>1)&1)*EXP2E(qubit2) + ((x&1)*EXP2E(qubit3)) )

//============== Unified 2-qubit Gate ================
//Perform communication optimization here
//Since this is for single-qubit density matrix gate,
//It is only possible qubit0 < qubit1 and qubit1 can be remote
__device__ __inline__ void C2V1_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag, const IdxType qubit0, const IdxType qubit1)
{
    assert (qubit0 != qubit1); //Non-cloning

    if (qubit1 < sim->lg2_m_gpu) 
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
        //BARRIER;
    }
}

//============== SWAP Gate ================
//This gate is for internal usage. It is used
//when two qubits of a C4 gate are remote qubis, we then
//swap one of them to a local qubit without noise,
//perform the C4 gate, and then swap back
//It is assumed qubit0 is local, qubit1 is remote
__device__ __inline__ void SWAP_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const IdxType qubit0, const IdxType qubit1)
{
    grid_group grid = this_grid(); 
    const int tid = blockDim.x * blockIdx.x + threadIdx.x; 
    assert (qubit0 < qubit1); //Non-cloning and qubit0<qubit1 assumption
    assert (qubit1 >= sim->lg2_m_gpu); //Otherwise don't need swap
    
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale+1));
    const IdxType per_pe_num = ((sim->dim)>>(sim->gpu_scale));

    const IdxType p = qubit0;
    const IdxType q = qubit1;

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
}

//============== Unified 4-qubit Gate ================
__device__ __inline__ void C4_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag, const IdxType qubit0, const IdxType qubit1,
        const IdxType qubit2, const IdxType qubit3)
{
    grid_group grid = this_grid(); 
    const int tid = blockDim.x * blockIdx.x + threadIdx.x; 
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale+4));
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
            PGAS_G(sv_real,term+SV16IDX(0)),  PGAS_G(sv_real,term+SV16IDX(1)),
            PGAS_G(sv_real,term+SV16IDX(2)),  PGAS_G(sv_real,term+SV16IDX(3)),
            PGAS_G(sv_real,term+SV16IDX(4)),  PGAS_G(sv_real,term+SV16IDX(5)),
            PGAS_G(sv_real,term+SV16IDX(6)),  PGAS_G(sv_real,term+SV16IDX(7)),
            PGAS_G(sv_real,term+SV16IDX(8)),  PGAS_G(sv_real,term+SV16IDX(9)),
            PGAS_G(sv_real,term+SV16IDX(10)), PGAS_G(sv_real,term+SV16IDX(11)),
            PGAS_G(sv_real,term+SV16IDX(12)), PGAS_G(sv_real,term+SV16IDX(13)),
            PGAS_G(sv_real,term+SV16IDX(14)), PGAS_G(sv_real,term+SV16IDX(15))
        };
        const ValType el_imag[16] = { 
            PGAS_G(sv_imag,term+SV16IDX(0)),  PGAS_G(sv_imag,term+SV16IDX(1)),
            PGAS_G(sv_imag,term+SV16IDX(2)),  PGAS_G(sv_imag,term+SV16IDX(3)),
            PGAS_G(sv_imag,term+SV16IDX(4)),  PGAS_G(sv_imag,term+SV16IDX(5)),
            PGAS_G(sv_imag,term+SV16IDX(6)),  PGAS_G(sv_imag,term+SV16IDX(7)),
            PGAS_G(sv_imag,term+SV16IDX(8)),  PGAS_G(sv_imag,term+SV16IDX(9)),
            PGAS_G(sv_imag,term+SV16IDX(10)), PGAS_G(sv_imag,term+SV16IDX(11)),
            PGAS_G(sv_imag,term+SV16IDX(12)), PGAS_G(sv_imag,term+SV16IDX(13)),
            PGAS_G(sv_imag,term+SV16IDX(14)), PGAS_G(sv_imag,term+SV16IDX(15))
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
            PGAS_P(sv_real, term+SV16IDX(j), res_real);
            PGAS_P(sv_imag, term+SV16IDX(j), res_imag);
        }
    }
    //BARRIER;
}

//============== Unified 4-qubit Gate with comm optimization ================
//This is with nvshmem merged communication optimization but without tensorcore
__device__ __inline__ void C4V1_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag, const IdxType qubit0, const IdxType qubit1,
        const IdxType qubit2, const IdxType qubit3)
{
    grid_group grid = this_grid(); 
    
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
        
    if (s < sim->lg2_m_gpu) //all 4 qubits are local
    {
        C4_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, qubit0, qubit1, qubit2, qubit3); 
    }
    else //s qubit is non-local
    {
        assert(r < sim->lg2_m_gpu);
        const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x; 
        const int laneid = (threadIdx.x & 31);
        const int warpid = (tid >> 5);
        const int wid = (threadIdx.x >> 5);
        const int nwarps = ((blockDim.x * gridDim.x)>>5);
        const IdxType per_pe_num = ((sim->dim)>>(sim->gpu_scale));
        const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale+3));

        //load data from pair GPU
        IdxType pair_gpu = (sim->i_gpu)^((IdxType)1<<(s-(sim->lg2_m_gpu)));
        if (sim->i_gpu>pair_gpu) return;

        ValType* sv_real_remote = sim->m_real;
        ValType* sv_imag_remote = sim->m_imag;

        if (tid == 0) nvshmem_double_get(sv_real_remote, sv_real, per_pe_num, pair_gpu );
        if (tid == 0) nvshmem_double_get(sv_imag_remote, sv_imag, per_pe_num, pair_gpu );
        grid.sync();

        for (IdxType i=(sim->i_gpu)*per_pe_work+tid; i<(sim->i_gpu+1)*per_pe_work;
                i+=blockDim.x*gridDim.x)
        {
            const IdxType term0 = MOD2E(i,p);
            const IdxType term1 = MOD2E(DIV2E(i,p),q-p-1)*EXP2E(p+1);
            const IdxType term2 = MOD2E(DIV2E(DIV2E(i,p),q-p-1),r-q-1)*EXP2E(q+1);
            const IdxType term3 = MOD2E(DIV2E(DIV2E(DIV2E(i,p),q-p-1),r-q-1),s-r-1)*EXP2E(r+1);
            const IdxType term4 = DIV2E(DIV2E(DIV2E(DIV2E(i,p),q-p-1),r-q-1),s-r-1)*EXP2E(s+1);
            const IdxType term = term4 + term3 + term2 + term1 + term0;
            ValType el_real[16];
            ValType el_imag[16];

            if (qubit3 == s) //qubit3 is remote qubit
            {
                el_real[0] = LOCAL_G(sv_real,term+SV16IDX(0));
                el_real[1] = LOCAL_G(sv_real_remote,term+SV16IDX(1));
                el_real[2] = LOCAL_G(sv_real,term+SV16IDX(2));
                el_real[3] = LOCAL_G(sv_real_remote,term+SV16IDX(3));
                el_real[4] = LOCAL_G(sv_real,term+SV16IDX(4));
                el_real[5] = LOCAL_G(sv_real_remote,term+SV16IDX(5));
                el_real[6] = LOCAL_G(sv_real,term+SV16IDX(6));
                el_real[7] = LOCAL_G(sv_real_remote,term+SV16IDX(7));
                el_real[8] = LOCAL_G(sv_real,term+SV16IDX(8));
                el_real[9] = LOCAL_G(sv_real_remote,term+SV16IDX(9));
                el_real[10] = LOCAL_G(sv_real,term+SV16IDX(10));
                el_real[11] = LOCAL_G(sv_real_remote,term+SV16IDX(11));
                el_real[12] = LOCAL_G(sv_real,term+SV16IDX(12));
                el_real[13] = LOCAL_G(sv_real_remote,term+SV16IDX(13));
                el_real[14] = LOCAL_G(sv_real,term+SV16IDX(14));
                el_real[15] = LOCAL_G(sv_real_remote,term+SV16IDX(15));

                el_imag[0] = LOCAL_G(sv_imag,term+SV16IDX(0));
                el_imag[1] = LOCAL_G(sv_imag_remote,term+SV16IDX(1));
                el_imag[2] = LOCAL_G(sv_imag,term+SV16IDX(2));
                el_imag[3] = LOCAL_G(sv_imag_remote,term+SV16IDX(3));
                el_imag[4] = LOCAL_G(sv_imag,term+SV16IDX(4));
                el_imag[5] = LOCAL_G(sv_imag_remote,term+SV16IDX(5));
                el_imag[6] = LOCAL_G(sv_imag,term+SV16IDX(6));
                el_imag[7] = LOCAL_G(sv_imag_remote,term+SV16IDX(7));
                el_imag[8] = LOCAL_G(sv_imag,term+SV16IDX(8));
                el_imag[9] = LOCAL_G(sv_imag_remote,term+SV16IDX(9));
                el_imag[10] = LOCAL_G(sv_imag,term+SV16IDX(10));
                el_imag[11] = LOCAL_G(sv_imag_remote,term+SV16IDX(11));
                el_imag[12] = LOCAL_G(sv_imag,term+SV16IDX(12));
                el_imag[13] = LOCAL_G(sv_imag_remote,term+SV16IDX(13));
                el_imag[14] = LOCAL_G(sv_imag,term+SV16IDX(14));
                el_imag[15] = LOCAL_G(sv_imag_remote,term+SV16IDX(15));
            }
            else //qubit2 is remote (not possible qubit0 or 1 is remote)
            {
                //if (( (laneid>>1)&1)!=0) el_real_s[j*16+laneid] = LOCAL_G(sv_real_remote,addr);
                //else el_real_s[j*16+laneid] = LOCAL_G(sv_real, addr);
                //laneid = 2,3,6,7,10,11,14,15;

                el_real[0] = LOCAL_G(sv_real,term+SV16IDX(0));
                el_real[1] = LOCAL_G(sv_real,term+SV16IDX(1));
                el_real[2] = LOCAL_G(sv_real_remote,term+SV16IDX(2));
                el_real[3] = LOCAL_G(sv_real_remote,term+SV16IDX(3));
                el_real[4] = LOCAL_G(sv_real,term+SV16IDX(4));
                el_real[5] = LOCAL_G(sv_real,term+SV16IDX(5));
                el_real[6] = LOCAL_G(sv_real_remote,term+SV16IDX(6));
                el_real[7] = LOCAL_G(sv_real_remote,term+SV16IDX(7));
                el_real[8] = LOCAL_G(sv_real,term+SV16IDX(8));
                el_real[9] = LOCAL_G(sv_real,term+SV16IDX(9));
                el_real[10] = LOCAL_G(sv_real_remote,term+SV16IDX(10));
                el_real[11] = LOCAL_G(sv_real_remote,term+SV16IDX(11));
                el_real[12] = LOCAL_G(sv_real,term+SV16IDX(12));
                el_real[13] = LOCAL_G(sv_real,term+SV16IDX(13));
                el_real[14] = LOCAL_G(sv_real_remote,term+SV16IDX(14));
                el_real[15] = LOCAL_G(sv_real_remote,term+SV16IDX(15));

                el_imag[0] = LOCAL_G(sv_imag,term+SV16IDX(0));
                el_imag[1] = LOCAL_G(sv_imag,term+SV16IDX(1));
                el_imag[2] = LOCAL_G(sv_imag_remote,term+SV16IDX(2));
                el_imag[3] = LOCAL_G(sv_imag_remote,term+SV16IDX(3));
                el_imag[4] = LOCAL_G(sv_imag,term+SV16IDX(4));
                el_imag[5] = LOCAL_G(sv_imag,term+SV16IDX(5));
                el_imag[6] = LOCAL_G(sv_imag_remote,term+SV16IDX(6));
                el_imag[7] = LOCAL_G(sv_imag_remote,term+SV16IDX(7));
                el_imag[8] = LOCAL_G(sv_imag,term+SV16IDX(8));
                el_imag[9] = LOCAL_G(sv_imag,term+SV16IDX(9));
                el_imag[10] = LOCAL_G(sv_imag_remote,term+SV16IDX(10));
                el_imag[11] = LOCAL_G(sv_imag_remote,term+SV16IDX(11));
                el_imag[12] = LOCAL_G(sv_imag,term+SV16IDX(12));
                el_imag[13] = LOCAL_G(sv_imag,term+SV16IDX(13));
                el_imag[14] = LOCAL_G(sv_imag_remote,term+SV16IDX(14));
                el_imag[15] = LOCAL_G(sv_imag_remote,term+SV16IDX(15));
            }
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
                //PGAS_P(sv_real, term+SV16IDX(j), res_real);
                //PGAS_P(sv_imag, term+SV16IDX(j), res_imag);
                if (qubit3 == s) //qubit3 is remote qubit
                {
                    if ((j&1)!=0)
                    {
                        LOCAL_P(sv_real_remote, term+SV16IDX(j), res_real);
                        LOCAL_P(sv_imag_remote, term+SV16IDX(j), res_imag);
                    }
                    else
                    {
                        LOCAL_P(sv_real, term+SV16IDX(j), res_real);
                        LOCAL_P(sv_imag, term+SV16IDX(j), res_imag);
                    }
                }
                else //qubit2 is remote (not possible qubit0 or 1 is remote)
                {
                    if (((j>>1)&1)!=0)
                    {
                        LOCAL_P(sv_real_remote, term+SV16IDX(j), res_real);
                        LOCAL_P(sv_imag_remote, term+SV16IDX(j), res_imag);
                    }
                    else
                    {
                        LOCAL_P(sv_real, term+SV16IDX(j), res_real);
                        LOCAL_P(sv_imag, term+SV16IDX(j), res_imag);
                    }
                }
            }
        }
        grid.sync();
        if (tid == 0) nvshmem_double_put(sv_real, sv_real_remote, per_pe_num, pair_gpu);
        if (tid == 0) nvshmem_double_put(sv_imag, sv_imag_remote, per_pe_num, pair_gpu);
        grid.sync();
    }
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
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale+4));
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

    for (IdxType i=(sim->i_gpu)*per_pe_work+(warpid<<3); i<(sim->i_gpu+1)*per_pe_work; i+=(nwarps<<3)) 
    {
        //if (laneid == 0) printf("=== warpid: %d ==== \n",warpid);
#pragma unroll
        for (int j=0; j<8; j++)
        {
            const IdxType term = get_term(i+j,p,q,r,s);
            if (laneid < 16) 
            {
                IdxType addr = term + SV16IDX(laneid);
                el_real_s[j*16+laneid] = PGAS_G(sv_real, addr);
            }
            else 
            {
                IdxType addr = term + SV16IDX(laneid-16);
                el_imag_s[j*16+laneid-16] = PGAS_G(sv_imag, addr);
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
        PGAS_P(sv_real, term0, c_frag_real_up.x[0]);
        PGAS_P(sv_imag, term0, c_frag_imag_up.x[0]);

        //Store second result per segment-C
        IdxType j1 = ((2*laneid+1)&7);
        IdxType k1 = ((2*laneid+1)>>3);
        const IdxType term1 = get_term(i+j1,p,q,r,s) + SV16IDX(k1);
        PGAS_P(sv_real, term1, c_frag_real_up.x[1]);
        PGAS_P(sv_imag, term1, c_frag_imag_up.x[1]);

        //Store first result per segment-C
        IdxType j2 = ((2*laneid)&7);
        IdxType k2 = ((2*laneid)>>3)+8;
        const IdxType term2 = get_term(i+j2,p,q,r,s) + SV16IDX(k2);
        PGAS_P(sv_real, term2, c_frag_real_dn.x[0]);
        PGAS_P(sv_imag, term2, c_frag_imag_dn.x[0]);

        //Store second result per segment-C
        IdxType j3 = ((2*laneid+1)&7);
        IdxType k3 = ((2*laneid+1)>>3)+8;
        const IdxType term3 = get_term(i+j3,p,q,r,s) + SV16IDX(k3);
        PGAS_P(sv_real, term3, c_frag_real_dn.x[1]);
        PGAS_P(sv_imag, term3, c_frag_imag_dn.x[1]);
    }
    //BARRIER;
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
    const IdxType per_pe_work = ((sim->dim)>>(sim->gpu_scale+4));
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

    //if (threadIdx.x == 0 && blockIdx.x == 0) 
    //printf("GPU-%lld: dim:%lld, per_pe_work:%lld, nwarps:%d, start:%lld, end:%lld, step:%lld\n",sim->i_gpu, sim->dim,per_pe_work, nwarps,
    //(sim->i_gpu)*per_pe_work+(warpid<<3), (sim->i_gpu+1)*per_pe_work, (nwarps<<3));

    for (IdxType i=(sim->i_gpu)*per_pe_work+(warpid<<3); i<(sim->i_gpu+1)*per_pe_work; i+=(nwarps<<3)) 
    {
        //load from matrix-B
        for (int j=0; j<8; j++)
        {
            const IdxType term = get_term(i+j,p,q,r,s);
            if (laneid < 16) el_real_s[j*16+laneid] = PGAS_G(sv_real, term+SV16IDX(laneid));
            else el_imag_s[j*16+laneid-16] = PGAS_G(sv_imag, term+SV16IDX(laneid-16));
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
        PGAS_P(sv_real, term0, c_frag_real.x[0]);
        PGAS_P(sv_imag, term0, c_frag_imag.x[0]);
        //Store second result per segment-C
        IdxType j1 = ((2*laneid+1)&7);
        IdxType k1 = ((2*laneid+1)>>3);
        const IdxType term1 = get_term(i+j1,p,q,r,s) + SV16IDX(k1);
        PGAS_P(sv_real, term1, c_frag_real.x[1]);
        PGAS_P(sv_imag, term1, c_frag_imag.x[1]);

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
        PGAS_P(sv_real, term2, c_frag_real.x[0]);
        PGAS_P(sv_imag, term2, c_frag_imag.x[0]);

        //Store second result per segment-C
        IdxType j3 = ((2*laneid+1)&7);
        IdxType k3 = ((2*laneid+1)>>3)+8;
        const IdxType term3 = get_term(i+j3,p,q,r,s) + SV16IDX(k3);
        PGAS_P(sv_real, term3, c_frag_real.x[1]);
        PGAS_P(sv_imag, term3, c_frag_imag.x[1]);
    }
    //BARRIER;
}

//============== Unified 4-qubit Gate with TC ================
//This is with tensorcore and nvshmem merged communication optimization
__device__ __inline__ void C4TCV3_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, 
        const ValType* gm_real, const ValType* gm_imag, const IdxType qubit0, const IdxType qubit1,
        const IdxType qubit2, const IdxType qubit3)
{
    grid_group grid = this_grid(); 
    
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

    if (s < sim->lg2_m_gpu) //all 4 qubits are local
    {
        C4TCV2_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, qubit0, qubit1, qubit2, qubit3); 
    }
    else //s qubit is non-local
    {
        const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x; 
        const int laneid = (threadIdx.x & 31);
        const int warpid = (tid >> 5);
        const int wid = (threadIdx.x >> 5);
        const int nwarps = ((blockDim.x * gridDim.x)>>5);
        const IdxType per_pe_num = ((sim->dim)>>(sim->gpu_scale));
        const IdxType whole_work = ((sim->dim)>>4);

        extern __shared__ ValType els[];
        ValType* el_real_s = &els[wid*8*16*2];//per warp 8*16 for real 
        ValType* el_imag_s = &els[wid*8*16*2+8*16]; //per warp 8*16 for imag

        //load data from pair GPU
        IdxType pair_gpu = (sim->i_gpu)^((IdxType)1<<(s-(sim->lg2_m_gpu)));
        if (sim->i_gpu>pair_gpu) return;

        ValType* sv_real_remote = sim->m_real;
        ValType* sv_imag_remote = sim->m_imag;

        if (tid == 0) nvshmem_double_get(sv_real_remote, sv_real, per_pe_num, pair_gpu );
        if (tid == 0) nvshmem_double_get(sv_imag_remote, sv_imag, per_pe_num, pair_gpu );
        grid.sync();

        wmma::fragment<wmma::matrix_a, 8, 8, 4, ValType, wmma::row_major> a_frag_real;
        wmma::fragment<wmma::matrix_a, 8, 8, 4, ValType, wmma::row_major> a_frag_imag;
        wmma::fragment<wmma::matrix_b, 8, 8, 4, ValType, wmma::col_major> b_frag_real;
        wmma::fragment<wmma::matrix_b, 8, 8, 4, ValType, wmma::col_major> b_frag_imag;
        wmma::fragment<wmma::accumulator, 8, 8, 4, ValType> c_frag_real;
        wmma::fragment<wmma::accumulator, 8, 8, 4, ValType> c_frag_imag;

        for (IdxType i=(warpid<<3); i<whole_work; i+=(nwarps<<3))
        {
            const IdxType base = get_term(i,p,q,r,s);
            if ( (per_pe_num*(sim->i_gpu))<= base && base < (per_pe_num*(sim->i_gpu+1)))
            {
                for (int j=0; j<8; j++)
                {
                    const IdxType term = get_term(i+j,p,q,r,s);
                    if (laneid < 16)
                    {
                        IdxType addr = term+SV16IDX(laneid);
                        //This is due to how SV16IDx is defined, qubit0 maps to MSB,
                        //qubit1 maps to MSB-1, qubit2 maps to MSB-2, qubit3 maps to LSB
                        if (qubit3 == s) //qubit3 is remote qubit
                        {
                            if ((laneid&1)!=0) el_real_s[j*16+laneid] = LOCAL_G(sv_real_remote, addr);
                            else el_real_s[j*16+laneid] = LOCAL_G(sv_real, addr);
                        }
                        else //qubit2 is remote (not possible qubit0 or 1 is remote)
                        {
                            if (( (laneid>>1)&1)!=0) el_real_s[j*16+laneid] = LOCAL_G(sv_real_remote,addr);
                            else el_real_s[j*16+laneid] = LOCAL_G(sv_real, addr);
                        }
                    }
                    else 
                    {
                        IdxType addr = term+SV16IDX(laneid-16);
                        if (qubit3 == s) //qubit3 is remote qubit
                        {
                            if ((laneid&1)!=0) el_imag_s[j*16+laneid] = LOCAL_G(sv_imag_remote, addr);
                            else el_imag_s[j*16+laneid] = LOCAL_G(sv_imag, addr);
                        }
                        else //qubit2 is remote (not possible qubit0 or 1 is remote)
                        {
                            if (( (laneid>>1)&1)!=0) el_imag_s[j*16+laneid] = LOCAL_G(sv_imag_remote, addr);
                            else el_imag_s[j*16+laneid] = LOCAL_G(sv_imag, addr);
                            //el_imag_s[j*16+laneid] = PGAS_G(sv_imag, addr);
                        }
                    }
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

                if (qubit3 == s) //qubit3 is remote qubit
                {
                    if ((k0&1)!=0)
                    {
                        LOCAL_P(sv_real_remote, term0, c_frag_real.x[0]); 
                        LOCAL_P(sv_imag_remote, term0, c_frag_imag.x[0]); 
                    }
                    else 
                    {
                        LOCAL_P(sv_real, term0, c_frag_real.x[0]);
                        LOCAL_P(sv_imag, term0, c_frag_imag.x[0]);
                    }
                }
                else //qubit2 is remote (not possible qubit0 or 1 is remote)
                {
                    if (((k0>>1)&1)!=0)
                    {
                        LOCAL_P(sv_real_remote, term0, c_frag_real.x[0]); 
                        LOCAL_P(sv_imag_remote, term0, c_frag_imag.x[0]); 
                    }
                    else 
                    {
                        LOCAL_P(sv_real, term0, c_frag_real.x[0]);
                        LOCAL_P(sv_imag, term0, c_frag_imag.x[0]);
                    }
                }
                //Store second result per segment-C
                IdxType j1 = ((2*laneid+1)&7);
                IdxType k1 = ((2*laneid+1)>>3);
                const IdxType term1 = get_term(i+j1,p,q,r,s) + SV16IDX(k1);

                if (qubit3 == s) //qubit3 is remote qubit
                {
                    if ((k1&1)!=0)
                    {
                        LOCAL_P(sv_real_remote, term1, c_frag_real.x[1]); 
                        LOCAL_P(sv_imag_remote, term1, c_frag_imag.x[1]); 
                    }
                    else 
                    {
                        LOCAL_P(sv_real, term1, c_frag_real.x[1]);
                        LOCAL_P(sv_imag, term1, c_frag_imag.x[1]);
                    }
                }
                else //qubit2 is remote (not possible qubit0 or 1 is remote)
                {
                    if (((k1>>1)&1)!=0)
                    {
                        LOCAL_P(sv_real_remote, term1, c_frag_real.x[1]); 
                        LOCAL_P(sv_imag_remote, term1, c_frag_imag.x[1]); 
                    }
                    else 
                    {
                        LOCAL_P(sv_real, term1, c_frag_real.x[1]);
                        LOCAL_P(sv_imag, term1, c_frag_imag.x[1]);
                    }
                }
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

                if (qubit3 == s) //qubit3 is remote qubit
                {
                    if ((k2&1)!=0)
                    {
                        LOCAL_P(sv_real_remote, term2, c_frag_real.x[0]); 
                        LOCAL_P(sv_imag_remote, term2, c_frag_imag.x[0]); 
                    }
                    else 
                    {
                        LOCAL_P(sv_real, term2, c_frag_real.x[0]);
                        LOCAL_P(sv_imag, term2, c_frag_imag.x[0]);
                    }
                }
                else //qubit2 is remote (not possible qubit0 or 1 is remote)
                {
                    if (((k2>>1)&1)!=0)
                    {
                        LOCAL_P(sv_real_remote, term2, c_frag_real.x[0]); 
                        LOCAL_P(sv_imag_remote, term2, c_frag_imag.x[0]); 
                    }
                    else 
                    {
                        LOCAL_P(sv_real, term2, c_frag_real.x[0]);
                        LOCAL_P(sv_imag, term2, c_frag_imag.x[0]);
                    }
                }
                //Store second result per segment-C
                IdxType j3 = ((2*laneid+1)&7);
                IdxType k3 = ((2*laneid+1)>>3)+8;
                const IdxType term3 = get_term(i+j3,p,q,r,s) + SV16IDX(k3);

                if (qubit3 == s) //qubit3 is remote qubit
                {
                    if ((k3&1)!=0)
                    {
                        LOCAL_P(sv_real_remote, term3, c_frag_real.x[1]); 
                        LOCAL_P(sv_imag_remote, term3, c_frag_imag.x[1]); 
                    }
                    else 
                    {
                        LOCAL_P(sv_real, term3, c_frag_real.x[1]);
                        LOCAL_P(sv_imag, term3, c_frag_imag.x[1]);
                    }
                }
                else //qubit2 is remote (not possible qubit0 or 1 is remote)
                {
                    if (((k3>>1)&1)!=0)
                    {
                        LOCAL_P(sv_real_remote, term3, c_frag_real.x[1]); 
                        LOCAL_P(sv_imag_remote, term3, c_frag_imag.x[1]); 
                    }
                    else 
                    {
                        LOCAL_P(sv_real, term3, c_frag_real.x[1]);
                        LOCAL_P(sv_imag, term3, c_frag_imag.x[1]);
                    }
                }
            }
        }
        grid.sync();
        if (tid == 0) nvshmem_double_put(sv_real, sv_real_remote, per_pe_num, pair_gpu);
        if (tid == 0) nvshmem_double_put(sv_imag, sv_imag_remote, per_pe_num, pair_gpu);
        //BARRIER;
    }
}

#endif



///*
__device__ __inline__ void M_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag, ValType* gm_real, ValType* gm_imag,
        const IdxType qubit, const ValType rand)
{
    grid_group grid = this_grid();
    assert(sim-> n_qubits < sim->lg2_m_gpu); 
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    const IdxType per_pe_work_dm = ((sim->dim)>>(sim->gpu_scale));
    ValType * m_real = sim->m_real;
    ValType * m_imag = sim->m_imag;
    IdxType mask = ((IdxType)1<<qubit);

    if (sim->i_gpu == 0)
    {
        for (IdxType i=tid; i<((IdxType)1<<(sim->n_qubits)); i+=blockDim.x*gridDim.x)
        {
            const ValType val = PGAS_G(sv_real, (i<<(sim->n_qubits))+i);
            if ( (i & mask) == 0) m_real[i] = 0;
            else m_real[i] = abs(val);
        }
        grid.sync();
        for (IdxType k=((IdxType)1<<(sim->n_qubits-1)); k>0; k>>=1)
        {
            for (IdxType i=tid; i<k; i+=blockDim.x*gridDim.x) m_real[i] += m_real[i+k];
            grid.sync();
        }
    }
    BARRIER;
    if (tid==0 && sim->i_gpu!=0) m_real[0] = PGAS_G(m_real,0);
    grid.sync();
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
    BARRIER;
    C2V1_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, qubit, qubit+sim->n_qubits);
    BARRIER;
    if (tid==0) sim->results_gpu[0] = (rand<=prob_of_one?1:0);
}
//*/


__device__ __inline__ void Normalization(const Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    grid_group grid = this_grid();
    assert(sim-> n_qubits < sim->lg2_m_gpu); 
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    const IdxType per_pe_work_dm = ((sim->dim)>>(sim->gpu_scale));
    ValType * m_real = sim->m_real;

    if (sim->i_gpu == 0)
    {
        for (IdxType i=tid; i<((IdxType)1<<(sim->n_qubits)); i+=blockDim.x*gridDim.x)
        {
            const ValType val = PGAS_G(sv_real, (i<<(sim->n_qubits))+i);
            m_real[i] = abs(val);
        }
        grid.sync();
        for (IdxType k=((IdxType)1<<(sim->n_qubits-1)); k>0; k>>=1)
        {
            for (IdxType i=tid; i<k; i+=blockDim.x*gridDim.x) m_real[i] += m_real[i+k];
            grid.sync();
        }
    }
    BARRIER;
    if (tid==0 && sim->i_gpu!=0) m_real[0] = PGAS_G(m_real,0);
    grid.sync();
    ValType purity = m_real[0];
    for (IdxType i=tid; i<per_pe_work_dm; i+=blockDim.x*gridDim.x)
    {
        sv_real[i] /= purity;
        sv_imag[i] /= purity;
    }
}




/*
__device__ __inline__ void M_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag,
        const IdxType qubit, const ValType rand)
{
    grid_group grid = this_grid();
    assert(sim-> n_qubits < sim->lg2_m_gpu); //ensure the diagonal can be fit into GPU-0's part of m_real
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    const IdxType per_pe_work_dm = ((sim->dim)>>(sim->gpu_scale));
    ValType * m_real = sim->m_real;
    ValType * m_imag = sim->m_imag;
    IdxType mask = ((IdxType)1<<qubit);

    //if (sim->i_gpu == 0 && threadIdx.x==0 && blockIdx.x==0) printf("per_pe_work:%lld, gpu_scale:%lld, qubit:%lld \n", per_pe_work_dm, sim->gpu_scale, qubit);
    

    if (sim->i_gpu == 0)
    {
        for (IdxType i=tid; i<((IdxType)1<<(sim->n_qubits)); i+=blockDim.x*gridDim.x)
        {
            const ValType val = PGAS_G(sv_real, (i<<(sim->n_qubits))+i);
            m_imag[i] = abs(val);
            if ( (i & mask) == 0) m_real[i] = 0;
            else m_real[i] = abs(val);
        }
        grid.sync();
        for (IdxType k=((IdxType)1<<(sim->n_qubits-1)); k>0; k>>=1)
        {
            for (IdxType i=tid; i<k; i+=blockDim.x*gridDim.x) m_real[i] += m_real[i+k];
            for (IdxType i=tid; i<k; i+=blockDim.x*gridDim.x) m_imag[i] += m_imag[i+k];
            grid.sync();
        }
    }
    BARRIER;
    if (tid==0 && sim->i_gpu!=0)
    { 
        m_real[0] = PGAS_G(m_real,0);
        m_imag[0] = PGAS_G(m_imag,0);
    }
    grid.sync();
    ValType prob_of_one = m_real[0];
    ValType overall = m_imag[0];

    //assert (prob_of_one <= 1.0);

    bool val = (rand <= (prob_of_one/overall));
    if (val) // we get 1, so we set all entires with (id&mask==0) to 0, and scale entires with (id&mask==1) by factor
    {
        ValType factor = 1./prob_of_one; //we compute 1/sqrt(prob), so other entries can times this val
        for (IdxType i=tid; i<per_pe_work_dm; i+=blockDim.x*gridDim.x)
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

                //if (sv_real[i] < ERROR_BAR) sv_real[i] = 0;
                //else sv_real[i] *= factor;
                //if (sv_imag[i] < ERROR_BAR) sv_imag[i] = 0;
                //else sv_imag[i] *= factor;
            }
        }
    }
    else // we get 0, so we set all entires with (id&mask!=0) to 0, and scale entires with (id&mask==0) by factor
    {
        //assert(prob_of_one != 1);
        ValType factor =  1./(overall-prob_of_one); //we compute 1/sqrt(prob), so other entries can times this val
        for (IdxType i=tid; i<per_pe_work_dm; i+=blockDim.x*gridDim.x)
        {
            if ( (i & mask) == 0)
            {
                sv_real[i] *= factor;
                sv_imag[i] *= factor;

                //if (sv_real[i] < ERROR_BAR) sv_real[i] = 0;
                //else sv_real[i] *= factor;
                //if (sv_imag[i] < ERROR_BAR) sv_imag[i] = 0;
                //else sv_imag[i] *= factor;
            }
            else
            {
                sv_real[i] = 0;
                sv_imag[i] = 0;
            }
        }
    }

    if (tid==0) sim->results_gpu[0] = (rand<=prob_of_one?1:0);
}

*/




__device__ __inline__ void MA_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag,
        const IdxType repetition)
{
    grid_group grid = this_grid();
    assert(sim-> n_qubits < sim->lg2_m_gpu); //ensure the diagonal can be fit into GPU-0's part of m_real
    const IdxType n_size = (IdxType)1<<(sim->n_qubits);
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    ValType * m_real = sim->m_real;
    ValType * m_imag = sim->m_imag;

    BARRIER;

    if (sim->i_gpu == 0)
    {    
        for (IdxType i=tid; i<n_size; i+=blockDim.x*gridDim.x)
            m_real[i] = abs(PGAS_G(sv_real, (i<<(sim->n_qubits))+i));
        grid.sync();

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
            m_real[n_size] = m_real[n_size-1];
            m_real[n_size-1] = 0;
            ValType purity = fabs(m_real[n_size]);
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
        grid.sync();

        for (IdxType i=0; i<repetition; i++)
        {
            ValType r = sim->randoms_gpu[i];
            for (IdxType j=tid; j<n_size; j+=blockDim.x*gridDim.x)
            {
                if (m_real[j]<=r && r<m_real[j+1]) PGAS_P(m_imag, i, j);
            }
            grid.sync();
        }
    }
    BARRIER;
    for (IdxType i=tid; i<repetition; i+=blockDim.x*gridDim.x)
    {
        sim->results_gpu[i] = PGAS_G(m_imag, i);
    }
    BARRIER;
}


/*
__device__ __inline__ void MA_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag,
        const IdxType repetition)
{
    grid_group grid = this_grid();
    assert(sim-> n_qubits < sim->lg2_m_gpu); //ensure the diagonal can be fit into GPU-0's part of m_real
    const IdxType n_size = (IdxType)1<<(sim->n_qubits);
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    ValType * m_real = sim->m_real;
    ValType * m_imag = sim->m_imag;

    BARRIER;

    if (sim->i_gpu == 0)
    {    
        for (IdxType i=tid; i<n_size; i+=blockDim.x*gridDim.x)
        {
            ValType val = abs(PGAS_G(sv_real, (i<<(sim->n_qubits))+i));
            PGAS_P(m_real, i, val);
        }
        grid.sync();

        //Parallel prefix sum
        for (IdxType d=0; d<(sim->n_qubits); d++)
        {
            IdxType step = (IdxType)1<<(d+1);
            for (IdxType k=tid*step; k<n_size; k+=step*blockDim.x*gridDim.x)
            {
                ValType val = PGAS_G(m_real, k+(1<<d)-1) + PGAS_G(m_real, k+(1<<(d+1))-1);
                PGAS_P(m_real, k+(1<<(d+1))-1, val);
            }
            grid.sync();
        }
        if (tid == 0)
        {
            PGAS_P(m_real, n_size, PGAS_G(m_real, n_size-1));
            PGAS_P(m_real, n_size-1, 0);
            ValType purity = fabs(PGAS_G(m_real, n_size));
            if ( abs(purity - 1.0) > ERROR_BAR )
                printf("MA: Purity Check fails with %lf\n", purity);
        }
        grid.sync();
        for (IdxType d=(sim->n_qubits)-1; d>=0; d--)
        {
            IdxType step = (IdxType)1<<(d+1);
            for (IdxType k=tid*step; k<n_size-1; k+=step*blockDim.x*gridDim.x)
            {
                ValType tmp = PGAS_G(m_real, k+((IdxType)1<<d)-1);
                PGAS_P(m_real, k+((IdxType)1<<d)-1, PGAS_G(m_real, k+((IdxType)1<<(d+1))-1));
                PGAS_P(m_real, k+((IdxType)1<<(d+1))-1, tmp+PGAS_G(m_real, k+((IdxType)1<<(d+1))-1)); 

            }
            grid.sync();
        }
        grid.sync();

        for (IdxType i=0; i<repetition; i++)
        {
            ValType r = sim->randoms_gpu[i];
            for (IdxType j=tid; j<n_size; j+=blockDim.x*gridDim.x)
            {
                if (PGAS_G(m_real,j)<=r && r<PGAS_G(m_real,j+1)) PGAS_P(m_imag, i, j);
            }
            grid.sync();
        }
    }
    BARRIER;
    for (IdxType i=tid; i<repetition; i+=blockDim.x*gridDim.x)
    {
        sim->results_gpu[i] = PGAS_G(m_imag, i);
    }
    BARRIER;
}
*/


__device__ __inline__ void RESET_GATE(const Simulation* sim, ValType* sv_real, ValType* sv_imag,
        const IdxType qubit)
{
    grid_group grid = this_grid();
    assert(sim-> n_qubits < sim->lg2_m_gpu); //ensure the diagonal can be fit into GPU-0's part of m_real
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    const IdxType per_pe_work_dm = ((sim->dim)>>(sim->gpu_scale));
    ValType * m_real = sim->m_real;
    IdxType mask = ((IdxType)1<<qubit);

    if (sim->i_gpu == 0)
    {
        for (IdxType i=tid; i<((IdxType)1<<(sim->n_qubits)); i+=blockDim.x*gridDim.x)
        {
            if ( (i & mask) == 0) //for all conditions with qubit=0, we set it to 0, so we sum up all prob that qubit=1
            {
                m_real[i] = 0;
            }
            else
            {
                const ValType val = PGAS_G(sv_real, (i<<(sim->n_qubits))+i);
                PGAS_P(m_real,i,abs(val));
            }
        }
        grid.sync();
        for (IdxType k=((IdxType)1<<(sim->n_qubits-1)); k>0; k>>=1)
        {
            for (IdxType i=tid; i<k; i+=blockDim.x*gridDim.x) m_real[i] += m_real[i+k];
            grid.sync();
        }
    }
    BARRIER;
    if (tid==0 && sim->i_gpu!=0) m_real[0] = PGAS_G(m_real,0);
    grid.sync();
    ValType prob_of_one = m_real[0];
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
    BARRIER;
}


__device__ __inline__ void Purity_Check(const Simulation* sim, const IdxType t, ValType* sv_real, ValType* sv_imag)
{
    grid_group grid = this_grid();
    assert(sim-> n_qubits < sim->lg2_m_gpu); //ensure the diagonal can be fit into GPU-0's part of m_real
    const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
    const IdxType per_pe_work_dm = ((sim->dim)>>(sim->gpu_scale));
    ValType * m_real = sim->m_real;

    BARRIER;

    if (sim->i_gpu == 0)
    {
        for (IdxType i=tid; i<((IdxType)1<<(sim->n_qubits)); i+=blockDim.x*gridDim.x)
        {
            const ValType val = PGAS_G(sv_real, (i<<(sim->n_qubits))+i);
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

    }
    BARRIER;
}




//=====================================================================================
//Per-gate execution function
__device__ void Gate::exe_op(Simulation* sim, ValType* sv_real, ValType* sv_imag)
{
    grid_group grid = this_grid(); 
    //only need sync when operating on remote qubits
    if (((ctrl+sim->n_qubits)>=sim->lg2_m_gpu) || ((qubit+sim->n_qubits)>=sim->lg2_m_gpu))
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
    }
    else if (op_name == C2) 
    {
        C2V1_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, qubit, qubit+(sim->n_qubits));
    }
    else if (op_name == C4) 
    {

        if (((ctrl+sim->n_qubits)>=sim->lg2_m_gpu) && ((qubit+sim->n_qubits)>=sim->lg2_m_gpu))
        {
            SWAP_GATE(sim, sv_real, sv_imag, 0, ctrl+(sim->n_qubits));
            BARRIER;
#ifdef ENABLE_TENSOR_CORES
            C4TCV3_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, ctrl, qubit, 0, qubit+(sim->n_qubits));
#else
            C4V1_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, ctrl, qubit, 0, qubit+(sim->n_qubits));
#endif
            BARRIER;
            SWAP_GATE(sim, sv_real, sv_imag, 0, ctrl+(sim->n_qubits));
        }
        else
        {
#ifdef ENABLE_TENSOR_CORES
            C4TCV3_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, ctrl, qubit, ctrl+(sim->n_qubits), qubit+(sim->n_qubits));
#else
            C4V1_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, ctrl, qubit, ctrl+(sim->n_qubits), qubit+(sim->n_qubits));
#endif
        }

    }
    //only need sync when operating on remote qubits
    if (((ctrl+sim->n_qubits)>=sim->lg2_m_gpu) || ((qubit+sim->n_qubits)>=sim->lg2_m_gpu))
        if( threadIdx.x==0 && blockIdx.x==0 ) nvshmem_barrier_all(); 
    grid.sync();
}

}; //namespace NWQSim

#endif
