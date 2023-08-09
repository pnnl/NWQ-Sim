// ---------------------------------------------------------------------------
// File: svsim_amdgpu_mpi.cuh
// MPI based implementation of the scale-out SV-Sim gates and
// simulation runtime using AMD GPU backend.
// ---------------------------------------------------------------------------

#pragma once

#include "../state.hpp"
#include "../nwq_util.hpp"
#include "../gate.hpp"
#include "../circuit.hpp"

#include "../private/config.hpp"
#include "../private/macros.hpp"
#include "../private/sim_gate.hpp"
#include "../private/hip_util.hpp"
#include "../circuit_pass/fusion.hpp"

#include <assert.h>
#include <random>
#include <complex.h>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <memory>

#include <mpi.h>
#include <hip/hip_runtime.h>
#include <hip/hip_cooperative_groups.h>

namespace NWQSim
{
    using namespace cooperative_groups;
    using namespace std;

    class SV_HIP_MPI;
    void Purity_Check(const SV_HIP_MPI *sim, const IdxType t, const OP gate_op, const IdxType ctrl, const IdxType qubit);

    // SV_HIP_MPI runtime
    void simulation_kernel(SV_HIP_MPI *, std::vector<SVGate> &);

    class SV_HIP_MPI : public QuantumState
    {
    public:
        SV_HIP_MPI(IdxType _n_qubits) : QuantumState(_n_qubits)
        {

            n_qubits = _n_qubits;
            dim = ((IdxType)1 << (n_qubits));
            half_dim = ((IdxType)1 << (n_qubits - 1));
            sv_size = dim * (IdxType)sizeof(ValType);

            int mpi_size;
            MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
            n_gpus = (IdxType)mpi_size;
            int mpi_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
            i_proc = (IdxType)mpi_rank;

            // always be 0 since 1-MPI maps to 1-GPU
            hipSafeCall(hipSetDevice(0));
            gpu_scale = floor(log((double)n_gpus + 0.5) / log(2.0));
            lg2_m_gpu = n_qubits - gpu_scale;
            m_gpu = ((IdxType)1 << (lg2_m_gpu));
            sv_size_per_gpu = sv_size / n_gpus;
            // CPU side initialization
            if (!is_power_of_2(n_gpus))
                throw std::logic_error("Error: Number of GPUs should be power of 2.");
            if (dim % n_gpus != 0)
                throw std::logic_error("Error: Number of GPUs is too large or too small.");
            if (lg2_m_gpu < 5)
                throw std::logic_error("Error: Each GPU should have at least 5 qubits for multi-node version. Please increase qubits or reduce GPUs");
            // CPU side initialization
            SAFE_ALOC_HOST_HIP(sv_real_cpu, sv_size_per_gpu);
            SAFE_ALOC_HOST_HIP(sv_imag_cpu, sv_size_per_gpu);
            memset(sv_real_cpu, 0, sv_size_per_gpu);
            memset(sv_imag_cpu, 0, sv_size_per_gpu);
            // State-vector initial state [0..0] = 1
            if (i_proc == 0)
                sv_real_cpu[0] = 1.;
            // GPU memory allocation
            SAFE_ALOC_GPU_HIP(sv_real, sv_size_per_gpu);
            SAFE_ALOC_GPU_HIP(sv_imag, sv_size_per_gpu);
            SAFE_ALOC_GPU_HIP(m_real, sv_size_per_gpu);
            SAFE_ALOC_GPU_HIP(m_imag, sv_size_per_gpu);
            hipCheckError();
            gpu_mem = sv_size_per_gpu * 4;

            // GPU memory initilization
            hipSafeCall(hipMemcpy(sv_real, sv_real_cpu, sv_size_per_gpu,
                                  hipMemcpyHostToDevice));
            hipSafeCall(hipMemcpy(sv_imag, sv_imag_cpu, sv_size_per_gpu,
                                  hipMemcpyHostToDevice));
            hipSafeCall(hipMemset(m_real, 0, sv_size_per_gpu));
            hipSafeCall(hipMemset(m_imag, 0, sv_size_per_gpu));
            SAFE_ALOC_GPU_HIP(sim_gpu, sizeof(SV_HIP_MPI));
            rng.seed(time(0));
        }

        ~SV_HIP_MPI()
        {
            // Release for GPU side
            SAFE_FREE_GPU_HIP(sv_real);
            SAFE_FREE_GPU_HIP(sv_imag);
            SAFE_FREE_GPU_HIP(m_real);
            SAFE_FREE_GPU_HIP(m_imag);
            SAFE_FREE_GPU_HIP(sim_gpu);

            // Release for CPU side
            SAFE_FREE_HOST_HIP(sv_real_cpu);
            SAFE_FREE_HOST_HIP(sv_imag_cpu);
            SAFE_FREE_HOST_HIP(randoms);
            SAFE_FREE_HOST_HIP(results);
            SAFE_FREE_HOST_HIP(results_local);
        }

        void reset_state() override
        {
            // Reset CPU input & output
            memset(sv_real_cpu, 0, sv_size_per_gpu);
            memset(sv_imag_cpu, 0, sv_size_per_gpu);
            // State Vector initial state [0..0] = 1
            if (i_proc == 0)
                sv_real_cpu[0] = 1.;
            // GPU side initialization
            hipSafeCall(hipMemcpy(sv_real, sv_real_cpu,
                                  sv_size_per_gpu, hipMemcpyHostToDevice));
            hipSafeCall(hipMemcpy(sv_imag, sv_imag_cpu,
                                  sv_size_per_gpu, hipMemcpyHostToDevice));
            hipSafeCall(hipMemset(m_real, 0, sv_size_per_gpu));
            hipSafeCall(hipMemset(m_imag, 0, sv_size_per_gpu));
        }

        void set_seed(IdxType seed) override
        {
            rng.seed(seed);
        }

        void sim(std::shared_ptr<NWQSim::Circuit> circuit) override
        {

            IdxType origional_gates = circuit->num_gates();
            IdxType n_measures = prepare_measure(circuit->get_gates());
            std::vector<SVGate> cpu_vec = fuse_circuit_sv(circuit);

            // Copy the circuit to GPU
            copy_gates_to_gpu(cpu_vec);
            IdxType n_gates = cpu_vec.size();

            hipSafeCall(hipSetDevice(0));

            hipSafeCall(hipMemcpy(sim_gpu, this,
                                  sizeof(SV_HIP_MPI), hipMemcpyHostToDevice));

            double *sim_times;
            double sim_time;
            cpu_timer sim_timer;

            if (Config::PRINT_SIM_TRACE)
            {
                if (i_proc == 0)
                {
                    SAFE_ALOC_HOST_HIP(sim_times, sizeof(double) * n_gpus);
                    memset(sim_times, 0, sizeof(double) * n_gpus);
                    printf("SVSim_GPU is running! Requesting %lld qubits.\n", n_qubits);
                }

                MPI_Barrier(MPI_COMM_WORLD);
                sim_timer.start_timer();
            }

            //=========================================
            simulation_kernel(this, cpu_vec);
            //=========================================

            if (Config::PRINT_SIM_TRACE)
            {
                sim_timer.stop_timer();
                sim_time = sim_timer.measure();

                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Gather(&sim_time, 1, MPI_DOUBLE,
                           &sim_times[i_proc], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                if (i_proc == 0)
                {
                    double avg_sim_time = 0;
                    for (unsigned d = 0; d < n_gpus; d++)
                    {
                        avg_sim_time += sim_times[d];
                    }
                    avg_sim_time /= (double)n_gpus;
                    printf("\n============== SV-Sim ===============\n");
                    printf("nqubits:%lld, ngates:%lld, sim_gates:%lld, ngpus:%lld, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms, mem:%.3lf MB, mem_per_gpu:%.3lf MB\n",
                           n_qubits, input_gates, n_gates, n_gpus, avg_sim_time, 0.,
                           avg_sim_time, gpu_mem / 1024 / 1024 * n_gpus, gpu_mem / 1024 / 1024);
                    printf("=====================================\n");
                    SAFE_FREE_HOST_HIP(sim_times);
                }
            }
        }

        IdxType *get_results() override
        {
            return results;
        }

        IdxType measure(IdxType qubit) override
        {
            std::shared_ptr<Circuit> circuit = std::make_shared<Circuit>(n_qubits);
            circuit->M(qubit);
            sim(circuit);
            return results[0];
        }

        IdxType *measure_all(IdxType repetition) override
        {
            std::shared_ptr<Circuit> circuit = std::make_shared<Circuit>(n_qubits);
            circuit->MA(repetition);
            sim(circuit);
            return results;
        }
        ValType get_exp_z(const std::vector<size_t> &in_bits) override
        {
            throw std::runtime_error("Get EXP Z Not implemented yet (SV_HIP_MPI)!");
        }
        ValType get_exp_z() override
        {
            throw std::runtime_error("Get EXP Z Not implemented yet (SV_HIP_MPI)!");
        }

        void print_res_state() override
        {
            hipSafeCall(hipMemcpy(sv_real_cpu, sv_real, sv_size_per_gpu, hipMemcpyDeviceToHost));
            hipSafeCall(hipMemcpy(sv_imag_cpu, sv_imag, sv_size_per_gpu, hipMemcpyDeviceToHost));

            ValType *sv_diag_real = NULL;
            ValType *sv_diag_imag = NULL;
            if (i_proc == 0)
                SAFE_ALOC_HOST_HIP(sv_diag_real, dim * sizeof(ValType));
            if (i_proc == 0)
                SAFE_ALOC_HOST_HIP(sv_diag_imag, dim * sizeof(ValType));

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Gather(sv_real_cpu, m_gpu, MPI_DOUBLE,
                       &sv_diag_real[i_proc * m_gpu], m_gpu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gather(sv_imag_cpu, m_gpu, MPI_DOUBLE,
                       &sv_diag_imag[i_proc * m_gpu], m_gpu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            if (i_proc == 0)
            {
                IdxType num = ((IdxType)1 << n_qubits);
                printf("----- SVSim ------\n");
                for (IdxType i = 0; i < num; i++)
                {
                    printf("(%.3lf,%.3lfj) ", sv_diag_real[i], sv_diag_imag[i]);
                    if ((i + 1) % 8 == 0)
                        printf("\n");
                }
                printf("\n");
                SAFE_FREE_HOST_HIP(sv_diag_real);
                SAFE_FREE_HOST_HIP(sv_diag_imag);
            }
        }

    public:
        // n_qubits is the number of qubits
        IdxType n_qubits;

        // which gpu
        IdxType gpu_scale;
        IdxType n_gpus;
        IdxType lg2_m_gpu;
        IdxType m_gpu;
        IdxType sv_size_per_gpu;

        // gpu_scale is 2^x of the number of GPUs, e.g., with 8 GPUs the gpu_scale is 3 (2^3=8)
        IdxType dim;
        IdxType half_dim;
        IdxType sv_size;
        // CPU arrays
        ValType *sv_real_cpu;
        ValType *sv_imag_cpu;
        // GPU arrays
        ValType *sv_real;
        ValType *sv_imag;
        // GPU Buffer
        ValType *m_real;
        ValType *m_imag;
        // For measurement randoms
        ValType *randoms = NULL;
        // For measurement result
        IdxType *results = NULL;
        // Local measurement result for MA
        IdxType *results_local = NULL;
        // Random
        std::mt19937 rng;
        std::uniform_real_distribution<ValType> uni_dist;
        // GPU memory usage
        ValType gpu_mem;

        // GPU-side gates instance
        SVGate *gates_gpu = NULL;

        void copy_gates_to_gpu(std::vector<SVGate> &cpu_vec)
        {
            // Allocate memory on CPU
            size_t vec_size = cpu_vec.size() * sizeof(SVGate);

            // Allocate memory on GPU
            SAFE_FREE_GPU_HIP(gates_gpu);
            SAFE_ALOC_GPU_HIP(gates_gpu, vec_size);
            hipSafeCall(hipMemcpy(gates_gpu, cpu_vec.data(), vec_size, hipMemcpyHostToDevice));
        }

        IdxType prepare_measure(std::vector<Gate> gates)
        {
            // Determine the total number of measurements
            IdxType n_slots = 0;
            for (auto g : gates)
            {
                if (g.op_name == OP::M)
                    n_slots++;
                else if (g.op_name == OP::MA)
                    n_slots += g.repetition;
            }

            // Prepare randoms and results memory
            SAFE_FREE_HOST_HIP(results);
            SAFE_ALOC_HOST_HIP(results, sizeof(IdxType) * n_slots);
            memset(results, 0, sizeof(IdxType) * n_slots);

            SAFE_FREE_GPU_HIP(results_gpu);
            SAFE_ALOC_GPU_HIP(results_gpu, sizeof(IdxType) * n_slots);
            hipSafeCall(hipMemset(results_gpu, 0, sizeof(IdxType) * n_slots));

            SAFE_FREE_HOST_HIP(randoms);
            SAFE_ALOC_HOST_HIP(randoms, sizeof(ValType) * n_slots);
            for (IdxType i = 0; i < n_slots; i++)
                randoms[i] = uni_dist(rng);

            SAFE_FREE_GPU_HIP(randoms_gpu);
            SAFE_ALOC_GPU_HIP(randoms_gpu, sizeof(ValType) * n_slots);
            hipSafeCall(hipMemcpy(randoms_gpu, randoms,
                                  sizeof(ValType) * n_slots, hipMemcpyHostToDevice));

            return n_slots;
        }
    };

#define LAUNCH_KERNEL(X)                                                                        \
    dim3 gridDim(1, 1, 1);                                                                      \
    hipDeviceProp_t deviceProp;                                                                 \
    unsigned smem_size = 0;                                                                     \
    hipSafeCall(hipGetDeviceProperties(&deviceProp, 0));                                        \
    int numBlocksPerSm;                                                                         \
    hipSafeCall(hipOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm,                   \
                                                             (X), THREADS_CTA_HIP, smem_size)); \
    gridDim.x = numBlocksPerSm * deviceProp.multiProcessorCount;                                \
    hipSafeCall(hipLaunchCooperativeKernel((void *)(X), gridDim,                                \
                                           THREADS_CTA_HIP, args, smem_size, 0));               \
    hipCheckError();                                                                            \
    hipSafeCall(hipDeviceSynchronize());

    //============== Local Unified 1-qubit Gate ================
    __global__ void C1LC_GATE(SV_HIP_MPI *sim, const IdxType t)
    {
        ValType *sv_real = sim->sv_real;
        ValType *sv_imag = sim->sv_imag;

        const ValType *gm_real = sim->gates_gpu[t].gm_real;
        const ValType *gm_imag = sim->gates_gpu[t].gm_imag;

        const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
        const IdxType per_pe_work = ((sim->half_dim) >> (sim->gpu_scale));
        IdxType qubit = sim->gates_gpu[t].qubit;

        for (IdxType i = (sim->i_proc) * per_pe_work + tid; i < (sim->i_proc + 1) * per_pe_work; i += blockDim.x * gridDim.x)
        {
            IdxType outer = (i >> qubit);
            IdxType inner = (i & (((IdxType)1 << qubit) - 1));
            IdxType offset = (outer << (qubit + 1));
            IdxType pos0 = ((offset + inner) & (sim->m_gpu - 1));
            IdxType pos1 = ((offset + inner + ((IdxType)1 << qubit)) & (sim->m_gpu - 1));
            const ValType el0_real = LOCAL_G_HIP_MPI(sv_real, pos0);
            const ValType el0_imag = LOCAL_G_HIP_MPI(sv_imag, pos0);
            const ValType el1_real = LOCAL_G_HIP_MPI(sv_real, pos1);
            const ValType el1_imag = LOCAL_G_HIP_MPI(sv_imag, pos1);

            ValType sv_real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag);
            ValType sv_imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real);
            ValType sv_real_pos1 = (gm_real[2] * el0_real) - (gm_imag[2] * el0_imag) + (gm_real[3] * el1_real) - (gm_imag[3] * el1_imag);
            ValType sv_imag_pos1 = (gm_real[2] * el0_imag) + (gm_imag[2] * el0_real) + (gm_real[3] * el1_imag) + (gm_imag[3] * el1_real);
            LOCAL_P_HIP_MPI(sv_real, pos0, sv_real_pos0);
            LOCAL_P_HIP_MPI(sv_imag, pos0, sv_imag_pos0);
            LOCAL_P_HIP_MPI(sv_real, pos1, sv_real_pos1);
            LOCAL_P_HIP_MPI(sv_imag, pos1, sv_imag_pos1);
        }
    }

    __global__ void C1RM_GATE(SV_HIP_MPI *sim, const IdxType t)
    {
        ValType *sv_real = sim->sv_real;
        ValType *sv_imag = sim->sv_imag;
        ValType *sv_real_remote = sim->m_real;
        ValType *sv_imag_remote = sim->m_imag;
        const ValType *gm_real = sim->gates_gpu[t].gm_real;
        const ValType *gm_imag = sim->gates_gpu[t].gm_imag;
        const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
        const IdxType per_pe_work = ((sim->dim) >> (sim->gpu_scale));
        for (IdxType i = tid; i < per_pe_work; i += blockDim.x * gridDim.x)
        {
            const ValType el0_real = LOCAL_G_HIP_MPI(sv_real, i);
            const ValType el0_imag = LOCAL_G_HIP_MPI(sv_imag, i);
            const ValType el1_real = LOCAL_G_HIP_MPI(sv_real_remote, i);
            const ValType el1_imag = LOCAL_G_HIP_MPI(sv_imag_remote, i);

            ValType real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag);
            ValType imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real);
            ValType real_pos1 = (gm_real[2] * el0_real) - (gm_imag[2] * el0_imag) + (gm_real[3] * el1_real) - (gm_imag[3] * el1_imag);
            ValType imag_pos1 = (gm_real[2] * el0_imag) + (gm_imag[2] * el0_real) + (gm_real[3] * el1_imag) + (gm_imag[3] * el1_real);
            LOCAL_P_HIP_MPI(sv_real, i, real_pos0);
            LOCAL_P_HIP_MPI(sv_imag, i, imag_pos0);
            LOCAL_P_HIP_MPI(sv_real_remote, i, real_pos1);
            LOCAL_P_HIP_MPI(sv_imag_remote, i, imag_pos1);
        }
    }

    void C1_GATE(const SV_HIP_MPI *sim, const IdxType t, const IdxType qubit)
    {
        ValType *sv_real = sim->sv_real;
        ValType *sv_imag = sim->sv_imag;
        ValType *m_real = sim->m_real;
        ValType *m_imag = sim->m_imag;

        if (qubit < sim->lg2_m_gpu)
        {
            void *args[] = {(void *)(&sim->sim_gpu), (void *)&t};
            LAUNCH_KERNEL(C1LC_GATE);
        }
        else
        {
            const IdxType per_pe_work = ((sim->dim) >> (sim->gpu_scale));
            IdxType pair_gpu = (sim->i_proc) ^ ((IdxType)1 << (qubit - (sim->lg2_m_gpu)));
            MPI_Barrier(MPI_COMM_WORLD);

            if (sim->i_proc > pair_gpu)
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

            if (sim->i_proc < pair_gpu)
            {
                void *args[] = {(void *)(&sim->sim_gpu), (void *)&t};
                LAUNCH_KERNEL(C1RM_GATE);
            }
            MPI_Barrier(MPI_COMM_WORLD);

            if (sim->i_proc > pair_gpu)
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
    __global__ void C2LC_GATE(const SV_HIP_MPI *sim, const IdxType qubit0,
                              const IdxType qubit1, const IdxType t)
    {
        ValType *sv_real = sim->sv_real;
        ValType *sv_imag = sim->sv_imag;
        const ValType *gm_real = sim->gates_gpu[t].gm_real;
        const ValType *gm_imag = sim->gates_gpu[t].gm_imag;
        const int tid = blockDim.x * blockIdx.x + threadIdx.x;
        const IdxType per_pe_work = ((sim->dim) >> (sim->gpu_scale + 2));
        assert(qubit0 != qubit1); // Non-cloning

        const IdxType q0dim = ((IdxType)1 << max(qubit0, qubit1));
        const IdxType q1dim = ((IdxType)1 << min(qubit0, qubit1));
        const IdxType outer_factor = ((sim->dim) + q0dim + q0dim - 1) >> (max(qubit0, qubit1) + 1);
        const IdxType mider_factor = (q0dim + q1dim + q1dim - 1) >> (min(qubit0, qubit1) + 1);
        const IdxType inner_factor = q1dim;
        const IdxType qubit0_dim = ((IdxType)1 << qubit0);
        const IdxType qubit1_dim = ((IdxType)1 << qubit1);

        for (IdxType i = (sim->i_proc) * per_pe_work + tid; i < (sim->i_proc + 1) * per_pe_work;
             i += blockDim.x * gridDim.x)
        {
            IdxType outer = ((i / inner_factor) / (mider_factor)) * (q0dim + q0dim);
            IdxType mider = ((i / inner_factor) % (mider_factor)) * (q1dim + q1dim);
            IdxType inner = i % inner_factor;
            IdxType pos0 = outer + mider + inner;
            IdxType pos1 = outer + mider + inner + qubit1_dim;
            IdxType pos2 = outer + mider + inner + qubit0_dim;
            IdxType pos3 = outer + mider + inner + q0dim + q1dim;

            const ValType el0_real = LOCAL_G_HIP_MPI(sv_real, pos0);
            const ValType el0_imag = LOCAL_G_HIP_MPI(sv_imag, pos0);
            const ValType el1_real = LOCAL_G_HIP_MPI(sv_real, pos1);
            const ValType el1_imag = LOCAL_G_HIP_MPI(sv_imag, pos1);
            const ValType el2_real = LOCAL_G_HIP_MPI(sv_real, pos2);
            const ValType el2_imag = LOCAL_G_HIP_MPI(sv_imag, pos2);
            const ValType el3_real = LOCAL_G_HIP_MPI(sv_real, pos3);
            const ValType el3_imag = LOCAL_G_HIP_MPI(sv_imag, pos3);

            // Real part
            ValType sv_real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag) + (gm_real[2] * el2_real) - (gm_imag[2] * el2_imag) + (gm_real[3] * el3_real) - (gm_imag[3] * el3_imag);
            ValType sv_real_pos1 = (gm_real[4] * el0_real) - (gm_imag[4] * el0_imag) + (gm_real[5] * el1_real) - (gm_imag[5] * el1_imag) + (gm_real[6] * el2_real) - (gm_imag[6] * el2_imag) + (gm_real[7] * el3_real) - (gm_imag[7] * el3_imag);
            ValType sv_real_pos2 = (gm_real[8] * el0_real) - (gm_imag[8] * el0_imag) + (gm_real[9] * el1_real) - (gm_imag[9] * el1_imag) + (gm_real[10] * el2_real) - (gm_imag[10] * el2_imag) + (gm_real[11] * el3_real) - (gm_imag[11] * el3_imag);
            ValType sv_real_pos3 = (gm_real[12] * el0_real) - (gm_imag[12] * el0_imag) + (gm_real[13] * el1_real) - (gm_imag[13] * el1_imag) + (gm_real[14] * el2_real) - (gm_imag[14] * el2_imag) + (gm_real[15] * el3_real) - (gm_imag[15] * el3_imag);

            // Imag part
            ValType sv_imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real) + (gm_real[2] * el2_imag) + (gm_imag[2] * el2_real) + (gm_real[3] * el3_imag) + (gm_imag[3] * el3_real);
            ValType sv_imag_pos1 = (gm_real[4] * el0_imag) + (gm_imag[4] * el0_real) + (gm_real[5] * el1_imag) + (gm_imag[5] * el1_real) + (gm_real[6] * el2_imag) + (gm_imag[6] * el2_real) + (gm_real[7] * el3_imag) + (gm_imag[7] * el3_real);
            ValType sv_imag_pos2 = (gm_real[8] * el0_imag) + (gm_imag[8] * el0_real) + (gm_real[9] * el1_imag) + (gm_imag[9] * el1_real) + (gm_real[10] * el2_imag) + (gm_imag[10] * el2_real) + (gm_real[11] * el3_imag) + (gm_imag[11] * el3_real);
            ValType sv_imag_pos3 = (gm_real[12] * el0_imag) + (gm_imag[12] * el0_real) + (gm_real[13] * el1_imag) + (gm_imag[13] * el1_real) + (gm_real[14] * el2_imag) + (gm_imag[14] * el2_real) + (gm_real[15] * el3_imag) + (gm_imag[15] * el3_real);

            LOCAL_P_HIP_MPI(sv_real, pos0, sv_real_pos0);
            LOCAL_P_HIP_MPI(sv_real, pos1, sv_real_pos1);
            LOCAL_P_HIP_MPI(sv_real, pos2, sv_real_pos2);
            LOCAL_P_HIP_MPI(sv_real, pos3, sv_real_pos3);

            LOCAL_P_HIP_MPI(sv_imag, pos0, sv_imag_pos0);
            LOCAL_P_HIP_MPI(sv_imag, pos1, sv_imag_pos1);
            LOCAL_P_HIP_MPI(sv_imag, pos2, sv_imag_pos2);
            LOCAL_P_HIP_MPI(sv_imag, pos3, sv_imag_pos3);
        }
    }

    //============== Unified 2-qubit Gate ================
    // Perform communication optimization here
    __global__ void C2RM_GATE(const SV_HIP_MPI *sim, const IdxType qubit0,
                              const IdxType qubit1, const IdxType t)
    {
        ValType *sv_real = sim->sv_real;
        ValType *sv_imag = sim->sv_imag;
        ValType *sv_real_remote = sim->m_real;
        ValType *sv_imag_remote = sim->m_imag;
        const ValType *gm_real = sim->gates_gpu[t].gm_real;
        const ValType *gm_imag = sim->gates_gpu[t].gm_imag;
        const int tid = blockDim.x * blockIdx.x + threadIdx.x;
        const IdxType per_pe_work = ((sim->dim) >> (sim->gpu_scale + 1));
        const IdxType per_pe_num = ((sim->dim) >> (sim->gpu_scale));
        const IdxType p = min(qubit0, qubit1);
        const IdxType q = max(qubit0, qubit1);

        for (IdxType i = (sim->i_proc) * per_pe_work + tid; i < (sim->i_proc + 1) * per_pe_work;
             i += blockDim.x * gridDim.x)
        {
            ValType el_real[4];
            ValType el_imag[4];
            ValType res_real[4] = {0};
            ValType res_imag[4] = {0};
            const IdxType term0 = MOD2E(i, p);
            const IdxType term1 = MOD2E(DIV2E(i, p), q - p - 1) * EXP2E(p + 1);
            const IdxType term2 = DIV2E(DIV2E(i, p), q - p - 1) * EXP2E(q + 1);
            const IdxType term = term2 + term1 + term0;
            el_real[0] = LOCAL_G_HIP_MPI(sv_real, term + SV4IDX(0));
            el_imag[0] = LOCAL_G_HIP_MPI(sv_imag, term + SV4IDX(0));
            el_real[3] = LOCAL_G_HIP_MPI(sv_real_remote, term + SV4IDX(3));
            el_imag[3] = LOCAL_G_HIP_MPI(sv_imag_remote, term + SV4IDX(3));

            if (qubit0 == q) // qubit0 is the remote qubit
            {
                el_real[1] = LOCAL_G_HIP_MPI(sv_real, term + SV4IDX(1));
                el_imag[1] = LOCAL_G_HIP_MPI(sv_imag, term + SV4IDX(1));
                el_real[2] = LOCAL_G_HIP_MPI(sv_real_remote, term + SV4IDX(2));
                el_imag[2] = LOCAL_G_HIP_MPI(sv_imag_remote, term + SV4IDX(2));
            }
            else // qubit1 is the remote qubit
            {
                el_real[1] = LOCAL_G_HIP_MPI(sv_real_remote, term + SV4IDX(1));
                el_imag[1] = LOCAL_G_HIP_MPI(sv_imag_remote, term + SV4IDX(1));
                el_real[2] = LOCAL_G_HIP_MPI(sv_real, term + SV4IDX(2));
                el_imag[2] = LOCAL_G_HIP_MPI(sv_imag, term + SV4IDX(2));
            }
#pragma unroll
            for (unsigned j = 0; j < 4; j++)
            {
#pragma unroll
                for (unsigned k = 0; k < 4; k++)
                {
                    res_real[j] += (el_real[k] * gm_real[j * 4 + k]) - (el_imag[k] * gm_imag[j * 4 + k]);
                    res_imag[j] += (el_real[k] * gm_imag[j * 4 + k]) + (el_imag[k] * gm_real[j * 4 + k]);
                }
            }
            LOCAL_P_HIP_MPI(sv_real, term + SV4IDX(0), res_real[0]);
            LOCAL_P_HIP_MPI(sv_imag, term + SV4IDX(0), res_imag[0]);
            LOCAL_P_HIP_MPI(sv_real_remote, term + SV4IDX(3), res_real[3]);
            LOCAL_P_HIP_MPI(sv_imag_remote, term + SV4IDX(3), res_imag[3]);

            if (qubit0 == q) // qubit0 is the remote qubit
            {
                LOCAL_P_HIP_MPI(sv_real, term + SV4IDX(1), res_real[1]);
                LOCAL_P_HIP_MPI(sv_imag, term + SV4IDX(1), res_imag[1]);
                LOCAL_P_HIP_MPI(sv_real_remote, term + SV4IDX(2), res_real[2]);
                LOCAL_P_HIP_MPI(sv_imag_remote, term + SV4IDX(2), res_imag[2]);
            }
            else // qubit1 is the remote qubit
            {
                LOCAL_P_HIP_MPI(sv_real_remote, term + SV4IDX(1), res_real[1]);
                LOCAL_P_HIP_MPI(sv_imag_remote, term + SV4IDX(1), res_imag[1]);
                LOCAL_P_HIP_MPI(sv_real, term + SV4IDX(2), res_real[2]);
                LOCAL_P_HIP_MPI(sv_imag, term + SV4IDX(2), res_imag[2]);
            }
        }
    }

    void C2RMM_GATE(const SV_HIP_MPI *sim, const IdxType qubit0, const IdxType qubit1, const IdxType t)
    {
        const IdxType per_pe_work = ((sim->dim) >> (sim->gpu_scale));
        const IdxType p = min(qubit0, qubit1);
        const IdxType q = max(qubit0, qubit1);
        // load data from pair GPU
        IdxType pair_gpu = (sim->i_proc) ^ ((IdxType)1 << (q - (sim->lg2_m_gpu)));
        MPI_Barrier(MPI_COMM_WORLD);

        if (sim->i_proc > pair_gpu)
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

        if (sim->i_proc < pair_gpu)
        {
            void *args[] = {(void *)(&sim->sim_gpu), (void *)&qubit0, (void *)&qubit1, (void *)&t};
            LAUNCH_KERNEL(C2RM_GATE);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (sim->i_proc > pair_gpu)
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
    // This gate is for internal usage. It is used
    // when ctrl and target qubits are remote qubis, we then
    // swap one of them to a local qubit,
    // perform the C2 gate, and then swap back
    // It is assumed qubit0 is local, qubit1 is remote

    __global__ void SWAPLC_GATE(const SV_HIP_MPI *sim, const IdxType qubit0, const IdxType qubit1)
    {
        assert(qubit0 != qubit1); // Non-cloning
        const int tid = blockDim.x * blockIdx.x + threadIdx.x;
        const IdxType per_pe_work = ((sim->dim) >> (sim->gpu_scale + 1));
        ValType *sv_real = sim->sv_real;
        ValType *sv_imag = sim->sv_imag;
        ValType *sv_real_remote = sim->m_real;
        ValType *sv_imag_remote = sim->m_imag;
        const IdxType p = min(qubit0, qubit1);
        const IdxType q = max(qubit0, qubit1);

        for (IdxType i = (sim->i_proc) * per_pe_work + tid; i < (sim->i_proc + 1) * per_pe_work;
             i += blockDim.x * gridDim.x)
        {
            ValType el_real[4];
            ValType el_imag[4];
            ValType res_real[4] = {0};
            ValType res_imag[4] = {0};
            const IdxType term0 = MOD2E(i, p);
            const IdxType term1 = MOD2E(DIV2E(i, p), q - p - 1) * EXP2E(p + 1);
            const IdxType term2 = DIV2E(DIV2E(i, p), q - p - 1) * EXP2E(q + 1);
            const IdxType term = term2 + term1 + term0;

            el_real[1] = LOCAL_G_HIP_MPI(sv_real_remote, term + SV4IDX(1));
            el_imag[1] = LOCAL_G_HIP_MPI(sv_imag_remote, term + SV4IDX(1));
            el_real[2] = LOCAL_G_HIP_MPI(sv_real, term + SV4IDX(2));
            el_imag[2] = LOCAL_G_HIP_MPI(sv_imag, term + SV4IDX(2));

            res_real[1] = el_real[2];
            res_imag[1] = el_imag[2];
            res_real[2] = el_real[1];
            res_imag[2] = el_imag[1];

            LOCAL_P_HIP_MPI(sv_real_remote, term + SV4IDX(1), res_real[1]);
            LOCAL_P_HIP_MPI(sv_imag_remote, term + SV4IDX(1), res_imag[1]);
            LOCAL_P_HIP_MPI(sv_real, term + SV4IDX(2), res_real[2]);
            LOCAL_P_HIP_MPI(sv_imag, term + SV4IDX(2), res_imag[2]);
        }
    }

    void SWAP_GATE(const SV_HIP_MPI *sim, const IdxType qubit0, const IdxType qubit1)
    {
        assert(qubit0 != qubit1); // Non-cloning
        const IdxType per_pe_num = ((sim->dim) >> (sim->gpu_scale));
        const IdxType p = min(qubit0, qubit1);
        const IdxType q = max(qubit0, qubit1);

        // load data from pair GPU
        IdxType pair_gpu = (sim->i_proc) ^ ((IdxType)1 << (q - (sim->lg2_m_gpu)));
        MPI_Barrier(MPI_COMM_WORLD);

        if (sim->i_proc > pair_gpu)
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

        if (sim->i_proc < pair_gpu)
        {
            void *args[] = {(void *)(&sim->sim_gpu), (void *)&qubit0, (void *)&qubit1};
            LAUNCH_KERNEL(SWAPLC_GATE);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (sim->i_proc > pair_gpu)
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

    void C2_GATE(const SV_HIP_MPI *sim, const IdxType t, const IdxType ctrl, const IdxType qubit)
    {

        assert(ctrl != qubit); // Non-cloning
        // Case-1: both ctrl and qubit are local

        if (ctrl < sim->lg2_m_gpu && qubit < sim->lg2_m_gpu)
        {
            void *args[] = {(void *)(&sim->sim_gpu), &ctrl, &qubit, (void *)&t};
            LAUNCH_KERNEL(C2LC_GATE);
        }

        // Case-2: both ctrl and qubit are global
        else if (ctrl >= sim->lg2_m_gpu && qubit >= sim->lg2_m_gpu)
        {
            SWAP_GATE(sim, 0, ctrl);
            MPI_Barrier(MPI_COMM_WORLD);
            C2RMM_GATE(sim, 0, qubit, t);
            MPI_Barrier(MPI_COMM_WORLD);
            SWAP_GATE(sim, 0, ctrl);
        }
        // Case-3: one of ctrl and qubit is global
        else
        {
            C2RMM_GATE(sim, ctrl, qubit, t);
        }
    }

    __global__ void MLC1_GATE(const SV_HIP_MPI *sim, const IdxType t)
    {
        grid_group grid = this_grid();
        const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
        const IdxType qubit = sim->gates_gpu[t].qubit;
        const IdxType per_pe_work = ((sim->dim) >> (sim->gpu_scale));
        IdxType mask = ((IdxType)1 << qubit);
        ValType *sv_real = sim->sv_real;
        ValType *sv_imag = sim->sv_imag;
        ValType *m_real = sim->m_real;
        for (IdxType i = tid; i < sim->m_gpu; i += blockDim.x * gridDim.x)
        {
            IdxType idx = (sim->i_proc) * per_pe_work + i;
            if ((idx & mask) == 0)
                m_real[i] = 0;
            else
                m_real[i] = sv_real[i] * sv_real[i] + sv_imag[i] * sv_imag[i];
        }
        grid.sync();
        // Parallel reduction
        for (IdxType k = (sim->m_gpu >> (IdxType)1); k > 0; k >>= 1)
        {
            for (IdxType i = tid; i < k; i += blockDim.x * gridDim.x)
            {
                m_real[i] += m_real[i + k];
            }
            grid.sync();
        }
    }

    __global__ void MLC2_GATE(const SV_HIP_MPI *sim, const ValType prob_of_one, const IdxType t, const IdxType cur_index)
    {
        grid_group grid = this_grid();
        const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
        const IdxType qubit = sim->gates_gpu[t].qubit;
        const IdxType per_pe_work = ((sim->dim) >> (sim->gpu_scale));
        ValType *sv_real = sim->sv_real;
        ValType *sv_imag = sim->sv_imag;
        ValType *m_real = sim->m_real;
        IdxType mask = ((IdxType)1 << qubit);
        const ValType rand = sim->randoms_gpu[cur_index];
        if (rand < prob_of_one)
        {
            ValType factor = 1. / sqrt(prob_of_one);
            for (IdxType i = tid; i < sim->m_gpu; i += blockDim.x * gridDim.x)
            {
                IdxType idx = (sim->i_proc) * per_pe_work + i;
                if ((idx & mask) == 0)
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
            ValType factor = 1. / sqrt(1. - prob_of_one);
            for (IdxType i = tid; i < sim->m_gpu; i += blockDim.x * gridDim.x)
            {
                IdxType idx = (sim->i_proc) * per_pe_work + i;
                if ((idx & mask) == 0)
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

    void M_GATE(const SV_HIP_MPI *sim, const IdxType t, const IdxType cur_index)
    {
        {
            void *args[] = {(void *)(&sim->sim_gpu), (void *)&t};
            LAUNCH_KERNEL(MLC1_GATE);
        }
        ValType sum = 0;
        ValType prob_of_one = 0;
        hipSafeCall(hipMemcpy(&sum, &(sim->m_real[0]), sizeof(ValType), hipMemcpyHostToDevice));
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&sum, &prob_of_one, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        {
            void *args[] = {(void *)(&sim->sim_gpu), &prob_of_one, (void *)&t, (void *)&cur_index};
            LAUNCH_KERNEL(MLC2_GATE);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    void MA_GATE(const SV_HIP_MPI *sim, const IdxType repetition, const IdxType cur_index)
    {
        const IdxType sv_size_per_gpu = sim->sv_size_per_gpu;

        ValType *sv_real_cpu = sim->sv_real_cpu;
        ValType *sv_imag_cpu = sim->sv_imag_cpu;
        hipSafeCall(hipMemcpy(sv_real_cpu, sim->sv_real, sv_size_per_gpu, hipMemcpyDeviceToHost));
        hipSafeCall(hipMemcpy(sv_imag_cpu, sim->sv_imag, sv_size_per_gpu, hipMemcpyDeviceToHost));

        ValType *m_buff = NULL;
        SAFE_ALOC_HOST_HIP(m_buff, sv_size_per_gpu);
        memset(m_buff, 0, sv_size_per_gpu);
        const IdxType n_size = (IdxType)1 << (sim->n_qubits);

        ValType reduce = 0;
        ValType partial = 0;
        for (IdxType i = 0; i < sim->m_gpu; i++)
            reduce += sv_real_cpu[i] * sv_real_cpu[i] + sv_imag_cpu[i] * sv_imag_cpu[i];

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Scan(&reduce, &partial, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        if (sim->i_proc == (sim->n_gpus - 1)) // last node
        {
            ValType purity = fabs(partial);
            if (fabs(purity - 1.0) > ERROR_BAR)
                printf("MA: Purity Check fails with %lf\n", purity);
        }

        m_buff[0] = (partial - reduce); // per-node incremental val
        for (IdxType i = 1; i < sim->m_gpu; i++)
            m_buff[i] = m_buff[i - 1] + ((sv_real_cpu[i - 1] * sv_real_cpu[i - 1]) + (sv_imag_cpu[i - 1] * sv_imag_cpu[i - 1]));

        for (IdxType j = 0; j < n_size; j++)
        {
            IdxType local_cpu = j >> (sim->lg2_m_gpu);
            IdxType local_j = j & (sim->m_gpu - 1);
            if (local_cpu == sim->i_proc)
            {
                ValType lower = LOCAL_G_HIP_MPI(m_buff, j);
                ValType upper = 0;
                if (j + 1 == n_size)
                    upper = 1.0; // last element
                else
                    upper = (local_j + 1 == sim->m_gpu) ? partial : m_buff[local_j + 1]; // last element per node

                for (IdxType i = 0; i < repetition; i++)
                {
                    ValType r = sim->randoms[cur_index + i];
                    // all nodes store partial results locally. since other entires are all
                    // zeros, we can do all-reduce to merge
                    if (lower <= r && r < upper)
                        sim->results_local[cur_index + i] = j;
                }
            }
        }
        SAFE_FREE_HOST_HIP(m_buff);
    }

    __global__ void RESETLC1_GATE(const SV_HIP_MPI *sim, const IdxType t)
    {
        grid_group grid = this_grid();
        const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
        const IdxType qubit = sim->gates_gpu[t].qubit;
        IdxType mask = ((IdxType)1 << qubit);
        const IdxType per_pe_work = ((sim->dim) >> (sim->gpu_scale));
        ValType *sv_real = sim->sv_real;
        ValType *sv_imag = sim->sv_imag;
        ValType *m_real = sim->m_real;
        for (IdxType i = tid; i < sim->m_gpu; i += blockDim.x * gridDim.x)
        {
            IdxType idx = (sim->i_proc) * per_pe_work + i;
            if ((idx & mask) == 0)
                m_real[i] = 0;
            else
                m_real[i] = sv_real[i] * sv_real[i] + sv_imag[i] * sv_imag[i];
        }
        grid.sync();
        // Parallel reduction
        for (IdxType k = (sim->m_gpu >> (IdxType)1); k > 0; k >>= 1)
        {
            for (IdxType i = tid; i < k; i += blockDim.x * gridDim.x)
            {
                m_real[i] += m_real[i + k];
            }
            grid.sync();
        }
    }

    __global__ void RESETLC2_GATE(const SV_HIP_MPI *sim, const ValType prob_of_one, const IdxType t)
    {
        const IdxType qubit = sim->gates_gpu[t].qubit;
        const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
        const IdxType per_pe_work = ((sim->dim) >> (sim->gpu_scale));
        ValType *sv_real = sim->sv_real;
        ValType *sv_imag = sim->sv_imag;
        ValType factor = 1.0 / sqrt(1.0 - prob_of_one);
        IdxType mask = ((IdxType)1 << qubit);
        for (IdxType i = tid; i < sim->m_gpu; i += blockDim.x * gridDim.x)
        {
            IdxType idx = (sim->i_proc) * per_pe_work + i;
            if ((idx & mask) == 0)
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

    __global__ void RESETLC3_GATE(const SV_HIP_MPI *sim, const IdxType t)
    {
        const IdxType qubit = sim->gates_gpu[t].qubit;
        const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
        const IdxType per_pe_work = ((sim->dim) >> (sim->gpu_scale));
        ValType *sv_real = sim->sv_real;
        ValType *sv_imag = sim->sv_imag;
        IdxType mask = ((IdxType)1 << qubit);

        for (IdxType i = tid; i < sim->m_gpu; i += blockDim.x * gridDim.x)
        {
            if ((i & mask) == 0)
            {
                IdxType dual_i = i ^ mask;
                sv_real[i] = sv_real[dual_i];
                sv_imag[i] = sv_imag[dual_i];
                sv_real[dual_i] = 0;
                sv_imag[dual_i] = 0;
            }
        }
    }

    void RESET_GATE(const SV_HIP_MPI *sim, const IdxType t, const IdxType qubit)
    {

        const IdxType per_pe_work = ((sim->dim) >> (sim->gpu_scale));
        ValType *sv_real = sim->sv_real;
        ValType *sv_imag = sim->sv_imag;

        void *args[] = {(void *)(&sim->sim_gpu), (void *)&t};
        LAUNCH_KERNEL(RESETLC1_GATE);

        ValType sum = 0;
        ValType prob_of_one = 0;
        hipSafeCall(hipMemcpy(&sum, &(sim->m_real[0]), sizeof(ValType), hipMemcpyHostToDevice));
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&sum, &prob_of_one, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        if (prob_of_one < 1.0) // still possible to normalize
        {
            void *args[] = {(void *)(&sim->sim_gpu), (void *)&prob_of_one, (void *)&t};
            LAUNCH_KERNEL(RESETLC2_GATE);
        }
        else
        {
            if (qubit >= sim->lg2_m_gpu) // remote qubit, need switch
            {
                IdxType pair_gpu = (sim->i_proc) ^ ((IdxType)1 << (qubit - (sim->lg2_m_gpu)));
                assert(pair_gpu != sim->i_proc);

                ValType *sv_real_remote = sim->m_real;
                ValType *sv_imag_remote = sim->m_imag;

                if (sim->i_proc > pair_gpu)
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

                hipSafeCall(hipMemcpy(sv_real, sv_real_remote, per_pe_work * sizeof(ValType), hipMemcpyDeviceToDevice));
                hipSafeCall(hipMemcpy(sv_imag, sv_imag_remote, per_pe_work * sizeof(ValType), hipMemcpyDeviceToDevice));
            }
            else
            {
                void *args[] = {(void *)(&sim->sim_gpu), (void *)&t};
                LAUNCH_KERNEL(RESETLC3_GATE);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    __global__ void PURITYLC_GATE(SV_HIP_MPI *sim, const IdxType t)
    {
        grid_group grid = this_grid();
        const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
        const IdxType per_pe_work = ((sim->dim) >> (sim->gpu_scale));
        ValType *sv_real = sim->sv_real;
        ValType *sv_imag = sim->sv_imag;
        ValType *m_real = sim->m_real;
        for (IdxType i = tid; i < sim->m_gpu; i += blockDim.x * gridDim.x)
        {
            m_real[i] = sv_real[i] * sv_real[i] + sv_imag[i] * sv_imag[i];
        }
        grid.sync();
        // Parallel reduction
        for (IdxType k = (sim->m_gpu >> (IdxType)1); k > 0; k >>= 1)
        {
            for (IdxType i = tid; i < k; i += blockDim.x * gridDim.x)
            {
                m_real[i] += m_real[i + k];
            }
            grid.sync();
        }
    }

    void Purity_Check(const SV_HIP_MPI *sim, const IdxType t, const OP gate_op, const IdxType ctrl, const IdxType qubit)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        void *args[] = {(void *)(&sim->sim_gpu), (void *)&t};
        LAUNCH_KERNEL(PURITYLC_GATE);
        ValType sum = 0;
        ValType purity = 0;
        hipSafeCall(hipMemcpy(&sum, &(sim->m_real[0]), sizeof(ValType), hipMemcpyDeviceToHost));
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Reduce(&sum, &purity, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        if (sim->i_proc == 0 && abs(purity - 1.0) > ERROR_BAR)
        {
            printf("Purity Check fails after Gate-%lld=>%s(ctrl:%lld,qubit:%lld) with %lf\n", t, OP_NAMES[gate_op], ctrl, qubit, purity);
        }
    }

    //=====================================================================================
    void simulation_kernel(SV_HIP_MPI *sim, std::vector<SVGate> &gates)
    {
        IdxType cur_index = 0;

        for (IdxType t = 0; t < gates.size(); t++)
        {
            auto op_name = gates[t].op_name;
            auto qubit = (sv_gpu->gates_gpu)[t].qubit;
            auto ctrl = (sv_gpu->gates_gpu)[t].ctrl;

            if (op_name == OP::RESET)
            {
                RESET_GATE(sim, t);
            }
            else if (op_name == OP::M)
            {
                M_GATE(sim, t, cur_index);
                cur_index++;
            }
            else if (op_name == OP::MA)
            {
                MA_GATE(sim, t, cur_index);
                cur_index += qubit;
            }
            else if (op_name == OP::C1)
            {
                C1_GATE(sim, t, qubit);
            }
            else if (op_name == OP::C2)
            {
                C2_GATE(sim, t, ctrl, qubit);
            }

#ifdef PURITY_CHECK
            Purity_Check(sim, t, op_name, ctrl, qubit);
#endif
        }

        // Reduce all the measurement results
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(sim->results_local, sim->results, cur_index, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }

}; // namespace NWQSim

#endif
