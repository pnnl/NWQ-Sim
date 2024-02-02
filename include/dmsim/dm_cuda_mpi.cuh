#pragma once

#include "../state.hpp"
#include "../nwq_util.hpp"
#include "../gate.hpp"
#include "../circuit.hpp"

#include "../private/config.hpp"
#include "../private/cuda_util.cuh"
#include "../private/macros.hpp"
#include "../private/sim_gate.hpp"

#include "../circuit_pass/fusion.hpp"

#include <assert.h>
#include <random>
#include <complex.h>
#include <cooperative_groups.h>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <cuda.h>
#include <mpi.h>

#include <nvshmem.h>
#include <nvshmemx.h>
//#include <nvshmemx_error.h>

#ifdef FP64_TC_AVAILABLE
#include <mma.h>
#endif

namespace NWQSim
{
    using namespace cooperative_groups;
    using namespace std;
    using namespace nvcuda;

    // Simulation kernel runtime
    class DM_CUDA_MPI;
    __global__ void dm_simulation_kernel_cuda_mpi(DM_CUDA_MPI *dm_gpu, IdxType n_gates, IdxType n_qubits, bool enable_tc);

    class DM_CUDA_MPI : public QuantumState
    {
    public:
        DM_CUDA_MPI(IdxType _n_qubits) : QuantumState(_n_qubits)
        {
            // Initialize the GPU
            n_qubits = _n_qubits;
            dim = (IdxType)1 << (2 * n_qubits);
            half_dim = (IdxType)1 << (2 * n_qubits - 1);
            dm_size = dim * (IdxType)sizeof(ValType);

            // set GPUs and communication

            comm_global = MPI_COMM_WORLD;
            nvshmemx_init_attr_t attr;
            MPI_Comm comm = comm_global;
            attr.mpi_comm = &comm;
            nvshmemx_init_attr(NVSHMEMX_INIT_WITH_MPI_COMM, &attr);
            n_gpus = nvshmem_n_pes();
            i_proc = nvshmem_my_pe();
            // always be 0 since 1-MPI maps to 1-GPU
            cudaSafeCall(cudaSetDevice(0));

            gpu_scale = floor(log((double)n_gpus + 0.5) / log(2.0));
            lg2_m_gpu = 2 * n_qubits - gpu_scale;
            m_gpu = ((IdxType)1 << (lg2_m_gpu));

            dm_size_per_gpu = dm_size / n_gpus;
            // CPU side initialization
            if (!is_power_of_2(n_gpus))
                throw std::logic_error("Error: Number of GPUs should be power of 2.");
            if (dim % n_gpus != 0)
                throw std::logic_error("Error: Number of GPUs is too large or too small.");
            if (lg2_m_gpu < 5)
                throw std::logic_error("Error: Each GPU should have at least 5 qubits for multi-node version. Please increase qubits or reduce GPUs");
            // CPU side initialization
            SAFE_ALOC_HOST_CUDA(dm_real_cpu, dm_size_per_gpu);
            SAFE_ALOC_HOST_CUDA(dm_imag_cpu, dm_size_per_gpu);
            memset(dm_real_cpu, 0, dm_size_per_gpu);
            memset(dm_imag_cpu, 0, dm_size_per_gpu);
            // State-vector initial state [0..0] = 1
            if (i_proc == 0)
                dm_real_cpu[0] = 1.;
            // NVSHMEM GPU memory allocation
            dm_real = (ValType *)nvshmem_malloc(dm_size_per_gpu);
            dm_imag = (ValType *)nvshmem_malloc(dm_size_per_gpu);
            m_real = (ValType *)nvshmem_malloc(dm_size_per_gpu);
            m_imag = (ValType *)nvshmem_malloc(dm_size_per_gpu);
            cudaCheckError();
            gpu_mem += dm_size_per_gpu * 4;

            // GPU memory initilization
            cudaSafeCall(cudaMemcpy(dm_real, dm_real_cpu, dm_size_per_gpu,
                                    cudaMemcpyHostToDevice));
            cudaSafeCall(cudaMemcpy(dm_imag, dm_imag_cpu, dm_size_per_gpu,
                                    cudaMemcpyHostToDevice));
            cudaSafeCall(cudaMemset(m_real, 0, dm_size_per_gpu));
            cudaSafeCall(cudaMemset(m_imag, 0, dm_size_per_gpu));
            rng.seed(time(0));
        }

        ~DM_CUDA_MPI()
        {
            // Release for CPU side
            SAFE_FREE_HOST_CUDA(dm_real_cpu);
            SAFE_FREE_HOST_CUDA(dm_imag_cpu);
            SAFE_FREE_HOST_CUDA(randoms);
            SAFE_FREE_HOST_CUDA(results);

            // Release for GPU side
            nvshmem_free(dm_real);
            nvshmem_free(dm_imag);
            nvshmem_free(m_real);
            nvshmem_free(m_imag);

            SAFE_FREE_GPU(randoms_gpu);
            SAFE_FREE_GPU(gates_gpu);

            nvshmem_free(results_gpu);

            nvshmem_finalize();
        }

        void reset_state() override
        {
            // Reset CPU input & output
            memset(dm_real_cpu, 0, dm_size_per_gpu);
            memset(dm_imag_cpu, 0, dm_size_per_gpu);
            // State Vector initial state [0..0] = 1
            if (i_proc == 0)
                dm_real_cpu[0] = 1.;
            // GPU side initialization
            cudaSafeCall(cudaMemcpy(dm_real, dm_real_cpu,
                                    dm_size_per_gpu, cudaMemcpyHostToDevice));
            cudaSafeCall(cudaMemcpy(dm_imag, dm_imag_cpu,
                                    dm_size_per_gpu, cudaMemcpyHostToDevice));
            cudaSafeCall(cudaMemset(m_real, 0, dm_size_per_gpu));
            cudaSafeCall(cudaMemset(m_imag, 0, dm_size_per_gpu));
        }

        void set_seed(IdxType seed) override
        {
            rng.seed(seed);
        }

        void sim(std::shared_ptr<NWQSim::Circuit> circuit) override
        {

            IdxType origional_gates = circuit->num_gates();

            IdxType n_measures = prepare_measure(circuit->get_gates());

            std::vector<DMGate> cpu_vec = fuse_circuit_dm(circuit);

            assert(circuit->num_qubits() == n_qubits);

            // Copy the circuit to GPU
            copy_gates_to_gpu(cpu_vec);

            IdxType n_gates = cpu_vec.size();

            bool is_tc = false;

#ifdef FP64_TC_AVAILABLE
            is_tc = Config::ENABLE_TENSOR_CORE;
#endif
            DM_CUDA_MPI *dm_gpu;
            SAFE_ALOC_GPU(dm_gpu, sizeof(DM_CUDA_MPI));
            // Copy the simulator instance to GPU
            cudaSafeCall(cudaMemcpy(dm_gpu, this,
                                    sizeof(DM_CUDA_MPI), cudaMemcpyHostToDevice));

            double *sim_times;
            double sim_time;
            gpu_timer sim_timer;

            if (Config::PRINT_SIM_TRACE && i_proc == 0)
            {
                SAFE_ALOC_HOST_CUDA(sim_times, sizeof(double) * n_gpus);
                memset(sim_times, 0, sizeof(double) * n_gpus);
                printf("DMSim Multi-GPU is running! Requesting %lld qubits.\n", circuit->num_qubits());
            }

            dim3 gridDim(1, 1, 1);
            cudaDeviceProp deviceProp;
            cudaSafeCall(cudaGetDeviceProperties(&deviceProp, 0));
            // 8*16 is per warp shared-memory usage for C4 TC, with real and imag
            unsigned smem_size;
            if (is_tc)
                smem_size = THREADS_CTA_CUDA / 32 * 8 * 16 * 2 * sizeof(ValType);
            else
                smem_size = 0;

            int numBlocksPerSm;
            cudaSafeCall(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm,
                                                                       dm_simulation_kernel_cuda_mpi, THREADS_CTA_CUDA, smem_size));
            gridDim.x = numBlocksPerSm * deviceProp.multiProcessorCount;
            void *args[] = {&dm_gpu, &n_gates, &n_qubits, &is_tc};
            cudaSafeCall(cudaDeviceSynchronize());

            if (Config::PRINT_SIM_TRACE)
                MPI_Barrier(MPI_COMM_WORLD);
            sim_timer.start_timer();

            NVSHMEM_CHECK(nvshmemx_collective_launch((const void *)dm_simulation_kernel_cuda_mpi, gridDim,
                                                     THREADS_CTA_CUDA, args, smem_size, 0));

            cudaSafeCall(cudaDeviceSynchronize());

            sim_timer.stop_timer();
            sim_time = sim_timer.measure();

            if (i_proc == 0)
                printf("GPU kernel time: %.3lf ms\n", sim_time);

            // Copy the results back to CPU
            cudaSafeCall(cudaMemcpy(results, results_gpu, n_measures * sizeof(IdxType), cudaMemcpyDeviceToHost));

            if (i_proc == 0)
                printf("results copied to cpu\n");

            cudaCheckError();

            if (Config::PRINT_SIM_TRACE)
            {
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
                    printf("\n============== DM-Sim ===============\n");
                    printf("n_qubits:%lld, n_gates:%lld, sim_gates:%lld, ngpus:%lld, TensorCore: %s, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms, mem:%.3lf MB, mem_per_gpu:%.3lf MB\n",
                           n_qubits, origional_gates, n_gates, n_gpus, (Config::ENABLE_TENSOR_CORE ? "Enabled" : "Disabled"), avg_sim_time, 0.,
                           avg_sim_time, gpu_mem / 1024 / 1024 * n_gpus, gpu_mem / 1024 / 1024);
                    printf("=====================================\n");
                    SAFE_FREE_HOST_CUDA(sim_times);
                }
            }

            SAFE_FREE_GPU(dm_gpu);
        }

        IdxType *get_results() override
        {
            return results;
        }

        ValType get_exp_z(const std::vector<size_t> &in_bits) override
        {
            throw std::logic_error("get_exp_z Not implemented (DM_CUDA_MPI)");
        }

        ValType get_exp_z() override
        {
            throw std::logic_error("get_exp_z Not implemented (DM_CUDA_MPI)");
        }

        IdxType measure(IdxType qubit) override
        {
            std::shared_ptr<NWQSim::Circuit> circuit = std::make_shared<Circuit>(n_qubits);
            circuit->M(qubit);
            sim(circuit);
            return results[0];
        }

        IdxType *measure_all(IdxType repetition) override
        {
            std::shared_ptr<NWQSim::Circuit> circuit = std::make_shared<Circuit>(n_qubits);
            circuit->MA(repetition);
            sim(circuit);
            return results;
        }

        void print_res_state() override
        {

            cudaCheckError();
            cudaSafeCall(cudaMemcpy(dm_real_cpu, dm_real, dm_size_per_gpu, cudaMemcpyDeviceToHost));
            cudaSafeCall(cudaMemcpy(dm_imag_cpu, dm_imag, dm_size_per_gpu, cudaMemcpyDeviceToHost));

            ValType *dm_diag_real = NULL;
            ValType *dm_diag_imag = NULL;
            if (i_proc == 0)
                SAFE_ALOC_HOST_CUDA(dm_diag_real, dim * sizeof(ValType));
            if (i_proc == 0)
                SAFE_ALOC_HOST_CUDA(dm_diag_imag, dim * sizeof(ValType));

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Gather(dm_real_cpu, m_gpu, MPI_DOUBLE,
                       &dm_diag_real[i_proc * m_gpu], m_gpu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gather(dm_imag_cpu, m_gpu, MPI_DOUBLE,
                       &dm_diag_imag[i_proc * m_gpu], m_gpu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            if (i_proc == 0)
            {
                IdxType num = ((IdxType)1 << n_qubits);
                printf("----- DM Diag ------\n");
                for (IdxType i = 0; i < num; i++)
                {
                    printf("(%.3lf,%.3lfj) ", dm_diag_real[i], dm_diag_imag[i]);
                    if ((i + 1) % 8 == 0)
                        printf("\n");
                }
                printf("\n");
                SAFE_FREE_HOST_CUDA(dm_diag_real);
                SAFE_FREE_HOST_CUDA(dm_diag_imag);
            }
        }

    public:
        // n_qubits is the number of qubits
        IdxType n_qubits;

        IdxType gpu_scale;
        IdxType n_gpus;
        IdxType lg2_m_gpu;
        IdxType m_gpu;
        IdxType dm_size_per_gpu;

        // gpu_scale is 2^x of the number of GPUs, e.g., with 8 GPUs the gpu_scale is 3 (2^3=8)
        IdxType dim;
        IdxType half_dim;
        IdxType dm_size;
        // CPU arrays
        ValType *dm_real_cpu;
        ValType *dm_imag_cpu;
        // GPU arrays
        ValType *dm_real;
        ValType *dm_imag;
        // For joint measurement
        ValType *m_real;
        ValType *m_imag;
        // For measurement randoms
        ValType *randoms = NULL;
        ValType *randoms_gpu = NULL;
        // For measurement result
        IdxType *results = NULL;
        IdxType *results_gpu = NULL;
        // Random
        std::mt19937 rng;
        std::uniform_real_distribution<ValType> uni_dist;
        // GPU memory usage
        ValType gpu_mem;

        // hold the GPU-side simulator instances
        MPI_Comm comm_global;

        // GPU-side simulator instance
        DMGate *gates_gpu = NULL;

        void copy_gates_to_gpu(std::vector<DMGate> &cpu_vec)
        {
            // Allocate memory on CPU
            size_t vec_size = cpu_vec.size() * sizeof(DMGate);

            // Allocate memory on GPU
            SAFE_FREE_GPU(gates_gpu);
            SAFE_ALOC_GPU(gates_gpu, vec_size);
            cudaSafeCall(cudaMemcpy(gates_gpu, cpu_vec.data(), vec_size, cudaMemcpyHostToDevice));
        }

        IdxType prepare_measure(std::vector<Gate> gates)
        {
            // Determine the total number of measurements
            IdxType n_measures = 0;
            for (auto g : gates)
            {
                if (g.op_name == OP::M)
                    n_measures++;
                else if (g.op_name == OP::MA)
                    n_measures += g.repetition;
            }

            // Prepare randoms and results memory
            SAFE_FREE_HOST_CUDA(results);
            SAFE_ALOC_HOST_CUDA(results, sizeof(IdxType) * n_measures);
            memset(results, 0, sizeof(IdxType) * n_measures);

            nvshmem_free(results_gpu);
            results_gpu = (IdxType *)nvshmem_malloc(sizeof(IdxType) * n_measures);

            cudaSafeCall(cudaMemset(results_gpu, 0, sizeof(IdxType) * n_measures));

            SAFE_FREE_HOST_CUDA(randoms);
            SAFE_ALOC_HOST_CUDA(randoms, sizeof(ValType) * n_measures);

            for (IdxType i = 0; i < n_measures; i++)
                randoms[i] = uni_dist(rng);
            SAFE_FREE_GPU(randoms_gpu);
            SAFE_ALOC_GPU(randoms_gpu, sizeof(ValType) * n_measures);
            cudaSafeCall(cudaMemcpy(randoms_gpu, randoms,
                                    sizeof(ValType) * n_measures, cudaMemcpyHostToDevice));

            return n_measures;
        }

        //================================= Gate Definition ========================================
        //============== Unified 2-qubit Gate without comm optimization ================
        __device__ __inline__ void C2_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit0, const IdxType qubit1)
        {
            grid_group grid = this_grid();
            const int tid = blockDim.x * blockIdx.x + threadIdx.x;
            const IdxType per_pe_work = ((dim) >> (gpu_scale + 2));
            assert(qubit0 != qubit1); // Non-cloning

            const IdxType q0dim = ((IdxType)1 << max(qubit0, qubit1));
            const IdxType q1dim = ((IdxType)1 << min(qubit0, qubit1));
            const IdxType outer_factor = ((dim) + q0dim + q0dim - 1) >> (max(qubit0, qubit1) + 1);
            const IdxType mider_factor = (q0dim + q1dim + q1dim - 1) >> (min(qubit0, qubit1) + 1);
            const IdxType inner_factor = q1dim;
            const IdxType qubit0_dim = ((IdxType)1 << qubit0);
            const IdxType qubit1_dim = ((IdxType)1 << qubit1);

            for (IdxType i = (i_proc)*per_pe_work + tid; i < (i_proc + 1) * per_pe_work;
                 i += blockDim.x * gridDim.x)
            {
                IdxType outer = ((i / inner_factor) / (mider_factor)) * (q0dim + q0dim);
                IdxType mider = ((i / inner_factor) % (mider_factor)) * (q1dim + q1dim);
                IdxType inner = i % inner_factor;
                IdxType pos0 = outer + mider + inner;
                IdxType pos1 = outer + mider + inner + qubit1_dim;
                IdxType pos2 = outer + mider + inner + qubit0_dim;
                IdxType pos3 = outer + mider + inner + q0dim + q1dim;

                const ValType el0_real = PGAS_G(dm_real, pos0);
                const ValType el0_imag = PGAS_G(dm_imag, pos0);
                const ValType el1_real = PGAS_G(dm_real, pos1);
                const ValType el1_imag = PGAS_G(dm_imag, pos1);
                const ValType el2_real = PGAS_G(dm_real, pos2);
                const ValType el2_imag = PGAS_G(dm_imag, pos2);
                const ValType el3_real = PGAS_G(dm_real, pos3);
                const ValType el3_imag = PGAS_G(dm_imag, pos3);

                // Real part
                ValType dm_real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag) + (gm_real[2] * el2_real) - (gm_imag[2] * el2_imag) + (gm_real[3] * el3_real) - (gm_imag[3] * el3_imag);
                ValType dm_real_pos1 = (gm_real[4] * el0_real) - (gm_imag[4] * el0_imag) + (gm_real[5] * el1_real) - (gm_imag[5] * el1_imag) + (gm_real[6] * el2_real) - (gm_imag[6] * el2_imag) + (gm_real[7] * el3_real) - (gm_imag[7] * el3_imag);
                ValType dm_real_pos2 = (gm_real[8] * el0_real) - (gm_imag[8] * el0_imag) + (gm_real[9] * el1_real) - (gm_imag[9] * el1_imag) + (gm_real[10] * el2_real) - (gm_imag[10] * el2_imag) + (gm_real[11] * el3_real) - (gm_imag[11] * el3_imag);
                ValType dm_real_pos3 = (gm_real[12] * el0_real) - (gm_imag[12] * el0_imag) + (gm_real[13] * el1_real) - (gm_imag[13] * el1_imag) + (gm_real[14] * el2_real) - (gm_imag[14] * el2_imag) + (gm_real[15] * el3_real) - (gm_imag[15] * el3_imag);

                // Imag part
                ValType dm_imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real) + (gm_real[2] * el2_imag) + (gm_imag[2] * el2_real) + (gm_real[3] * el3_imag) + (gm_imag[3] * el3_real);
                ValType dm_imag_pos1 = (gm_real[4] * el0_imag) + (gm_imag[4] * el0_real) + (gm_real[5] * el1_imag) + (gm_imag[5] * el1_real) + (gm_real[6] * el2_imag) + (gm_imag[6] * el2_real) + (gm_real[7] * el3_imag) + (gm_imag[7] * el3_real);
                ValType dm_imag_pos2 = (gm_real[8] * el0_imag) + (gm_imag[8] * el0_real) + (gm_real[9] * el1_imag) + (gm_imag[9] * el1_real) + (gm_real[10] * el2_imag) + (gm_imag[10] * el2_real) + (gm_real[11] * el3_imag) + (gm_imag[11] * el3_real);
                ValType dm_imag_pos3 = (gm_real[12] * el0_imag) + (gm_imag[12] * el0_real) + (gm_real[13] * el1_imag) + (gm_imag[13] * el1_real) + (gm_real[14] * el2_imag) + (gm_imag[14] * el2_real) + (gm_real[15] * el3_imag) + (gm_imag[15] * el3_real);

                PGAS_P(dm_real, pos0, dm_real_pos0);
                PGAS_P(dm_real, pos1, dm_real_pos1);
                PGAS_P(dm_real, pos2, dm_real_pos2);
                PGAS_P(dm_real, pos3, dm_real_pos3);

                PGAS_P(dm_imag, pos0, dm_imag_pos0);
                PGAS_P(dm_imag, pos1, dm_imag_pos1);
                PGAS_P(dm_imag, pos2, dm_imag_pos2);
                PGAS_P(dm_imag, pos3, dm_imag_pos3);
            }
            // BARR;
        }
        //============== Unified 2-qubit Gate ================
        // Perform communication optimization here
        // Since this is for single-qubit density matrix gate,
        // It is only possible qubit0 < qubit1 and qubit1 can be remote
        __device__ __inline__ void C2V1_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit0, const IdxType qubit1)
        {
            assert(qubit0 != qubit1); // Non-cloning

            if (qubit1 < lg2_m_gpu)
            {
                C2_GATE(gm_real, gm_imag, qubit0, qubit1);
            }
            else
            {
                grid_group grid = this_grid();
                const int tid = blockDim.x * blockIdx.x + threadIdx.x;
                // For processing a remote qubit, half of the GPUs will be idle for
                // better communication efficiency. Depending on qubit position,
                // for a GPU pair, we only use GPUs with smaller ids.
                // Consequently, each GPU should take double of the workload than before
                // Therefore, here it is gpu_scale+1 not gpu_scale+2
                const IdxType per_pe_work = ((dim) >> (gpu_scale + 1));
                const IdxType per_pe_num = ((dim) >> (gpu_scale));
                const IdxType p = min(qubit0, qubit1);
                const IdxType q = max(qubit0, qubit1);

                // load data from pair GPU
                IdxType pair_gpu = (i_proc) ^ ((IdxType)1 << (q - (lg2_m_gpu)));
                if (i_proc > pair_gpu)
                    return;

                ValType *dm_real_remote = m_real;
                ValType *dm_imag_remote = m_imag;

                if (tid == 0)
                    nvshmem_double_get(dm_real_remote, dm_real, per_pe_num, pair_gpu);
                if (tid == 0)
                    nvshmem_double_get(dm_imag_remote, dm_imag, per_pe_num, pair_gpu);
                grid.sync();

                for (IdxType i = (i_proc)*per_pe_work + tid; i < (i_proc + 1) * per_pe_work;
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
                    el_real[0] = LOCAL_G_CUDA_MPI(dm_real, term + SV4IDX(0));
                    el_imag[0] = LOCAL_G_CUDA_MPI(dm_imag, term + SV4IDX(0));
                    el_real[3] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV4IDX(3));
                    el_imag[3] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV4IDX(3));

                    if (qubit0 == q) // qubit0 is the remote qubit
                    {
                        el_real[1] = LOCAL_G_CUDA_MPI(dm_real, term + SV4IDX(1));
                        el_imag[1] = LOCAL_G_CUDA_MPI(dm_imag, term + SV4IDX(1));
                        el_real[2] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV4IDX(2));
                        el_imag[2] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV4IDX(2));
                    }
                    else // qubit1 is the remote qubit
                    {
                        el_real[1] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV4IDX(1));
                        el_imag[1] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV4IDX(1));
                        el_real[2] = LOCAL_G_CUDA_MPI(dm_real, term + SV4IDX(2));
                        el_imag[2] = LOCAL_G_CUDA_MPI(dm_imag, term + SV4IDX(2));
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
                    LOCAL_P_CUDA_MPI(dm_real, term + SV4IDX(0), res_real[0]);
                    LOCAL_P_CUDA_MPI(dm_imag, term + SV4IDX(0), res_imag[0]);
                    LOCAL_P_CUDA_MPI(dm_real_remote, term + SV4IDX(3), res_real[3]);
                    LOCAL_P_CUDA_MPI(dm_imag_remote, term + SV4IDX(3), res_imag[3]);

                    if (qubit0 == q) // qubit0 is the remote qubit
                    {
                        LOCAL_P_CUDA_MPI(dm_real, term + SV4IDX(1), res_real[1]);
                        LOCAL_P_CUDA_MPI(dm_imag, term + SV4IDX(1), res_imag[1]);
                        LOCAL_P_CUDA_MPI(dm_real_remote, term + SV4IDX(2), res_real[2]);
                        LOCAL_P_CUDA_MPI(dm_imag_remote, term + SV4IDX(2), res_imag[2]);
                    }
                    else // qubit1 is the remote qubit
                    {
                        LOCAL_P_CUDA_MPI(dm_real_remote, term + SV4IDX(1), res_real[1]);
                        LOCAL_P_CUDA_MPI(dm_imag_remote, term + SV4IDX(1), res_imag[1]);
                        LOCAL_P_CUDA_MPI(dm_real, term + SV4IDX(2), res_real[2]);
                        LOCAL_P_CUDA_MPI(dm_imag, term + SV4IDX(2), res_imag[2]);
                    }
                }
                grid.sync();
                if (tid == 0)
                    nvshmem_double_put(dm_real, dm_real_remote, per_pe_num, pair_gpu);
                if (tid == 0)
                    nvshmem_double_put(dm_imag, dm_imag_remote, per_pe_num, pair_gpu);
                // BARR_NVSHMEM;
            }
        }

        //============== SWAP Gate ================
        // This gate is for internal usage. It is used
        // when two qubits of a C4 gate are remote qubis, we then
        // swap one of them to a local qubit without noise,
        // perform the C4 gate, and then swap back
        // It is assumed qubit0 is local, qubit1 is remote
        __device__ __inline__ void SWAP_GATE(
            const IdxType qubit0, const IdxType qubit1)
        {
            grid_group grid = this_grid();
            const int tid = blockDim.x * blockIdx.x + threadIdx.x;
            assert(qubit0 < qubit1);     // Non-cloning and qubit0<qubit1 assumption
            assert(qubit1 >= lg2_m_gpu); // Otherwise don't need swap

            const IdxType per_pe_work = ((dim) >> (gpu_scale + 1));
            const IdxType per_pe_num = ((dim) >> (gpu_scale));

            const IdxType p = qubit0;
            const IdxType q = qubit1;

            // load data from pair GPU
            IdxType pair_gpu = (i_proc) ^ ((IdxType)1 << (q - (lg2_m_gpu)));
            if (i_proc > pair_gpu)
                return;

            ValType *dm_real_remote = m_real;
            ValType *dm_imag_remote = m_imag;

            if (tid == 0)
                nvshmem_double_get(dm_real_remote, dm_real, per_pe_num, pair_gpu);
            if (tid == 0)
                nvshmem_double_get(dm_imag_remote, dm_imag, per_pe_num, pair_gpu);
            grid.sync();

            for (IdxType i = (i_proc)*per_pe_work + tid; i < (i_proc + 1) * per_pe_work;
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

                el_real[1] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV4IDX(1));
                el_imag[1] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV4IDX(1));
                el_real[2] = LOCAL_G_CUDA_MPI(dm_real, term + SV4IDX(2));
                el_imag[2] = LOCAL_G_CUDA_MPI(dm_imag, term + SV4IDX(2));

                res_real[1] = el_real[2];
                res_imag[1] = el_imag[2];
                res_real[2] = el_real[1];
                res_imag[2] = el_imag[1];

                LOCAL_P_CUDA_MPI(dm_real_remote, term + SV4IDX(1), res_real[1]);
                LOCAL_P_CUDA_MPI(dm_imag_remote, term + SV4IDX(1), res_imag[1]);
                LOCAL_P_CUDA_MPI(dm_real, term + SV4IDX(2), res_real[2]);
                LOCAL_P_CUDA_MPI(dm_imag, term + SV4IDX(2), res_imag[2]);
            }
            grid.sync();
            if (tid == 0)
                nvshmem_double_put(dm_real, dm_real_remote, per_pe_num, pair_gpu);
            if (tid == 0)
                nvshmem_double_put(dm_imag, dm_imag_remote, per_pe_num, pair_gpu);
        }

        //============== Unified 4-qubit Gate ================
        __device__ __inline__ void C4_GATE(const ValType *gm_real, const ValType *gm_imag,
                                           const IdxType qubit0, const IdxType qubit1,
                                           const IdxType qubit2, const IdxType qubit3)
        {
            grid_group grid = this_grid();
            const int tid = blockDim.x * blockIdx.x + threadIdx.x;
            const IdxType per_pe_work = ((dim) >> (gpu_scale + 4));
            assert(qubit0 != qubit1); // Non-cloning
            assert(qubit0 != qubit2); // Non-cloning
            assert(qubit0 != qubit3); // Non-cloning
            assert(qubit1 != qubit2); // Non-cloning
            assert(qubit1 != qubit3); // Non-cloning
            assert(qubit2 != qubit3); // Non-cloning

            // need to sort qubits: min->max: p, q, r, s
            const IdxType v0 = min(qubit0, qubit1);
            const IdxType v1 = min(qubit2, qubit3);
            const IdxType v2 = max(qubit0, qubit1);
            const IdxType v3 = max(qubit2, qubit3);
            const IdxType p = min(v0, v1);
            const IdxType q = min(min(v2, v3), max(v0, v1));
            const IdxType r = max(min(v2, v3), max(v0, v1));
            const IdxType s = max(v2, v3);

            for (IdxType i = (i_proc)*per_pe_work + tid; i < (i_proc + 1) * per_pe_work;
                 i += blockDim.x * gridDim.x)
            {
                const IdxType term0 = MOD2E(i, p);
                const IdxType term1 = MOD2E(DIV2E(i, p), q - p - 1) * EXP2E(p + 1);
                const IdxType term2 = MOD2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1) * EXP2E(q + 1);
                const IdxType term3 = MOD2E(DIV2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(r + 1);
                const IdxType term4 = DIV2E(DIV2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(s + 1);
                const IdxType term = term4 + term3 + term2 + term1 + term0;

                const ValType el_real[16] = {
                    PGAS_G(dm_real, term + SV16IDX(0)), PGAS_G(dm_real, term + SV16IDX(1)),
                    PGAS_G(dm_real, term + SV16IDX(2)), PGAS_G(dm_real, term + SV16IDX(3)),
                    PGAS_G(dm_real, term + SV16IDX(4)), PGAS_G(dm_real, term + SV16IDX(5)),
                    PGAS_G(dm_real, term + SV16IDX(6)), PGAS_G(dm_real, term + SV16IDX(7)),
                    PGAS_G(dm_real, term + SV16IDX(8)), PGAS_G(dm_real, term + SV16IDX(9)),
                    PGAS_G(dm_real, term + SV16IDX(10)), PGAS_G(dm_real, term + SV16IDX(11)),
                    PGAS_G(dm_real, term + SV16IDX(12)), PGAS_G(dm_real, term + SV16IDX(13)),
                    PGAS_G(dm_real, term + SV16IDX(14)), PGAS_G(dm_real, term + SV16IDX(15))};
                const ValType el_imag[16] = {
                    PGAS_G(dm_imag, term + SV16IDX(0)), PGAS_G(dm_imag, term + SV16IDX(1)),
                    PGAS_G(dm_imag, term + SV16IDX(2)), PGAS_G(dm_imag, term + SV16IDX(3)),
                    PGAS_G(dm_imag, term + SV16IDX(4)), PGAS_G(dm_imag, term + SV16IDX(5)),
                    PGAS_G(dm_imag, term + SV16IDX(6)), PGAS_G(dm_imag, term + SV16IDX(7)),
                    PGAS_G(dm_imag, term + SV16IDX(8)), PGAS_G(dm_imag, term + SV16IDX(9)),
                    PGAS_G(dm_imag, term + SV16IDX(10)), PGAS_G(dm_imag, term + SV16IDX(11)),
                    PGAS_G(dm_imag, term + SV16IDX(12)), PGAS_G(dm_imag, term + SV16IDX(13)),
                    PGAS_G(dm_imag, term + SV16IDX(14)), PGAS_G(dm_imag, term + SV16IDX(15))};
#pragma unroll
                for (unsigned j = 0; j < 16; j++)
                {
                    ValType res_real = 0;
                    ValType res_imag = 0;
#pragma unroll
                    for (unsigned k = 0; k < 16; k++)
                    {
                        res_real += (el_real[k] * gm_real[j * 16 + k]) - (el_imag[k] * gm_imag[j * 16 + k]);
                        res_imag += (el_real[k] * gm_imag[j * 16 + k]) + (el_imag[k] * gm_real[j * 16 + k]);
                    }
                    PGAS_P(dm_real, term + SV16IDX(j), res_real);
                    PGAS_P(dm_imag, term + SV16IDX(j), res_imag);
                }
            }
            // BARR_NVSHMEM;
        }

        //============== Unified 4-qubit Gate with comm optimization ================
        // This is with nvshmem merged communication optimization but without tensorcore
        __device__ __inline__ void C4V1_GATE(const ValType *gm_real, const ValType *gm_imag,
                                             const IdxType qubit0, const IdxType qubit1,
                                             const IdxType qubit2, const IdxType qubit3)
        {
            grid_group grid = this_grid();

            assert(qubit0 != qubit1); // Non-cloning
            assert(qubit0 != qubit2); // Non-cloning
            assert(qubit0 != qubit3); // Non-cloning
            assert(qubit1 != qubit2); // Non-cloning
            assert(qubit1 != qubit3); // Non-cloning
            assert(qubit2 != qubit3); // Non-cloning

            // need to sort qubits: min->max: p, q, r, s
            const IdxType v0 = min(qubit0, qubit1);
            const IdxType v1 = min(qubit2, qubit3);
            const IdxType v2 = max(qubit0, qubit1);
            const IdxType v3 = max(qubit2, qubit3);
            const IdxType p = min(v0, v1);
            const IdxType q = min(min(v2, v3), max(v0, v1));
            const IdxType r = max(min(v2, v3), max(v0, v1));
            const IdxType s = max(v2, v3);

            if (s < lg2_m_gpu) // all 4 qubits are local
            {
                C4_GATE(gm_real, gm_imag, qubit0, qubit1, qubit2, qubit3);
            }
            else // s qubit is non-local
            {
                assert(r < lg2_m_gpu);
                const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
                const IdxType per_pe_num = ((dim) >> (gpu_scale));
                const IdxType per_pe_work = ((dim) >> (gpu_scale + 3));

                // load data from pair GPU
                IdxType pair_gpu = (i_proc) ^ ((IdxType)1 << (s - (lg2_m_gpu)));
                if (i_proc > pair_gpu)
                    return;

                ValType *dm_real_remote = m_real;
                ValType *dm_imag_remote = m_imag;

                if (tid == 0)
                    nvshmem_double_get(dm_real_remote, dm_real, per_pe_num, pair_gpu);
                if (tid == 0)
                    nvshmem_double_get(dm_imag_remote, dm_imag, per_pe_num, pair_gpu);
                grid.sync();

                for (IdxType i = (i_proc)*per_pe_work + tid; i < (i_proc + 1) * per_pe_work;
                     i += blockDim.x * gridDim.x)
                {
                    const IdxType term0 = MOD2E(i, p);
                    const IdxType term1 = MOD2E(DIV2E(i, p), q - p - 1) * EXP2E(p + 1);
                    const IdxType term2 = MOD2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1) * EXP2E(q + 1);
                    const IdxType term3 = MOD2E(DIV2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(r + 1);
                    const IdxType term4 = DIV2E(DIV2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(s + 1);
                    const IdxType term = term4 + term3 + term2 + term1 + term0;
                    ValType el_real[16];
                    ValType el_imag[16];

                    if (qubit3 == s) // qubit3 is remote qubit
                    {
                        el_real[0] = LOCAL_G_CUDA_MPI(dm_real, term + SV16IDX(0));
                        el_real[1] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV16IDX(1));
                        el_real[2] = LOCAL_G_CUDA_MPI(dm_real, term + SV16IDX(2));
                        el_real[3] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV16IDX(3));
                        el_real[4] = LOCAL_G_CUDA_MPI(dm_real, term + SV16IDX(4));
                        el_real[5] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV16IDX(5));
                        el_real[6] = LOCAL_G_CUDA_MPI(dm_real, term + SV16IDX(6));
                        el_real[7] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV16IDX(7));
                        el_real[8] = LOCAL_G_CUDA_MPI(dm_real, term + SV16IDX(8));
                        el_real[9] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV16IDX(9));
                        el_real[10] = LOCAL_G_CUDA_MPI(dm_real, term + SV16IDX(10));
                        el_real[11] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV16IDX(11));
                        el_real[12] = LOCAL_G_CUDA_MPI(dm_real, term + SV16IDX(12));
                        el_real[13] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV16IDX(13));
                        el_real[14] = LOCAL_G_CUDA_MPI(dm_real, term + SV16IDX(14));
                        el_real[15] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV16IDX(15));

                        el_imag[0] = LOCAL_G_CUDA_MPI(dm_imag, term + SV16IDX(0));
                        el_imag[1] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV16IDX(1));
                        el_imag[2] = LOCAL_G_CUDA_MPI(dm_imag, term + SV16IDX(2));
                        el_imag[3] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV16IDX(3));
                        el_imag[4] = LOCAL_G_CUDA_MPI(dm_imag, term + SV16IDX(4));
                        el_imag[5] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV16IDX(5));
                        el_imag[6] = LOCAL_G_CUDA_MPI(dm_imag, term + SV16IDX(6));
                        el_imag[7] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV16IDX(7));
                        el_imag[8] = LOCAL_G_CUDA_MPI(dm_imag, term + SV16IDX(8));
                        el_imag[9] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV16IDX(9));
                        el_imag[10] = LOCAL_G_CUDA_MPI(dm_imag, term + SV16IDX(10));
                        el_imag[11] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV16IDX(11));
                        el_imag[12] = LOCAL_G_CUDA_MPI(dm_imag, term + SV16IDX(12));
                        el_imag[13] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV16IDX(13));
                        el_imag[14] = LOCAL_G_CUDA_MPI(dm_imag, term + SV16IDX(14));
                        el_imag[15] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV16IDX(15));
                    }
                    else // qubit2 is remote (not possible qubit0 or 1 is remote)
                    {
                        // if (( (laneid>>1)&1)!=0) el_real_s[j*16+laneid] = LOCAL_G_CUDA_MPI(dm_real_remote,addr);
                        // else el_real_s[j*16+laneid] = LOCAL_G_CUDA_MPI(dm_real, addr);
                        // laneid = 2,3,6,7,10,11,14,15;

                        el_real[0] = LOCAL_G_CUDA_MPI(dm_real, term + SV16IDX(0));
                        el_real[1] = LOCAL_G_CUDA_MPI(dm_real, term + SV16IDX(1));
                        el_real[2] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV16IDX(2));
                        el_real[3] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV16IDX(3));
                        el_real[4] = LOCAL_G_CUDA_MPI(dm_real, term + SV16IDX(4));
                        el_real[5] = LOCAL_G_CUDA_MPI(dm_real, term + SV16IDX(5));
                        el_real[6] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV16IDX(6));
                        el_real[7] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV16IDX(7));
                        el_real[8] = LOCAL_G_CUDA_MPI(dm_real, term + SV16IDX(8));
                        el_real[9] = LOCAL_G_CUDA_MPI(dm_real, term + SV16IDX(9));
                        el_real[10] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV16IDX(10));
                        el_real[11] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV16IDX(11));
                        el_real[12] = LOCAL_G_CUDA_MPI(dm_real, term + SV16IDX(12));
                        el_real[13] = LOCAL_G_CUDA_MPI(dm_real, term + SV16IDX(13));
                        el_real[14] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV16IDX(14));
                        el_real[15] = LOCAL_G_CUDA_MPI(dm_real_remote, term + SV16IDX(15));

                        el_imag[0] = LOCAL_G_CUDA_MPI(dm_imag, term + SV16IDX(0));
                        el_imag[1] = LOCAL_G_CUDA_MPI(dm_imag, term + SV16IDX(1));
                        el_imag[2] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV16IDX(2));
                        el_imag[3] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV16IDX(3));
                        el_imag[4] = LOCAL_G_CUDA_MPI(dm_imag, term + SV16IDX(4));
                        el_imag[5] = LOCAL_G_CUDA_MPI(dm_imag, term + SV16IDX(5));
                        el_imag[6] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV16IDX(6));
                        el_imag[7] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV16IDX(7));
                        el_imag[8] = LOCAL_G_CUDA_MPI(dm_imag, term + SV16IDX(8));
                        el_imag[9] = LOCAL_G_CUDA_MPI(dm_imag, term + SV16IDX(9));
                        el_imag[10] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV16IDX(10));
                        el_imag[11] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV16IDX(11));
                        el_imag[12] = LOCAL_G_CUDA_MPI(dm_imag, term + SV16IDX(12));
                        el_imag[13] = LOCAL_G_CUDA_MPI(dm_imag, term + SV16IDX(13));
                        el_imag[14] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV16IDX(14));
                        el_imag[15] = LOCAL_G_CUDA_MPI(dm_imag_remote, term + SV16IDX(15));
                    }
#pragma unroll
                    for (unsigned j = 0; j < 16; j++)
                    {
                        ValType res_real = 0;
                        ValType res_imag = 0;
#pragma unroll
                        for (unsigned k = 0; k < 16; k++)
                        {
                            res_real += (el_real[k] * gm_real[j * 16 + k]) - (el_imag[k] * gm_imag[j * 16 + k]);
                            res_imag += (el_real[k] * gm_imag[j * 16 + k]) + (el_imag[k] * gm_real[j * 16 + k]);
                        }
                        // PGAS_P(dm_real, term+SV16IDX(j), res_real);
                        // PGAS_P(dm_imag, term+SV16IDX(j), res_imag);
                        if (qubit3 == s) // qubit3 is remote qubit
                        {
                            if ((j & 1) != 0)
                            {
                                LOCAL_P_CUDA_MPI(dm_real_remote, term + SV16IDX(j), res_real);
                                LOCAL_P_CUDA_MPI(dm_imag_remote, term + SV16IDX(j), res_imag);
                            }
                            else
                            {
                                LOCAL_P_CUDA_MPI(dm_real, term + SV16IDX(j), res_real);
                                LOCAL_P_CUDA_MPI(dm_imag, term + SV16IDX(j), res_imag);
                            }
                        }
                        else // qubit2 is remote (not possible qubit0 or 1 is remote)
                        {
                            if (((j >> 1) & 1) != 0)
                            {
                                LOCAL_P_CUDA_MPI(dm_real_remote, term + SV16IDX(j), res_real);
                                LOCAL_P_CUDA_MPI(dm_imag_remote, term + SV16IDX(j), res_imag);
                            }
                            else
                            {
                                LOCAL_P_CUDA_MPI(dm_real, term + SV16IDX(j), res_real);
                                LOCAL_P_CUDA_MPI(dm_imag, term + SV16IDX(j), res_imag);
                            }
                        }
                    }
                }
                grid.sync();
                if (tid == 0)
                    nvshmem_double_put(dm_real, dm_real_remote, per_pe_num, pair_gpu);
                if (tid == 0)
                    nvshmem_double_put(dm_imag, dm_imag_remote, per_pe_num, pair_gpu);
                grid.sync();
            }
        }

        __device__ __inline__ IdxType get_term(IdxType idx, IdxType p, IdxType q, IdxType r, IdxType s)
        {
            const IdxType term0 = MOD2E(idx, p);
            const IdxType term1 = MOD2E(DIV2E(idx, p), q - p - 1) * EXP2E(p + 1);
            const IdxType term2 = MOD2E(DIV2E(DIV2E(idx, p), q - p - 1), r - q - 1) * EXP2E(q + 1);
            const IdxType term3 = MOD2E(DIV2E(DIV2E(DIV2E(idx, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(r + 1);
            const IdxType term4 = DIV2E(DIV2E(DIV2E(DIV2E(idx, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(s + 1);
            const IdxType term = term4 + term3 + term2 + term1 + term0;
            return term;
        }

#ifdef FP64_TC_AVAILABLE
        //============== Unified 4-qubit Gate with TC V1 ================
        // This is with tensorcore optimization
        __device__ __inline__ void C4TCV1_GATE(const ValType *gm_real, const ValType *gm_imag,
                                               const IdxType qubit0, const IdxType qubit1,
                                               const IdxType qubit2, const IdxType qubit3)
        {
            grid_group grid = this_grid();
            const int tid = blockDim.x * blockIdx.x + threadIdx.x;
            const int laneid = (threadIdx.x & 31);
            const int warpid = (tid >> 5);
            const int wid = (threadIdx.x >> 5);
            const int nwarps = ((blockDim.x * gridDim.x) >> 5);
            const IdxType per_pe_work = ((dim) >> (gpu_scale + 4));
            // Our design require at least (gpu_scale + 4 + log(8)) qubits, because
            // tensorcore demands [8,4]*[4,8], we need at least 8 on B's row dim.
            // This implies for 2 GPUs, we need at least 8 qubit for SV, thus 4 qubit for DM.
            // The smaller condition is difficult to handle due to the complex get_term function,
            // We need to figure out, after the get_term, which row/col it maps to and filter out.
            // But the row and col will be quite difficult to obtain.

            extern __shared__ ValType els[];
            ValType *el_real_s = &els[wid * 8 * 16 * 2];          // per warp 8*16 for real
            ValType *el_imag_s = &els[wid * 8 * 16 * 2 + 8 * 16]; // per warp 8*16 for imag

            assert(qubit0 != qubit1); // Non-cloning
            assert(qubit0 != qubit2); // Non-cloning
            assert(qubit0 != qubit3); // Non-cloning
            assert(qubit1 != qubit2); // Non-cloning
            assert(qubit1 != qubit3); // Non-cloning
            assert(qubit2 != qubit3); // Non-cloning

            // need to sort qubits: min->max: p, q, r, s
            const IdxType v0 = min(qubit0, qubit1);
            const IdxType v1 = min(qubit2, qubit3);
            const IdxType v2 = max(qubit0, qubit1);
            const IdxType v3 = max(qubit2, qubit3);
            const IdxType p = min(v0, v1);
            const IdxType q = min(min(v2, v3), max(v0, v1));
            const IdxType r = max(min(v2, v3), max(v0, v1));
            const IdxType s = max(v2, v3);

            wmma::fragment<wmma::matrix_a, 8, 8, 4, ValType, wmma::row_major> a_frag_real_up;
            wmma::fragment<wmma::matrix_a, 8, 8, 4, ValType, wmma::row_major> a_frag_imag_up;
            wmma::fragment<wmma::matrix_a, 8, 8, 4, ValType, wmma::row_major> a_frag_real_dn;
            wmma::fragment<wmma::matrix_a, 8, 8, 4, ValType, wmma::row_major> a_frag_imag_dn;
            wmma::fragment<wmma::matrix_b, 8, 8, 4, ValType, wmma::col_major> b_frag_real;
            wmma::fragment<wmma::matrix_b, 8, 8, 4, ValType, wmma::col_major> b_frag_imag;
            wmma::fragment<wmma::accumulator, 8, 8, 4, ValType> c_frag_real_up;
            wmma::fragment<wmma::accumulator, 8, 8, 4, ValType> c_frag_imag_up;
            ;
            wmma::fragment<wmma::accumulator, 8, 8, 4, ValType> c_frag_real_dn;
            wmma::fragment<wmma::accumulator, 8, 8, 4, ValType> c_frag_imag_dn;
            ;

            // if (threadIdx.x == 0 && blockIdx.x == 0)
            // printf("GPU-%lld: dim:%lld, per_pe_work:%lld, nwarps:%d, start:%lld, end:%lld, step:%lld\n",i_proc, dim,per_pe_work, nwarps,
            //(i_proc)*per_pe_work+(warpid<<3), (i_proc+1)*per_pe_work, (nwarps<<3));

            for (IdxType i = (i_proc)*per_pe_work + (warpid << 3); i < (i_proc + 1) * per_pe_work; i += (nwarps << 3))
            {
                // if (laneid == 0) printf("=== warpid: %d ==== \n",warpid);
#pragma unroll
                for (int j = 0; j < 8; j++)
                {
                    const IdxType term = get_term(i + j, p, q, r, s);
                    if (laneid < 16)
                    {
                        IdxType addr = term + SV16IDX(laneid);
                        el_real_s[j * 16 + laneid] = PGAS_G(dm_real, addr);
                    }
                    else
                    {
                        IdxType addr = term + SV16IDX(laneid - 16);
                        el_imag_s[j * 16 + laneid - 16] = PGAS_G(dm_imag, addr);
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
                    // load A from const mem
                    load_matrix_sync(a_frag_real_up, &gm_real[c * 4], 16);
                    load_matrix_sync(a_frag_imag_up, &gm_imag[c * 4], 16);
                    load_matrix_sync(a_frag_real_dn, &gm_real[16 * 8 + c * 4], 16);
                    load_matrix_sync(a_frag_imag_dn, &gm_imag[16 * 8 + c * 4], 16);

                    // load B from shared mem
                    load_matrix_sync(b_frag_real, &el_real_s[c * 4], 16);
                    load_matrix_sync(b_frag_imag, &el_imag_s[c * 4], 16);
                    // complex multiplication
                    mma_sync(c_frag_imag_up, a_frag_real_up, b_frag_imag, c_frag_imag_up);
                    mma_sync(c_frag_imag_up, a_frag_imag_up, b_frag_real, c_frag_imag_up);
                    mma_sync(c_frag_real_up, a_frag_real_up, b_frag_real, c_frag_real_up);
                    a_frag_imag_up.x[0] = -a_frag_imag_up.x[0];
                    mma_sync(c_frag_real_up, a_frag_imag_up, b_frag_imag, c_frag_real_up);

                    mma_sync(c_frag_imag_dn, a_frag_real_dn, b_frag_imag, c_frag_imag_dn);
                    mma_sync(c_frag_imag_dn, a_frag_imag_dn, b_frag_real, c_frag_imag_dn);
                    mma_sync(c_frag_real_dn, a_frag_real_dn, b_frag_real, c_frag_real_dn);
                    a_frag_imag_dn.x[0] = -a_frag_imag_dn.x[0];
                    mma_sync(c_frag_real_dn, a_frag_imag_dn, b_frag_imag, c_frag_real_dn);
                }
                // Store first result per segment-C
                IdxType j0 = ((2 * laneid) & 7);
                IdxType k0 = ((2 * laneid) >> 3);
                const IdxType term0 = get_term(i + j0, p, q, r, s) + SV16IDX(k0);
                PGAS_P(dm_real, term0, c_frag_real_up.x[0]);
                PGAS_P(dm_imag, term0, c_frag_imag_up.x[0]);

                // Store second result per segment-C
                IdxType j1 = ((2 * laneid + 1) & 7);
                IdxType k1 = ((2 * laneid + 1) >> 3);
                const IdxType term1 = get_term(i + j1, p, q, r, s) + SV16IDX(k1);
                PGAS_P(dm_real, term1, c_frag_real_up.x[1]);
                PGAS_P(dm_imag, term1, c_frag_imag_up.x[1]);

                // Store first result per segment-C
                IdxType j2 = ((2 * laneid) & 7);
                IdxType k2 = ((2 * laneid) >> 3) + 8;
                const IdxType term2 = get_term(i + j2, p, q, r, s) + SV16IDX(k2);
                PGAS_P(dm_real, term2, c_frag_real_dn.x[0]);
                PGAS_P(dm_imag, term2, c_frag_imag_dn.x[0]);

                // Store second result per segment-C
                IdxType j3 = ((2 * laneid + 1) & 7);
                IdxType k3 = ((2 * laneid + 1) >> 3) + 8;
                const IdxType term3 = get_term(i + j3, p, q, r, s) + SV16IDX(k3);
                PGAS_P(dm_real, term3, c_frag_real_dn.x[1]);
                PGAS_P(dm_imag, term3, c_frag_imag_dn.x[1]);
            }
            // BARR_NVSHMEM;
        }

        //============== Unified 4-qubit Gate with TC V2 ================
        // This is with further tensorcore optimization
        __device__ __inline__ void C4TCV2_GATE(const ValType *gm_real, const ValType *gm_imag,
                                               const IdxType qubit0, const IdxType qubit1,
                                               const IdxType qubit2, const IdxType qubit3)
        {
            grid_group grid = this_grid();
            const int tid = blockDim.x * blockIdx.x + threadIdx.x;
            const int laneid = (threadIdx.x & 31);
            const int warpid = (tid >> 5);
            const int wid = (threadIdx.x >> 5);
            const int nwarps = ((blockDim.x * gridDim.x) >> 5);
            const IdxType per_pe_work = ((dim) >> (gpu_scale + 4));
            // Our design require at least (gpu_scale + 4 + log(8)) qubits, because
            // tensorcore demands [8,4]*[4,8], we need at least 8 on B's row dim.
            // This implies for 2 GPUs, we need at least 8 qubit for SV, thus 4 qubit for DM.
            // The smaller condition is difficult to handle due to the complex get_term function,
            // We need to figure out, after the get_term, which row/col it maps to and filter out.
            // But the row and col will be quite difficult to obtain.

            extern __shared__ ValType els[];
            ValType *el_real_s = &els[wid * 8 * 16 * 2];          // per warp 8*16 for real
            ValType *el_imag_s = &els[wid * 8 * 16 * 2 + 8 * 16]; // per warp 8*16 for imag

            assert(qubit0 != qubit1); // Non-cloning
            assert(qubit0 != qubit2); // Non-cloning
            assert(qubit0 != qubit3); // Non-cloning
            assert(qubit1 != qubit2); // Non-cloning
            assert(qubit1 != qubit3); // Non-cloning
            assert(qubit2 != qubit3); // Non-cloning

            // need to sort qubits: min->max: p, q, r, s
            const IdxType v0 = min(qubit0, qubit1);
            const IdxType v1 = min(qubit2, qubit3);
            const IdxType v2 = max(qubit0, qubit1);
            const IdxType v3 = max(qubit2, qubit3);
            const IdxType p = min(v0, v1);
            const IdxType q = min(min(v2, v3), max(v0, v1));
            const IdxType r = max(min(v2, v3), max(v0, v1));
            const IdxType s = max(v2, v3);

            wmma::fragment<wmma::matrix_a, 8, 8, 4, ValType, wmma::row_major> a_frag_real;
            wmma::fragment<wmma::matrix_a, 8, 8, 4, ValType, wmma::row_major> a_frag_imag;
            wmma::fragment<wmma::matrix_b, 8, 8, 4, ValType, wmma::col_major> b_frag_real;
            wmma::fragment<wmma::matrix_b, 8, 8, 4, ValType, wmma::col_major> b_frag_imag;
            wmma::fragment<wmma::accumulator, 8, 8, 4, ValType> c_frag_real;
            wmma::fragment<wmma::accumulator, 8, 8, 4, ValType> c_frag_imag;

            // if (threadIdx.x == 0 && blockIdx.x == 0)
            // printf("GPU-%lld: dim:%lld, per_pe_work:%lld, nwarps:%d, start:%lld, end:%lld, step:%lld\n",i_proc, dim,per_pe_work, nwarps,
            //(i_proc)*per_pe_work+(warpid<<3), (i_proc+1)*per_pe_work, (nwarps<<3));

            for (IdxType i = (i_proc)*per_pe_work + (warpid << 3); i < (i_proc + 1) * per_pe_work; i += (nwarps << 3))
            {
                // load from matrix-B
                for (int j = 0; j < 8; j++)
                {
                    const IdxType term = get_term(i + j, p, q, r, s);
                    if (laneid < 16)
                        el_real_s[j * 16 + laneid] = PGAS_G(dm_real, term + SV16IDX(laneid));
                    else
                        el_imag_s[j * 16 + laneid - 16] = PGAS_G(dm_imag, term + SV16IDX(laneid - 16));
                }
                __syncwarp();

                //============= For matrix-A upper-part ============
                wmma::fill_fragment(c_frag_real, 0.0);
                wmma::fill_fragment(c_frag_imag, 0.0);
                for (unsigned c = 0; c < 4; c++)
                {
                    // load A from const mem
                    load_matrix_sync(a_frag_real, &gm_real[c * 4], 16);
                    load_matrix_sync(a_frag_imag, &gm_imag[c * 4], 16);
                    // load B from shared mem
                    load_matrix_sync(b_frag_real, &el_real_s[c * 4], 16);
                    load_matrix_sync(b_frag_imag, &el_imag_s[c * 4], 16);
                    // complex multiplication
                    mma_sync(c_frag_imag, a_frag_real, b_frag_imag, c_frag_imag);
                    mma_sync(c_frag_imag, a_frag_imag, b_frag_real, c_frag_imag);
                    mma_sync(c_frag_real, a_frag_real, b_frag_real, c_frag_real);
                    a_frag_imag.x[0] = -a_frag_imag.x[0];
                    mma_sync(c_frag_real, a_frag_imag, b_frag_imag, c_frag_real);
                }
                // Store first result per segment-C
                IdxType j0 = ((2 * laneid) & 7);
                IdxType k0 = ((2 * laneid) >> 3);
                const IdxType term0 = get_term(i + j0, p, q, r, s) + SV16IDX(k0);
                PGAS_P(dm_real, term0, c_frag_real.x[0]);
                PGAS_P(dm_imag, term0, c_frag_imag.x[0]);
                // Store second result per segment-C
                IdxType j1 = ((2 * laneid + 1) & 7);
                IdxType k1 = ((2 * laneid + 1) >> 3);
                const IdxType term1 = get_term(i + j1, p, q, r, s) + SV16IDX(k1);
                PGAS_P(dm_real, term1, c_frag_real.x[1]);
                PGAS_P(dm_imag, term1, c_frag_imag.x[1]);

                //============= For matrix-A lower-part ============
                wmma::fill_fragment(c_frag_real, 0.0);
                wmma::fill_fragment(c_frag_imag, 0.0);
                for (unsigned c = 0; c < 4; c++)
                {
                    // load A from const mem
                    load_matrix_sync(a_frag_real, &gm_real[16 * 8 + c * 4], 16);
                    load_matrix_sync(a_frag_imag, &gm_imag[16 * 8 + c * 4], 16);
                    // load B from shared mem
                    load_matrix_sync(b_frag_real, &el_real_s[c * 4], 16);
                    load_matrix_sync(b_frag_imag, &el_imag_s[c * 4], 16);
                    // complex multiplication
                    mma_sync(c_frag_imag, a_frag_real, b_frag_imag, c_frag_imag);
                    mma_sync(c_frag_imag, a_frag_imag, b_frag_real, c_frag_imag);
                    mma_sync(c_frag_real, a_frag_real, b_frag_real, c_frag_real);
                    a_frag_imag.x[0] = -a_frag_imag.x[0];
                    mma_sync(c_frag_real, a_frag_imag, b_frag_imag, c_frag_real);
                }
                // Store first result per segment-C
                IdxType j2 = ((2 * laneid) & 7);
                IdxType k2 = ((2 * laneid) >> 3) + 8;
                const IdxType term2 = get_term(i + j2, p, q, r, s) + SV16IDX(k2);
                PGAS_P(dm_real, term2, c_frag_real.x[0]);
                PGAS_P(dm_imag, term2, c_frag_imag.x[0]);

                // Store second result per segment-C
                IdxType j3 = ((2 * laneid + 1) & 7);
                IdxType k3 = ((2 * laneid + 1) >> 3) + 8;
                const IdxType term3 = get_term(i + j3, p, q, r, s) + SV16IDX(k3);
                PGAS_P(dm_real, term3, c_frag_real.x[1]);
                PGAS_P(dm_imag, term3, c_frag_imag.x[1]);
            }
            // BARR_NVSHMEM;
        }

        //============== Unified 4-qubit Gate with TC ================
        // This is with tensorcore and nvshmem merged communication optimization
        __device__ __inline__ void C4TCV3_GATE(const ValType *gm_real, const ValType *gm_imag,
                                               const IdxType qubit0, const IdxType qubit1,
                                               const IdxType qubit2, const IdxType qubit3)
        {
            grid_group grid = this_grid();

            assert(qubit0 != qubit1); // Non-cloning
            assert(qubit0 != qubit2); // Non-cloning
            assert(qubit0 != qubit3); // Non-cloning
            assert(qubit1 != qubit2); // Non-cloning
            assert(qubit1 != qubit3); // Non-cloning
            assert(qubit2 != qubit3); // Non-cloning

            // need to sort qubits: min->max: p, q, r, s
            const IdxType v0 = min(qubit0, qubit1);
            const IdxType v1 = min(qubit2, qubit3);
            const IdxType v2 = max(qubit0, qubit1);
            const IdxType v3 = max(qubit2, qubit3);
            const IdxType p = min(v0, v1);
            const IdxType q = min(min(v2, v3), max(v0, v1));
            const IdxType r = max(min(v2, v3), max(v0, v1));
            const IdxType s = max(v2, v3);

            if (s < lg2_m_gpu) // all 4 qubits are local
            {
                C4TCV2_GATE(gm_real, gm_imag, qubit0, qubit1, qubit2, qubit3);
            }
            else // s qubit is non-local
            {
                const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
                const int laneid = (threadIdx.x & 31);
                const int warpid = (tid >> 5);
                const int wid = (threadIdx.x >> 5);
                const int nwarps = ((blockDim.x * gridDim.x) >> 5);
                const IdxType per_pe_num = ((dim) >> (gpu_scale));
                const IdxType whole_work = ((dim) >> 4);

                extern __shared__ ValType els[];
                ValType *el_real_s = &els[wid * 8 * 16 * 2];          // per warp 8*16 for real
                ValType *el_imag_s = &els[wid * 8 * 16 * 2 + 8 * 16]; // per warp 8*16 for imag

                // load data from pair GPU
                IdxType pair_gpu = (i_proc) ^ ((IdxType)1 << (s - (lg2_m_gpu)));
                if (i_proc > pair_gpu)
                    return;

                ValType *dm_real_remote = m_real;
                ValType *dm_imag_remote = m_imag;

                if (tid == 0)
                    nvshmem_double_get(dm_real_remote, dm_real, per_pe_num, pair_gpu);
                if (tid == 0)
                    nvshmem_double_get(dm_imag_remote, dm_imag, per_pe_num, pair_gpu);
                grid.sync();

                wmma::fragment<wmma::matrix_a, 8, 8, 4, ValType, wmma::row_major> a_frag_real;
                wmma::fragment<wmma::matrix_a, 8, 8, 4, ValType, wmma::row_major> a_frag_imag;
                wmma::fragment<wmma::matrix_b, 8, 8, 4, ValType, wmma::col_major> b_frag_real;
                wmma::fragment<wmma::matrix_b, 8, 8, 4, ValType, wmma::col_major> b_frag_imag;
                wmma::fragment<wmma::accumulator, 8, 8, 4, ValType> c_frag_real;
                wmma::fragment<wmma::accumulator, 8, 8, 4, ValType> c_frag_imag;

                for (IdxType i = (warpid << 3); i < whole_work; i += (nwarps << 3))
                {
                    const IdxType base = get_term(i, p, q, r, s);
                    if ((per_pe_num * (i_proc)) <= base && base < (per_pe_num * (i_proc + 1)))
                    {
                        for (int j = 0; j < 8; j++)
                        {
                            const IdxType term = get_term(i + j, p, q, r, s);
                            if (laneid < 16)
                            {
                                IdxType addr = term + SV16IDX(laneid);
                                // This is due to how SV16IDx is defined, qubit0 maps to MSB,
                                // qubit1 maps to MSB-1, qubit2 maps to MSB-2, qubit3 maps to LSB
                                if (qubit3 == s) // qubit3 is remote qubit
                                {
                                    if ((laneid & 1) != 0)
                                        el_real_s[j * 16 + laneid] = LOCAL_G_CUDA_MPI(dm_real_remote, addr);
                                    else
                                        el_real_s[j * 16 + laneid] = LOCAL_G_CUDA_MPI(dm_real, addr);
                                }
                                else // qubit2 is remote (not possible qubit0 or 1 is remote)
                                {
                                    if (((laneid >> 1) & 1) != 0)
                                        el_real_s[j * 16 + laneid] = LOCAL_G_CUDA_MPI(dm_real_remote, addr);
                                    else
                                        el_real_s[j * 16 + laneid] = LOCAL_G_CUDA_MPI(dm_real, addr);
                                }
                            }
                            else
                            {
                                IdxType addr = term + SV16IDX(laneid - 16);
                                if (qubit3 == s) // qubit3 is remote qubit
                                {
                                    if ((laneid & 1) != 0)
                                        el_imag_s[j * 16 + laneid] = LOCAL_G_CUDA_MPI(dm_imag_remote, addr);
                                    else
                                        el_imag_s[j * 16 + laneid] = LOCAL_G_CUDA_MPI(dm_imag, addr);
                                }
                                else // qubit2 is remote (not possible qubit0 or 1 is remote)
                                {
                                    if (((laneid >> 1) & 1) != 0)
                                        el_imag_s[j * 16 + laneid] = LOCAL_G_CUDA_MPI(dm_imag_remote, addr);
                                    else
                                        el_imag_s[j * 16 + laneid] = LOCAL_G_CUDA_MPI(dm_imag, addr);
                                    // el_imag_s[j*16+laneid] = PGAS_G(dm_imag, addr);
                                }
                            }
                        }
                        __syncwarp();

                        //============= For matrix-A upper-part ============
                        wmma::fill_fragment(c_frag_real, 0.0);
                        wmma::fill_fragment(c_frag_imag, 0.0);

                        for (unsigned c = 0; c < 4; c++)
                        {
                            // load A from const mem
                            load_matrix_sync(a_frag_real, &gm_real[c * 4], 16);
                            load_matrix_sync(a_frag_imag, &gm_imag[c * 4], 16);
                            // load B from shared mem
                            load_matrix_sync(b_frag_real, &el_real_s[c * 4], 16);
                            load_matrix_sync(b_frag_imag, &el_imag_s[c * 4], 16);
                            // complex multiplication
                            mma_sync(c_frag_imag, a_frag_real, b_frag_imag, c_frag_imag);
                            mma_sync(c_frag_imag, a_frag_imag, b_frag_real, c_frag_imag);
                            mma_sync(c_frag_real, a_frag_real, b_frag_real, c_frag_real);
                            a_frag_imag.x[0] = -a_frag_imag.x[0];
                            mma_sync(c_frag_real, a_frag_imag, b_frag_imag, c_frag_real);
                        }
                        // Store first result per segment-C
                        IdxType j0 = ((2 * laneid) & 7);
                        IdxType k0 = ((2 * laneid) >> 3);
                        const IdxType term0 = get_term(i + j0, p, q, r, s) + SV16IDX(k0);

                        if (qubit3 == s) // qubit3 is remote qubit
                        {
                            if ((k0 & 1) != 0)
                            {
                                LOCAL_P_CUDA_MPI(dm_real_remote, term0, c_frag_real.x[0]);
                                LOCAL_P_CUDA_MPI(dm_imag_remote, term0, c_frag_imag.x[0]);
                            }
                            else
                            {
                                LOCAL_P_CUDA_MPI(dm_real, term0, c_frag_real.x[0]);
                                LOCAL_P_CUDA_MPI(dm_imag, term0, c_frag_imag.x[0]);
                            }
                        }
                        else // qubit2 is remote (not possible qubit0 or 1 is remote)
                        {
                            if (((k0 >> 1) & 1) != 0)
                            {
                                LOCAL_P_CUDA_MPI(dm_real_remote, term0, c_frag_real.x[0]);
                                LOCAL_P_CUDA_MPI(dm_imag_remote, term0, c_frag_imag.x[0]);
                            }
                            else
                            {
                                LOCAL_P_CUDA_MPI(dm_real, term0, c_frag_real.x[0]);
                                LOCAL_P_CUDA_MPI(dm_imag, term0, c_frag_imag.x[0]);
                            }
                        }
                        // Store second result per segment-C
                        IdxType j1 = ((2 * laneid + 1) & 7);
                        IdxType k1 = ((2 * laneid + 1) >> 3);
                        const IdxType term1 = get_term(i + j1, p, q, r, s) + SV16IDX(k1);

                        if (qubit3 == s) // qubit3 is remote qubit
                        {
                            if ((k1 & 1) != 0)
                            {
                                LOCAL_P_CUDA_MPI(dm_real_remote, term1, c_frag_real.x[1]);
                                LOCAL_P_CUDA_MPI(dm_imag_remote, term1, c_frag_imag.x[1]);
                            }
                            else
                            {
                                LOCAL_P_CUDA_MPI(dm_real, term1, c_frag_real.x[1]);
                                LOCAL_P_CUDA_MPI(dm_imag, term1, c_frag_imag.x[1]);
                            }
                        }
                        else // qubit2 is remote (not possible qubit0 or 1 is remote)
                        {
                            if (((k1 >> 1) & 1) != 0)
                            {
                                LOCAL_P_CUDA_MPI(dm_real_remote, term1, c_frag_real.x[1]);
                                LOCAL_P_CUDA_MPI(dm_imag_remote, term1, c_frag_imag.x[1]);
                            }
                            else
                            {
                                LOCAL_P_CUDA_MPI(dm_real, term1, c_frag_real.x[1]);
                                LOCAL_P_CUDA_MPI(dm_imag, term1, c_frag_imag.x[1]);
                            }
                        }
                        //============= For matrix-A lower-part ============
                        wmma::fill_fragment(c_frag_real, 0.0);
                        wmma::fill_fragment(c_frag_imag, 0.0);
                        for (unsigned c = 0; c < 4; c++)
                        {
                            // load A from const mem
                            load_matrix_sync(a_frag_real, &gm_real[16 * 8 + c * 4], 16);
                            load_matrix_sync(a_frag_imag, &gm_imag[16 * 8 + c * 4], 16);
                            // load B from shared mem
                            load_matrix_sync(b_frag_real, &el_real_s[c * 4], 16);
                            load_matrix_sync(b_frag_imag, &el_imag_s[c * 4], 16);
                            // complex multiplication
                            mma_sync(c_frag_imag, a_frag_real, b_frag_imag, c_frag_imag);
                            mma_sync(c_frag_imag, a_frag_imag, b_frag_real, c_frag_imag);
                            mma_sync(c_frag_real, a_frag_real, b_frag_real, c_frag_real);
                            a_frag_imag.x[0] = -a_frag_imag.x[0];
                            mma_sync(c_frag_real, a_frag_imag, b_frag_imag, c_frag_real);
                        }
                        // Store first result per segment-C
                        IdxType j2 = ((2 * laneid) & 7);
                        IdxType k2 = ((2 * laneid) >> 3) + 8;
                        const IdxType term2 = get_term(i + j2, p, q, r, s) + SV16IDX(k2);

                        if (qubit3 == s) // qubit3 is remote qubit
                        {
                            if ((k2 & 1) != 0)
                            {
                                LOCAL_P_CUDA_MPI(dm_real_remote, term2, c_frag_real.x[0]);
                                LOCAL_P_CUDA_MPI(dm_imag_remote, term2, c_frag_imag.x[0]);
                            }
                            else
                            {
                                LOCAL_P_CUDA_MPI(dm_real, term2, c_frag_real.x[0]);
                                LOCAL_P_CUDA_MPI(dm_imag, term2, c_frag_imag.x[0]);
                            }
                        }
                        else // qubit2 is remote (not possible qubit0 or 1 is remote)
                        {
                            if (((k2 >> 1) & 1) != 0)
                            {
                                LOCAL_P_CUDA_MPI(dm_real_remote, term2, c_frag_real.x[0]);
                                LOCAL_P_CUDA_MPI(dm_imag_remote, term2, c_frag_imag.x[0]);
                            }
                            else
                            {
                                LOCAL_P_CUDA_MPI(dm_real, term2, c_frag_real.x[0]);
                                LOCAL_P_CUDA_MPI(dm_imag, term2, c_frag_imag.x[0]);
                            }
                        }
                        // Store second result per segment-C
                        IdxType j3 = ((2 * laneid + 1) & 7);
                        IdxType k3 = ((2 * laneid + 1) >> 3) + 8;
                        const IdxType term3 = get_term(i + j3, p, q, r, s) + SV16IDX(k3);

                        if (qubit3 == s) // qubit3 is remote qubit
                        {
                            if ((k3 & 1) != 0)
                            {
                                LOCAL_P_CUDA_MPI(dm_real_remote, term3, c_frag_real.x[1]);
                                LOCAL_P_CUDA_MPI(dm_imag_remote, term3, c_frag_imag.x[1]);
                            }
                            else
                            {
                                LOCAL_P_CUDA_MPI(dm_real, term3, c_frag_real.x[1]);
                                LOCAL_P_CUDA_MPI(dm_imag, term3, c_frag_imag.x[1]);
                            }
                        }
                        else // qubit2 is remote (not possible qubit0 or 1 is remote)
                        {
                            if (((k3 >> 1) & 1) != 0)
                            {
                                LOCAL_P_CUDA_MPI(dm_real_remote, term3, c_frag_real.x[1]);
                                LOCAL_P_CUDA_MPI(dm_imag_remote, term3, c_frag_imag.x[1]);
                            }
                            else
                            {
                                LOCAL_P_CUDA_MPI(dm_real, term3, c_frag_real.x[1]);
                                LOCAL_P_CUDA_MPI(dm_imag, term3, c_frag_imag.x[1]);
                            }
                        }
                    }
                }
                grid.sync();
                if (tid == 0)
                    nvshmem_double_put(dm_real, dm_real_remote, per_pe_num, pair_gpu);
                if (tid == 0)
                    nvshmem_double_put(dm_imag, dm_imag_remote, per_pe_num, pair_gpu);
                // BARR_NVSHMEM;
            }
        }
#endif
        ///*
        __device__ __inline__ void M_GATE(ValType *gm_real, ValType *gm_imag,
                                          const IdxType qubit, const IdxType cur_index)
        {
            grid_group grid = this_grid();
            assert(n_qubits < lg2_m_gpu);
            const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;

            IdxType mask = ((IdxType)1 << qubit);

            if (i_proc == 0)
            {
                for (IdxType i = tid; i < ((IdxType)1 << (n_qubits)); i += blockDim.x * gridDim.x)
                {
                    const ValType val = PGAS_G(dm_real, (i << (n_qubits)) + i);
                    if ((i & mask) == 0)
                        m_real[i] = 0;
                    else
                        m_real[i] = abs(val);
                }
                grid.sync();
                for (IdxType k = ((IdxType)1 << (n_qubits - 1)); k > 0; k >>= 1)
                {
                    for (IdxType i = tid; i < k; i += blockDim.x * gridDim.x)
                        m_real[i] += m_real[i + k];
                    grid.sync();
                }
            }
            BARR_NVSHMEM;
            if (tid == 0 && i_proc != 0)
                m_real[0] = PGAS_G(m_real, 0);
            grid.sync();
            ValType prob_of_one = m_real[0];
            grid.sync();

            ValType rand = randoms_gpu[cur_index];
            if (tid == 0)
            {
                if (rand <= prob_of_one)
                {
                    gm_real[15] = 1.0 / prob_of_one;
                }
                else
                {
                    gm_real[0] = 1.0 / (1.0 - prob_of_one);
                }
            }
            BARR_NVSHMEM;
            C2V1_GATE(gm_real, gm_imag, qubit, qubit + n_qubits);
            BARR_NVSHMEM;
            if (tid == 0)
                results_gpu[cur_index] = (rand <= prob_of_one ? 1 : 0);
        }
        //*/

        __device__ __inline__ void Normalization(ValType *dm_real, ValType *dm_imag)
        {
            grid_group grid = this_grid();
            assert(n_qubits < lg2_m_gpu);
            const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
            const IdxType per_pe_work_dm = ((dim) >> (gpu_scale));

            if (i_proc == 0)
            {
                for (IdxType i = tid; i < ((IdxType)1 << (n_qubits)); i += blockDim.x * gridDim.x)
                {
                    const ValType val = PGAS_G(dm_real, (i << (n_qubits)) + i);
                    m_real[i] = abs(val);
                }
                grid.sync();
                for (IdxType k = ((IdxType)1 << (n_qubits - 1)); k > 0; k >>= 1)
                {
                    for (IdxType i = tid; i < k; i += blockDim.x * gridDim.x)
                        m_real[i] += m_real[i + k];
                    grid.sync();
                }
            }
            BARR_NVSHMEM;
            if (tid == 0 && i_proc != 0)
                m_real[0] = PGAS_G(m_real, 0);
            grid.sync();
            ValType purity = m_real[0];
            for (IdxType i = tid; i < per_pe_work_dm; i += blockDim.x * gridDim.x)
            {
                dm_real[i] /= purity;
                dm_imag[i] /= purity;
            }
        }

        __device__ __inline__ void MA_GATE(const IdxType repetition, const IdxType cur_index)
        {
            grid_group grid = this_grid();
            assert(n_qubits < lg2_m_gpu); // ensure the diagonal can be fit into GPU-0's part of m_real
            const IdxType n_size = (IdxType)1 << (n_qubits);
            const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;

            BARR_NVSHMEM;
            if (i_proc == 0)
            {
                for (IdxType i = tid; i < n_size; i += blockDim.x * gridDim.x)
                    m_real[i] = abs(PGAS_G(dm_real, (i << (n_qubits)) + i));
                grid.sync();

                // Parallel prefix sum
                for (IdxType d = 0; d < (n_qubits); d++)
                {
                    IdxType step = (IdxType)1 << (d + 1);
                    for (IdxType k = tid * step; k < n_size; k += step * blockDim.x * gridDim.x)
                        m_real[k + (1 << (d + 1)) - 1] = m_real[k + (1 << d) - 1] + m_real[k + (1 << (d + 1)) - 1];
                    grid.sync();
                }
                if (tid == 0)
                {
                    ValType purity = fabs(m_real[n_size - 1]);
                    m_real[n_size - 1] = 0;
                    if (fabs(purity - 1.0) > ERROR_BAR)
                        printf("MA: Purity Check fails with %lf\n", purity);
                }
                grid.sync();
                for (IdxType d = (n_qubits)-1; d >= 0; d--)
                {
                    IdxType step = (IdxType)1 << (d + 1);
                    for (IdxType k = tid * step; k < n_size - 1; k += step * blockDim.x * gridDim.x)
                    {
                        ValType tmp = m_real[k + ((IdxType)1 << d) - 1];
                        m_real[k + ((IdxType)1 << d) - 1] = m_real[k + ((IdxType)1 << (d + 1)) - 1];
                        m_real[k + ((IdxType)1 << (d + 1)) - 1] = tmp + m_real[k + ((IdxType)1 << (d + 1)) - 1];
                    }
                    grid.sync();
                }
                grid.sync();

                for (IdxType j = tid; j < n_size; j += blockDim.x * gridDim.x)
                {
                    ValType lower = m_real[j];
                    ValType upper = (j + 1 == n_size) ? 1 : m_real[j + 1];
                    for (IdxType i = 0; i < repetition; i++)
                    {
                        ValType r = randoms_gpu[cur_index + i];
                        if (lower <= r && r < upper)
                            results_gpu[cur_index + i] = j;
                    }
                }
            }
            BARR_NVSHMEM;
            if (i_proc != 0 && tid == 0)
                nvshmem_longlong_get(results_gpu, results_gpu, repetition, 0);
            BARR_NVSHMEM;
        }

        __device__ __inline__ void RESET_GATE(const IdxType qubit)
        {
            grid_group grid = this_grid();
            assert(n_qubits < lg2_m_gpu); // ensure the diagonal can be fit into GPU-0's part of m_real
            const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
            const IdxType per_pe_work_dm = ((dim) >> (gpu_scale));

            IdxType mask = ((IdxType)1 << qubit);

            if (i_proc == 0)
            {
                for (IdxType i = tid; i < ((IdxType)1 << (n_qubits)); i += blockDim.x * gridDim.x)
                {
                    if ((i & mask) == 0) // for all conditions with qubit=0, we set it to 0, so we sum up all prob that qubit=1
                    {
                        m_real[i] = 0;
                    }
                    else
                    {
                        const ValType val = PGAS_G(dm_real, (i << (n_qubits)) + i);
                        PGAS_P(m_real, i, abs(val));
                    }
                }
                grid.sync();
                for (IdxType k = ((IdxType)1 << (n_qubits - 1)); k > 0; k >>= 1)
                {
                    for (IdxType i = tid; i < k; i += blockDim.x * gridDim.x)
                        m_real[i] += m_real[i + k];
                    grid.sync();
                }
            }
            BARR_NVSHMEM;
            if (tid == 0 && i_proc != 0)
                m_real[0] = PGAS_G(m_real, 0);
            grid.sync();
            ValType prob_of_one = m_real[0];

            if (prob_of_one < 1.0) // still possible to normalize
            {
                ValType factor = 1.0 / (1.0 - prob_of_one);
                for (IdxType i = tid; i < per_pe_work_dm; i += blockDim.x * gridDim.x)
                {
                    if ((i & mask) == 0)
                    {
                        dm_real[i] *= factor;
                        dm_imag[i] *= factor;
                    }
                    else
                    {
                        dm_real[i] = 0;
                        dm_imag[i] = 0;
                    }
                }
            }
            else
            {
                for (IdxType i = tid; i < per_pe_work_dm; i += blockDim.x * gridDim.x)
                {
                    if ((i & mask) == 0)
                    {
                        IdxType dual_i = i ^ mask;
                        dm_real[i] = dm_real[dual_i];
                        dm_imag[i] = dm_imag[dual_i];
                        dm_real[dual_i] = 0;
                        dm_imag[dual_i] = 0;
                    }
                }
            }
            BARR_NVSHMEM;
        }
    };

    __global__ void dm_simulation_kernel_cuda_mpi(DM_CUDA_MPI *dm_gpu, IdxType n_gates, IdxType n_qubits, bool enable_tc)
    {
        IdxType cur_index = 0;
        IdxType lg2_m_gpu = dm_gpu->lg2_m_gpu;
        grid_group grid = this_grid();

        for (IdxType t = 0; t < n_gates; t++)
        {
            OP op_name = (dm_gpu->gates_gpu)[t].op_name;
            IdxType qubit = (dm_gpu->gates_gpu)[t].qubit;

            IdxType ctrl = (dm_gpu->gates_gpu)[t].ctrl;
            ValType *gm_real = (dm_gpu->gates_gpu)[t].gm_real;
            ValType *gm_imag = (dm_gpu->gates_gpu)[t].gm_imag;

            IdxType repetition = (dm_gpu->gates_gpu)[t].qubit;

            grid_group grid = this_grid();
            // only need sync when operating on remote qubits
            if (((ctrl + n_qubits) >= lg2_m_gpu) || ((qubit + n_qubits) >= lg2_m_gpu))
            {
                grid.sync();
                if (threadIdx.x == 0 && blockIdx.x == 0)
                    nvshmem_barrier_all();
            }

            if (op_name == OP::RESET)
            {
                dm_gpu->RESET_GATE(qubit);
            }
            else if (op_name == OP::M)
            {
                dm_gpu->M_GATE(gm_real, gm_imag, qubit, cur_index);
                cur_index++;
            }
            else if (op_name == OP::MA)
            {
                dm_gpu->MA_GATE(repetition, cur_index);
                cur_index += repetition;
            }
            else if (op_name == OP::C2)
            {
                dm_gpu->C2V1_GATE(gm_real, gm_imag, qubit, qubit + (n_qubits));
            }
            else if (op_name == OP::C4)
            {

                if (((ctrl + n_qubits) >= lg2_m_gpu) && ((qubit + n_qubits) >= lg2_m_gpu))
                {
                    dm_gpu->SWAP_GATE(0, ctrl + (n_qubits));
                    BARR_NVSHMEM;
#ifdef FP64_TC_AVAILABLE
                    if (enable_tc)
                        dm_gpu->C4TCV3_GATE(gm_real, gm_imag, ctrl, qubit, 0, qubit + (n_qubits));
                    else
#endif
                        dm_gpu->C4V1_GATE(gm_real, gm_imag, ctrl, qubit, 0, qubit + (n_qubits));

                    BARR_NVSHMEM;
                    dm_gpu->SWAP_GATE(0, ctrl + (n_qubits));
                }
                else
                {
#ifdef FP64_TC_AVAILABLE
                    if (enable_tc)
                        dm_gpu->C4TCV3_GATE(gm_real, gm_imag, ctrl, qubit, ctrl + (n_qubits), qubit + (n_qubits));
                    else
#endif
                        dm_gpu->C4V1_GATE(gm_real, gm_imag, ctrl, qubit, ctrl + (n_qubits), qubit + (n_qubits));
                }
            }
            // only need sync when operating on remote qubits
            if (((ctrl + n_qubits) >= lg2_m_gpu) || ((qubit + n_qubits) >= lg2_m_gpu))
                if (threadIdx.x == 0 && blockIdx.x == 0)
                    nvshmem_barrier_all();
            grid.sync();
        }
    }
}
