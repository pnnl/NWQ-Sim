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
#include <sstream>

#include <hip/hip_runtime.h>
#include <hip/hip_cooperative_groups.h>

namespace NWQSim
{
    using namespace cooperative_groups;
    using namespace std;

    // Simulation kernel runtime
    class DM_HIP;
    __global__ void dm_simulation_kernel_hip(DM_HIP *dm_gpu, IdxType n_gates, IdxType n_qubits, bool tensor_core);

    class DM_HIP : public QuantumState
    {
    public:
        DM_HIP(IdxType _n_qubits) : QuantumState(_n_qubits)
        {
            // Initialize the GPU
            n_qubits = _n_qubits;
            dim = (IdxType)1 << (2 * n_qubits);
            half_dim = (IdxType)1 << (2 * n_qubits - 1);
            dm_size = dim * (IdxType)sizeof(ValType);

            // always be 0 since 1-MPI maps to 1-GPU
            i_proc = 0;
            hipSafeCall(hipSetDevice(i_proc));

            // CPU side initialization
            SAFE_ALOC_HOST_HIP(dm_real_cpu, dm_size);
            SAFE_ALOC_HOST_HIP(dm_imag_cpu, dm_size);
            memset(dm_real_cpu, 0, dm_size);
            memset(dm_imag_cpu, 0, dm_size);
            // State-vector initial state [0..0] = 1
            dm_real_cpu[0] = 1.;

            // NVSHMEM GPU memory allocation
            SAFE_ALOC_GPU_HIP(dm_real, dm_size);
            SAFE_ALOC_GPU_HIP(dm_imag, dm_size);
            SAFE_ALOC_GPU_HIP(m_real, dm_size + sizeof(ValType));
            SAFE_ALOC_GPU_HIP(m_imag, dm_size + sizeof(ValType));

            hipCheckError();
            gpu_mem = dm_size * 4;

            // GPU memory initilization
            hipSafeCall(hipMemcpy(dm_real, dm_real_cpu, dm_size,
                                  hipMemcpyHostToDevice));
            hipSafeCall(hipMemcpy(dm_imag, dm_imag_cpu, dm_size,
                                  hipMemcpyHostToDevice));
            hipSafeCall(hipMemset(m_real, 0, dm_size + sizeof(ValType)));
            hipSafeCall(hipMemset(m_imag, 0, dm_size + sizeof(ValType)));

            rng.seed(time(0));
        }

        ~DM_HIP()
        {
            // Release for CPU side
            SAFE_FREE_HOST_HIP(dm_real_cpu);
            SAFE_FREE_HOST_HIP(dm_imag_cpu);
            SAFE_FREE_HOST_HIP(randoms);
            SAFE_FREE_HOST_HIP(results);

            // Release for GPU side
            SAFE_FREE_GPU_HIP(dm_real);
            SAFE_FREE_GPU_HIP(dm_imag);
            SAFE_FREE_GPU_HIP(m_real);
            SAFE_FREE_GPU_HIP(m_imag);

            SAFE_FREE_GPU_HIP(gates_gpu);

            SAFE_FREE_GPU_HIP(randoms_gpu);
            SAFE_FREE_GPU_HIP(results_gpu);
        }

        void reset_state() override
        {
            // Reset CPU input & output
            memset(dm_real_cpu, 0, dm_size);
            memset(dm_imag_cpu, 0, dm_size);
            // State Vector initial state [0..0] = 1
            dm_real_cpu[0] = 1.;
            // GPU side initialization
            hipSafeCall(hipMemcpy(dm_real, dm_real_cpu,
                                  dm_size, hipMemcpyHostToDevice));
            hipSafeCall(hipMemcpy(dm_imag, dm_imag_cpu,
                                  dm_size, hipMemcpyHostToDevice));
            hipSafeCall(hipMemset(m_real, 0, dm_size + sizeof(ValType)));
            hipSafeCall(hipMemset(m_imag, 0, dm_size + sizeof(ValType)));
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

#ifdef HIP_TC_AVAILABLE
            is_tc = Config::ENABLE_TENSOR_CORE;
#endif

            DM_HIP *dm_gpu;
            SAFE_ALOC_GPU_HIP(dm_gpu, sizeof(DM_HIP));
            // Copy the simulator instance to GPU
            hipSafeCall(hipMemcpy(dm_gpu, this,
                                  sizeof(DM_HIP), hipMemcpyHostToDevice));

            double sim_time;
            gpu_timer sim_timer;

            if (Config::PRINT_SIM_TRACE)
            {
                printf("DMSIM_HIP is running! Requesting %lld qubits.\n", n_qubits);
            }

            dim3 gridDim(1, 1, 1);
            hipDeviceProp_t deviceProp;
            hipSafeCall(hipGetDeviceProperties(&deviceProp, 0));
            // 8*16 is per warp shared-memory usage for C4 TC, with real and imag
            unsigned smem_size;
            if (is_tc)
                smem_size = THREADS_CTA_HIP / 32 * 8 * 16 * 2 * sizeof(ValType);
            else
                smem_size = 0;
            int numBlocksPerSm;
            hipSafeCall(hipOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm,
                                                                     dm_simulation_kernel_hip, THREADS_CTA_HIP, smem_size));
            gridDim.x = numBlocksPerSm * deviceProp.multiProcessorCount;
            void *args[] = {&dm_gpu, &n_gates, &n_qubits, &is_tc};
            hipSafeCall(hipDeviceSynchronize());

            sim_timer.start_timer();

            hipLaunchCooperativeKernel((void *)dm_simulation_kernel_hip, gridDim,
                                       THREADS_CTA_HIP, args, smem_size, 0);
            hipSafeCall(hipDeviceSynchronize());

            sim_timer.stop_timer();
            sim_time = sim_timer.measure();

            // Copy the results back to CPU
            hipSafeCall(hipMemcpy(results, results_gpu, n_measures * sizeof(IdxType), hipMemcpyDeviceToHost));
            hipCheckError();

            if (Config::PRINT_SIM_TRACE)
            {
                printf("\n============== DM-Sim (HIP) ===============\n");
                printf("n_qubits:%lld, n_gates:%lld, sim_gates:%lld, ngpus:%d, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms, mem:%.3lf MB, mem_per_gpu:%.3lf MB\n",
                       n_qubits, origional_gates, n_gates, 1, sim_time, 0.,
                       sim_time, gpu_mem / 1024 / 1024, gpu_mem / 1024 / 1024);
                printf("=====================================\n");
            }
            SAFE_FREE_GPU_HIP(dm_gpu);
        }

        IdxType *get_results() override
        {
            return results;
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

        ValType get_exp_z(const std::vector<size_t> &in_bits) override
        {
            throw std::logic_error("get_exp_z Not implemented (DM_HIP)");
        }

        ValType get_exp_z() override
        {
            throw std::logic_error("get_exp_z Not implemented (DM_HIP)");
        }

        void print_res_state() override
        {
            hipSafeCall(hipMemcpy(dm_real_cpu, dm_real, dm_size, hipMemcpyDeviceToHost));
            hipSafeCall(hipMemcpy(dm_imag_cpu, dm_imag, dm_size, hipMemcpyDeviceToHost));

            IdxType num = ((IdxType)1 << n_qubits);
            printf("----- DMSim ------\n");
            for (IdxType i = 0; i < num; i++)
            {
                printf("(%.3lf,%.3lfj) ", dm_real_cpu[i], dm_imag_cpu[i]);
                if ((i + 1) % 8 == 0)
                    printf("\n");
            }
            printf("\n");
        }

    public:
        // n_qubits is the number of qubits
        IdxType n_qubits;

        IdxType dm_size;

        IdxType dim;
        IdxType half_dim;

        // GPU memory usage
        ValType gpu_mem;
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

        // GPU-side simulator instance
        DMGate *gates_gpu = NULL;

        void copy_gates_to_gpu(std::vector<DMGate> &cpu_vec)
        {
            // Allocate memory on CPU
            size_t vec_size = cpu_vec.size() * sizeof(DMGate);

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

        //================================= Gate Definition ========================================
        //============== Unified 2-qubit Gate without comm optimization ================
        __device__ __inline__ void C2_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit0, const IdxType qubit1)
        {
            grid_group grid = this_grid();
            const int tid = blockDim.x * blockIdx.x + threadIdx.x;
            const IdxType per_pe_work = ((dim) >> 2);
            assert(qubit0 != qubit1); // Non-cloning

            const IdxType q0dim = ((IdxType)1 << max(qubit0, qubit1));
            const IdxType q1dim = ((IdxType)1 << min(qubit0, qubit1));
            const IdxType outer_factor = ((dim) + q0dim + q0dim - 1) >> (max(qubit0, qubit1) + 1);
            const IdxType mider_factor = (q0dim + q1dim + q1dim - 1) >> (min(qubit0, qubit1) + 1);
            const IdxType inner_factor = q1dim;
            const IdxType qubit0_dim = ((IdxType)1 << qubit0);
            const IdxType qubit1_dim = ((IdxType)1 << qubit1);

            for (IdxType i = tid; i < per_pe_work; i += blockDim.x * gridDim.x)
            {
                IdxType outer = ((i / inner_factor) / (mider_factor)) * (q0dim + q0dim);
                IdxType mider = ((i / inner_factor) % (mider_factor)) * (q1dim + q1dim);
                IdxType inner = i % inner_factor;
                IdxType pos0 = outer + mider + inner;
                IdxType pos1 = outer + mider + inner + qubit1_dim;
                IdxType pos2 = outer + mider + inner + qubit0_dim;
                IdxType pos3 = outer + mider + inner + q0dim + q1dim;

                const ValType el0_real = LOCAL_G_HIP(dm_real, pos0);
                const ValType el0_imag = LOCAL_G_HIP(dm_imag, pos0);
                const ValType el1_real = LOCAL_G_HIP(dm_real, pos1);
                const ValType el1_imag = LOCAL_G_HIP(dm_imag, pos1);
                const ValType el2_real = LOCAL_G_HIP(dm_real, pos2);
                const ValType el2_imag = LOCAL_G_HIP(dm_imag, pos2);
                const ValType el3_real = LOCAL_G_HIP(dm_real, pos3);
                const ValType el3_imag = LOCAL_G_HIP(dm_imag, pos3);

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

                LOCAL_P_HIP(dm_real, pos0, dm_real_pos0);
                LOCAL_P_HIP(dm_real, pos1, dm_real_pos1);
                LOCAL_P_HIP(dm_real, pos2, dm_real_pos2);
                LOCAL_P_HIP(dm_real, pos3, dm_real_pos3);

                LOCAL_P_HIP(dm_imag, pos0, dm_imag_pos0);
                LOCAL_P_HIP(dm_imag, pos1, dm_imag_pos1);
                LOCAL_P_HIP(dm_imag, pos2, dm_imag_pos2);
                LOCAL_P_HIP(dm_imag, pos3, dm_imag_pos3);
            }
        }

        //============== Unified 4-qubit Gate ================
        __device__ __inline__ void C4_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit0, const IdxType qubit1,
                                           const IdxType qubit2, const IdxType qubit3)
        {
            grid_group grid = this_grid();
            const int tid = blockDim.x * blockIdx.x + threadIdx.x;
            const IdxType per_pe_work = ((dim) >> 4);
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
                    LOCAL_G_HIP(dm_real, term + SV16IDX(0)), LOCAL_G_HIP(dm_real, term + SV16IDX(1)),
                    LOCAL_G_HIP(dm_real, term + SV16IDX(2)), LOCAL_G_HIP(dm_real, term + SV16IDX(3)),
                    LOCAL_G_HIP(dm_real, term + SV16IDX(4)), LOCAL_G_HIP(dm_real, term + SV16IDX(5)),
                    LOCAL_G_HIP(dm_real, term + SV16IDX(6)), LOCAL_G_HIP(dm_real, term + SV16IDX(7)),
                    LOCAL_G_HIP(dm_real, term + SV16IDX(8)), LOCAL_G_HIP(dm_real, term + SV16IDX(9)),
                    LOCAL_G_HIP(dm_real, term + SV16IDX(10)), LOCAL_G_HIP(dm_real, term + SV16IDX(11)),
                    LOCAL_G_HIP(dm_real, term + SV16IDX(12)), LOCAL_G_HIP(dm_real, term + SV16IDX(13)),
                    LOCAL_G_HIP(dm_real, term + SV16IDX(14)), LOCAL_G_HIP(dm_real, term + SV16IDX(15))};
                const ValType el_imag[16] = {
                    LOCAL_G_HIP(dm_imag, term + SV16IDX(0)), LOCAL_G_HIP(dm_imag, term + SV16IDX(1)),
                    LOCAL_G_HIP(dm_imag, term + SV16IDX(2)), LOCAL_G_HIP(dm_imag, term + SV16IDX(3)),
                    LOCAL_G_HIP(dm_imag, term + SV16IDX(4)), LOCAL_G_HIP(dm_imag, term + SV16IDX(5)),
                    LOCAL_G_HIP(dm_imag, term + SV16IDX(6)), LOCAL_G_HIP(dm_imag, term + SV16IDX(7)),
                    LOCAL_G_HIP(dm_imag, term + SV16IDX(8)), LOCAL_G_HIP(dm_imag, term + SV16IDX(9)),
                    LOCAL_G_HIP(dm_imag, term + SV16IDX(10)), LOCAL_G_HIP(dm_imag, term + SV16IDX(11)),
                    LOCAL_G_HIP(dm_imag, term + SV16IDX(12)), LOCAL_G_HIP(dm_imag, term + SV16IDX(13)),
                    LOCAL_G_HIP(dm_imag, term + SV16IDX(14)), LOCAL_G_HIP(dm_imag, term + SV16IDX(15))};
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
                    LOCAL_P_HIP(dm_real, term + SV16IDX(j), res_real);
                    LOCAL_P_HIP(dm_imag, term + SV16IDX(j), res_imag);
                }
            }
            // BARR_HIP;
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

#ifdef HIP_TC_AVAILABLE

        // TC GATE OPERATIONS
#endif
        __device__ __inline__ void M_GATE(ValType *gm_real, ValType *gm_imag,
                                          const IdxType qubit, const IdxType cur_index)
        {
            grid_group grid = this_grid();
            const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;

            IdxType mask = ((IdxType)1 << qubit);

            for (IdxType i = tid; i < ((IdxType)1 << (n_qubits)); i += blockDim.x * gridDim.x)
            {
                const ValType val = LOCAL_G_HIP(dm_real, (i << (n_qubits)) + i);
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
            BARR_HIP;
            C2_GATE(gm_real, gm_imag, qubit, qubit + n_qubits);
            BARR_HIP;
            if (tid == 0)
                results_gpu[cur_index] = (rand <= prob_of_one ? 1 : 0);
        }

        __device__ __inline__ void MA_GATE(const IdxType repetition, const IdxType cur_index)
        {
            grid_group grid = this_grid();
            const IdxType n_size = (IdxType)1 << (n_qubits);
            const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;

            BARR_HIP;

            for (IdxType i = tid; i < n_size; i += blockDim.x * gridDim.x)
                m_real[i] = abs(LOCAL_G_HIP(dm_real, (i << (n_qubits)) + i));

            BARR_HIP;

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
                if (abs(purity - 1.0) > ERROR_BAR)
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
            BARR_HIP;
        }

        __device__ __inline__ void RESET_GATE(const IdxType qubit)
        {
            grid_group grid = this_grid();
            const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
            const IdxType per_pe_work_dm = (dim);

            IdxType mask = ((IdxType)1 << qubit);
            mask = (mask<<n_qubits) + mask;

            for (IdxType i = tid; i < per_pe_work_dm; i += blockDim.x * gridDim.x)
            {
                if ((i & mask) == 0)
                {
                    IdxType dual_i = i ^ mask;
                    dm_real[i] += dm_real[dual_i];
                    dm_imag[i] += dm_imag[dual_i];
                }
                else
                {
                    dm_real[i] = 0;
                    dm_imag[i] = 0;
                }
            }
            BARR_HIP;
        }
    };

    __global__ void dm_simulation_kernel_hip(DM_HIP *dm_gpu, IdxType n_gates, IdxType n_qubits, bool tensor_core)
    {
        IdxType cur_index = 0;
        grid_group grid = this_grid();

        for (IdxType t = 0; t < n_gates; t++)
        {
            auto op_name = (dm_gpu->gates_gpu)[t].op_name;
            auto qubit = (dm_gpu->gates_gpu)[t].qubit;

            auto ctrl = (dm_gpu->gates_gpu)[t].ctrl;
            auto gm_real = (dm_gpu->gates_gpu)[t].gm_real;
            auto gm_imag = (dm_gpu->gates_gpu)[t].gm_imag;

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
                dm_gpu->MA_GATE(qubit, cur_index);
                cur_index += qubit;
            }
            else if (op_name == OP::C2)
            {
                dm_gpu->C2_GATE(gm_real, gm_imag, qubit, qubit + (n_qubits));
            }
            else if (op_name == OP::C4)
            {
#ifdef HIP_TC_AVAILABLE
                // TC GATE OPERATIONS
                if (tensor_core)
                    dm_gpu->C4TCV1_GATE(gm_real, gm_imag, ctrl, qubit, ctrl + (n_qubits), qubit + (n_qubits));
                else
#endif
                    dm_gpu->C4_GATE(gm_real, gm_imag, ctrl, qubit, ctrl + (n_qubits), qubit + (n_qubits));
            }
            grid.sync();
        }
    }
} // namespace NWQSim
