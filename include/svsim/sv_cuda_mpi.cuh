#pragma once

#include "../state.hpp"
#include "../nwq_util.hpp"
#include "../gate.hpp"
#include "../circuit.hpp"

#include "../config.hpp"
#include "../private/cuda_util.cuh"
#include "../private/macros.hpp"
#include "../private/sim_gate.hpp"
#include "../private/exp_gate_declarations.cuh"

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
#include <mma.h>
#include <mpi.h>

#include <nvshmem.h>
#include <nvshmemx.h>
// #include <nvshmemx_error.h>

namespace NWQSim
{
    using namespace cooperative_groups;
    using namespace std;

    // Simulation kernel runtime
    class SV_CUDA_MPI;
    __global__ void simulation_kernel_cuda_mpi(SV_CUDA_MPI *sv_gpu, IdxType n_gates);

    class SV_CUDA_MPI : public QuantumState
    {
    public:
        SV_CUDA_MPI(IdxType _n_qubits) : QuantumState(SimType::SV)
        {
            // Initialize the GPU
            n_qubits = _n_qubits;
            dim = (IdxType)1 << (n_qubits);
            half_dim = (IdxType)1 << (n_qubits - 1);
            sv_size = dim * (IdxType)sizeof(ValType);

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
            lg2_m_gpu = n_qubits - gpu_scale;
            m_gpu = ((IdxType)1 << (lg2_m_gpu));

            sv_size_per_gpu = sv_size / n_gpus;
            sv_offset = i_proc * dim / n_gpus;
            dim_per_gpu = dim / n_gpus;

            // CPU side initialization
            if (!is_power_of_2(n_gpus))
                throw std::logic_error("Error: Number of GPUs should be power of 2.");
            if (dim % n_gpus != 0)
                throw std::logic_error("Error: Number of GPUs is too large or too small.");
            if (lg2_m_gpu < 5)
                throw std::logic_error("Error: Each GPU should have at least 5 qubits for multi-node version. Please increase qubits or reduce GPUs");
            // CPU side initialization
            SAFE_ALOC_HOST_CUDA(sv_real_cpu, sv_size_per_gpu);
            SAFE_ALOC_HOST_CUDA(sv_imag_cpu, sv_size_per_gpu);
            memset(sv_real_cpu, 0, sv_size_per_gpu);
            memset(sv_imag_cpu, 0, sv_size_per_gpu);
            // State-vector initial state [0..0] = 1
            if (i_proc == 0)
                sv_real_cpu[0] = 1.;
            // NVSHMEM GPU memory allocation
            sv_real = (ValType *)nvshmem_malloc(sv_size_per_gpu);
            sv_imag = (ValType *)nvshmem_malloc(sv_size_per_gpu);
            m_real = (ValType *)nvshmem_malloc(sv_size_per_gpu);
            m_imag = (ValType *)nvshmem_malloc(sv_size_per_gpu);
            cudaCheckError();
            gpu_mem += sv_size_per_gpu * 4;

            // GPU memory initilization
            cudaSafeCall(cudaMemcpy(sv_real, sv_real_cpu, sv_size_per_gpu,
                                    cudaMemcpyHostToDevice));
            cudaSafeCall(cudaMemcpy(sv_imag, sv_imag_cpu, sv_size_per_gpu,
                                    cudaMemcpyHostToDevice));
            cudaSafeCall(cudaMemset(m_real, 0, sv_size_per_gpu));
            cudaSafeCall(cudaMemset(m_imag, 0, sv_size_per_gpu));
            rng.seed(Config::RANDOM_SEED);
        }

        ~SV_CUDA_MPI()
        {
            // Release for CPU side
            SAFE_FREE_HOST_CUDA(sv_real_cpu);
            SAFE_FREE_HOST_CUDA(sv_imag_cpu);

            // Release for GPU side
            nvshmem_free(sv_real);
            nvshmem_free(sv_imag);
            nvshmem_free(m_real);
            nvshmem_free(m_imag);

            SAFE_FREE_HOST_CUDA(results);

            // SAFE_FREE_GPU(randoms_gpu);
            // SAFE_FREE_GPU(gates_gpu);
            // SAFE_FREE_HOST_CUDA(randoms);
            // SAFE_FREE_HOST_CUDA(results);
            nvshmem_finalize();
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
            cudaSafeCall(cudaMemcpy(sv_real, sv_real_cpu,
                                    sv_size_per_gpu, cudaMemcpyHostToDevice));
            cudaSafeCall(cudaMemcpy(sv_imag, sv_imag_cpu,
                                    sv_size_per_gpu, cudaMemcpyHostToDevice));
            cudaSafeCall(cudaMemset(m_real, 0, sv_size_per_gpu));
            cudaSafeCall(cudaMemset(m_imag, 0, sv_size_per_gpu));
        }

        void set_seed(IdxType seed) override
        {
            rng.seed(seed);
        }
        /**
         * @brief Read an initial quantum state (TODO: Check for state purity/validity)
         * @note
         * @param fpath: Path to the input binary file
         * @param format: Format of the input (either "sv" or "dm")
         */
        virtual void set_initial(std::string fpath, std::string format) override
        {
            std::ifstream instream;
            if (format != "sv")
            {
                throw std::runtime_error("SV-Sim only supports statevector input states\n");
            }
            instream.open(fpath, std::ios::in | std::ios::binary);
            if (instream.is_open())
            {
                instream.seekg(sv_size_per_gpu * i_proc);
                instream.read((char *)sv_real_cpu, sv_size_per_gpu);
                instream.seekg(dim * sizeof(ValType) + sv_size_per_gpu * i_proc);
                instream.read((char *)sv_imag_cpu, sv_size_per_gpu);
                cudaSafeCall(cudaMemcpy(sv_real, sv_real_cpu,
                                        sv_size_per_gpu, cudaMemcpyHostToDevice));
                cudaSafeCall(cudaMemcpy(sv_imag, sv_imag_cpu,
                                        sv_size_per_gpu, cudaMemcpyHostToDevice));
                instream.close();
            }
            std::cout << i_proc << std::endl;
        }

        /**
         * @brief Dump the current quantum state to a binary output file
         * @note
         * @param outpath: Output path
         */
        virtual void dump_res_state(std::string outpath) override
        {
            std::ofstream outstream;

            IdxType ticket = 1;
            // synchronize the file writes with a basic point-point ticket lock
            if (i_proc != 0)
            {
                MPI_Recv(&ticket, 1, MPI_INT64_T, i_proc - 1, i_proc, comm_global, MPI_STATUS_IGNORE);
                outstream.open(outpath, std::ios::app | std::ios::binary);
            }
            else
            {
                outstream.open(outpath, std::ios::trunc | std::ios::binary); // remove existing file
            }
            if (!outstream.is_open())
            {

                MPI_Send(&ticket, 1, MPI_INT64_T, i_proc + 1, i_proc + 1, comm_global);
                if (i_proc == 0)
                    std::cout << "Could not open file " << outpath << std::endl;
                return;
            }
            cudaSafeCall(cudaMemcpy(sv_real_cpu, sv_real, sv_size_per_gpu, cudaMemcpyDeviceToHost));
            // append to the end of the file
            outstream.write((char *)sv_real_cpu, sv_size_per_gpu);
            // outstream.write((char*)dm_imag_cpu, sizeof(ValType) * dm_size_per_gpu);
            // now write the imaginary part
            outstream.flush();
            outstream.close(); // close to flush the stream
            if (i_proc != n_gpus - 1)
            {
                MPI_Send(&ticket, 1, MPI_INT64_T, i_proc + 1, i_proc + 1, comm_global);
            }
            // synchronize the file writes with a basic point-point ticket lock
            if (i_proc != 0)
            {
                MPI_Recv(&ticket, 1, MPI_INT64_T, i_proc - 1, i_proc, comm_global, MPI_STATUS_IGNORE);
            }

            // reopen the file with the changes from the other threads
            outstream.open(outpath, std::ios::app | std::ios::binary);
            cudaSafeCall(cudaMemcpy(sv_imag_cpu, sv_imag, sv_size_per_gpu, cudaMemcpyDeviceToHost));
            // outstream.write((char*)dm_real_cpu, sizeof(ValType) * dm_size_per_gpu);
            outstream.write((char *)sv_imag_cpu, sv_size_per_gpu);

            outstream.flush();
            outstream.close();
            if (i_proc != n_gpus - 1)
            {
                MPI_Send(&ticket, 1, MPI_INT64_T, i_proc + 1, i_proc + 1, comm_global);
            }
        };
        void sim(std::shared_ptr<NWQSim::Circuit> circuit) override
        {
            IdxType origional_gates = circuit->num_gates();

            IdxType n_measures = prepare_measure(circuit->get_gates());

            std::vector<SVGate> cpu_vec = fuse_circuit_sv(circuit);

            // Copy the circuit to GPU
            copy_gates_to_gpu(cpu_vec);

            IdxType n_gates = cpu_vec.size();

            SV_CUDA_MPI *sv_gpu;
            SAFE_ALOC_GPU(sv_gpu, sizeof(SV_CUDA_MPI));

            // Copy the simulator instance to GPU
            cudaSafeCall(cudaMemcpy(sv_gpu, this,
                                    sizeof(SV_CUDA_MPI), cudaMemcpyHostToDevice));

            double *sim_times;
            double sim_time;
            gpu_timer sim_timer;

            if (Config::PRINT_SIM_TRACE && i_proc == 0)
            {
                SAFE_ALOC_HOST_CUDA(sim_times, sizeof(double) * n_gpus);
                memset(sim_times, 0, sizeof(double) * n_gpus);
                printf("SVSim_gpu is running! Requesting %lld qubits.\n", circuit->num_qubits());
            }

            dim3 gridDim(1, 1, 1);
            cudaDeviceProp deviceProp;
            cudaSafeCall(cudaGetDeviceProperties(&deviceProp, 0));
            // 8*16 is per warp shared-memory usage for C4 TC, with real and imag
            // unsigned smem_size = THREADS_CTA_CUDA/32*8*16*2*sizeof(ValType);
            unsigned smem_size = 0 * sizeof(ValType);
            int numBlocksPerSm;
            cudaSafeCall(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm,
                                                                       simulation_kernel_cuda_mpi, THREADS_CTA_CUDA, smem_size));
            gridDim.x = numBlocksPerSm * deviceProp.multiProcessorCount;
            void *args[] = {&sv_gpu, &n_gates};
            cudaSafeCall(cudaDeviceSynchronize());
            if (Config::PRINT_SIM_TRACE)
                MPI_Barrier(MPI_COMM_WORLD);
            sim_timer.start_timer();

            NVSHMEM_CHECK(nvshmemx_collective_launch((const void *)simulation_kernel_cuda_mpi, gridDim,
                                                     THREADS_CTA_CUDA, args, smem_size, 0));
            cudaSafeCall(cudaDeviceSynchronize());

            sim_timer.stop_timer();
            sim_time = sim_timer.measure();

            // Copy the results back to CPU

            if (n_measures > 0)
            {
                cudaSafeCall(cudaMemcpy(results, results_gpu, n_measures * sizeof(IdxType), cudaMemcpyDeviceToHost));

                nvshmem_free(results_gpu);

                SAFE_FREE_HOST_CUDA(randoms);
                SAFE_FREE_GPU(randoms_gpu);
            }

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
                    printf("\n============== SV-Sim ===============\n");
                    printf("n_qubits:%lld, n_gates:%lld, sim_gates:%lld, ngpus:%lld, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms, mem:%.3lf MB, mem_per_gpu:%.3lf MB\n",
                           n_qubits, origional_gates, n_gates, n_gpus, avg_sim_time, 0.,
                           avg_sim_time, gpu_mem / 1024 / 1024 * n_gpus, gpu_mem / 1024 / 1024);
                    printf("=====================================\n");
                    SAFE_FREE_HOST_CUDA(sim_times);
                }
            }

            SAFE_FREE_GPU(sv_gpu);
        }

        IdxType *get_results() override
        {
            return results;
        }

        ValType get_exp_z(const std::vector<size_t> &in_bits) override
        {
            // copy these vectors to device
            size_t *in_bits_gpu;

            SAFE_ALOC_GPU(in_bits_gpu, in_bits.size() * sizeof(size_t));
            cudaSafeCall(cudaMemcpy(in_bits_gpu, in_bits.data(), in_bits.size() * sizeof(size_t), cudaMemcpyHostToDevice));

            // result
            double result = 0.0;
            double *result_gpu;
            SAFE_ALOC_GPU(result_gpu, sizeof(ValType));
            cudaSafeCall(cudaMemcpy(result_gpu, &result, sizeof(double), cudaMemcpyHostToDevice));

            // call kernel
            int threadsPerBlock = THREADS_CTA_CUDA; // Choose based on your GPU's capability
            int blocksPerGrid = (dim_per_gpu + threadsPerBlock - 1) / threadsPerBlock;
            int sharedMemorySize = threadsPerBlock * sizeof(double); // Size of shared memory

            // Launch the kernel
            gpu_exp_z_bits<<<blocksPerGrid, threadsPerBlock, sharedMemorySize>>>(in_bits_gpu, in_bits.size(), sv_real, sv_imag, result_gpu, dim_per_gpu, i_proc * dim_per_gpu);

            // copy result back
            cudaSafeCall(cudaMemcpy(&result, result_gpu, sizeof(ValType), cudaMemcpyDeviceToHost));

            SAFE_FREE_GPU(in_bits_gpu);
            SAFE_FREE_GPU(result_gpu);

            double final_result;
            MPI_Allreduce(&result, &final_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            return final_result;
        }

        ValType get_exp_z() override
        {
            // result
            double result = 0.0;
            double *result_gpu;
            SAFE_ALOC_GPU(result_gpu, sizeof(double));
            cudaSafeCall(cudaMemcpy(result_gpu, &result, sizeof(double), cudaMemcpyHostToDevice));

            // call kernel
            int threadsPerBlock = THREADS_CTA_CUDA; // Choose based on your GPU's capability
            int blocksPerGrid = (dim_per_gpu + threadsPerBlock - 1) / threadsPerBlock;
            int sharedMemorySize = threadsPerBlock * sizeof(double); // Size of shared memory

            // Launch the kernel
            gpu_exp_z<<<blocksPerGrid, threadsPerBlock, sharedMemorySize>>>(sv_real, sv_imag, result_gpu, dim_per_gpu, i_proc * dim_per_gpu);

            // copy result back
            cudaSafeCall(cudaMemcpy(&result, result_gpu, sizeof(double), cudaMemcpyDeviceToHost));

            SAFE_FREE_GPU(result_gpu);

            double final_result;
            MPI_Allreduce(&result, &final_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            return final_result;
        }

        void save_state() override
        {
            // Save current state to CPU memory
            cudaSafeCall(cudaMemcpy(sv_real_cpu, sv_real, sv_size_per_gpu, cudaMemcpyDeviceToHost));
            cudaSafeCall(cudaMemcpy(sv_imag_cpu, sv_imag, sv_size_per_gpu, cudaMemcpyDeviceToHost));
        }

        void load_state() override
        {
            // Load state from CPU memory
            cudaSafeCall(cudaMemcpy(sv_real, sv_real_cpu, sv_size_per_gpu, cudaMemcpyHostToDevice));
            cudaSafeCall(cudaMemcpy(sv_imag, sv_imag_cpu, sv_size_per_gpu, cudaMemcpyHostToDevice));
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
            cudaSafeCall(cudaMemcpy(sv_real_cpu, sv_real, sv_size_per_gpu, cudaMemcpyDeviceToHost));
            cudaSafeCall(cudaMemcpy(sv_imag_cpu, sv_imag, sv_size_per_gpu, cudaMemcpyDeviceToHost));

            ValType *sv_diag_real = NULL;
            ValType *sv_diag_imag = NULL;
            if (i_proc == 0)
                SAFE_ALOC_HOST_CUDA(sv_diag_real, dim * sizeof(ValType));
            if (i_proc == 0)
                SAFE_ALOC_HOST_CUDA(sv_diag_imag, dim * sizeof(ValType));

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
                SAFE_FREE_HOST_CUDA(sv_diag_real);
                SAFE_FREE_HOST_CUDA(sv_diag_imag);
            }
        }

    public:
        // n_qubits is the number of qubits
        IdxType n_qubits;

        IdxType gpu_scale;
        IdxType n_gpus;
        IdxType lg2_m_gpu;
        IdxType m_gpu;
        IdxType sv_size_per_gpu;

        // gpu_scale is 2^x of the number of GPUs, e.g., with 8 GPUs the gpu_scale is 3 (2^3=8)
        IdxType dim;
        IdxType half_dim;
        IdxType dim_per_gpu;
        IdxType sv_size;
        IdxType sv_offset;
        // CPU arrays
        ValType *sv_real_cpu;
        ValType *sv_imag_cpu;
        // GPU arrays
        ValType *sv_real;
        ValType *sv_imag;
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
        SVGate *gates_gpu = NULL;

        void copy_gates_to_gpu(std::vector<SVGate> &cpu_vec)
        {
            // Allocate memory on CPU
            size_t vec_size = cpu_vec.size() * sizeof(SVGate);

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

            if (n_measures > 0)
            {
                SAFE_FREE_HOST_CUDA(results);

                // Prepare randoms and results memory
                SAFE_ALOC_HOST_CUDA(results, sizeof(IdxType) * n_measures);
                memset(results, 0, sizeof(IdxType) * n_measures);

                results_gpu = (IdxType *)nvshmem_malloc(sizeof(IdxType) * n_measures);
                cudaSafeCall(cudaMemset(results_gpu, 0, sizeof(IdxType) * n_measures));

                SAFE_ALOC_HOST_CUDA(randoms, sizeof(ValType) * n_measures);

                for (IdxType i = 0; i < n_measures; i++)
                    randoms[i] = uni_dist(rng);

                SAFE_ALOC_GPU(randoms_gpu, sizeof(ValType) * n_measures);
                cudaSafeCall(cudaMemcpy(randoms_gpu, randoms,
                                        sizeof(ValType) * n_measures, cudaMemcpyHostToDevice));
            }

            return n_measures;
        }

        //============== Local Unified 1-qubit Gate ================
        __device__ __inline__ void C1_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit)
        {
            grid_group grid = this_grid();
            const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
            const IdxType per_pe_work = ((half_dim) >> (gpu_scale));

            for (IdxType i = (i_proc)*per_pe_work + tid; i < (i_proc + 1) * per_pe_work; i += blockDim.x * gridDim.x)
            {
                IdxType outer = (i >> qubit);
                IdxType inner = (i & (((IdxType)1 << qubit) - 1));
                IdxType offset = (outer << (qubit + 1));
                IdxType pos0 = ((offset + inner) & (m_gpu - 1));
                IdxType pos1 = ((offset + inner + ((IdxType)1 << qubit)) & (m_gpu - 1));
                const ValType el0_real = LOCAL_G_CUDA_MPI(sv_real, pos0);
                const ValType el0_imag = LOCAL_G_CUDA_MPI(sv_imag, pos0);
                const ValType el1_real = LOCAL_G_CUDA_MPI(sv_real, pos1);
                const ValType el1_imag = LOCAL_G_CUDA_MPI(sv_imag, pos1);
                ValType sv_real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag);
                ValType sv_imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real);
                ValType sv_real_pos1 = (gm_real[2] * el0_real) - (gm_imag[2] * el0_imag) + (gm_real[3] * el1_real) - (gm_imag[3] * el1_imag);
                ValType sv_imag_pos1 = (gm_real[2] * el0_imag) + (gm_imag[2] * el0_real) + (gm_real[3] * el1_imag) + (gm_imag[3] * el1_real);
                LOCAL_P_CUDA_MPI(sv_real, pos0, sv_real_pos0);
                LOCAL_P_CUDA_MPI(sv_imag, pos0, sv_imag_pos0);
                LOCAL_P_CUDA_MPI(sv_real, pos1, sv_real_pos1);
                LOCAL_P_CUDA_MPI(sv_imag, pos1, sv_imag_pos1);
            }
            grid.sync();
        }

        /** The C1V1 version relies on bidirectional communication between each GPU pair. However, on Perlmutter ss11 with the recent slingshot-V2 upgrade, bidirectional communication
          encounters a bug. If that is the case, we use C1V2.
        */
        __device__ __inline__ void C1V1_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit)
        {
            grid_group grid = this_grid();
            const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
            const IdxType per_pe_work = ((half_dim) >> (gpu_scale));

            if (qubit < lg2_m_gpu)
            {
                C1_GATE(gm_real, gm_imag, qubit);
            }
            else
            {
                IdxType pair_gpu = (i_proc) ^ ((IdxType)1 << (qubit - (lg2_m_gpu)));
                IdxType flag = (i_proc < pair_gpu) ? 0 : 1;
                ValType *sv_real_pos0 = NULL;
                ValType *sv_imag_pos0 = NULL;
                ValType *sv_real_pos1 = NULL;
                ValType *sv_imag_pos1 = NULL;
                if (flag == 0) // copy right to local
                {
                    sv_real_pos0 = sv_real;
                    sv_imag_pos0 = sv_imag;
                    sv_real_pos1 = m_real;
                    sv_imag_pos1 = m_imag;
                    if (tid == 0)
                        nvshmem_double_get(sv_real_pos1, sv_real, per_pe_work, pair_gpu);
                    if (tid == 0)
                        nvshmem_double_get(sv_imag_pos1, sv_imag, per_pe_work, pair_gpu);
                }
                else // copy left to local
                {
                    sv_real_pos0 = m_real;
                    sv_imag_pos0 = m_imag;
                    sv_real_pos1 = &sv_real[per_pe_work];
                    sv_imag_pos1 = &sv_imag[per_pe_work];
                    if (tid == 0)
                        nvshmem_double_get(sv_real_pos0, &sv_real[per_pe_work], per_pe_work, pair_gpu);
                    if (tid == 0)
                        nvshmem_double_get(sv_imag_pos0, &sv_imag[per_pe_work], per_pe_work, pair_gpu);
                }
                grid.sync();
                for (IdxType i = tid; i < per_pe_work; i += blockDim.x * gridDim.x)
                {
                    const ValType el0_real = LOCAL_G_CUDA_MPI(sv_real_pos0, i);
                    const ValType el0_imag = LOCAL_G_CUDA_MPI(sv_imag_pos0, i);
                    const ValType el1_real = LOCAL_G_CUDA_MPI(sv_real_pos1, i);
                    const ValType el1_imag = LOCAL_G_CUDA_MPI(sv_imag_pos1, i);
                    ValType real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag);
                    ValType imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real);
                    ValType real_pos1 = (gm_real[2] * el0_real) - (gm_imag[2] * el0_imag) + (gm_real[3] * el1_real) - (gm_imag[3] * el1_imag);
                    ValType imag_pos1 = (gm_real[2] * el0_imag) + (gm_imag[2] * el0_real) + (gm_real[3] * el1_imag) + (gm_imag[3] * el1_real);
                    LOCAL_P_CUDA_MPI(sv_real_pos0, i, real_pos0);
                    LOCAL_P_CUDA_MPI(sv_imag_pos0, i, imag_pos0);
                    LOCAL_P_CUDA_MPI(sv_real_pos1, i, real_pos1);
                    LOCAL_P_CUDA_MPI(sv_imag_pos1, i, imag_pos1);
                }
                grid.sync();
                if (flag == 0) // copy local to right
                {
                    if (tid == 0)
                        nvshmem_double_put(sv_real, sv_real_pos1, per_pe_work, pair_gpu);
                    if (tid == 0)
                        nvshmem_double_put(sv_imag, sv_imag_pos1, per_pe_work, pair_gpu);
                }
                else // copy local to left
                {
                    if (tid == 0)
                        nvshmem_double_put(&sv_real[per_pe_work], sv_real_pos0, per_pe_work, pair_gpu);
                    if (tid == 0)
                        nvshmem_double_put(&sv_imag[per_pe_work], sv_imag_pos0, per_pe_work, pair_gpu);
                }
            }
        }

        /** This is a less optimized version. Half of the GPUs stay idle for better communication efficiency. Use this version only when bidirectional NVSHMEM communicaiton is not well-supported.
         */
        __device__ __inline__ void C1V2_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit)
        {
            grid_group grid = this_grid();
            const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
            const IdxType per_pe_work = ((dim) >> (gpu_scale));

            if (qubit < lg2_m_gpu)
            {
                C1_GATE(gm_real, gm_imag, qubit);
            }
            else
            {
                IdxType pair_gpu = (i_proc) ^ ((IdxType)1 << (qubit - (lg2_m_gpu)));
                if (i_proc > pair_gpu)
                    return;
                ValType *sv_real_remote = m_real;
                ValType *sv_imag_remote = m_imag;
                if (tid == 0)
                    nvshmem_double_get(sv_real_remote, sv_real, per_pe_work, pair_gpu);
                if (tid == 0)
                    nvshmem_double_get(sv_imag_remote, sv_imag, per_pe_work, pair_gpu);
                grid.sync();

                for (IdxType i = tid; i < per_pe_work; i += blockDim.x * gridDim.x)
                {
                    const ValType el0_real = LOCAL_G_CUDA_MPI(sv_real, i);
                    const ValType el0_imag = LOCAL_G_CUDA_MPI(sv_imag, i);
                    const ValType el1_real = LOCAL_G_CUDA_MPI(sv_real_remote, i);
                    const ValType el1_imag = LOCAL_G_CUDA_MPI(sv_imag_remote, i);

                    ValType real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag);
                    ValType imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real);
                    ValType real_pos1 = (gm_real[2] * el0_real) - (gm_imag[2] * el0_imag) + (gm_real[3] * el1_real) - (gm_imag[3] * el1_imag);
                    ValType imag_pos1 = (gm_real[2] * el0_imag) + (gm_imag[2] * el0_real) + (gm_real[3] * el1_imag) + (gm_imag[3] * el1_real);
                    LOCAL_P_CUDA_MPI(sv_real, i, real_pos0);
                    LOCAL_P_CUDA_MPI(sv_imag, i, imag_pos0);
                    LOCAL_P_CUDA_MPI(sv_real_remote, i, real_pos1);
                    LOCAL_P_CUDA_MPI(sv_imag_remote, i, imag_pos1);
                }
                grid.sync();
                if (tid == 0)
                    nvshmem_double_put(sv_real, sv_real_remote, per_pe_work, pair_gpu);
                if (tid == 0)
                    nvshmem_double_put(sv_imag, sv_imag_remote, per_pe_work, pair_gpu);
            }
        }

        //============== Local 2-qubit Gate  ================
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

                const ValType el0_real = LOCAL_G_CUDA_MPI(sv_real, pos0);
                const ValType el0_imag = LOCAL_G_CUDA_MPI(sv_imag, pos0);
                const ValType el1_real = LOCAL_G_CUDA_MPI(sv_real, pos1);
                const ValType el1_imag = LOCAL_G_CUDA_MPI(sv_imag, pos1);
                const ValType el2_real = LOCAL_G_CUDA_MPI(sv_real, pos2);
                const ValType el2_imag = LOCAL_G_CUDA_MPI(sv_imag, pos2);
                const ValType el3_real = LOCAL_G_CUDA_MPI(sv_real, pos3);
                const ValType el3_imag = LOCAL_G_CUDA_MPI(sv_imag, pos3);

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

                LOCAL_P_CUDA_MPI(sv_real, pos0, sv_real_pos0);
                LOCAL_P_CUDA_MPI(sv_real, pos1, sv_real_pos1);
                LOCAL_P_CUDA_MPI(sv_real, pos2, sv_real_pos2);
                LOCAL_P_CUDA_MPI(sv_real, pos3, sv_real_pos3);

                LOCAL_P_CUDA_MPI(sv_imag, pos0, sv_imag_pos0);
                LOCAL_P_CUDA_MPI(sv_imag, pos1, sv_imag_pos1);
                LOCAL_P_CUDA_MPI(sv_imag, pos2, sv_imag_pos2);
                LOCAL_P_CUDA_MPI(sv_imag, pos3, sv_imag_pos3);
            }
            // BARR_NVSHMEM;
        }

        //============== Unified 2-qubit Gate ================
        // Perform communication optimization here
        __device__ __inline__ void C2V1_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit0, const IdxType qubit1)
        {
            assert(qubit0 != qubit1); // Non-cloning

            if (qubit0 < lg2_m_gpu && qubit1 < lg2_m_gpu)
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

                ValType *sv_real_remote = m_real;
                ValType *sv_imag_remote = m_imag;

                if (tid == 0)
                    nvshmem_double_get(sv_real_remote, sv_real, per_pe_num, pair_gpu);
                if (tid == 0)
                    nvshmem_double_get(sv_imag_remote, sv_imag, per_pe_num, pair_gpu);
                grid.sync();

                IdxType index = (i_proc >> (q - (lg2_m_gpu) + 1)) << q - (lg2_m_gpu);
                index |= i_proc & (((IdxType)1 << (q - (lg2_m_gpu))) - 1);
                for (IdxType i = (index)*per_pe_work + tid; i < (index + 1) * per_pe_work;
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
                    el_real[0] = LOCAL_G_CUDA_MPI(sv_real, term + SV4IDX(0));
                    el_imag[0] = LOCAL_G_CUDA_MPI(sv_imag, term + SV4IDX(0));
                    el_real[3] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV4IDX(3));
                    el_imag[3] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV4IDX(3));

                    if (qubit0 == q) // qubit0 is the remote qubit
                    {
                        el_real[1] = LOCAL_G_CUDA_MPI(sv_real, term + SV4IDX(1));
                        el_imag[1] = LOCAL_G_CUDA_MPI(sv_imag, term + SV4IDX(1));
                        el_real[2] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV4IDX(2));
                        el_imag[2] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV4IDX(2));
                    }
                    else // qubit1 is the remote qubit
                    {
                        el_real[1] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV4IDX(1));
                        el_imag[1] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV4IDX(1));
                        el_real[2] = LOCAL_G_CUDA_MPI(sv_real, term + SV4IDX(2));
                        el_imag[2] = LOCAL_G_CUDA_MPI(sv_imag, term + SV4IDX(2));
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
                    LOCAL_P_CUDA_MPI(sv_real, term + SV4IDX(0), res_real[0]);
                    LOCAL_P_CUDA_MPI(sv_imag, term + SV4IDX(0), res_imag[0]);
                    LOCAL_P_CUDA_MPI(sv_real_remote, term + SV4IDX(3), res_real[3]);
                    LOCAL_P_CUDA_MPI(sv_imag_remote, term + SV4IDX(3), res_imag[3]);

                    if (qubit0 == q) // qubit0 is the remote qubit
                    {
                        LOCAL_P_CUDA_MPI(sv_real, term + SV4IDX(1), res_real[1]);
                        LOCAL_P_CUDA_MPI(sv_imag, term + SV4IDX(1), res_imag[1]);
                        LOCAL_P_CUDA_MPI(sv_real_remote, term + SV4IDX(2), res_real[2]);
                        LOCAL_P_CUDA_MPI(sv_imag_remote, term + SV4IDX(2), res_imag[2]);
                    }
                    else // qubit1 is the remote qubit
                    {
                        LOCAL_P_CUDA_MPI(sv_real_remote, term + SV4IDX(1), res_real[1]);
                        LOCAL_P_CUDA_MPI(sv_imag_remote, term + SV4IDX(1), res_imag[1]);
                        LOCAL_P_CUDA_MPI(sv_real, term + SV4IDX(2), res_real[2]);
                        LOCAL_P_CUDA_MPI(sv_imag, term + SV4IDX(2), res_imag[2]);
                    }
                }
                grid.sync();
                if (tid == 0)
                    nvshmem_double_put(sv_real, sv_real_remote, per_pe_num, pair_gpu);
                if (tid == 0)
                    nvshmem_double_put(sv_imag, sv_imag_remote, per_pe_num, pair_gpu);
                // BARR_NVSHMEM;
            }
        }
        //============== Local 2-qubit Expectation Gate  ================
        __device__ inline void EXPECT_C2_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit0, const IdxType qubit1, const IdxType mask)
        {
            // Functions very similarly to the two-qubit gate. Computes diagonalization and expectation value in-place for a weight-2 Pauli operator
            grid_group grid = this_grid();
            const int tid = blockDim.x * blockIdx.x + threadIdx.x;
            const IdxType per_pe_work = ((dim) >> (gpu_scale + 2));
            assert(qubit0 != qubit1); // Non-cloning
            const IdxType q0dim = ((IdxType)1 << max(qubit0, qubit1));
            const IdxType q1dim = ((IdxType)1 << min(qubit0, qubit1));
            // const IdxType outer_factor = ((dim) + q0dim + q0dim - 1) >> (max(qubit0, qubit1) + 1);
            const IdxType mider_factor = (q0dim + q1dim + q1dim - 1) >> (min(qubit0, qubit1) + 1);
            const IdxType inner_factor = q1dim;
            const IdxType qubit0_dim = ((IdxType)1 << qubit0);
            const IdxType qubit1_dim = ((IdxType)1 << qubit1);
            ValType exp_val = 0.0;
            for (IdxType i = (i_proc)*per_pe_work + tid; i < (i_proc + 1) * per_pe_work; i += blockDim.x * gridDim.x)
            {
                IdxType outer = ((i / inner_factor) / (mider_factor)) * (q0dim + q0dim);
                IdxType mider = ((i / inner_factor) % (mider_factor)) * (q1dim + q1dim);
                IdxType inner = i % inner_factor;
                IdxType pos0 = outer + mider + inner;
                IdxType pos1 = outer + mider + inner + qubit1_dim;
                IdxType pos2 = outer + mider + inner + qubit0_dim;
                IdxType pos3 = outer + mider + inner + q0dim + q1dim;

                const ValType el0_real = LOCAL_G_CUDA_MPI(sv_real, pos0);
                const ValType el0_imag = LOCAL_G_CUDA_MPI(sv_imag, pos0);
                const ValType el1_real = LOCAL_G_CUDA_MPI(sv_real, pos1);
                const ValType el1_imag = LOCAL_G_CUDA_MPI(sv_imag, pos1);
                const ValType el2_real = LOCAL_G_CUDA_MPI(sv_real, pos2);
                const ValType el2_imag = LOCAL_G_CUDA_MPI(sv_imag, pos2);
                const ValType el3_real = LOCAL_G_CUDA_MPI(sv_real, pos3);
                const ValType el3_imag = LOCAL_G_CUDA_MPI(sv_imag, pos3);

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

                ValType v0 = sv_real_pos0 * sv_real_pos0 + sv_imag_pos0 * sv_imag_pos0;
                ValType v1 = sv_real_pos1 * sv_real_pos1 + sv_imag_pos1 * sv_imag_pos1;
                ValType v2 = sv_real_pos2 * sv_real_pos2 + sv_imag_pos2 * sv_imag_pos2;
                ValType v3 = sv_real_pos3 * sv_real_pos3 + sv_imag_pos3 * sv_imag_pos3;

                exp_val += hasEvenParity_cu(pos0 & mask, n_qubits) ? v0 : -v0;
                exp_val += hasEvenParity_cu(pos1 & mask, n_qubits) ? v1 : -v1;
                exp_val += hasEvenParity_cu(pos2 & mask, n_qubits) ? v2 : -v2;
                exp_val += hasEvenParity_cu(pos3 & mask, n_qubits) ? v3 : -v3;
            }
            // BARR_NVSHMEM;
            if (tid < per_pe_work)
            {
                m_real[tid] = exp_val;
            }
        }

        //============== Unified 2-qubit Gate ================
        // Perform communication optimization here
        __device__ inline void EXPECT_C2V1_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit0, const IdxType qubit1, IdxType mask)
        {
            // Functions very similarly to the two-qubit gate with remote optimization.
            assert(qubit0 != qubit1); // Non-cloning
            grid_group grid = this_grid();
            const int tid = blockDim.x * blockIdx.x + threadIdx.x;
            if (qubit0 < lg2_m_gpu && qubit1 < lg2_m_gpu)
            {
                EXPECT_C2_GATE(gm_real, gm_imag, qubit0, qubit1, mask);
            }
            else
            {
                ValType exp_val = 0.0;
                const IdxType per_pe_work = ((dim) >> (gpu_scale + 1));
                const IdxType per_pe_num = ((dim) >> (gpu_scale));
                const IdxType p = min(qubit0, qubit1);
                const IdxType q = max(qubit0, qubit1);

                // load data from pair node
                IdxType pair_gpu = (i_proc) ^ ((IdxType)1 << (q - (lg2_m_gpu)));
                assert(pair_gpu != i_proc);

                if (i_proc > pair_gpu)
                {
                    if (tid < per_pe_work)
                    {
                        m_real[tid] = exp_val;
                    }
                    return;
                }

                ValType *sv_real_remote = m_real;
                ValType *sv_imag_remote = m_imag;
                if (tid == 0)
                    nvshmem_double_get(sv_real_remote, sv_real, per_pe_num, pair_gpu);
                if (tid == 0)
                    nvshmem_double_get(sv_imag_remote, sv_imag, per_pe_num, pair_gpu);
                grid.sync();

                IdxType index = (i_proc >> (q - (lg2_m_gpu) + 1)) << q - (lg2_m_gpu);
                index |= i_proc & (((IdxType)1 << (q - (lg2_m_gpu))) - 1);

                for (IdxType i = (index)*per_pe_work + tid; i < (index + 1) * per_pe_work; i += blockDim.x * gridDim.x)
                {
                    ValType el_real[4];
                    ValType el_imag[4];
                    const IdxType term0 = MOD2E(i, p);
                    const IdxType term1 = MOD2E(DIV2E(i, p), q - p - 1) * EXP2E(p + 1);
                    const IdxType term2 = DIV2E(DIV2E(i, p), q - p - 1) * EXP2E(q + 1);
                    const IdxType term = term2 + term1 + term0;
                    IdxType indices[4] = {
                        term + SV4IDX(0),
                        term + SV4IDX(1),
                        term + SV4IDX(2),
                        term + SV4IDX(3)};
                    el_real[0] = LOCAL_G_CUDA_MPI(sv_real, term + SV4IDX(0));
                    el_imag[0] = LOCAL_G_CUDA_MPI(sv_imag, term + SV4IDX(0));
                    el_real[3] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV4IDX(3));
                    el_imag[3] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV4IDX(3));

                    if (qubit0 == q) // qubit0 is the remote qubit
                    {
                        el_real[1] = LOCAL_G_CUDA_MPI(sv_real, term + SV4IDX(1));
                        el_imag[1] = LOCAL_G_CUDA_MPI(sv_imag, term + SV4IDX(1));
                        el_real[2] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV4IDX(2));
                        el_imag[2] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV4IDX(2));
                    }
                    else // qubit1 is the remote qubit
                    {
                        el_real[1] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV4IDX(1));
                        el_imag[1] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV4IDX(1));
                        el_real[2] = LOCAL_G_CUDA_MPI(sv_real, term + SV4IDX(2));
                        el_imag[2] = LOCAL_G_CUDA_MPI(sv_imag, term + SV4IDX(2));
                    }
                    // #pragma unroll
                    for (unsigned j = 0; j < 4; j++)
                    {
                        ValType r_real = 0.0;
                        ValType r_imag = 0.0;
                        // #pragma unroll
                        for (unsigned k = 0; k < 4; k++)
                        {
                            r_real += (el_real[k] * gm_real[j * 4 + k]) - (el_imag[k] * gm_imag[j * 4 + k]);
                            r_imag += (el_real[k] * gm_imag[j * 4 + k]) + (el_imag[k] * gm_real[j * 4 + k]);
                        }
                        ValType val = r_real * r_real + r_imag * r_imag;
                        exp_val += hasEvenParity_cu(indices[j] & mask, n_qubits) ? val : -val;
                    }
                }
                if (tid < per_pe_work)
                {
                    m_real[tid] = exp_val;
                }
            }
        }

        // Compute the expectation value for a weight-four Pauli string
        __device__ inline void EXPECT_C4_GATE(const ValType *gm_real, const ValType *gm_imag, IdxType qubit0, IdxType qubit1, IdxType qubit2, IdxType qubit3, IdxType mask)
        {
            grid_group grid = this_grid();
            const int tid = blockDim.x * blockIdx.x + threadIdx.x;

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
            ValType exp_val = 0.0;
            const IdxType per_pe_work = ((dim) >> (gpu_scale + 4));

            for (IdxType i = (i_proc)*per_pe_work + tid; i < (i_proc + 1) * per_pe_work; i += blockDim.x * gridDim.x)
            {
                const IdxType term0 = MOD2E(i, p);
                const IdxType term1 = MOD2E(DIV2E(i, p), q - p - 1) * EXP2E(p + 1);
                const IdxType term2 = MOD2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1) * EXP2E(q + 1);
                const IdxType term3 = MOD2E(DIV2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(r + 1);
                const IdxType term4 = DIV2E(DIV2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(s + 1);
                const IdxType term = term4 + term3 + term2 + term1 + term0;
                const ValType el_real[16] = {
                    LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(0)), LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(1)),
                    LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(2)), LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(3)),
                    LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(4)), LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(5)),
                    LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(6)), LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(7)),
                    LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(8)), LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(9)),
                    LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(10)), LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(11)),
                    LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(12)), LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(13)),
                    LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(14)), LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(15))};
                const ValType el_imag[16] = {
                    LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(0)), LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(1)),
                    LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(2)), LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(3)),
                    LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(4)), LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(5)),
                    LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(6)), LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(7)),
                    LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(8)), LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(9)),
                    LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(10)), LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(11)),
                    LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(12)), LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(13)),
                    LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(14)), LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(15))};
                // #pragma unroll
                for (unsigned j = 0; j < 16; j++)
                {
                    ValType res_real = 0;
                    ValType res_imag = 0;
                    // #pragma unroll
                    for (unsigned k = 0; k < 16; k++)
                    {
                        res_real += (el_real[k] * gm_real[j * 16 + k]) - (el_imag[k] * gm_imag[j * 16 + k]);
                        res_imag += (el_real[k] * gm_imag[j * 16 + k]) + (el_imag[k] * gm_real[j * 16 + k]);
                    }

                    ValType val = res_real * res_real + res_imag * res_imag;
                    exp_val += hasEvenParity_cu((term + SV16IDX(j)) & mask, n_qubits) ? val : -val;
                }
            }
            if (tid < per_pe_work)
            {
                m_real[tid] = exp_val;
            }
        }

        // Compute the expectation value of a weight-four Pauli operator (this time with communication optimizations, one qubit is remote)
        __device__ void EXPECT_C4V1_GATE(const ValType *gm_real, const ValType *gm_imag, IdxType qubit0, IdxType qubit1, IdxType qubit2, IdxType qubit3, IdxType mask)
        {
            grid_group grid = this_grid();
            const int tid = blockDim.x * blockIdx.x + threadIdx.x;
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
            const IdxType per_pe_work = ((dim) >> (gpu_scale + 3));
            const IdxType per_pe_num = ((dim) >> (gpu_scale));
            if (s < lg2_m_gpu)
            {
                EXPECT_C4_GATE(gm_real, gm_imag, qubit0, qubit1, qubit2, qubit3, mask);
            }
            else
            {
                ValType exp_val = 0.0;

                // load data from pair node
                IdxType pair_gpu = (i_proc) ^ ((IdxType)1 << (s - (lg2_m_gpu)));
                assert(pair_gpu != i_proc);
                if (i_proc > pair_gpu)
                {

                    if (tid < per_pe_work)
                    {
                        m_real[tid] = exp_val;
                    }
                    return;
                }
                // make the indices continguous (e.g. map from i_proc = 0, 2, 4, 6 - > 0, 1, 2, 3)
                // Basically, we're just ``deleting'' bit s - (lg2_m_gpu). If s - (lg2_m_gpu) is 1 and i_proc = 5 = 0b101, then we get 0b011 (3)
                IdxType index = (i_proc >> (s - (lg2_m_gpu) + 1)) << (s - (lg2_m_gpu));
                index |= i_proc & (((IdxType)1 << (s - (lg2_m_gpu))) - 1);
                ValType *sv_real_remote = m_real;
                ValType *sv_imag_remote = m_imag;
                if (tid == 0)
                    nvshmem_double_get(sv_real_remote, sv_real, per_pe_num, pair_gpu);
                if (tid == 0)
                    nvshmem_double_get(sv_imag_remote, sv_imag, per_pe_num, pair_gpu);
                grid.sync();
                for (IdxType i = (index)*per_pe_work + tid; i < (index + 1) * per_pe_work; i += blockDim.x * gridDim.x)
                {
                    ValType el_real[16];
                    ValType el_imag[16];
                    const IdxType term0 = MOD2E(i, p);
                    const IdxType term1 = MOD2E(DIV2E(i, p), q - p - 1) * EXP2E(p + 1);
                    const IdxType term2 = MOD2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1) * EXP2E(q + 1);
                    const IdxType term3 = MOD2E(DIV2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(r + 1);
                    const IdxType term4 = DIV2E(DIV2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(s + 1);
                    const IdxType term = term4 + term3 + term2 + term1 + term0;
                    if (qubit3 == s) // qubit3 is remote qubit
                    {
                        // printf("Qubit 3");
                        el_real[0] = LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(0));
                        el_real[1] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV16IDX(1));
                        el_real[2] = LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(2));
                        el_real[3] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV16IDX(3));
                        el_real[4] = LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(4));
                        el_real[5] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV16IDX(5));
                        el_real[6] = LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(6));
                        el_real[7] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV16IDX(7));
                        el_real[8] = LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(8));
                        el_real[9] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV16IDX(9));
                        el_real[10] = LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(10));
                        el_real[11] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV16IDX(11));
                        el_real[12] = LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(12));
                        el_real[13] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV16IDX(13));
                        el_real[14] = LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(14));
                        el_real[15] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV16IDX(15));

                        el_imag[0] = LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(0));
                        el_imag[1] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV16IDX(1));
                        el_imag[2] = LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(2));
                        el_imag[3] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV16IDX(3));
                        el_imag[4] = LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(4));
                        el_imag[5] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV16IDX(5));
                        el_imag[6] = LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(6));
                        el_imag[7] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV16IDX(7));
                        el_imag[8] = LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(8));
                        el_imag[9] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV16IDX(9));
                        el_imag[10] = LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(10));
                        el_imag[11] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV16IDX(11));
                        el_imag[12] = LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(12));
                        el_imag[13] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV16IDX(13));
                        el_imag[14] = LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(14));
                        el_imag[15] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV16IDX(15));
                    }
                    else // qubit2 is remote (not possible qubit0 or 1 is remote)
                    {
                        // if (( (laneid>>1)&1)!=0) el_real_s[j*16+laneid] = LOCAL_G_CUDA_MPI(sv_real_remote,addr);
                        // else el_real_s[j*16+laneid] = LOCAL_G_CUDA_MPI(sv_real, addr);
                        // laneid = 2,3,6,7,10,11,14,15;

                        el_real[0] = LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(0));
                        el_real[1] = LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(1));
                        el_real[2] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV16IDX(2));
                        el_real[3] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV16IDX(3));
                        el_real[4] = LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(4));
                        el_real[5] = LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(5));
                        el_real[6] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV16IDX(6));
                        el_real[7] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV16IDX(7));
                        el_real[8] = LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(8));
                        el_real[9] = LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(9));
                        el_real[10] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV16IDX(10));
                        el_real[11] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV16IDX(11));
                        el_real[12] = LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(12));
                        el_real[13] = LOCAL_G_CUDA_MPI(sv_real, term + SV16IDX(13));
                        el_real[14] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV16IDX(14));
                        el_real[15] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV16IDX(15));

                        el_imag[0] = LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(0));
                        el_imag[1] = LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(1));
                        el_imag[2] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV16IDX(2));
                        el_imag[3] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV16IDX(3));
                        el_imag[4] = LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(4));
                        el_imag[5] = LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(5));
                        el_imag[6] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV16IDX(6));
                        el_imag[7] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV16IDX(7));
                        el_imag[8] = LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(8));
                        el_imag[9] = LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(9));
                        el_imag[10] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV16IDX(10));
                        el_imag[11] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV16IDX(11));
                        el_imag[12] = LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(12));
                        el_imag[13] = LOCAL_G_CUDA_MPI(sv_imag, term + SV16IDX(13));
                        el_imag[14] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV16IDX(14));
                        el_imag[15] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV16IDX(15));
                    }

                    // #pragma unroll
                    for (unsigned j = 0; j < 16; j++)
                    {
                        ValType res_real = 0.0;
                        ValType res_imag = 0.0;
                        // #pragma unroll
                        for (unsigned k = 0; k < 16; k++)
                        {
                            res_real += (el_real[k] * gm_real[j * 16 + k]) - (el_imag[k] * gm_imag[j * 16 + k]);
                            res_imag += (el_real[k] * gm_imag[j * 16 + k]) + (el_imag[k] * gm_real[j * 16 + k]);
                        }
                        ValType val = res_real * res_real + res_imag * res_imag;

                        exp_val += hasEvenParity_cu((term + SV16IDX(j)) & mask, n_qubits) ? val : -val;
                        // if (abs(val) > 1e-3) {
                        //     printf("%lld %lld %lld %lld %lld %lld %lld %f %lld %f\n", tid, i_proc, index, j, term, i, term + SV16IDX(j), val, hasEvenParity_cu((term + SV16IDX(j)) & mask, n_qubits), exp_val);
                        // }
                    }
                }
                if (tid < per_pe_work)
                {
                    m_real[tid] = exp_val;
                }
            }
        }

        // Compute the expectation value of a diagonal-basis Pauli operator. No communication needed, just need the local sum
        __device__ inline void EXPECT_C0_GATE(IdxType mask, ValType coeff)
        {
            grid_group grid = this_grid();
            const int tid = blockDim.x * blockIdx.x + threadIdx.x;
            double exp_val = 0.0;
            const IdxType per_pe_work = m_gpu;
            // loop over local quantum state
            for (IdxType i = (i_proc)*per_pe_work + tid; i < (i_proc + 1) * per_pe_work; i += blockDim.x * gridDim.x)
            {
                // compute the probability amplitude
                double val = LOCAL_G_CUDA_MPI(sv_real, i) * LOCAL_G_CUDA_MPI(sv_real, i) + LOCAL_G_CUDA_MPI(sv_imag, i) * LOCAL_G_CUDA_MPI(sv_imag, i);
                // add/subtract depending on the Pauli parity relative to the state index
                exp_val += hasEvenParity_cu(i & mask, n_qubits) ? val : -val;
            }
            if (tid < per_pe_work)
            {
                // store accumulated result in m_real
                m_real[tid] += coeff * exp_val;
            }
        }

        // Perform a parallel reduction over the locally cached expectation sums
        __device__ inline void EXPECT_REDUCE(ValType *output)
        {
            grid_group grid = this_grid();
            const int tid = blockDim.x * blockIdx.x + threadIdx.x;
            // ensure the reduction dimension is a power of 2
            IdxType gridlog2 = 63 - __clzll(blockDim.x * gridDim.x);
            if (blockDim.x * gridDim.x & (((IdxType)1 << gridlog2) - 1))
            {
                gridlog2 += 1;
            }
            IdxType reduce_limit = (IdxType)1 << gridlog2;
            reduce_limit = min(reduce_limit, dim >> gpu_scale);
            // Parallel reduction
            for (IdxType k = (reduce_limit >> 1); k > 0; k >>= 1)
            {
                if (tid < k && tid + k < blockDim.x * gridDim.x)
                {
                    m_real[tid] += m_real[tid + k];
                }
                grid.sync();
            }
            grid.sync();
            if (tid == 0)
            {
                // write the output to the assigned pointer
                *output = m_real[0];
            }
        }
        // Compute the expectation value for a QWC group
        __device__ __inline__ void EXPECT_GATE(ObservableList *o)
        {
            grid_group grid = this_grid();
            const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
            if (tid < m_gpu)
            {
                m_real[tid] = 0;
            }
            size_t nterms = o->numterms;
            IdxType *zmasks = o->zmasks;
            ValType *coeffs = o->coeffs;
            for (IdxType obs_ind = 0; obs_ind < nterms; obs_ind++)
            {
                EXPECT_C0_GATE(zmasks[obs_ind], coeffs[obs_ind]);
            }
            BARR_NVSHMEM;
            EXPECT_REDUCE(&o->exp_output);
        }

        __device__ __inline__ void M_GATE(const IdxType qubit, const IdxType cur_index)
        {
            ValType rand = randoms_gpu[cur_index];

            grid_group grid = this_grid();
            const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
            const IdxType per_pe_work = ((dim) >> (gpu_scale));
            // ValType *m_real = m_real;
            IdxType mask = ((IdxType)1 << qubit);

            for (IdxType i = tid; i < m_gpu; i += blockDim.x * gridDim.x)
            {
                IdxType idx = (i_proc)*per_pe_work + i;
                if ((idx & mask) == 0)
                    m_real[i] = 0;
                else
                    m_real[i] = sv_real[i] * sv_real[i] + sv_imag[i] * sv_imag[i];
            }
            BARR_NVSHMEM;

            // Parallel reduction
            for (IdxType k = ((IdxType)1 << (n_qubits - 1)); k > 0; k >>= 1)
            {
                for (IdxType i = tid; i < k; i += blockDim.x * gridDim.x)
                {
                    IdxType local_gpu = i >> (lg2_m_gpu);
                    if (local_gpu == i_proc)
                    {
                        IdxType local_idx = i & (m_gpu - 1);
                        m_real[local_idx] += PGAS_G(m_real, i + k);
                    }
                }
                BARR_NVSHMEM;
            }

            if (tid == 0 && i_proc != 0)
                m_real[0] = PGAS_G(m_real, 0);
            grid.sync();
            ValType prob_of_one = m_real[0];
            grid.sync();

            if (rand < prob_of_one)
            {
                ValType factor = 1. / sqrt(prob_of_one);
                for (IdxType i = tid; i < m_gpu; i += blockDim.x * gridDim.x)
                {
                    IdxType idx = (i_proc)*per_pe_work + i;
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
                for (IdxType i = tid; i < m_gpu; i += blockDim.x * gridDim.x)
                {
                    IdxType idx = (i_proc)*per_pe_work + i;
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
            if (tid == 0)
                results_gpu[cur_index] = (rand <= prob_of_one ? 1 : 0);
            BARR_NVSHMEM;
        }

        __device__ __inline__ void MA_GATE(const IdxType repetition, const IdxType cur_index)
        {
            grid_group grid = this_grid();
            const IdxType n_size = (IdxType)1 << (n_qubits);
            const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
            // ValType *m_real = m_real;
            for (IdxType i = tid; i < m_gpu; i += blockDim.x * gridDim.x)
            {
                m_real[i] = sv_real[i] * sv_real[i] + sv_imag[i] * sv_imag[i];
            }
            BARR_NVSHMEM;

            // Parallel prefix sum
            for (IdxType d = 0; d < (n_qubits); d++)
            {
                IdxType step = (IdxType)1 << (d + 1);
                for (IdxType k = tid * step; k < n_size; k += step * blockDim.x * gridDim.x)
                {
                    IdxType idx = k + ((IdxType)1 << (d + 1)) - 1;
                    IdxType local_gpu = idx >> (lg2_m_gpu);
                    if (local_gpu == i_proc)
                    {
                        IdxType local_idx = idx & ((m_gpu)-1);
                        m_real[local_idx] += PGAS_G(m_real, k + ((IdxType)1 << d) - 1);
                    }
                }
                BARR_NVSHMEM;
            }

            if (i_proc == (n_gpus - 1) && tid == 0) // last GPU
            {
                ValType val = LOCAL_G_CUDA_MPI(m_real, n_size - 1);
                LOCAL_P_CUDA_MPI(m_real, n_size - 1, 0);
                ValType purity = fabs(val);
                if (fabs(purity - 1.0) > ERROR_BAR)
                    printf("MA: Purity Check fails with %lf\n", purity);
            }
            BARR_NVSHMEM;

            for (IdxType d = (n_qubits)-1; d >= 0; d--)
            {
                IdxType step = (IdxType)1 << (d + 1);
                for (IdxType k = tid * step; k < n_size - 1; k += step * blockDim.x * gridDim.x)
                {
                    IdxType idx = k + ((IdxType)1 << (d + 1)) - 1;
                    IdxType local_gpu = idx >> (lg2_m_gpu);
                    if (local_gpu == i_proc)
                    {
                        IdxType local_idx = idx & ((m_gpu)-1);
                        IdxType remote_idx = k + ((IdxType)1 << d) - 1;
                        ValType tmp = PGAS_G(m_real, remote_idx);
                        ValType tmp2 = m_real[local_idx];
                        PGAS_P(m_real, remote_idx, tmp2);
                        m_real[local_idx] = tmp + tmp2;
                    }
                }
                BARR_NVSHMEM;
            }

            for (IdxType j = tid; j < n_size; j += blockDim.x * gridDim.x)
            {
                IdxType local_gpu = j >> (lg2_m_gpu);
                if (local_gpu == i_proc)
                {
                    ValType lower = LOCAL_G_CUDA_MPI(m_real, j);
                    ValType upper = (j + 1 == n_size) ? 1 : PGAS_G(m_real, j + 1);
                    for (IdxType i = 0; i < repetition; i++)
                    {
                        ValType r = randoms_gpu[cur_index + i];
                        if (lower <= r && r < upper)
                            nvshmem_longlong_p(&results_gpu[cur_index + i], j, 0);
                    }
                }
            }

            BARR_NVSHMEM;

            if (i_proc != 0 && tid == 0)
                nvshmem_longlong_get(results_gpu, results_gpu, repetition, 0);
            BARR_NVSHMEM;
        }

        __device__ __inline__ void MAV1_GATE(const IdxType repetition, const IdxType cur_index)
        {
            grid_group grid = this_grid();
            const IdxType n_size = (IdxType)1 << (n_qubits);
            const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
            // ValType *m_real = m_real;
            // ValType *m_imag = m_imag;
            for (IdxType i = tid; i < m_gpu; i += blockDim.x * gridDim.x)
            {
                m_real[i] = sv_real[i] * sv_real[i] + sv_imag[i] * sv_imag[i];
            }
            grid.sync();

            // local parallel prefix sum
            for (IdxType d = 0; d < (lg2_m_gpu); d++)
            {
                IdxType step = (IdxType)1 << (d + 1);
                for (IdxType k = tid * step; k < m_gpu; k += step * blockDim.x * gridDim.x)
                {
                    m_real[(k + ((IdxType)1 << (d + 1)) - 1)] += m_real[k + ((IdxType)1 << d) - 1];
                }
                grid.sync();
            }
            if (tid == 0)
            {
                m_imag[0] = m_real[m_gpu - 1]; // local sum
                m_real[m_gpu - 1] = 0;
            }
            grid.sync();
            for (IdxType d = (lg2_m_gpu)-1; d >= 0; d--)
            {
                IdxType step = (IdxType)1 << (d + 1);
                for (IdxType k = tid * step; k < m_gpu - 1; k += step * blockDim.x * gridDim.x)
                {
                    ValType tmp = m_real[k + ((IdxType)1 << d) - 1];
                    ValType tmp2 = m_real[(k + ((IdxType)1 << (d + 1)) - 1)];
                    m_real[k + ((IdxType)1 << d) - 1] = tmp2;
                    m_real[(k + ((IdxType)1 << (d + 1)) - 1)] = tmp + tmp2;
                }
                grid.sync();
            }
            BARR_NVSHMEM;

            if (i_proc == 0 && tid == 0) // first GPU
            {
                ValType partial = 0;
                for (IdxType g = 0; g < n_gpus; g++)
                {
                    nvshmem_double_p(&m_imag[1], partial, g);
                    ValType inc = nvshmem_double_g(&m_imag[0], g);
                    partial += inc;
                }
                ValType purity = fabs(partial);
                if (fabs(purity - 1.0) > ERROR_BAR)
                    printf("MA: Purity Check fails with %lf\n", purity);
            }

            BARR_NVSHMEM;
            for (IdxType i = tid; i < m_gpu; i += blockDim.x * gridDim.x)
            {
                m_real[i] += m_imag[1];
            }

            BARR_NVSHMEM;
            for (IdxType j = tid; j < n_size; j += blockDim.x * gridDim.x)
            {
                IdxType local_gpu = j >> (lg2_m_gpu);
                if (local_gpu == i_proc)
                {
                    ValType lower = LOCAL_G_CUDA_MPI(m_real, j);
                    ValType upper = (j + 1 == n_size) ? 1 : PGAS_G(m_real, j + 1);
                    for (IdxType i = 0; i < repetition; i++)
                    {
                        ValType r = randoms_gpu[cur_index + i];
                        if (lower <= r && r < upper)
                            nvshmem_longlong_p(&results_gpu[cur_index + i], j, 0);
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
            const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
            const IdxType per_pe_work = ((dim) >> (gpu_scale));
            // ValType *m_real = m_real;
            IdxType mask = ((IdxType)1 << qubit);

            for (IdxType i = tid; i < m_gpu; i += blockDim.x * gridDim.x)
            {
                IdxType idx = (i_proc)*per_pe_work + i;
                if ((idx & mask) == 0)
                    m_real[i] = 0;
                else
                    m_real[i] = sv_real[i] * sv_real[i] + sv_imag[i] * sv_imag[i];
            }
            BARR_NVSHMEM;

            // Parallel reduction
            for (IdxType k = ((IdxType)1 << (n_qubits - 1)); k > 0; k >>= 1)
            {
                for (IdxType i = tid; i < k; i += blockDim.x * gridDim.x)
                {
                    IdxType local_gpu = i >> (lg2_m_gpu);
                    if (local_gpu == i_proc)
                    {
                        IdxType local_idx = i & (m_gpu - 1);
                        m_real[local_idx] += PGAS_G(m_real, i + k);
                    }
                }
                BARR_NVSHMEM;
            }

            if (tid == 0 && i_proc != 0)
                m_real[0] = PGAS_G(m_real, 0);
            grid.sync();
            ValType prob_of_one = m_real[0];
            grid.sync();

            if (prob_of_one < 1.0) // still possible to normalize
            {
                ValType factor = 1.0 / sqrt(1.0 - prob_of_one);
                for (IdxType i = tid; i < m_gpu; i += blockDim.x * gridDim.x)
                {
                    IdxType idx = (i_proc)*per_pe_work + i;
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
            else
            {
                if (qubit >= lg2_m_gpu) // remote qubit, need switch
                {
                    IdxType pair_gpu = (i_proc) ^ ((IdxType)1 << (qubit - (lg2_m_gpu)));
                    assert(pair_gpu != i_proc);
                    ValType *sv_real_remote = m_real;
                    ValType *sv_imag_remote = m_imag;
                    if (tid == 0)
                        nvshmem_double_get(sv_real_remote, sv_real, per_pe_work, pair_gpu);
                    if (tid == 0)
                        nvshmem_double_get(sv_imag_remote, sv_imag, per_pe_work, pair_gpu);
                    BARR_NVSHMEM;
                    for (IdxType i = tid; i < m_gpu; i += blockDim.x * gridDim.x)
                    {
                        sv_real[i] = sv_real_remote[i];
                        sv_imag[i] = sv_imag_remote[i];
                    }
                }
                else
                {
                    for (IdxType i = tid; i < m_gpu; i += blockDim.x * gridDim.x)
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
            }
            BARR_NVSHMEM;
        }

        //============== SWAP Gate ================
        // This gate is for internal usage. It is used
        // when ctrl and target qubits are remote qubis, we then
        // swap one of them to a local qubit,
        // perform the C2 gate, and then swap back
        // It is assumed qubit0 is local, qubit1 is remote
        __device__ __inline__ void SWAP_GATE(const IdxType qubit0, const IdxType qubit1)
        {
            assert(qubit0 != qubit1); // Non-cloning
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

            ValType *sv_real_remote = m_real;
            ValType *sv_imag_remote = m_imag;

            if (tid == 0)
                nvshmem_double_get(sv_real_remote, sv_real, per_pe_num, pair_gpu);
            if (tid == 0)
                nvshmem_double_get(sv_imag_remote, sv_imag, per_pe_num, pair_gpu);
            grid.sync();
            IdxType index = (i_proc >> (q - (lg2_m_gpu) + 1)) << (q - (lg2_m_gpu));
            index |= i_proc & (((IdxType)1 << (q - (lg2_m_gpu))) - 1);
            for (IdxType i = (index)*per_pe_work + tid; i < (index + 1) * per_pe_work;
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

                el_real[1] = LOCAL_G_CUDA_MPI(sv_real_remote, term + SV4IDX(1));
                el_imag[1] = LOCAL_G_CUDA_MPI(sv_imag_remote, term + SV4IDX(1));
                el_real[2] = LOCAL_G_CUDA_MPI(sv_real, term + SV4IDX(2));
                el_imag[2] = LOCAL_G_CUDA_MPI(sv_imag, term + SV4IDX(2));

                res_real[1] = el_real[2];
                res_imag[1] = el_imag[2];
                res_real[2] = el_real[1];
                res_imag[2] = el_imag[1];

                LOCAL_P_CUDA_MPI(sv_real_remote, term + SV4IDX(1), res_real[1]);
                LOCAL_P_CUDA_MPI(sv_imag_remote, term + SV4IDX(1), res_imag[1]);
                LOCAL_P_CUDA_MPI(sv_real, term + SV4IDX(2), res_real[2]);
                LOCAL_P_CUDA_MPI(sv_imag, term + SV4IDX(2), res_imag[2]);
            }
            grid.sync();
            if (tid == 0)
                nvshmem_double_put(sv_real, sv_real_remote, per_pe_num, pair_gpu);
            if (tid == 0)
                nvshmem_double_put(sv_imag, sv_imag_remote, per_pe_num, pair_gpu);
            // BARR_NVSHMEM;
        }

        __device__ __inline__ void Purity_Check(const IdxType t)
        {
            grid_group grid = this_grid();
            const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
            const IdxType per_pe_work = ((dim) >> (gpu_scale));
            // ValType *m_real = m_real;

            for (IdxType i = (i_proc)*per_pe_work + tid; i < (i_proc + 1) * per_pe_work; i += blockDim.x * gridDim.x)
            {
                ValType val_real = PGAS_G(sv_real, i);
                ValType val_imag = PGAS_G(sv_imag, i);
                PGAS_P(m_real, i, (val_real * val_real) + (val_imag * val_imag));
            }
            BARR_NVSHMEM;

            if (i_proc == 0)
            {
                for (IdxType k = ((IdxType)1 << (n_qubits - 1)); k > 0; k >>= 1)
                {
                    for (IdxType i = tid; i < k; i += blockDim.x * gridDim.x)
                    {
                        ValType a = PGAS_G(m_real, i);
                        ValType b = PGAS_G(m_real, i + k);
                        PGAS_P(m_real, i, a + b);
                    }
                    grid.sync();
                }
                if (threadIdx.x == 0 && blockIdx.x == 0)
                {
                    ValType purity = m_real[0];
                    if (abs(purity - 1.0) > ERROR_BAR)
                    {
                        printf("Purity Check fails after Gate-%lld with %lf\n", t, purity);
                        // Gate *g = &circuit_handle_gpu[t];
                        // printf("Purity Check fails after Gate-%lld=>%s(ctrl:%lld,qubit:%lld,theta:%lf) with %lf\n", t, OP_NAMES_CUDA[g->op_name], g->ctrl, g->qubit, g->theta, purity);
                    }
                }
            }
            BARR_NVSHMEM;
        }
        virtual ValType *get_real() const override { return sv_real; };
        virtual ValType *get_imag() const override { return sv_imag; };
    };

    __global__ void simulation_kernel_cuda_mpi(SV_CUDA_MPI *sv_gpu, IdxType n_gates)
    {
        IdxType cur_index = 0;
        IdxType lg2_m_gpu = sv_gpu->lg2_m_gpu;
        grid_group grid = this_grid();
        bool already_sync = false;
        for (IdxType t = 0; t < n_gates; t++)
        {
            OP op_name = (sv_gpu->gates_gpu)[t].op_name;
            IdxType qubit = (sv_gpu->gates_gpu)[t].qubit;

            IdxType ctrl = (sv_gpu->gates_gpu)[t].ctrl;
            ValType *gm_real = (sv_gpu->gates_gpu)[t].gm_real;
            ValType *gm_imag = (sv_gpu->gates_gpu)[t].gm_imag;

            IdxType repetition = (sv_gpu->gates_gpu)[t].qubit;

            // only need sync when operating on remote qubits
            if ((ctrl >= lg2_m_gpu) || (qubit >= lg2_m_gpu))
            {
                if (!already_sync) // do not need repeated sync
                {
                    if (threadIdx.x == 0 && blockIdx.x == 0)
                        nvshmem_barrier_all();
                    grid.sync();
                }
            }
            already_sync = false;

            // only need sync when operating on remote qubits
            if (op_name == OP::C1)
            {
                sv_gpu->C1V2_GATE(gm_real, gm_imag, qubit);
                // sv_gpu->C1V1_GATE(gm_real, gm_imag, qubit);
            }
            else if (op_name == OP::C2)
            {
                if ((ctrl >= lg2_m_gpu) && (qubit >= lg2_m_gpu))
                {
                    sv_gpu->SWAP_GATE(0, ctrl);
                    BARR_NVSHMEM;
                    sv_gpu->C2V1_GATE(gm_real, gm_imag, 0, qubit);
                    BARR_NVSHMEM;
                    sv_gpu->SWAP_GATE(0, ctrl);
                }
                else
                {
                    sv_gpu->C2V1_GATE(gm_real, gm_imag, ctrl, qubit);
                }
            }
            else if (op_name == OP::RESET)
            {
                sv_gpu->RESET_GATE(qubit);
            }
            else if (op_name == OP::M)
            {
                sv_gpu->M_GATE(qubit, cur_index);
                cur_index++;
            }
            else if (op_name == OP::MA)
            {
                sv_gpu->MA_GATE(repetition, cur_index);
                cur_index += repetition;
            }
            else if (op_name == OP::EXPECT)
            {

                ObservableList *o = (ObservableList *)((sv_gpu->gates_gpu)[t].data);
                BARR_NVSHMEM;
                sv_gpu->EXPECT_GATE(o);

                BARR_NVSHMEM;
            }
            // only need sync when operating on remote qubits
            if ((ctrl >= lg2_m_gpu) || (qubit >= lg2_m_gpu))
            {
                if (threadIdx.x == 0 && blockIdx.x == 0)
                    nvshmem_barrier_all();
                already_sync = true;
            }
            grid.sync();
        }
    }
}
