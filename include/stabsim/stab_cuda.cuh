#pragma once

#include "state.hpp"
#include "nwq_util.hpp"
#include "gate.hpp"
#include "circuit.hpp"

#include "config.hpp"
#include "private/cuda_util.cuh"
#include "private/macros.hpp"
#include "private/sim_gate.hpp"
#include "utils/prng.hpp"

#include "circuit_pass/fusion.hpp"
#include "private/exp_gate_declarations.cuh"

#include <random>
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <set>
#include <assert.h>
#include <cuda_runtime_api.h>
#include <cooperative_groups.h>
#include <curand_kernel.h>
#include <iostream>
#include <cuda.h>
#include <memory>

#ifdef FP64_TC_AVAILABLE
#include <mma.h>
#endif

namespace NWQSim
{
    using namespace cooperative_groups;
    namespace cg = cooperative_groups;
    using namespace std;

    // Simulation kernel runtime
    class STAB_CUDA;
    __global__ void simulation_kernel_cuda(STAB_CUDA *stab_gpu, Gate *gates_gpu, IdxType n_gates, int32_t seed);
    __global__ void simulation_kernel_cuda2D(STAB_CUDA *stab_gpu, Gate *gates_gpu, IdxType gate_chunk);

    class STAB_CUDA : public QuantumState
    {
    public:
        // Default identity constructor
        STAB_CUDA(IdxType _n_qubits) : QuantumState(SimType::STAB)
        {
            std::cout << "CUDA constructor" << std::endl;
            /*Initialize basic data*/
            n = _n_qubits;
            rows = 2 * n + 1;
            cols = n;
            stabCounts = n;
            packed_bits = 32;
            packed_rows = (rows + packed_bits - 1) / packed_bits;
            packed_r_size = packed_rows * sizeof(int32_t);
            packed_matrix_size = (packed_rows * cols) * sizeof(int32_t);
            bit_r_size = rows * sizeof(int32_t);
            bit_matrix_size = rows * cols * sizeof(int32_t);

            i_proc = 0;
            cudaSafeCall(cudaSetDevice(i_proc));

            // Space out the 2D vectors for x z and r
            reset_state();

            seed = 0;

            /*End initialization*/

            // Allocate the packed data to the GPU side using NVSHMEM and a tempRow for row swapping
            //  SAFE_ALOC_GPU(x_packed_gpu, packed_matrix_size);
            //  SAFE_ALOC_GPU(z_packed_gpu, packed_matrix_size);
            //  SAFE_ALOC_GPU(r_packed_gpu, packed_r_size);

            SAFE_ALOC_GPU(x_bit_gpu, bit_matrix_size);
            SAFE_ALOC_GPU(z_bit_gpu, bit_matrix_size);
            SAFE_ALOC_GPU(r_bit_gpu, bit_r_size);

            cudaCheckError();
        }

        ~STAB_CUDA()
        {
            SAFE_FREE_HOST_CUDA(x_packed_cpu);
            SAFE_FREE_HOST_CUDA(z_packed_cpu);
            SAFE_FREE_HOST_CUDA(r_packed_cpu);

            // Allocate the packed data to the GPU side using NVSHMEM
            SAFE_FREE_GPU(x_packed_gpu);
            SAFE_FREE_GPU(z_packed_gpu);
            SAFE_FREE_GPU(r_packed_gpu);
            SAFE_FREE_GPU(x_bit_gpu);
            SAFE_FREE_GPU(z_bit_gpu);
            SAFE_FREE_GPU(r_bit_gpu);
            SAFE_FREE_GPU(gates_gpu);
            SAFE_FREE_GPU(d_local_sums);
        }

        // Packs down 32 rows in each column and flattens
        void pack_tableau()
        {
            // print_res_state();
            // Allocate memory for x_packed_cpu (2D array), but as a flattened 1D array to ensure memory is contiguous
            x_packed_cpu = new int32_t[packed_rows * cols]();
            z_packed_cpu = new int32_t[packed_rows * cols]();
            r_packed_cpu = new int32_t[packed_rows]();
            int32_t index;

            for (int32_t i = 0; i < rows; i++)
            {
                mask = i % packed_bits;
                r_packed_cpu[i / packed_bits] |= r[i] << (mask);

                for (int32_t j = 0; j < cols; j++)
                {
                    // Index conversion
                    index = ((i / packed_bits) * cols) + j;

                    x_packed_cpu[index] |= (x[i][j] << (mask));
                    z_packed_cpu[index] |= (z[i][j] << (mask));
                }
            }

            std::cout << "Pack tableau done!" << std::endl;

            // print_packed_state();
        }

        // Unpacks packed CPU arrays back to bit values
        void unpack_tableau()
        {
            std::cout << "Unpacking tableau!" << std::endl;
            int32_t index;

            for (int32_t i = 0; i < rows; i++)
            {
                mask = 1 << (i % packed_bits);
                r[i] = (r_packed_cpu[i / packed_bits] & mask) ? 1 : 0;
                for (int32_t j = 0; j < cols; j++)
                {
                    index = (i / packed_bits * cols) + j;

                    x[i][j] = (x_packed_cpu[index] & mask) ? 1 : 0;
                    z[i][j] = (z_packed_cpu[index] & mask) ? 1 : 0;
                }
            }
        }

        void copy_bits_to_gpu()
        {
            IdxType index;

            int32_t *temp_r_bit = new int32_t[rows];
            int32_t *temp_x_bit = new int32_t[rows * cols];
            int32_t *temp_z_bit = new int32_t[rows * cols];

            for (int32_t i = 0; i < rows; i++)
            {
                temp_r_bit[i] = r[i];

                for (int32_t j = 0; j < cols; j++)
                {
                    index = (i * cols) + j;
                    temp_x_bit[index] = x[i][j];
                    temp_z_bit[index] = z[i][j];
                }
            }

            cudaMemcpy(r_bit_gpu, temp_r_bit, rows * sizeof(int32_t), cudaMemcpyHostToDevice);
            cudaMemcpy(x_bit_gpu, temp_x_bit, rows * cols * sizeof(int32_t), cudaMemcpyHostToDevice);
            cudaMemcpy(z_bit_gpu, temp_z_bit, rows * cols * sizeof(int32_t), cudaMemcpyHostToDevice);

            delete[] temp_r_bit;
            delete[] temp_x_bit;
            delete[] temp_z_bit;
        }

        void copy_bits_from_gpu()
        {
            IdxType index;

            int32_t *temp_r_bit = new int32_t[rows];
            int32_t *temp_x_bit = new int32_t[rows * cols];
            int32_t *temp_z_bit = new int32_t[rows * cols];

            cudaMemcpy(temp_r_bit, r_bit_gpu, rows * sizeof(int32_t), cudaMemcpyDeviceToHost);
            cudaMemcpy(temp_x_bit, x_bit_gpu, rows * cols * sizeof(int32_t), cudaMemcpyDeviceToHost);
            cudaMemcpy(temp_z_bit, z_bit_gpu, rows * cols * sizeof(int32_t), cudaMemcpyDeviceToHost);

            for (int32_t i = 0; i < rows; i++)
            {
                r[i] = temp_r_bit[i];

                for (int32_t j = 0; j < cols; j++)
                {
                    index = (i * cols) + j;
                    x[i][j] = temp_x_bit[index];
                    z[i][j] = temp_z_bit[index];
                }
            }

            delete[] temp_r_bit;
            delete[] temp_x_bit;
            delete[] temp_z_bit;
        }
        // Prints the tableau including phase
        void print_packed_state()
        {
            int32_t index;
            std::cout << "---------- Packed Tableau: ----------" << std::endl;
            for (int32_t i = 0; i < rows; i++)
            {
                mask = 1 << (i % packed_bits);

                if (((i == (rows / 2)) || (i == rows - 1)))
                {
                    for (int32_t j = -5; j < (n * 2); j++)
                    {
                        std::cout << "-";
                    }
                    std::cout << std::endl;
                }
                for (int32_t j = 0; j < cols; j++)
                {
                    // Index conversion
                    index = (i / packed_bits * cols) + j;
                    std::cout << ((x_packed_cpu[index] & mask) ? 1 : 0);
                }
                std::cout << " | ";
                for (int32_t j = 0; j < cols; j++)
                {
                    // Index conversion
                    index = (i / packed_bits * cols) + j;
                    std::cout << ((z_packed_cpu[index] & mask) ? 1 : 0);
                }
                std::cout << "|" << ((r_packed_cpu[i / packed_bits] & mask) ? 1 : 0) << std::endl;
            }
            std::cout << "---------- ----------" << std::endl;
        }

        // resets the tableau to a full identity
        void reset_state() override
        {
            x.resize(rows, std::vector<int>(cols, 0)); // first 2n+1 x n block. first n represents destabilizers
                                                       // second n represents stabilizers + 1 extra row
            z.resize(rows, std::vector<int>(cols, 0)); // second 2n+1 x n block to form the 2n+1 x 2n sized tableau
            r.resize(rows, 0);                         // column on the right with 2n+1 rows
            // The 2n+1 th row is scratch space

            // Intialize the identity tableau
            for (int32_t i = 0; i < n; i++)
            {
                x[i][i] = 1;
                z[i + n][i] = 1;
            }
            std::cout << "Reset state done" << std::endl;
            /*End initialization*/
        }

        void set_seed(IdxType s) override
        {
            seed = s;
        }

        // Prints the tableau including phase
        void print_res_state() override
        {
            for (int32_t i = 0; i < rows; i++)
            {
                if (((i == (rows / 2)) || (i == rows - 1)))
                {
                    for (int32_t j = -5; j < (n * 2); j++)
                    {
                        std::cout << "-";
                    }
                    std::cout << std::endl;
                }
                for (int32_t j = 0; j < cols; j++)
                {
                    std::cout << x[i][j];
                }
                std::cout << " | ";
                for (int32_t j = 0; j < cols; j++)
                {
                    std::cout << z[i][j];
                }
                std::cout << "|" << r[i] << std::endl;
            }
            std::cout << std::endl;
        }

        // Prints a 2n x m tableau matrix w/o phase bits
        void print_table(std::vector<std::vector<int>> &M)
        {
            for (int32_t i = 0; i < M[0].size() + 1; i++)
            {
                std::cout << "- ";
            }
            std::cout << std::endl;
            for (int32_t i = 0; i < M.size(); i++)
            {
                for (int32_t j = 0; j < M[i].size() / 2; j++)
                {
                    std::cout << M[i][j] << " ";
                }
                std::cout << "| ";
                for (int32_t j = M[i].size() / 2; j < M[i].size(); j++)
                {
                    std::cout << M[i][j] << " ";
                }
                std::cout << std::endl;
            }
        }

        std::vector<int32_t> get_measurement_results() override
        {
            std::cout << "Getting measurement results: " << measurement_results.size() << " measurements" << std::endl;
            return measurement_results;
        }

        // Get the stabilizers in the tableau, i.e. all of the Pauli strings that stabilize a certain circuit
        std::vector<std::string> get_stabilizers() override
        {
            int32_t x_val;
            int32_t z_val;
            std::string stabilizers;
            std::vector<std::string> pauliStrings;
            for (int32_t i = (rows >> 1); i < rows - 1; i++) // rows of stabilizers
            {
                stabilizers.clear();               // clear the temporary stabilizers string
                for (int32_t j = 0; j < cols; j++) // qubits/cols
                {
                    x_val = x[i][j];
                    z_val = z[i][j];
                    assert((x_val < 2) && (z_val < 2));
                    if (x_val)
                    {
                        if (z_val)
                            stabilizers += 'Y';
                        else
                            stabilizers += 'X';
                    }
                    else
                    {
                        if (z_val)
                            stabilizers += 'Z';
                        else
                            stabilizers += 'I';
                    }
                } // For columns(qubits)
                pauliStrings.push_back(stabilizers);
            } // For rows(pauli strings)
            return pauliStrings;
        }

        // Simulate the gates from a circuit in the tableau
        void sim2D(std::shared_ptr<Circuit> circuit2D, std::vector<int> gate_chunks, double &sim_time) override
        {
            STAB_CUDA *stab_gpu;
            SAFE_ALOC_GPU(stab_gpu, sizeof(STAB_CUDA));
            // Copy the simulator instance to GPU
            cudaSafeCall(cudaMemcpy(stab_gpu, this,
                                    sizeof(STAB_CUDA), cudaMemcpyHostToDevice));

            gpu_timer sim_timer;

            // Pack tableau and copy to GPU
            pack_tableau();
            cudaSafeCall(cudaMemcpy(x_packed_gpu, x_packed_cpu, packed_matrix_size,
                                    cudaMemcpyHostToDevice));
            cudaSafeCall(cudaMemcpy(z_packed_gpu, z_packed_cpu, packed_matrix_size,
                                    cudaMemcpyHostToDevice));
            cudaSafeCall(cudaMemcpy(r_packed_gpu, r_packed_cpu, packed_r_size, cudaMemcpyHostToDevice));

            std::cout << "Data copied" << std::endl;

            int32_t minGridSize, blockSize;
            cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, simulation_kernel_cuda2D, 0, 0);
            int32_t threadsPerBlockX = 32;
            int32_t threadsPerBlockY = blockSize / threadsPerBlockX;

            if (threadsPerBlockY > 16)
                threadsPerBlockY = 16;

            dim3 threadsPerBlock(threadsPerBlockX, threadsPerBlockY);
            dim3 blocksPerGrid((rows - 1 + threadsPerBlockX - 1) / threadsPerBlockX,
                               (n + threadsPerBlockY - 1) / threadsPerBlockY);

            std::cout << "Blocks calculated" << std::endl;

            std::vector<Gate> gates2D = circuit2D->get_gates();
            int32_t num_gates = gates2D.size();
            copy_gates_to_gpu(gates2D);
            std::cout << "2D circuit parsed" << std::endl;

            /*Simulate*/
            if (Config::PRINT_SIM_TRACE)
            {
                printf("STABSim_gpu is running! Using %lld qubits.\n", n);
            }

            // Call the kernel for each set of gates
            IdxType gate_index_sum = 0;
            int32_t gate_chunk_size = 0;
            sim_timer.start_timer();
            for (int32_t i = 0; i < gate_chunks.size(); i++)
            {
                gate_index_sum += gate_chunk_size;
                gate_chunk_size = gate_chunks[i];
                // std::cout << "Launching 2D kernel " << num_gates << std::endl;
                // std::cout << "Gate chunk sie " << gate_chunk_size << std::endl;

                simulation_kernel_cuda2D<<<blocksPerGrid, threadsPerBlock>>>(stab_gpu, &gates_gpu[gate_index_sum], gate_chunk_size);
                cudaSafeCall(cudaDeviceSynchronize());
            }
            sim_timer.stop_timer();
            /*End simulate*/
            sim_time += sim_timer.measure();

            cudaCheckError();
            // Copy data to the CPU side and unpack
            cudaSafeCall(cudaMemcpy(x_packed_cpu, x_packed_gpu, packed_matrix_size,
                                    cudaMemcpyDeviceToHost));
            cudaSafeCall(cudaMemcpy(z_packed_cpu, z_packed_gpu, packed_matrix_size,
                                    cudaMemcpyDeviceToHost));
            cudaSafeCall(cudaMemcpy(r_packed_cpu, r_packed_gpu, packed_r_size, cudaMemcpyDeviceToHost));
            unpack_tableau();

            if (Config::PRINT_SIM_TRACE)
            {
                printf("\n============== STAB-Sim ===============\n");
                printf("n_qubits:%lld, n_gates:%lld, ncpus:%d, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms\n",
                       n, num_gates, 1, sim_time, 0.,
                       sim_time);
                printf("=====================================\n");
            }

            SAFE_FREE_GPU(gates_gpu);

            //=========================================
        }

        // Simulate the gates from a circuit in the tableau
        void sim(std::shared_ptr<Circuit> circuit, double &sim_time) override
        {
            SAFE_ALOC_GPU(d_local_sums, cols * sizeof(int32_t));
            cudaMemset(d_local_sums, 0, cols * sizeof(int32_t));

            STAB_CUDA *stab_gpu;
            SAFE_ALOC_GPU(stab_gpu, sizeof(STAB_CUDA));
            this->d_local_sums = d_local_sums; // Set the pointer before copying

            gpu_timer sim_timer;

            // Copy bit data to GPU
            copy_bits_to_gpu();

            std::cout << "Data copied" << std::endl;

            std::vector<Gate> gates = circuit->get_gates();

            // Copy gates to the gpu side
            copy_gates_to_gpu(gates);
            IdxType n_gates = gates.size();

            cudaDeviceProp prop;
            cudaGetDeviceProperties(&prop, 0);

            printf("Max grid size: (%d, %d, %d)\n", prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
            printf("Max threads per block: %d\n", prop.maxThreadsPerBlock);
            printf("Device Name: %s\n", prop.name);
            printf("Max Blocks per Multiprocessor: %d\n", prop.maxThreadsPerMultiProcessor);
            printf("Number of SMs: %d\n", prop.multiProcessorCount);

            int32_t threadsPerBlock = prop.maxThreadsPerBlock; // 1024
            int32_t maxBlocksPerSM;
            cudaOccupancyMaxActiveBlocksPerMultiprocessor(
                &maxBlocksPerSM,                /* out: max active blocks */
                (void *)simulation_kernel_cuda, /* kernel */
                threadsPerBlock,                /* threads per block */
                0                               /* shared memory per block */
            );

            printf("Max active blocks per SM: %d\n", maxBlocksPerSM);
            int32_t maxBlocks = prop.multiProcessorCount * maxBlocksPerSM;
            printf("Total max cooperative blocks: %d\n", maxBlocks);
            cudaFuncAttributes attr;
            cudaFuncGetAttributes(&attr, simulation_kernel_cuda);
            printf("Registers per thread: %d\n", attr.numRegs);
            if (attr.numRegs > 64)
            {
                std::cerr << "Warning: High register usage per thread: " << attr.numRegs << std::endl;
            }

            int32_t minGridSize, blockSize;
            cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, simulation_kernel_cuda, 0, 0);
            printf("Max block size: %d\n", blockSize);
            int32_t blocksPerGrid = maxBlocks;

            // In sim() before creating stab_gpu:
            // Allocate and initialize d_local_sums
            // Create and copy stab_gpu AFTER allocation

            /*Simulate*/
            if (Config::PRINT_SIM_TRACE)
            {
                printf("STABSim_gpu is running! Using %lld qubits.\n", n);
            }

            void *args[] = {&stab_gpu, &gates_gpu, &n_gates, &seed};

            // Initialize the measurement counter on the device to 0
            if (d_measurement_idx_counter)
            {
                int32_t initial_count = 0;
                cudaMemcpy(d_measurement_idx_counter, &initial_count, sizeof(int32_t), cudaMemcpyHostToDevice);
            }
            // Copy the simulator instance to GPU
            cudaSafeCall(cudaMemcpy(stab_gpu, this,
                                    sizeof(STAB_CUDA), cudaMemcpyHostToDevice));

            int32_t warpsPerBlock = threadsPerBlock / 32;
            size_t sharedMemSize = warpsPerBlock * sizeof(int32_t);

            cudaDeviceSynchronize(); // Ensure the initialization is complete
            /*Simulate*/
            std::cout << "\n -------------------- \n Simulation starting! \n -------------------- \n"
                      << std::endl;
            sim_timer.start_timer();

            // Launch with cooperative kernel
            cudaLaunchCooperativeKernel((void *)simulation_kernel_cuda, blocksPerGrid, threadsPerBlock, args, sharedMemSize);

            sim_timer.stop_timer();

            cudaSafeCall(cudaDeviceSynchronize());

            cudaCheckError();

            if (cudaGetLastError() != cudaSuccess)
            {
                cudaError_t err = cudaGetLastError();
                if (err != cudaSuccess)
                {
                    std::cerr << "CUDA error after kernel: " << cudaGetErrorString(err) << std::endl;
                    throw std::runtime_error("CUDA kernel had an issue during execution.");
                }
            }

            /*End simulate*/
            sim_time += sim_timer.measure();

            // Add to free_measurement_buffers():

            // Copy measurement results from device to host
            if (d_measurement_idx_counter)
            {
                int32_t num_measurements;
                cudaMemcpy(&num_measurements, d_measurement_idx_counter, sizeof(int32_t), cudaMemcpyDeviceToHost);
                // std::cout << "Number of measurements: " << num_measurements << std::endl;
                if (num_measurements > 0)
                {
                    cudaMemcpy(h_measurement_results, d_measurement_results, num_measurements * sizeof(int32_t), cudaMemcpyDeviceToHost);
                    measurement_results.resize(num_measurements);
                    std::copy(h_measurement_results, h_measurement_results + num_measurements, measurement_results.begin());
                }
            }

            cudaCheckError();
            // Copy data to the CPU side and unpack
            copy_bits_from_gpu();

            if (Config::PRINT_SIM_TRACE)
            {
                printf("\n============== STAB-Sim ===============\n");
                printf("n_qubits:%lld, n_gates:%lld, ncpus:%d, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms\n",
                       n, n_gates, 1, sim_time, 0.,
                       sim_time);
                printf("=====================================\n\n\n");
            }

            SAFE_FREE_GPU(gates_gpu);
            SAFE_FREE_GPU(d_measurement_idx_counter);
            SAFE_FREE_GPU(d_measurement_results);

            SAFE_FREE_GPU(d_local_sums);
            SAFE_FREE_GPU(stab_gpu);

            //=========================================
        }

        IdxType measure(IdxType qubit) override
        {
            throw std::logic_error("measure Not implemented (STAB_GPU)");
        }
        IdxType *measure_all(IdxType repetition) override
        {
            throw std::logic_error("measure_all Not implemented (STAB_GPU)");
        }
        void set_initial(std::string fpath, std::string format) override
        {
            throw std::logic_error("set_initial Not implemented (STAB_GPU)");
        }
        void dump_res_state(std::string outfile) override
        {
            throw std::logic_error("dump_res_state Not implemented (STAB_GPU)");
        }
        ValType *get_real() const override
        {
            throw std::logic_error("get_real Not implemented (STAB_GPU)");
        }
        ValType *get_imag() const override
        {
            throw std::logic_error("get_imag Not implemented (STAB_GPU)");
        }
        ValType get_exp_z(const std::vector<size_t> &in_bits) override
        {
            throw std::logic_error("get_exp_z Not implemented (STAB_GPU)");
        }
        ValType get_exp_z() override
        {
            throw std::logic_error("get_exp_z Not implemented (STAB_GPU)");
        }
        IdxType *get_results() override
        {
            throw std::logic_error("get_results Not implemented (STAB_CPU)");
        }

    public:
        IdxType n;
        IdxType n_gates;
        IdxType stabCounts;
        IdxType mask;
        IdxType rows;
        IdxType packed_rows;
        IdxType cols;
        int32_t packed_bits;
        IdxType packed_r_size;
        IdxType packed_matrix_size;
        IdxType bit_matrix_size;
        IdxType bit_r_size;
        std::vector<std::vector<int>> x;
        std::vector<std::vector<int>> z;
        std::vector<int> r;
        std::vector<int> measurement_results; // Store measurement results from the device

        std::vector<std::vector<Gate>> layered_gates;
        IdxType num_layers;
        // CPU Arrays
        int32_t *x_packed_cpu = nullptr;
        int32_t *z_packed_cpu = nullptr;
        int32_t *r_packed_cpu = nullptr;
        // GPU Arrays
        int32_t *x_packed_gpu = nullptr;
        int32_t *z_packed_gpu = nullptr;
        int32_t *r_packed_gpu = nullptr;
        int32_t *x_bit_gpu = nullptr;
        int32_t *z_bit_gpu = nullptr;
        int32_t *r_bit_gpu = nullptr;

        int32_t *d_local_sums = nullptr;

        IdxType seed;

        int32_t *d_measurement_results = nullptr;
        int32_t *h_measurement_results = nullptr;
        int32_t *d_measurement_idx_counter = nullptr;
        IdxType measurement_count = 0;

        void allocate_measurement_buffers(IdxType max_measurements) override
        {
            if (measurement_count == max_measurements && d_measurement_results != nullptr)
                return;

            free_measurement_buffers();
            measurement_count = max_measurements;

            if (measurement_count > 0)
            {
                SAFE_ALOC_GPU(d_measurement_results, measurement_count * sizeof(IdxType));
                h_measurement_results = new IdxType[measurement_count];
                SAFE_ALOC_GPU(d_measurement_idx_counter, sizeof(IdxType));
            }
        }

        void free_measurement_buffers() override
        {
            SAFE_FREE_GPU(d_measurement_results);
            delete[] h_measurement_results;
            SAFE_FREE_GPU(d_measurement_idx_counter);
            d_measurement_results = nullptr;
            h_measurement_results = nullptr;
            d_measurement_idx_counter = nullptr;
        }

        // Random
        std::mt19937 rng;
        std::uniform_int_distribution<int> dist;

        // GPU-side gate instance
        Gate *gates_gpu = nullptr;
        std::vector<std::vector<Gate>> gates2D;

        void copy_gates_to_gpu(std::vector<Gate> &cpu_vec)
        {
            // Allocate memory on CPU
            size_t vec_size = cpu_vec.size() * sizeof(Gate);

            // Allocate memory on GPU
            SAFE_FREE_GPU(gates_gpu);
            SAFE_ALOC_GPU(gates_gpu, vec_size);
            cudaSafeCall(cudaMemcpy(gates_gpu, cpu_vec.data(), vec_size, cudaMemcpyHostToDevice));
        }
    }; // End tableau class

    __device__ int32_t global_p;

    __inline__ __device__ void global_32_reduction(int32_t *input, int32_t &lane_id, IdxType &i)
    {
        // Warp-level reduction of a large array smaller than the threads available, down to an array of size/32

        int32_t local = input[i];

        int32_t total = __reduce_add_sync(__activemask(), local);

        if (lane_id == 0)
            input[i / 5] = total; // Overwrite the first size/32 values
    }

    __inline__ __device__ int32_t singular_reduction(int32_t *input, int32_t &lane_id, IdxType max_index)
    {
        int32_t local_sum = 0;

        // Each thread in warp processes multiple elements with stride 32
        for (int32_t j = lane_id; j < max_index; j += 32)
        {
            local_sum += input[j];
        }

        // Warp-level reduction of the partial sums
        int32_t total = __reduce_add_sync(__activemask(), local_sum);

        // All threads return the result, but only the first should be used
        if ((lane_id) == 0)
            return total;
        return 0;
    }

    // if (i==0) printf("====== rows:%llu, cols:%llu, blockDim.x:%u, gridDim.x:%u\n", rows, cols, blockDim.x, gridDim.x);

    __device__ int32_t smuggle_scratch;
    __launch_bounds__(1024, 1)
        __global__ void simulation_kernel_cuda(STAB_CUDA *stab_gpu, Gate *gates_gpu, IdxType n_gates, int32_t seed)
    {
        cg::grid_group grid = cg::this_grid();
        IdxType i = blockIdx.x * blockDim.x + threadIdx.x;
        IdxType rows = stab_gpu->rows;
        IdxType cols = stab_gpu->cols;

        int32_t lane_id = i % 32;
        int32_t warp_id = i / 32;
        int32_t n_warps = blockDim.x * gridDim.x / 32;
        extern __shared__ int32_t shared_mem[];
        int32_t block_warp_id = threadIdx.x / 32;
        int32_t block_lane_id = threadIdx.x % 32;
        IdxType n_qubits = stab_gpu->n;
        int32_t *x_arr = stab_gpu->x_bit_gpu;
        int32_t *z_arr = stab_gpu->z_bit_gpu;
        int32_t *r_arr = stab_gpu->r_bit_gpu;
        int32_t *global_sums = stab_gpu->d_local_sums;

        int *d_measurement_idx_counter = stab_gpu->d_measurement_idx_counter;
        int32_t x, z, x_ctrl, z_ctrl;
        OP op_name;
        int32_t a, b;
        int32_t p;

        // Initialize shared memory (per block)
        __shared__ int32_t block_p_shared;
        int32_t scratch_row = rows - 1;

        for (int32_t k = 0; k < n_gates; k++)
        {
            op_name = gates_gpu[k].op_name;
            a = gates_gpu[k].qubit;
            b = gates_gpu[k].ctrl;
            IdxType index = i * cols + a;
            IdxType ctrl_index = i * cols + b;

            switch (op_name)
            {
            case OP::H:
            {
                if (i < scratch_row)
                {
                    x = x_arr[index];
                    z = z_arr[index];
                    // Phase
                    r_arr[i] ^= (x & z);
                    // Entry -- swap x and z bits
                    x_arr[index] = z;
                    z_arr[index] = x;
                }
                break;
            }

            case OP::S:
                if (i < scratch_row)
                {
                    x = x_arr[index];
                    z = z_arr[index];
                    // Phase
                    r_arr[i] ^= (x & z);
                    // Entry
                    z_arr[index] ^= x;
                    // if(i == 0) printf("S\n\n");
                }
                break;

            case OP::SDG:
                if (i < scratch_row)
                {
                    x = x_arr[index];
                    z = z_arr[index];
                    // Phase
                    r_arr[i] ^= (x ^ (x & z));
                    // Entry
                    z_arr[index] = z ^ x;
                }
                break;

            case OP::CX:
            {
                if (i < scratch_row)
                {
                    x = x_arr[index];
                    z = z_arr[index];
                    x_ctrl = x_arr[ctrl_index];
                    z_ctrl = z_arr[ctrl_index];
                    // Phase
                    r_arr[i] ^= (x_ctrl & z & (x ^ z_ctrl ^ 1));
                    // Entry
                    x_arr[index] ^= x_ctrl;
                    z_arr[ctrl_index] ^= z;
                }
                break;
            }

            case OP::M:
            {
                // --- Find p (minimum row `i` in [cols, scratch_row) with x_arr[i*cols+a] = 1) ---

                // 1. Initialization
                // The first thread of each block initializes the block's local minimum.
                if (threadIdx.x == 0)
                {
                    block_p_shared = rows;
                }
                // The first thread of the entire grid initializes the global minimum.
                if (i == 0)
                {
                    global_p = rows;
                    smuggle_scratch = 0; // Also reset the accumulator here
                }
                // All threads must wait for initialization to complete.
                // A grid.sync() is necessary to ensure global_p is visible to all.
                grid.sync();

                // 2. Stage 1: In-Block Reduction
                // Each thread in the valid range checks its row.
                if (i >= cols && i < scratch_row)
                {
                    if (x_arr[index])
                    {
                        // Update the block's local minimum atomically.
                        atomicMin(&block_p_shared, i);
                    }
                }
                // Synchronize within the block to ensure all threads have contributed to block_p_shared.
                __syncthreads();

                // 3. Stage 2: In-Grid Reduction
                // The first thread of each block contributes its block's minimum to the global minimum.
                if (threadIdx.x == 0)
                {
                    atomicMin(&global_p, block_p_shared);
                }
                // Synchronize the entire grid to ensure global_p has the final result.
                grid.sync();

                // 4. Finalization
                // All threads read the final value.
                p = global_p;
                // A final sync ensures all threads have the correct `p` before branching.
                grid.sync();

                if (p != rows)
                {
                    // if(i==0) printf("rand p=%d\n", p);

                    //     Random measurement
                    //     If the stabilizer anticommutes, update it (aside from p)

                    for (IdxType local_i = warp_id; local_i < scratch_row; local_i += n_warps)
                    {
                        if (local_i != p && x_arr[local_i * cols + a])
                        {
                            // Many rowsums
                            int32_t local_sum = 0;
                            for (int32_t j = lane_id; j < n_qubits; j += 32)
                            {
                                IdxType h_idx = local_i * cols + j;
                                IdxType p_idx = p * cols + j;
                                int32_t x_p = x_arr[p_idx];
                                int32_t z_p = z_arr[p_idx];
                                int32_t x_h = x_arr[h_idx];
                                int32_t z_h = z_arr[h_idx];

                                local_sum += x_p * (z_p * (z_h - x_h) + (1 - z_p) * z_h * (2 * x_h - 1)) + (1 - x_p) * z_p * x_h * (1 - 2 * z_h);

                                atomicXor(&x_arr[h_idx], x_p);
                                atomicXor(&z_arr[h_idx], z_p);
                            }
                            // __syncthreads();
                            // merge warp-wise local_sum
                            for (int32_t offset = 16; offset > 0; offset /= 2)
                            {
                                local_sum += __shfl_down_sync(0xffffffff, local_sum, offset);
                            }
                            // __syncthreads();
                            // per head-lane owns the merged local_sum and performs adjustment
                            if (lane_id == 0)
                            {
                                // if((local_sum % 4 != 0) && (abs(local_sum % 2) != 0))
                                //     printf("random m%d local = %d\n", *d_measurement_idx_counter, local_sum);
                                local_sum += 2 * r_arr[p] + 2 * r_arr[local_i];
                                r_arr[local_i] = (local_sum % 4 == 0) ? 0 : 1;
                                // printf("rand r%d row%lld\n", r_arr[local_i], local_i);
                            }
                        }
                    }

                    //                     for (int k = 0; k < scratch_row; k++)
                    //                     {
                    //                         if (k != p && x_arr[k * cols + a])
                    //                         {
                    //                             extern __shared__ int32_t shared_mem[];

                    //                             if(i < cols)
                    //                             {
                    //                                 // All threads check for anticommutation. This is a safe, redundant read.
                    //                                 IdxType h_idx = k * cols + i;
                    //                                 IdxType p_idx = p * cols + i;

                    //                                 int32_t thread_sum = 0;

                    //                                 // --- STAGE 1: Parallel Phase Calculation ---
                    //                                 // ALL threads (including i=0) participate in calculating the phase contribution.
                    //                                 // for (int32_t j = i; j < n_qubits; j += (gridDim.x * blockDim.x)) {
                    //                                 int32_t x_p = x_arr[p_idx];
                    //                                 int32_t z_p = z_arr[p_idx];
                    //                                 int32_t x_h = x_arr[h_idx];
                    //                                 int32_t z_h = z_arr[h_idx];
                    // // ... existing code ...
                    //                                 thread_sum += x_p * ( z_p * (z_h - x_h) + (1 - z_p) * z_h * (2 * x_h - 1) )
                    //                                             + (1 - x_p) * z_p * x_h * (1 - 2 * z_h);

                    //                                 x_arr[h_idx]^= x_p;
                    //                                 z_arr[h_idx]^= z_p;

                    //                                 // --- FASTER WARP-BASED REDUCTION ---

                    //                                 // --- CORRECTED: Use block-local threadIdx.x for shared memory operations ---

                    //                                 // STAGE 1: Intra-Warp Reduction
                    //                                 // Each warp reduces its 32 thread_sum values. This is divergence-free.
                    //                                 int32_t warp_sum = __reduce_add_sync(__activemask(), thread_sum);

                    //                                 // The first thread of each warp writes its partial sum to shared memory.
                    //                                 if (block_lane_id == 0)
                    //                                 {
                    //                                     shared_mem[block_warp_id] = warp_sum;
                    //                                 }
                    //                             }

                    //                             // Synchronize to ensure all warp sums are written to shared memory.
                    //                             __syncthreads();

                    //                             if(i < cols)
                    //                             {
                    //                                 // STAGE 2: Inter-Warp Reduction
                    //                                 // The first warp (32 threads) now reduces the 32 partial sums from shared memory.
                    //                                 // All other warps are idle, but this is a very small amount of work.
                    //                                 thread_sum = shared_mem[block_lane_id];
                    //                                 int32_t block_sum = __reduce_add_sync(__activemask(), thread_sum);

                    //                                 // The first thread of the block (thread 0) has the final sum and adds it to the global accumulator.
                    //                                 if(threadIdx.x == 0) {
                    //                                     atomicAdd(&smuggle_scratch, block_sum);
                    //                                 }//The first thread of the grid now has a global accumulated value to use
                    //                             }

                    //                             grid.sync(); //Sync for global sum, and update the scratch row while we're here

                    //                             if(i == 0)
                    //                             {
                    //                                 r_arr[k] =
                    //                                     (((smuggle_scratch + 2 * r_arr[p] + 2 * r_arr[k])%4)==0) ? 0 : 1;
                    //                                 smuggle_scratch = 0;
                    //                             }

                    //                             // Wait for all blocks to finish their atomicAdds.
                    //                             grid.sync();
                    // // ... existing code ...
                    //                         }
                    //                     }

                    grid.sync();

                    if (i < cols)
                    {
                        IdxType p_index = (p * cols) + i;
                        IdxType row_col_index = ((p - cols) * cols) + i;
                        x_arr[row_col_index] = x_arr[p_index];
                        z_arr[row_col_index] = z_arr[p_index];
                        x_arr[p_index] = 0;

                        if (i == a)
                            z_arr[p_index] = 1;
                        if (i != a)
                            z_arr[p_index] = 0;
                    }

                    if (i == 0)
                    {
                        // Draw bit from seed and measurement index
                        int32_t meas_idx = atomicAdd(d_measurement_idx_counter, 1);
                        // printf("Seed for measurement %d: %d\n", meas_idx, seed);
                        int32_t bit = prng_bit(seed, meas_idx);
                        // printf("rand m%d, res%d\n",  meas_idx, bit);

                        r_arr[p] = bit;
                        stab_gpu->d_measurement_results[meas_idx] = bit;
                    }
                }
                else
                { // Deterministic measurement
                    if (i < cols)
                    {
                        IdxType scratch_index = (scratch_row * cols) + i;
                        x_arr[scratch_index] = 0;
                        z_arr[scratch_index] = 0;
                    }
                    if (i == 0)
                    {
                        r_arr[scratch_row] = 0;
                        smuggle_scratch = 0;
                    }

                    // ... existing code ...
                    // This sync ensures the scratch row is zeroed out before we begin.
                    grid.sync();

                    // --- CORRECT HYBRID SERIAL-PARALLEL LOOP ---
                    // The loop over stabilizers 'k' is serial to respect the recursive phase dependency.
                    // ALL threads execute this loop. Thread 0 has the special role of updating the
                    // tableau state after each parallel computation is complete.
                    IdxType h_idx = scratch_row * cols + i;

                    for (int k = 0; k < cols; k++)
                    {
                        // All threads check for anticommutation. This is a safe, redundant read.
                        bool anticommutes = (x_arr[k * cols + a]);
                        if (anticommutes)
                        {
                            int32_t p_row = k + cols;
                            IdxType p_idx = p_row * cols + i;
                            int32_t thread_sum = 0;

                            if (i < cols)
                            {
                                // --- STAGE 1: Parallel Phase Calculation ---
                                // ALL threads (including i=0) participate in calculating the phase contribution.
                                // for (int32_t j = i; j < n_qubits; j += (gridDim.x * blockDim.x)) {
                                int32_t x_p = x_arr[p_idx];
                                int32_t z_p = z_arr[p_idx];
                                int32_t x_h = x_arr[h_idx];
                                int32_t z_h = z_arr[h_idx];
                                // ... existing code ...
                                thread_sum += x_p * (z_p * (z_h - x_h) + (1 - z_p) * z_h * (2 * x_h - 1)) + (1 - x_p) * z_p * x_h * (1 - 2 * z_h);
                                // thread_sum += (x_h & z_p) - (x_p & z_h) +
                                //     (((x_h & x_p & z_h) - (x_h & x_p & z_p) - (x_h & z_h & z_p) +
                                //     (x_p & z_h & z_p))<<1);

                                x_arr[h_idx] ^= x_p;
                                z_arr[h_idx] ^= z_p;

                                // --- FASTER WARP-BASED REDUCTION ---

                                // --- CORRECTED: Use block-local threadIdx.x for shared memory operations ---

                                // STAGE 1: Intra-Warp Reduction
                                // Each warp reduces its 32 thread_sum values. This is divergence-free.
                                int32_t warp_sum = __reduce_add_sync(__activemask(), thread_sum);

                                // The first thread of each warp writes its partial sum to shared memory.
                                if (block_lane_id == 0)
                                {
                                    shared_mem[block_warp_id] = warp_sum;
                                }
                            }

                            // Synchronize to ensure all warp sums are written to shared memory.
                            // There are 32 shared memory slots per 1024 thread block
                            __syncthreads();

                            int32_t block_sum = 0;
                            if (threadIdx.x < 32)
                            {
                                // Only threads in the first warp of each block load data.
                                //  The first warp in each block reduces its partial sums
                                block_sum = __reduce_add_sync(__activemask(), shared_mem[threadIdx.x]);
                                // Every block adds to the final sum as it gets its value
                                if (threadIdx.x == 0 && block_sum != 0)
                                {
                                    atomicAdd(&smuggle_scratch, block_sum);
                                }
                            }

                            grid.sync(); // Sync for global sum, and update the scratch row

                            if (i == 0) // Update the scratch bit and reset the global accumulator
                            {
                                // printf("r[row]%d ", r_arr[p_row]);

                                r_arr[scratch_row] =
                                    (((smuggle_scratch + 2 * r_arr[p_row] + 2 * r_arr[scratch_row]) % 4) == 0) ? 0 : 1;

                                // printf("determ r%d row%d\n", r_arr[scratch_row], p_row);

                                smuggle_scratch = 0;
                            }

                            // No particular reason
                            grid.sync();
                            // ... existing code ...
                        }
                    }

                    // After the loop, thread 0 stores the final measurement result.
                    if (i == 0)
                    {
                        // printf("determ m%d, res%d\n", *d_measurement_idx_counter, r_arr[scratch_row]);
                        int32_t meas_idx = atomicAdd(d_measurement_idx_counter, 1);
                        stab_gpu->d_measurement_results[meas_idx] = r_arr[scratch_row];
                    }
                    // ... existing code ...

                    /*
                    Unsafe version
                    */

                    // for (IdxType local_i = warp_id; local_i < cols; local_i += n_warps) //capture only the stabilizer rows (num stabilizer rows == num cols)
                    // {
                    //     IdxType local_index = local_i * cols + a;
                    //     if (x_arr[local_index])
                    //     {
                    //         //Parallel rowsum
                    //         int32_t local_sum = 0;
                    //         IdxType p_row = local_i + cols;

                    //         for (int32_t j = lane_id; j < n_qubits; j+=32) //Capture all the columns in a row
                    //         {
                    //             IdxType h_idx = scratch_row * cols + j;
                    //             IdxType p_idx = p_row * cols + j;
                    //             int32_t x_p = x_arr[p_idx];
                    //             int32_t z_p = z_arr[p_idx];
                    //             int32_t x_h = x_arr[h_idx];
                    //             int32_t z_h = z_arr[h_idx];

                    //             local_sum += x_p * ( z_p * (z_h - x_h) + (1 - z_p) * z_h * (2 * x_h - 1) )
                    //                 + (1 - x_p) * z_p * x_h * (1 - 2 * z_h);

                    //             atomicXor(&x_arr[h_idx], x_p);
                    //             atomicXor(&z_arr[h_idx], z_p);
                    //         }
                    //         //Merge warp-wise local_sum
                    //         for (int32_t offset = 16; offset > 0; offset /= 2)
                    //         {
                    //             local_sum += __shfl_down_sync(0xffffffff, local_sum, offset);
                    //         }
                    //         //Per head-lane owns the merged local_sum and contributes to the final sum
                    //         if (lane_id == 0)
                    //         {
                    //             int32_t row_term = 2*((((local_sum + 2*r_arr[p_row])%4)==0) ? 0:1);

                    //             //The length of the queue for this atomicAdd is \Theta(n/32) (or O(n/32) if we don't let 0's contribute)
                    //             atomicAdd(&smuggle_scratch, row_term); //Only one add per warp.

                    //             // printf("global_sums[%lld] = %d + 2 * %d\n", local_i, local_sum, r_arr[p_row]);
                    //             // global_sums[local_i] = local_sum + 2 * r_arr[p_row];
                    //             // printf("m%d sum: %d\n", *d_measurement_idx_counter, global_sums[local_i]%4);

                    //             // if((((global_sums[local_i]+2*smuggle_scratch)%4) == 0) ? false:true) atomicXor(&smuggle_scratch, 1);
                    //         }
                    //     }
                    //     // else
                    //     // {
                    //     //     if (lane_id == 0)
                    //     //     {
                    //     //         printf("global_sums[%lld] = 0\n", local_i);
                    //     //         global_sums[local_i] = 0;
                    //     //     }
                    //     // }
                    // }

                    // grid.sync();

                    // // if(i < 32) //First warp does the final reduction
                    // // {
                    // // int32_t total = singular_reduction(global_sums, lane_id, cols);
                    // if(i == 0)
                    // {
                    //     // if((total % 4 != 0) && (abs(total % 2) != 0))
                    //     //     printf("determ m%d total = %d\n", *d_measurement_idx_counter, total);
                    //     // printf("m%d sum: %d\n", *d_measurement_idx_counter, total);

                    //     // total += 2 * smuggle_scratch;

                    //     // printf("r[h]%d: %d\n", *d_measurement_idx_counter, r_arr[scratch_row]);
                    //     r_arr[scratch_row] = ((smuggle_scratch%4) == 0) ? 0: 1;
                    //     // printf("m%d r_scratch = %d\n", *d_measurement_idx_counter, r_arr[scratch_row]);

                    //     int32_t meas_idx = atomicAdd(d_measurement_idx_counter, 1);
                    //     stab_gpu->d_measurement_results[meas_idx] = r_arr[scratch_row];
                    // }

                    /*
                        End unsafe version
                    */

                    /*
                    Direct match version
                    This directly matches stab-cpu as a sanity check
                    */

                    // if(i == 0)
                    // {
                    //     for(int k = 0; k < cols; k++)
                    //     {
                    //         if(x_arr[k * cols + a])
                    //         {
                    //             IdxType p_row = k + cols;
                    //             int32_t local_sum = 0;
                    //             for(int j = 0; j < cols; j++)
                    //             {
                    //                 IdxType h_idx = scratch_row * cols + j;
                    //                 IdxType p_idx = p_row * cols + j;
                    //                 int32_t x_p = x_arr[p_idx];
                    //                 int32_t z_p = z_arr[p_idx];
                    //                 int32_t x_h = x_arr[h_idx];
                    //                 int32_t z_h = z_arr[h_idx];
                    //                 local_sum += x_p * ( z_p * (z_h - x_h) + (1 - z_p) * z_h * (2 * x_h - 1) )
                    //                     + (1 - x_p) * z_p * x_h * (1 - 2 * z_h); //AL version
                    //                 x_arr[h_idx] ^= x_p;
                    //                 z_arr[h_idx] ^= z_p;
                    //             }
                    //             // if((total % 4 != 0) && (abs(total % 2) != 0))
                    //             //     printf("determ m%d tot: %d\n", *d_measurement_idx_counter, local_sum);
                    //             local_sum = local_sum + 2 * r_arr[scratch_row] + 2 * r_arr[p_row];
                    //             r_arr[scratch_row] = (local_sum % 4 == 0) ? 0 : 1;
                    //         }
                    //     }
                    //     int32_t meas_idx = *d_measurement_idx_counter;
                    //     stab_gpu->d_measurement_results[meas_idx] = r_arr[scratch_row];
                    //     *d_measurement_idx_counter = meas_idx + 1;
                    // }

                    /*End Sanity Check Version*/
                } // End deterministic logic

                grid.sync(); // Make sure the tableau is up to date before resuming operations
                break;
            } // End M

                /*
                End sanity check verison
                */

            case OP::RESET:
            {
                if (i < scratch_row)
                {
                    if (threadIdx.x == 0)
                    {
                        block_p_shared = rows;
                        if (blockIdx.x == 0)
                        {
                            global_p = rows;
                        }
                    }
                }
                grid.sync();

                // In-block
                if (i < scratch_row)
                {
                    if ((i >= cols) && (x_arr[index]))
                    {
                        atomicMin(&block_p_shared, i);
                    }
                }
                grid.sync(); // block sync

                // In-grid
                if ((i < scratch_row) && threadIdx.x == 0)
                {
                    atomicMin(&global_p, block_p_shared);
                }
                grid.sync();

                p = global_p;

                grid.sync();

                if (p != rows)
                {
                    // Random measurement
                    // If the stabilizer anticommutes, update it (aside from p)

                    for (IdxType local_i = warp_id; local_i < scratch_row; local_i += n_warps)
                    {
                        if (local_i != p && x_arr[local_i * cols + a])
                        {
                            // Many rowsums
                            int32_t local_sum = 0;
                            for (int32_t j = lane_id; j < n_qubits; j += 32)
                            {
                                IdxType h_idx = local_i * cols + j;
                                IdxType p_idx = p * cols + j;
                                int32_t x_p = x_arr[p_idx];
                                int32_t z_p = z_arr[p_idx];
                                int32_t x_h = x_arr[h_idx];
                                int32_t z_h = z_arr[h_idx];

                                local_sum += x_p * (z_p * (z_h - x_h) + (1 - z_p) * z_h * (2 * x_h - 1)) + (1 - x_p) * z_p * x_h * (1 - 2 * z_h);

                                atomicXor(&x_arr[h_idx], x_p);
                                atomicXor(&z_arr[h_idx], z_p);
                            }
                            // __syncthreads();
                            // merge warp-wise local_sum
                            for (int32_t offset = 16; offset > 0; offset /= 2)
                            {
                                local_sum += __shfl_down_sync(0xffffffff, local_sum, offset);
                            }
                            // __syncthreads();
                            // per head-lane owns the merged local_sum and performs adjustment
                            if (lane_id == 0)
                            {
                                // if((local_sum % 4 != 0) && (abs(local_sum % 2) != 0))
                                //     printf("random m%d local = %d\n", *d_measurement_idx_counter, local_sum);
                                local_sum += 2 * r_arr[p] + 2 * r_arr[local_i];
                                r_arr[local_i] = (local_sum % 4 == 0) ? 0 : 1;
                            }
                        }
                    }

                    grid.sync();

                    if (i < cols)
                    {
                        IdxType p_index = (p * cols) + i;
                        IdxType row_col_index = ((p - cols) * cols) + i;
                        x_arr[row_col_index] = x_arr[p_index];
                        z_arr[row_col_index] = z_arr[p_index];
                        x_arr[p_index] = 0;

                        if (i == (p - cols))
                            z_arr[p_index] = 1;
                        if (i != (p - cols))
                            z_arr[p_index] = 0;
                    }

                    if (i == 0)
                    {
                        // Draw bit from seed and measurement index
                        int32_t meas_idx = atomicAdd(d_measurement_idx_counter, 1);
                        // printf("Seed for measurement %d: %d\n", meas_idx, seed);
                        int32_t bit = prng_bit(seed, meas_idx);
                        r_arr[p] = bit;
                        stab_gpu->d_measurement_results[meas_idx] = bit;
                    }

                    grid.sync();

                    if ((i < scratch_row) && (r_arr[p] == 1))
                    {
                        r_arr[i] ^= z_arr[index];
                    }

                    grid.sync();
                }
                else
                { // Deterministic measurement
                    if ((i < scratch_row) && (i < cols))
                    {
                        IdxType scratch_index = (scratch_row * cols) + i;
                        x_arr[scratch_index] = 0;
                        z_arr[scratch_index] = 0;
                    }
                    if (i == 0)
                    {
                        r_arr[scratch_row] = 0;
                        smuggle_scratch = 0;
                    }

                    // ... existing code ...
                    // This sync ensures the scratch row is zeroed out before we begin.
                    grid.sync();

                    // --- CORRECT HYBRID SERIAL-PARALLEL LOOP ---
                    // The loop over stabilizers 'k' is serial to respect the recursive phase dependency.
                    // ALL threads execute this loop. Thread 0 has the special role of updating the
                    // tableau state after each parallel computation is complete.
                    IdxType h_idx = scratch_row * cols + i;

                    for (int k = 0; k < cols; k++)
                    {

                        // All threads check for anticommutation. This is a safe, redundant read.
                        bool anticommutes = (x_arr[k * cols + a]);
                        if (anticommutes)
                        {
                            int32_t p_row = k + cols;
                            IdxType p_idx = p_row * cols + i;

                            int32_t thread_sum = 0;

                            extern __shared__ int32_t shared_mem[];
                            int32_t block_warp_id = threadIdx.x / 32;
                            int32_t block_lane_id = threadIdx.x % 32;

                            if (i < cols)
                            {
                                // --- STAGE 1: Parallel Phase Calculation ---
                                // ALL threads (including i=0) participate in calculating the phase contribution.
                                // for (int32_t j = i; j < n_qubits; j += (gridDim.x * blockDim.x)) {
                                int32_t x_p = x_arr[p_idx];
                                int32_t z_p = z_arr[p_idx];
                                int32_t x_h = x_arr[h_idx];
                                int32_t z_h = z_arr[h_idx];
                                // ... existing code ...
                                thread_sum += x_p * (z_p * (z_h - x_h) + (1 - z_p) * z_h * (2 * x_h - 1)) + (1 - x_p) * z_p * x_h * (1 - 2 * z_h);

                                x_arr[h_idx] ^= x_p;
                                z_arr[h_idx] ^= z_p;

                                // --- FASTER WARP-BASED REDUCTION ---

                                // --- CORRECTED: Use block-local threadIdx.x for shared memory operations ---

                                // STAGE 1: Intra-Warp Reduction
                                // Each warp reduces its 32 thread_sum values. This is divergence-free.
                                int32_t warp_sum = __reduce_add_sync(__activemask(), thread_sum);

                                // The first thread of each warp writes its partial sum to shared memory.
                                if (block_lane_id == 0)
                                {
                                    shared_mem[block_warp_id] = warp_sum;
                                }
                            }

                            // Synchronize to ensure all warp sums are written to shared memory.
                            __syncthreads();

                            if (i < cols)
                            {
                                // STAGE 2: Inter-Warp Reduction
                                // The first warp (32 threads) now reduces the 32 partial sums from shared memory.
                                // All other warps are idle, but this is a very small amount of work.
                                thread_sum = (threadIdx.x < 32) ? shared_mem[block_lane_id] : 0;
                                int32_t block_sum = __reduce_add_sync(__activemask(), thread_sum);

                                // The first thread of the block (thread 0) has the final sum and adds it to the global accumulator.
                                if (threadIdx.x == 0)
                                {
                                    atomicAdd(&smuggle_scratch, block_sum);
                                } // The first thread of the grid now has a global accumulated value to use
                            }

                            grid.sync(); // Sync for global sum, and update the scratch row while we're here

                            if (i == 0)
                            {
                                r_arr[scratch_row] =
                                    (((smuggle_scratch + 2 * r_arr[p_row] + 2 * r_arr[scratch_row]) % 4) == 0) ? 0 : 1;
                                smuggle_scratch = 0;
                            }

                            // Wait for all blocks to finish their atomicAdds.
                            grid.sync();
                            // ... existing code ...
                        }
                    }

                    grid.sync();

                    if ((i < scratch_row) && (r_arr[scratch_row] == 1))
                    {
                        r_arr[i] ^= z_arr[index];
                    }
                } // End deterministic logic
                break;
            } // End RESET

            default:
                printf("Non-Clifford or unrecognized gate: %d\n", op_name);
            }
        }

        // if (threadIdx.x == 0 && blockIdx.x == 0) printf("Kernel is done!\n");
        //  printf("Kernel is done!\n");
    } // end kernel

    __global__ void simulation_kernel_cuda2D(STAB_CUDA *stab_gpu, Gate *gates_gpu, IdxType gate_chunk)
    {
        int32_t row = blockIdx.x * blockDim.x + threadIdx.x; // Index for stabilizers (rows)
        int32_t col = blockIdx.y * blockDim.y + threadIdx.y; // Index for gates (columns)

        if (row >= stab_gpu->packed_rows)
            return; // Check out-of-bounds for qubits

        if (col >= gate_chunk)
            return; // Check out-of-bounds for gate

        // printf("Inside 2D kernel %d, gate chunk %lld gate qubit %lld \n", col, gate_chunk, gates_gpu[col].qubit);
        // printf("Gates gpu size %lld \n", (sizeof(gates_gpu)));

        int32_t target = gates_gpu[col].qubit; // Qubit target
        OP op_name = gates_gpu[col].op_name;   // Operation to perform
        int32_t *x_arr = stab_gpu->x_packed_gpu;
        int32_t *z_arr = stab_gpu->z_packed_gpu;

        // Calculate the index for this qubit in the packed arrays
        IdxType index = row * stab_gpu->cols + target;

        // Perform operations for all possible gates, but mask non-relevant ones
        // Start with the common operation - calculate phase and entry for all gates
        int32_t x = x_arr[index];
        int32_t z = z_arr[index];

        // Gate operations
        if (op_name == OP::H)
        {
            // printf("Inside h gate \n");
            // H Gate: Swap x and z
            stab_gpu->r_packed_gpu[row] ^= (x & z);
            // Entry -- swap x and z bits
            x_arr[index] = z;
            z_arr[index] = x;
            return;
        }
        if (op_name == OP::S)
        {
            // S Gate: Entry (z_arr[index] ^= x_arr[index])
            stab_gpu->r_packed_gpu[row] ^= (x & z);
            z_arr[index] = x ^ z;
            return;
        }
        if (op_name == OP::SDG)
        {
            // SDG Gate: Phase (x_arr[index] ^ (x_arr[index] & z_arr[index]))
            stab_gpu->r_packed_gpu[row] ^= (x ^ (x & z)); // SDG Phase operation
            z_arr[index] = x ^ z;                         // Entry (same as S gate)
            return;
        }
        if (op_name == OP::CX)
        {
            int32_t ctrl_index = row * stab_gpu->cols + gates_gpu[col].ctrl;

            int32_t x_ctrl = x_arr[ctrl_index];
            int32_t z_ctrl = z_arr[ctrl_index];

            // Phase
            stab_gpu->r_packed_gpu[row] ^= ((x_ctrl & z) & (x ^ z_ctrl ^ 1));

            // Entry
            x_arr[index] = x ^ x_ctrl;
            z_arr[ctrl_index] = z ^ z_ctrl;
            return;
        }
    }

} // namespace NWQSim

// #endif
