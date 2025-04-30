#pragma once

#include "../state.hpp"
#include "../nwq_util.hpp"
#include "../gate.hpp"
#include "../circuit.hpp"

#include "../config.hpp"
#include "../private/cuda_util.cuh"
#include "../private/macros.hpp"
#include "../private/sim_gate.hpp"

#include "../circuit_pass/fusion.hpp"
#include "../private/exp_gate_declarations.cuh"

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



#define CHECK_CUDA_CALL(call) \
    { \
        cudaError_t err = call; \
        if(err != cudaSuccess) { \
            printf("CUDA Error: %s (code %d) at %s:%d\n", cudaGetErrorString(err), err, __FILE__, __LINE__); \
            exit(EXIT_FAILURE); \
        } \
    }

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
    __global__ void simulation_kernel_cuda(STAB_CUDA* stab_gpu, Gate* gates_gpu, IdxType n_gates, int seed = 0);
    __global__ void simulation_kernel_cuda2D(STAB_CUDA* stab_gpu, Gate* gates_gpu, IdxType gate_chunk);
    __global__ void pauli_frame_sim(int* m_results, Gate* gates_gpu, IdxType n_gates, int n, int shots, int* x_frame, int* z_frame);

    class STAB_CUDA : public QuantumState
    {
    public:
        //Default identity constructor
        STAB_CUDA(IdxType _n_qubits) : QuantumState(SimType::STAB)
        {
            std::cout << "CUDA constructor" << std::endl;
            /*Initialize basic data*/
            n = _n_qubits;
            rows = 2*n+1;
            cols = n;
            stabCounts = n;
            packed_bits = 32;
            packed_rows = (rows + packed_bits - 1) / packed_bits;
            packed_r_size = packed_rows * sizeof(uint32_t);
            packed_matrix_size = (packed_rows * cols) * sizeof(uint32_t);
            bit_r_size = rows * sizeof(uint32_t);
            bit_matrix_size = rows * cols * sizeof(uint32_t);

            i_proc = 0;
            cudaSafeCall(cudaSetDevice(i_proc));

            //Space out the 2D vectors for x z and r
            reset_state();

            rng.seed(Config::RANDOM_SEED);
            dist = std::uniform_int_distribution<int>(0,1);

            SAFE_ALOC_HOST(singleResultHost, n*sizeof(int));
            memset(singleResultHost, 0, n*sizeof(int));
            SAFE_ALOC_GPU(singleResultGPU, n*sizeof(int));
            /*End initialization*/
            
            //Allocate the packed data to the GPU side using NVSHMEM and a tempRow for row swapping
            SAFE_ALOC_GPU(x_packed_gpu, packed_matrix_size);
            SAFE_ALOC_GPU(z_packed_gpu, packed_matrix_size);
            SAFE_ALOC_GPU(r_packed_gpu, packed_r_size);

            SAFE_ALOC_GPU(x_bit_gpu, bit_matrix_size);
            SAFE_ALOC_GPU(z_bit_gpu, bit_matrix_size);
            SAFE_ALOC_GPU(r_bit_gpu, bit_r_size);


            cudaCheckError();
        }

        ~STAB_CUDA()
        {
            // Release for CPU side
            SAFE_FREE_HOST(totalResults);

            SAFE_FREE_HOST_CUDA(x_packed_cpu);
            SAFE_FREE_HOST_CUDA(z_packed_cpu);
            SAFE_FREE_HOST_CUDA(r_packed_cpu);

            //Allocate the packed data to the GPU side using NVSHMEM
            SAFE_FREE_GPU(x_packed_gpu);
            SAFE_FREE_GPU(z_packed_gpu);
            SAFE_FREE_GPU(r_packed_gpu);
            SAFE_FREE_GPU(x_bit_gpu);
            SAFE_FREE_GPU(z_bit_gpu);
            SAFE_FREE_GPU(r_bit_gpu);
        }

        //Packs down 32 rows in each column and flattens
        void pack_tableau()
        {
            // print_res_state();
            //Allocate memory for x_packed_cpu (2D array), but as a flattened 1D array to ensure memory is contiguous
            x_packed_cpu = new uint32_t[packed_rows * cols]();
            z_packed_cpu = new uint32_t[packed_rows * cols]();
            r_packed_cpu = new uint32_t[packed_rows]();
            int index;
            
            for(int i = 0; i < rows; i++)
            {
                mask = i % packed_bits;
                r_packed_cpu[i/packed_bits] |= r[i] << (mask);

                for(int j = 0; j < cols; j++)
                {
                    //Index conversion
                    index = ((i/packed_bits) * cols) + j;

                    x_packed_cpu[index] |= (x[i][j] << (mask));
                    z_packed_cpu[index] |= (z[i][j] << (mask));
                }
            }

            std::cout << "Pack tableau done!" << std::endl;

            //print_packed_state();
        }

        //Unpacks packed CPU arrays back to bit values
        void unpack_tableau()
        {
            std::cout << "Unpacking tableau!" << std::endl;
            int index;

            for(int i = 0; i < rows; i++)
            {
                mask = 1 << (i % packed_bits);
                r[i] = (r_packed_cpu[i/packed_bits] & mask) ? 1 : 0;
                for(int j = 0; j < cols; j++)
                {
                    index = (i/packed_bits * cols) + j;

                    x[i][j] = (x_packed_cpu[index] & mask) ? 1 : 0;
                    z[i][j] = (z_packed_cpu[index] & mask) ? 1 : 0;
                }
            }
        }
        
        void copy_bits_to_gpu()
        {
            IdxType index;

            uint32_t* temp_r_bit = new uint32_t[rows];       
            uint32_t* temp_x_bit = new uint32_t[rows * cols];    
            uint32_t* temp_z_bit = new uint32_t[rows * cols];   
            
            for (int i = 0; i < rows; i++)
            {
                temp_r_bit[i] = r[i];

                for (int j = 0; j < cols; j++)
                {
                    index = (i * cols) + j;
                    temp_x_bit[index] = x[i][j];
                    temp_z_bit[index] = z[i][j];  
                }
            }

            cudaMemcpy(r_bit_gpu, temp_r_bit, rows * sizeof(uint32_t), cudaMemcpyHostToDevice);
            cudaMemcpy(x_bit_gpu, temp_x_bit, rows * cols * sizeof(uint32_t), cudaMemcpyHostToDevice);
            cudaMemcpy(z_bit_gpu, temp_z_bit, rows * cols * sizeof(uint32_t), cudaMemcpyHostToDevice);

            delete[] temp_r_bit;
            delete[] temp_x_bit;
            delete[] temp_z_bit;
        }

        void copy_bits_from_gpu()
        {
            IdxType index;

            uint32_t* temp_r_bit = new uint32_t[rows];
            uint32_t* temp_x_bit = new uint32_t[rows * cols];
            uint32_t* temp_z_bit = new uint32_t[rows * cols];

            cudaMemcpy(temp_r_bit, r_bit_gpu, rows * sizeof(uint32_t), cudaMemcpyDeviceToHost);
            cudaMemcpy(temp_x_bit, x_bit_gpu, rows * cols * sizeof(uint32_t), cudaMemcpyDeviceToHost);
            cudaMemcpy(temp_z_bit, z_bit_gpu, rows * cols * sizeof(uint32_t), cudaMemcpyDeviceToHost);
            

            for (int i = 0; i < rows; i++)
            {
                r[i] = temp_r_bit[i];

                for (int j = 0; j < cols; j++)
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
        //Prints the tableau including phase
        void print_packed_state()
        {
            int index;
            std::cout << "---------- Packed Tableau: ----------" << std::endl;
            for(int i = 0; i < rows; i++)
            {
                mask = 1 << (i % packed_bits);

                if(((i == (rows/2)) || (i == rows-1)))
                {
                    for(int j = -5; j < (n*2); j++)
                    {
                        std::cout << "-";
                    }
                    std::cout << std::endl;
                }
                for(int j = 0; j < cols; j++)
                {
                    //Index conversion
                    index = (i/packed_bits * cols) + j;
                    std::cout << ((x_packed_cpu[index] & mask) ? 1 : 0);
                }
                std::cout << " | ";
                for(int j = 0; j < cols; j++)
                {
                    //Index conversion
                    index = (i/packed_bits * cols) + j;
                    std::cout << ((z_packed_cpu[index] & mask) ? 1 : 0);
                }
                std::cout << "|" << ((r_packed_cpu[i/packed_bits] & mask) ? 1 : 0) << std::endl;
            }
            std::cout << "---------- ----------" << std::endl;
        }

        // Function to check CUDA errors
        void checkCudaError(cudaError_t error, const char* message) {
            if(error != cudaSuccess) {
                fprintf(stderr, "CUDA Error: %s: %s\n", message, cudaGetErrorString(error));
                exit(EXIT_FAILURE);
            }
        }

        // template <typename T>
        // __global__ void PackTo32Col(const T* X, unsigned* B, int height, int width)
        // {
        //     unsigned Bval, laneid;//fetch lane−id
        //     asm("mov.u32 %0, %%laneid;":"=r"(laneid));
        //     #pragma unroll
        //     for(int i = 0; i < WARP_SIZE; i++)
        //     {
        //         T f0 = X[(blockIdx.x * WARP_SIZE + i) * width + blockIdx.y * WARP_SIZE + laneid];
        //         //rotate anticlockwise for Col packing
        //         unsigned r0 = __brev(__ballot(f0> = 0));
        //         if(laneid == i ) 
        //             Bval = r0;
        //     }
        //     B[blockIdx.y * height + blockIdx.x * WARP_SIZE + laneid] = Bval;
        // }
        // __global__ void PackTo32Row(const int* input, unsigned* output, int width, int height) 
        // {
        //     int laneId = threadIdx.x; // Warp lane ID
        //     int blockRow = blockIdx.x * 32; // Block row start
        //     int blockCol = blockIdx.y * 32; // Block column start
        //     if(blockRow + laneId < height) {
        //         unsigned packed = 0;
        //         for (int i = 0; i < 32; ++i) {
        //             int colIdx = blockCol + i;
        //             if(colIdx < width) {
        //                 int bit = input[(blockRow + laneId) * width + colIdx];
        //                 packed |= (bit > 0 ? 1 : 0) << (31 - i);
        //             }
        //         }
        //         output[(blockRow + laneId) * (width / 32) + blockCol / 32] = packed;
        //     }
        // }

        

        //resets the tableau to a full identity
        void reset_state() override
        {   
            x.resize(rows, std::vector<int>(cols,0)); //first 2n+1 x n block. first n represents destabilizers
                                                       //second n represents stabilizers + 1 extra row
            z.resize(rows, std::vector<int>(cols,0)); //second 2n+1 x n block to form the 2n+1 x 2n sized tableau
            r.resize(rows, 0); //column on the right with 2n+1 rows
            //The 2n+1 th row is scratch space

            //Intialize the identity tableau
            for(int i = 0; i < n; i++)
            {
                x[i][i] = 1;
                z[i+n][i] = 1;
            }
            std::cout << "Reset state done" << std::endl;
            /*End initialization*/
        }
        
        void set_seed(IdxType s) override
        {
            rng.seed(s);
            seed = s;
        }

        //Prints the tableau including phase
        void print_res_state() override
        {
            for(int i = 0; i < rows; i++)
            {
                if(((i == (rows/2)) || (i == rows-1)))
                {
                    for(int j = -5; j < (n*2); j++)
                    {
                        std::cout << "-";
                    }
                    std::cout << std::endl;
                }
                for(int j = 0; j < cols; j++)
                {
                    std::cout << x[i][j];
                }
                std::cout << " | ";
                for(int j = 0; j < cols; j++)
                {
                    std::cout << z[i][j];
                }
                std::cout << "|" << r[i] << std::endl;
            }
        }

        //Prints a 2n x m tableau matrix w/o phase bits
        void print_table(std::vector<std::vector<int>>& M)
        {
            for(int i = 0; i < M[0].size()+1; i++)
            {
                std::cout << "- ";
            }
            std::cout << std::endl;
            for(int i = 0; i < M.size(); i++)
            {
                for(int j = 0; j < M[i].size()/2; j++)
                {
                    std::cout << M[i][j] << " ";
                }
                std::cout << "| ";
                for(int j = M[i].size()/2; j < M[i].size(); j++)
                {
                    std::cout << M[i][j] << " ";
                }
                std::cout << std::endl;
            }
        }

        //Convert a vector of 1's and 0's to an IdxType decimal number
        IdxType *get_results() override
        {
            return totalResults;  //Return a pointer to totalResults
        }

        //Get the stabilizers in the tableau, i.e. all of the Pauli strings that stabilize a certain circuit
        std::vector<std::string> get_stabilizers() override
        {
            int x_val;
            int z_val;
            std::string stabilizers;
            std::vector<std::string> pauliStrings;
            for(int i = (rows>>1); i < rows-1; i++) //rows of stabilizers
            {
                stabilizers.clear(); //clear the temporary stabilizers string
                for(int j = 0; j < cols; j++) //qubits/cols
                {
                    x_val = x[i][j];
                    z_val = z[i][j];
                    assert((x_val < 2) && (z_val < 2));
                    if(x_val)
                    {
                        if(z_val)
                            stabilizers += 'Y';
                        else
                            stabilizers += 'X';
                    }
                    else
                    {
                        if(z_val)
                            stabilizers += 'Z';
                        else
                            stabilizers += 'I';
                    }
                }//For columns(qubits)
                pauliStrings.push_back(stabilizers);
            }//For rows(pauli strings)
            return pauliStrings;
        }

        //Sub-process in measurement gates
        void rowsum(int h, int i)
        {
            int sum = 0;
            for(int j = 0; j < n; j++)
            {
                //Sum every column in the row
                if(x[i][j])
                {
                    if(z[i][j])
                        sum += z[h][j] - x[h][j];
                    else
                        sum += z[h][j] * (2*x[h][j]-1);
                }
                else if(z[i][j])
                    sum += x[h][j] * (1-2*z[h][j]);

                //XOR x's and z's
                x[h][j] = x[i][j] ^ x[h][j];
                z[h][j] = z[i][j] ^ z[h][j];
            }
            sum += 2*r[h] + 2*r[i];

            if(sum % 4 == 0)
                r[h] = 0;
            else
                r[h] = 1;
        } //End rowsum

        //For shot based measurement
        //h and i are row indices
        void tempRowsum(int h, int i, std::vector<std::vector<int>>& temp_x, std::vector<std::vector<int>>& temp_z, std::vector<int>& temp_r)
        {
            int sum = 0;
            for(int j = 0; j < n; j++)
            {
                //Sum every column in the row
                if(temp_x[i][j])
                {
                    if(z[i][j])
                        sum += temp_z[h][j] - temp_x[h][j];
                    else
                        sum += temp_z[h][j] * (2*temp_x[h][j]-1);
                }
                else if(temp_z[i][j])
                    sum += temp_x[h][j] * (1-2*temp_z[h][j]);

                //XOR x's and z's
                temp_x[h][j] = temp_x[i][j] ^ temp_x[h][j];
                temp_z[h][j] = temp_z[i][j] ^ temp_z[h][j];
            }
            sum += 2*temp_r[h] + 2*temp_r[i];

            if(sum % 4 == 0)
                temp_r[h] = 0;
            else
                temp_r[h] = 1;
        } //End rowsum

        int *getSingleResult() override{
            return singleResultHost;
        }

        //Provides a bit/shot measurement output without affecting the original tableau
        IdxType *measure_all(IdxType shots = 2048) override
        {
            //Each index of shotResults is a possible result from a full measurement, ex. 01100101 in an 8 qubit system
            //The integer at that i is the number of times that result occured
            SAFE_ALOC_HOST(totalResults, sizeof(IdxType) * shots);
            memset(totalResults, 0, sizeof(IdxType) * shots);

            int half_rows = rows >> 1;
            
            std::vector<std::vector<int>> temp_x;
            std::vector<std::vector<int>> temp_z;
            std::vector<int> temp_r;
            for(int shot = 0; shot < shots; shot++)
            {
                //Make a copy of the class being measured so many shots can be performed
                temp_x = x;
                temp_z = z;
                temp_r = r;
                for(int a = 0; a < n; a++)
                {  
                    int p = -1;
                    for(int p_index = half_rows; p_index < rows-1; p_index++)
                    {  
                        if(temp_x[p_index][a])
                        {
                            p = p_index;
                            break;
                        }
                    }
                    //A p such that x[p][a] = 1 exists
                    if(p > -1) //Random
                    {
                        for(int i = 0; i < rows-1; i++)
                        {
                            if((x[i][a]) && (i != p))
                            {
                                tempRowsum(i, p, temp_x, temp_z, temp_r);
                            }
                        }
                        temp_x[p-half_rows] = temp_x[p];
                        temp_z[p-half_rows] = temp_z[p];
                        //Change all the columns in row p to be 0
                        for(int i = 0; i < n; i++)
                        {
                            temp_x[p][i] = 0;
                            temp_z[p][i] = 0;                        
                        }

                        int randomBit = dist(rng);
                        
                        if(randomBit)
                        {
                            //std::cout << "Random result of 1" << std::endl;
                            temp_r[p] = 1;
                        }
                        else
                        {
                            //std::cout << "Random result of 0" << std::endl;
                            temp_r[p] = 0;
                        }
                        temp_z[p][a] = 1;

                        totalResults[shot] |= (temp_r[p] << a);
                        // std::cout << "Random measurement at qubit " << a << " value: " << (temp_r[p] << a) << std::endl;
                    }
                    else //Deterministic
                    {
                        //Set the scratch space row to be 0
                        //i is the column indexer in this case
                        for(int i = 0; i < n; i++)
                        {
                            temp_x[rows-1][i] = 0;
                            temp_z[rows-1][i] = 0;
                        }
                        temp_r[rows-1] = 0;

                        //Run rowsum subroutine
                        for(int i = 0; i < half_rows; i++)
                        {
                            if(temp_x[i][a] == 1)
                            {
                                //std::cout << "Perform rowsum at " << i << " + n" << std::endl;
                                tempRowsum(rows-1, i+half_rows, temp_x, temp_z, temp_r);
                            }
                        }

                        // std::cout << "Deterministc measurement at qubit " << a << " value: " << (temp_r[rows-1] << a) << std::endl;
                        totalResults[shot] |=  (temp_r[rows-1] << a);
                    } //End if else
                } //End single shot for all qubits
            }//End shots
            return totalResults;
        }//End measure_all

        //Provides a bit/shot measurement output without affecting the original tableau
        IdxType **measure_all_long(IdxType shots = 2048) override
        {
            //Each index of shotResults is a possible result from a full measurement, ex. 01100101 in an 8 qubit system
            //The integer at that i is the number of times that result occured
            totalResultsLong = new IdxType*[shots];
            for (size_t i = 0; i < shots; i++) {
                totalResultsLong[i] = new IdxType[((n+63)/64)]();  //Zero-initialize
            }

            int half_rows = rows >> 1;
            
            std::vector<std::vector<int>> temp_x;
            std::vector<std::vector<int>> temp_z;
            std::vector<int> temp_r;
            for(int shot = 0; shot < shots; shot++)
            {
                //Make a copy of the class being measured so many shots can be performed
                temp_x = x;
                temp_z = z;
                temp_r = r;
                for(int a = 0; a < n; a++)
                {  
                    int p = -1;
                    for(int p_index = half_rows; p_index < rows-1; p_index++)
                    {  
                        if(temp_x[p_index][a])
                        {
                            p = p_index;
                            break;
                        }
                    }
                    //A p such that x[p][a] = 1 exists
                    if(p > -1) //Random
                    {
                        for(int i = 0; i < rows-1; i++)
                        {
                            if((x[i][a]) && (i != p))
                            {
                                tempRowsum(i, p, temp_x, temp_z, temp_r);
                            }
                        }
                        temp_x[p-half_rows] = temp_x[p];
                        temp_z[p-half_rows] = temp_z[p];
                        //Change all the columns in row p to be 0
                        for(int i = 0; i < n; i++)
                        {
                            temp_x[p][i] = 0;
                            temp_z[p][i] = 0;                        
                        }

                        int randomBit = dist(rng);
                        
                        if(randomBit)
                        {
                            //std::cout << "Random result of 1" << std::endl;
                            temp_r[p] = 1;
                        }
                        else
                        {
                            //std::cout << "Random result of 0" << std::endl;
                            temp_r[p] = 0;
                        }
                        temp_z[p][a] = 1;

                        totalResultsLong[shot][a/64] |=  (temp_r[rows-1] << (a%64));
                        // std::cout << "Random measurement at qubit " << a << " value: " << (temp_r[p] << a) << std::endl;
                    }
                    else //Deterministic
                    {
                        //Set the scratch space row to be 0
                        //i is the column indexer in this case
                        for(int i = 0; i < n; i++)
                        {
                            temp_x[rows-1][i] = 0;
                            temp_z[rows-1][i] = 0;
                        }
                        temp_r[rows-1] = 0;

                        //Run rowsum subroutine
                        for(int i = 0; i < half_rows; i++)
                        {
                            if(temp_x[i][a] == 1)
                            {
                                //std::cout << "Perform rowsum at " << i << " + n" << std::endl;
                                tempRowsum(rows-1, i+half_rows, temp_x, temp_z, temp_r);
                            }
                        }

                        // std::cout << "Deterministc measurement at qubit " << a << " value: " << (temp_r[rows-1] << a) << std::endl;
                        totalResultsLong[shot][a/64] |=  (temp_r[rows-1] << (a%64));
                    } //End if else
                } //End single shot for all qubits
            }//End shots
            return totalResultsLong;
        }//End measure_all

        //Simulate the gates from a circuit in the tableau
        void sim2D(std::shared_ptr<Circuit> circuit2D, std::vector<int> gate_chunks, double &sim_time) override
        {
            STAB_CUDA *stab_gpu;
            SAFE_ALOC_GPU(stab_gpu, sizeof(STAB_CUDA));
            // Copy the simulator instance to GPU
            cudaSafeCall(cudaMemcpy(stab_gpu, this,
                                    sizeof(STAB_CUDA), cudaMemcpyHostToDevice));

            gpu_timer sim_timer;

            //Pack tableau and copy to GPU
            pack_tableau();
            cudaSafeCall(cudaMemcpy(x_packed_gpu, x_packed_cpu, packed_matrix_size,
                                    cudaMemcpyHostToDevice));
            cudaSafeCall(cudaMemcpy(z_packed_gpu, z_packed_cpu, packed_matrix_size,
                                    cudaMemcpyHostToDevice));
            cudaSafeCall(cudaMemcpy(r_packed_gpu, r_packed_cpu, packed_r_size, cudaMemcpyHostToDevice));

            std::cout << "Data copied" << std::endl;

            int minGridSize, blockSize;
            cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, simulation_kernel_cuda2D, 0, 0);
            int threadsPerBlockX = 32;
            int threadsPerBlockY = blockSize / threadsPerBlockX; 

            if(threadsPerBlockY > 16) threadsPerBlockY = 16;

            dim3 threadsPerBlock(threadsPerBlockX, threadsPerBlockY);
            dim3 blocksPerGrid((rows - 1 + threadsPerBlockX - 1) / threadsPerBlockX,
                            (n + threadsPerBlockY - 1) / threadsPerBlockY);

            std::cout << "Blocks calculated" << std::endl;

           

            std::vector<Gate> gates2D = circuit2D->get_gates();
            int num_gates = gates2D.size();
            copy_gates_to_gpu(gates2D);
            std::cout << "2D circuit parsed" << std::endl;

            /*Simulate*/
            if(Config::PRINT_SIM_TRACE)
            {
                printf("STABSim_gpu is running! Using %lld qubits.\n", n);
            }

            //Call the kernel for each set of gates
            IdxType gate_index_sum = 0;
            int gate_chunk_size = 0;
            sim_timer.start_timer();
            for(int i = 0; i < gate_chunks.size(); i++)
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
            //Copy data to the CPU side and unpack
            cudaSafeCall(cudaMemcpy(x_packed_cpu, x_packed_gpu, packed_matrix_size,
                                    cudaMemcpyDeviceToHost));
            cudaSafeCall(cudaMemcpy(z_packed_cpu, z_packed_gpu, packed_matrix_size,
                                    cudaMemcpyDeviceToHost));
            cudaSafeCall(cudaMemcpy(r_packed_cpu, r_packed_gpu, packed_r_size, cudaMemcpyDeviceToHost));
            unpack_tableau();

            if(Config::PRINT_SIM_TRACE)
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

        //Simulate the gates from a circuit in the tableau
        void sim(std::shared_ptr<Circuit> circuit, double &sim_time) override
        {
            STAB_CUDA *stab_gpu;
            SAFE_ALOC_GPU(stab_gpu, sizeof(STAB_CUDA));
            //Copy the simulator instance to GPU
            cudaSafeCall(cudaMemcpy(stab_gpu, this,
                                    sizeof(STAB_CUDA), cudaMemcpyHostToDevice));
            
            gpu_timer sim_timer;

            //Copy bit data to GPU
            copy_bits_to_gpu();
            singleResultHost = new int[n]();
            cudaMemcpy(singleResultHost, singleResultGPU, n*sizeof(int), cudaMemcpyHostToDevice);

            std::cout << "Data copied" << std::endl;

            std::vector<Gate> gates = circuit->get_gates();
            //Copy gates to the gpu side
            copy_gates_to_gpu(gates);
            IdxType n_gates = gates.size();

            cudaDeviceProp prop;
            cudaGetDeviceProperties(&prop, 0);
            printf("Max grid size: (%d, %d, %d)\n", prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
            printf("Max threads per block: %d\n", prop.maxThreadsPerBlock);
            printf("Device Name: %s\n", prop.name);
            printf("Max Blocks per Multiprocessor: %d\n", prop.maxThreadsPerMultiProcessor);
            printf("Number of SMs: %d\n", prop.multiProcessorCount);

            int maxBlocks;
            cudaOccupancyMaxActiveBlocksPerMultiprocessor(
                &maxBlocks, /* out: max active blocks */
                (void*)simulation_kernel_cuda, /* kernel */
                1024, /* threads per block */
                3*sizeof(int) /* shared memory per block */
            );

            printf("Max active blocks per SM: %d\n", maxBlocks);
            maxBlocks = prop.multiProcessorCount;
            printf("Total max cooperative blocks: %d\n", maxBlocks);
            cudaFuncAttributes attr;
            cudaFuncGetAttributes(&attr, simulation_kernel_cuda);
            printf("Registers per thread: %d\n", attr.numRegs);


            int minGridSize, blockSize;
            cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, simulation_kernel_cuda, 3* sizeof(int), 0);
            printf("Max block size: %d\n", blockSize);
            int threadsPerBlockX = 1024;

            int threadsPerBlock = (threadsPerBlockX);
            int blocksPerGrid = ((rows-1 + threadsPerBlockX - 1) / threadsPerBlockX);
            // int blocksPerGrid = 132;


            if(blocksPerGrid > maxBlocks)
            {
                fprintf(stderr, "Error: Attempting to launch more blocks than can fit on the GPU. Cooperative kernel requires 1 block per SM.\n", maxBlocks);
                exit(EXIT_FAILURE);  // Exit the program if more than 1 block per SM is attempted
            }

            std::cout << "Blocks calculated" << std::endl;
            std::cout << "X blocks = "  << blocksPerGrid << std::endl;

            /*Simulate*/
            if(Config::PRINT_SIM_TRACE)
            {
                printf("STABSim_gpu is running! Using %lld qubits.\n", n);
            }

            void *args[] = {&stab_gpu, &gates_gpu, &n_gates, &seed};

            /*Simulate*/
            std::cout << "\n -------------------- \n Simulation starting! \n -------------------- \n" << std::endl;
            sim_timer.start_timer();

            //Launch with cooperative kernel
            cudaLaunchCooperativeKernel((void*)simulation_kernel_cuda, blocksPerGrid, threadsPerBlock, args, 3* sizeof(int));
            // simulation_kernel_cuda_bitwise<<<blocksPerGrid, threadsPerBlock>>>(stab_gpu, gates_gpu, n_gates);

            sim_timer.stop_timer();
            // CHECK_CUDA_CALL(cudaPeekAtLastError());
            cudaSafeCall(cudaDeviceSynchronize());
            /*End simulate*/
            sim_time += sim_timer.measure();

            cudaCheckError();
            //Copy data to the CPU side and unpack
            copy_bits_from_gpu();
            cudaMemcpy(singleResultHost, singleResultGPU, n*sizeof(int), cudaMemcpyDeviceToHost);


            if(Config::PRINT_SIM_TRACE)
            {
                printf("\n============== STAB-Sim ===============\n");
                printf("n_qubits:%lld, n_gates:%lld, ncpus:%d, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms\n",
                       n, n_gates, 1, sim_time, 0.,
                       sim_time);
                printf("=====================================\n\n\n");
            }

            SAFE_FREE_GPU(gates_gpu);

            //=========================================
        }

        // void sample()
        // {
        //     int* d_x_frame;
        //     int* d_z_frame;
        //     cudaMalloc(&d_x_frame, n * sizeof(int));
        //     cudaMalloc(&d_z_frame, n * sizeof(int));

        //     void* args[] = {&m_results, &gates_gpu, &n_gates, &n, &d_x_frame, &d_z_frame};
        //     cudaLaunchCooperativeKernel((void*)pauli_frame_sim, gridDim, blockDim, args);
        // }

        IdxType measure(IdxType qubit) override
        {
            throw std::logic_error("measure Not implemented (STAB_GPU)");
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

    public:
        IdxType n;
        IdxType n_gates;
        IdxType stabCounts;
        IdxType mask;
        IdxType rows;
        IdxType packed_rows;
        IdxType cols;
        int packed_bits;
        IdxType packed_r_size;
        IdxType packed_matrix_size;
        IdxType bit_matrix_size;
        IdxType bit_r_size;
        std::vector<std::vector<int>> x;
        std::vector<std::vector<int>> z;
        std::vector<int> r;

        std::vector<std::vector<Gate>> layered_gates;
        IdxType num_layers;
        //CPU Arrays
        uint32_t* x_packed_cpu = nullptr;
        uint32_t* z_packed_cpu = nullptr;
        uint32_t* r_packed_cpu = nullptr;
        //GPU Arrays
        uint32_t* x_packed_gpu = nullptr;
        uint32_t* z_packed_gpu = nullptr;
        uint32_t* r_packed_gpu = nullptr;
        uint32_t* x_bit_gpu = nullptr;
        uint32_t* z_bit_gpu = nullptr;
        uint32_t* r_bit_gpu = nullptr;

        uint32_t* row_sum_gpu = nullptr;
        // int* p_shared;

        int* singleResultHost;
        int* singleResultGPU;
        int seed;

        IdxType* totalResults = nullptr;
        IdxType** totalResultsLong = nullptr;

        //Random
        std::mt19937 rng;
        std::uniform_int_distribution<int> dist;   

        //GPU-side gate instance
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
    }; //End tableau class

    __device__ int p_shared;
    __device__ uint32_t row_sum;
    __global__ void simulation_kernel_cuda(STAB_CUDA* stab_gpu, Gate* gates_gpu, IdxType n_gates, int seed)
    {
        cg::grid_group grid = cg::this_grid();
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        IdxType rows = stab_gpu->rows;
        IdxType cols = stab_gpu->cols;
        if(i >= rows-1) return;

        int n_qubits = stab_gpu->n;
        uint32_t* x_arr = stab_gpu->x_bit_gpu;
        uint32_t* z_arr = stab_gpu->z_bit_gpu;
        uint32_t* r_arr = stab_gpu->r_bit_gpu;
        uint32_t x, z;
        OP op_name;
        IdxType index;
        IdxType row_col_index;
        int a;
        int p;
        uint32_t local_sum;

        //Initialize shared memory (per block)
        __shared__ int local_p_shared;
        __shared__ int half_row;
        __shared__ int scratch_row;
        if(threadIdx.x == 0)
        {
            half_row = rows/2;
            scratch_row = rows-1;
            local_p_shared = rows;
        }
        __shared__ curandState state;
        if(i == rows/2)
        {
            curand_init(1234, i, 0, &state);
        }
        
        // printf("Starting for loop!! \n\n");
        for (int k = 0; k < n_gates; k++) 
        {
            op_name = gates_gpu[k].op_name;
            a = gates_gpu[k].qubit;
            index = i * cols + a;
            
            grid.sync();

            

            switch (op_name) 
            {
                case OP::H:
                    x = x_arr[index];
                    z = z_arr[index];
                    //Phase
                    r_arr[i] ^= (x & z);

                    //Entry -- swap x and z bits
                    x_arr[index] = z;
                    z_arr[index] = x;
                    break;

                case OP::S:
                    x = x_arr[index];
                    z = z_arr[index];

                    //Phase
                    r_arr[i] ^= (x & z);

                    //Entry
                    z_arr[index] = z ^ x;
                    // if(i == 0) printf("S\n\n");
                    break;

                case OP::SDG:
                    x = x_arr[index];
                    z = z_arr[index];

                    //Phase
                    r_arr[i] ^= (x ^ (x & z));

                    //Entry
                    z_arr[index] = z ^ x;
                    break;

                case OP::CX:
                    row_col_index = i * cols + gates_gpu[k].ctrl;

                    x = x_arr[index];
                    z = z_arr[index];

                    uint32_t x_ctrl = x_arr[row_col_index];
                    uint32_t z_ctrl = z_arr[row_col_index];

                    //Phase
                    r_arr[i] ^= ((x_ctrl & z) & (x^z_ctrl^1));

                    //Entry
                    x_arr[index] = x ^ x_ctrl;
                    z_arr[row_col_index] = z ^ z_ctrl;

                    break;

                case OP::M:
                {
                    //Reset p_shared in the first thread of the first block
                    if(threadIdx.x == 0)
                    {
                        //Initialize shared memory
                        local_p_shared = rows;
                        if(blockIdx.x == 0) 
                        {
                            // printf("M\n\n");
                            //Initialize global memory
                            p_shared = rows;
                        }
                    }
                    

                    //Initialize thread memory
                    p = rows;
                    if((i >= (half_row)) && (x_arr[index]))
                        p = i;
                        
                    
                    //Reduce within a block
                    atomicMin(&local_p_shared, p);
                    // if(i == 0 && j == 0)
                    //     printf("Shared reduced\n");
                    
                    //Reduce across all blocks (all blocks need to be caught up)
                    grid.sync();
                    // __syncthreads();
                    if(threadIdx.x == 0)
                        atomicMin(&p_shared, local_p_shared);
                    // if(i == 0 && j == 0)
                    //     printf("Global reduced\n");
                    // __syncthreads();
                    grid.sync();

                    //Debugging output
                    // if(i == 0)
                    //     printf("p_shared = %d\n", p_shared);

                    //If no p among the stabilizers is found, the measurement will be random
                    if(p_shared != rows)
                    {
                        //Rowsum for all rows < 2n
                        for(int k = 0; k < scratch_row; k++)
                        {
                            //printf("Entering row/2 for loop\n");
                            if(x_arr[k * cols + a] && k != p_shared)
                            {
                                //Initialize the sums from all rows we're interested in to 0
                                //Using i as columns
                                if(i == 0)
                                    row_sum = 0;
                                grid.sync();
                                if(i < half_row)
                                {
                                    index = (k * cols) + i;
                                    row_col_index = (p_shared * cols) + i;
                                    x = x_arr[row_col_index];
                                    z = z_arr[row_col_index];

                                    if(x && z) 
                                    {
                                        atomicAdd(&row_sum, z_arr[index] - x_arr[index]);
                                        x_arr[index] ^= x;
                                        z_arr[index] ^= z;
                                    }
                                    if(x && !z) 
                                    {
                                        atomicAdd(&row_sum, z_arr[index] * (2 * x_arr[index] - 1));
                                        x_arr[index] ^= x;
                                    }
                                    if(!x && z) 
                                    {
                                        atomicAdd(&row_sum, x_arr[index] * (1 - 2 * z_arr[index]));
                                        z_arr[index] ^= z;
                                    }
                                }

                                grid.sync();
                                
                                if(i == 0)
                                {
                                    //Add the stabilizer r value to the corresponding row sum
                                    row_sum += 2 * r_arr[k] + 2 * r_arr[p_shared];

                                    if(row_sum % 4)
                                        r_arr[k] = 1;
                                    else
                                        r_arr[k] = 0;
                                    // printf("row_sum[%d] = %d\n", i, row_sum);
                                }
                                //End Rowsum
                            }
                        }

                        grid.sync();

                        if(i < cols)
                        {
                            //Set every column of the destab of row p to the stab of row p
                            index = (p_shared * cols) + i;
                            row_col_index = ((p_shared-(half_row)) * cols) + i;
                            x_arr[row_col_index] = x_arr[index];    
                            z_arr[row_col_index] = z_arr[index];
                            x_arr[index] = 0;
                            z_arr[index] = 0; 
                        }
                        if(i == cols)
                        {
                            //Generate a random bit (0 or 1)
                            r_arr[p_shared] = curand(&state) & 1;
                            // printf("Random measurement at qubit %d value: %d\n", a, r_arr[p_shared]);

                            //Update z to reflect a z measurement
                            z_arr[(p_shared * cols) + a] = 1; 

                            stab_gpu->singleResultGPU[a] = r_arr[p_shared];
                            // printf("Random done %d\n ", stab_gpu->singleResultGPU[a]);
                        }
                    }
                    else
                    {
                        // if(i == 0)
                        // {
                        //     printf("Deterministic\n");
                        // }

                        //Set the scratch row to 0
                        if(i < cols) //using i as cols
                        {
                            index = (scratch_row * cols) + i;
                            x_arr[index] = 0;
                            z_arr[index] = 0;
                        }
                        if(i == cols)
                        {
                            r_arr[scratch_row] = 0;
                        }

                        //Wait for the scratch row to be reset before proceeding
                        grid.sync();
                    
                        //Rowsum for all rows < n
                        for(int k = 0; k < half_row; k++)
                        {
                            grid.sync();
                            //printf("Entering row/2 for loop\n");
                            if(x_arr[k * cols + a])
                            {
                                //printf("x_arr[i] = %d \n", x_arr[(i * cols) + a]);
                                //Start Rowsum

                                // Over every column (j)
                                // printf("rows = %d \n", rows);
                                // printf("i = %d \n", i);
                                // printf("j = %d \n", j);
                                // printf("row_sum[i] = %d \n", row_sum[i]);

                                //Initialize the sums from all rows we're interested in to 0
                                //Using i as columns
                                if(threadIdx.x == 0)
                                {
                                    row_sum = 0;
                                }
                                if(i < cols)
                                {
                                    row_col_index = (scratch_row * cols) + i;
                                    index = ((k+(n_qubits)) * cols) + i;
                                    local_sum = 0;            
                                    x = x_arr[index];
                                    z = z_arr[index];

                                    if(x && z) 
                                    {
                                        local_sum = z_arr[row_col_index] - x_arr[row_col_index];
                                        x_arr[row_col_index] ^= x;
                                        z_arr[row_col_index] ^= z;

                                        // printf("Col_val in %d = %d \n", i, col_val);
                                    }
                                    if(x && !z) 
                                    {
                                        local_sum = z_arr[row_col_index] * (2 * x_arr[row_col_index] - 1);
                                        x_arr[row_col_index] ^= x;
                                        // printf("Col_val in %d = %d \n", i, col_val);
                                    }
                                    if(!x && z) 
                                    {
                                        local_sum = x_arr[row_col_index] * (1 - 2 * z_arr[row_col_index]);
                                        z_arr[row_col_index] ^= z;
                                        // printf("Col_val in %d = %d \n", i, col_val);
                                    }
                                
                                    //Add all of the columns together for a given row
                                    atomicAdd(&row_sum, local_sum);
                                }
                                grid.sync();
                                if(i == 0)
                                {
                                    //Add the stabilizer r value to the corresponding row sum
                                    row_sum += 2 * r_arr[k+(half_row)] + 2 * r_arr[scratch_row];

                                    if(row_sum % 4)
                                        r_arr[scratch_row] = 1;
                                    else
                                        r_arr[scratch_row] = 0;
                                    // printf("row_sum[%d] = %d\n", i, row_sum);
                                }
                            }
                        } //End Rowsum

                        if(i == 0)
                        {
                            stab_gpu->singleResultGPU[a] = r_arr[scratch_row];
                            // printf("Deterministic measurement at qubit %d value: %d\n", 
                            //     a, r_arr[scratch_row]);
                        }
                    }
                    break;
                }

                case OP::RESET:
                {
                    //Reset p_shared in the first thread of the first block
                    if(threadIdx.x == 0)
                    {
                        //Initialize shared memory
                        local_p_shared = rows;
                        if(blockIdx.x == 0) 
                        {
                            // printf("M\n\n");
                            //Initialize global memory
                            p_shared = rows;
                        }
                    }
                    

                    //Initialize thread memory
                    p = rows;
                    if((i >= (half_row)) && (x_arr[index]))
                        p = i;
                        
                    
                    //Reduce within a block
                    atomicMin(&local_p_shared, p);
                    // if(i == 0 && j == 0)
                    //     printf("Shared reduced\n");
                    
                    //Reduce across all blocks (all blocks need to be caught up)
                    grid.sync();
                    // __syncthreads();
                    if(threadIdx.x == 0)
                        atomicMin(&p_shared, local_p_shared);
                    // if(i == 0 && j == 0)
                    //     printf("Global reduced\n");
                    // __syncthreads();
                    grid.sync();

                    //Debugging output
                    // if(i == 0)
                    //     printf("p_shared = %d\n", p_shared);

                    //If no p among the stabilizers is found, the measurement will be random
                    if(p_shared != rows)
                    {
                        if(i < cols)
                        {
                            //Set every column of the destab of row p to the stab of row p
                            index = (p_shared * cols) + i;
                            row_col_index = ((p_shared-(half_row)) * cols) + i;
                            x_arr[row_col_index] = x_arr[index];    
                            z_arr[row_col_index] = z_arr[index];
                            x_arr[index] = 0;
                            z_arr[index] = 0; 
                        }
                        
                        grid.sync();

                        //Rowsum for all rows < 2n
                        for(int k = 0; k < scratch_row; k++)
                        {
                            //printf("Entering row/2 for loop\n");
                            if(x_arr[k * cols + a] && k != p_shared)
                            {
                                //Initialize the sums from all rows we're interested in to 0
                                //Using i as columns
                                if(i == 0)
                                    row_sum = 0;
                                grid.sync();
                                if(i < half_row)
                                {
                                    index = (k * cols) + i;
                                    row_col_index = (p_shared * cols) + i;
                                    x = x_arr[row_col_index];
                                    z = z_arr[row_col_index];

                                    if(x && z) 
                                    {
                                        atomicAdd(&row_sum, z_arr[index] - x_arr[index]);
                                        x_arr[index] ^= x;
                                        z_arr[index] ^= z;
                                    }
                                    if(x && !z) 
                                    {
                                        atomicAdd(&row_sum, z_arr[index] * (2 * x_arr[index] - 1));
                                        x_arr[index] ^= x;
                                    }
                                    if(!x && z) 
                                    {
                                        atomicAdd(&row_sum, x_arr[index] * (1 - 2 * z_arr[index]));
                                        z_arr[index] ^= z;
                                    }
                                }

                                grid.sync();
                                
                                if(i == 0)
                                {
                                    //Add the stabilizer r value to the corresponding row sum
                                    row_sum += 2 * r_arr[k] + 2 * r_arr[p_shared];

                                    if(row_sum % 4)
                                        r_arr[k] = 1;
                                    else
                                        r_arr[k] = 0;
                                }
                                //End Rowsum
                            }
                        }
                    }
                    else
                    {
                        //Set the first stabilizer to Z
                        if(i < cols) //using i as cols
                        {
                            x_arr[i] = 0;
                            z_arr[i] = 0;
                            index = (n_qubits * cols) + i;
                            x_arr[index] = 0;
                            z_arr[index] = 0;
                        }
                        if(i == cols)
                        {
                            z_arr[(n_qubits * cols) + a] = 1;
                            x_arr[a] = 1;
                            r_arr[scratch_row] = 0;
                        }
                    }
                    break;
                }

                default:
                    printf("Non-Clifford or unrecognized gate: %d\n", op_name);
            }
        }
        // printf("Kernel is done!\n");
    }//end kernel

    __global__ void simulation_kernel_cuda2D(STAB_CUDA* stab_gpu, Gate* gates_gpu, IdxType gate_chunk) 
    {
        int row = blockIdx.x * blockDim.x + threadIdx.x;  //Index for stabilizers (rows)
        int col = blockIdx.y * blockDim.y + threadIdx.y;  //Index for gates (columns)

        if(row >= stab_gpu->packed_rows) return;  //Check out-of-bounds for qubits

        if(col >= gate_chunk) return;  //Check out-of-bounds for gate

        // printf("Inside 2D kernel %d, gate chunk %lld gate qubit %lld \n", col, gate_chunk, gates_gpu[col].qubit);
        // printf("Gates gpu size %lld \n", (sizeof(gates_gpu)));

        int target = gates_gpu[col].qubit; //Qubit target
        OP op_name = gates_gpu[col].op_name;  //Operation to perform
        uint32_t* x_arr = stab_gpu->x_packed_gpu;
        uint32_t* z_arr = stab_gpu->z_packed_gpu;

        //Calculate the index for this qubit in the packed arrays
        IdxType index = row * stab_gpu->cols + target;

        //Perform operations for all possible gates, but mask non-relevant ones
        //Start with the common operation - calculate phase and entry for all gates
        uint32_t x = x_arr[index];
        uint32_t z = z_arr[index];

        //Gate operations
        if(op_name == OP::H) {
            // printf("Inside h gate \n");
            // H Gate: Swap x and z
            stab_gpu->r_packed_gpu[row] ^= (x & z);
            //Entry -- swap x and z bits
            x_arr[index] = z;
            z_arr[index] = x;
            return;
        } 
        if(op_name == OP::S) {
            //S Gate: Entry (z_arr[index] ^= x_arr[index])
            stab_gpu->r_packed_gpu[row] ^= (x & z);
            z_arr[index] = x ^ z;
            return;
        }
        if(op_name == OP::SDG) {
            // SDG Gate: Phase (x_arr[index] ^ (x_arr[index] & z_arr[index]))
            stab_gpu->r_packed_gpu[row] ^= (x ^ (x & z));  // SDG Phase operation
            z_arr[index] = x ^ z;  // Entry (same as S gate)
            return;
        }
        if(op_name == OP::CX) {
            int ctrl_index = row * stab_gpu->cols + gates_gpu[col].ctrl;

            uint32_t x_ctrl = x_arr[ctrl_index];
            uint32_t z_ctrl = z_arr[ctrl_index];

            //Phase
            stab_gpu->r_packed_gpu[row] ^= ((x_ctrl & z) & (x^z_ctrl^1));

            //Entry
            x_arr[index] = x ^ x_ctrl;
            z_arr[ctrl_index] = z ^ z_ctrl;
            return;
        }
        if(op_name == OP::RESET) {
            x_arr[index] = 0;
            z_arr[index] = 0;
            if(row == (target/32))
            {   
                x_arr[index] ^= (1U << (target%32));
                stab_gpu->r_packed_gpu[row] &= ~(1U << (target % 32));
            }
            if(row == ((2*target)/32))
                z_arr[index] ^= (1U << ((2*target)%32));


            return;
        }
    }

    // __global__ void pauli_frame_sim(int* m_results, Gate* gates_gpu, IdxType n_gates, int n, int shots)
    // {
    //     cg::grid_group grid = cg::this_grid();
    //     int i = blockIdx.x * blockDim.x + threadIdx.x;
    //     IdxType rows = stab_gpu->rows;
    //     IdxType cols = stab_gpu->cols;
    //     if(i >= shots) return;
    //     int x_frame[n](0);
    //     int z_frame[n](0);
    //     int r = 0;

    //     OP op_name;

    //     curandState state;
    //     curand_init(1234, i, 0, &state);
        
    //     for (int k = 0; k < n_gates; k++) 
    //     {
    //         op_name = gates_gpu[k].op_name;
    //         a = gates_gpu[k].qubit;
            
    //         grid.sync();

    //         switch (op_name) 
    //         {
    //             case OP::H:
    //                 int p = (x_frame[a] << 1) | z_frame[a];

    //                 r ^= (p == 3);
    //                 x_frame[a] = (p == 1);
    //                 z_frame[a] = (p == 2); 

    //                 break;

    //             case OP::S:
    //                 int p = (x_frame[a] << 1) | z_frame[a];

    //                 r ^= (p == 3);
    //                 z_frame[a] = !(p == 3);

                    
    //                 break;

    //             case OP::CX:
    //                 int cntrl = gates_gpu[k].cntrl;
    //                 z_frame[cntrl] ^= z_frame[a];
    //                 x_frame[a] ^= x_frame[cntrl];

    //                 break;

    //             case OP::M:
    //         }
    //     }
    // }

} //namespace NWQSim

//#endif