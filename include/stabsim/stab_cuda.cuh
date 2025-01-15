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
#include <iostream>
#include <cuda.h>
#include <memory>

#ifdef FP64_TC_AVAILABLE
#include <mma.h>
#endif

namespace NWQSim
{
    using namespace cooperative_groups;
    using namespace std;

    // Simulation kernel runtime
    class STAB_CUDA;
    __global__ void stab_simulation_kernel_cuda(Gate* gates, IdxType g, uint32_t** x, uint32_t** z, uint32_t* r, IdxType rows, IdxType cols);

    class STAB_CUDA : public QuantumState
    {
    public:
        //Default identity constructor
        STAB_CUDA(IdxType _n_qubits) : QuantumState(SimType::STAB)
        {
            /*Initialize basic bit data*/
            n = _n_qubits;
            rows = 2*n+1;
            cols = n;
            stabCounts = n;
            packed_bits = 32;
            packed_rows = (rows/packed_bits)+1;
            packed_r_size = packed_rows * sizeof(uint32_t);
            packed_matrix_size = (packed_rows * cols) * sizeof(uint32_t);

            //Space out the vectors
            x.resize(rows, std::vector<int>(cols,0)); //first 2n+1 x n block. first n represents destabilizers
                                                       //second n represents stabilizers + 1 extra row
            z.resize(rows, std::vector<int>(cols,0)); //second 2n+1 x n block to form the 2n+1 x 2n sized tableau
            r.resize(rows, 0); //column on the right with 2n+1 rows
            //The 2n+1 th row is scratch space

            rng.seed(Config::RANDOM_SEED);
            dist = std::uniform_int_distribution<int>(0,1);
            SAFE_ALOC_HOST(totalResults, sizeof(IdxType));
            memset(totalResults, 0, sizeof(IdxType));
            /*End initialization*/

            //Initialize packed matrices
            SAFE_ALOC_HOST_CUDA(x_packed_cpu, packed_matrix_size);
            SAFE_ALOC_HOST_CUDA(z_packed_cpu, packed_matrix_size);
            SAFE_ALOC_HOST_CUDA(r_packed_cpu, packed_r_size);
            memset(x_packed_cpu, 0, packed_matrix_size);
            memset(z_packed_cpu, 0, packed_matrix_size);
            memset(r_packed_cpu, 0, packed_r_size);

            //Transpose and initialize identity tableau in the packed data
            for(int i = 0; i < n; i++)
            {
                x_packed_cpu[i/packed_bits][i] |= (1 << (i % packed_bits));
                z_packed_cpu[(i/packed_bits)+(n/packed_bits)][i] |= (1 << (i % packed_bits));
            }

            //Allocate the packed data to the GPU side using NVSHMEM
            SAFE_ALOC_GPU(x_packed_gpu, packed_matrix_size);
            SAFE_ALOC_GPU(z_packed_gpu, packed_matrix_size);
            SAFE_ALOC_GPU(r_packed_gpu, packed_r_size);

            cudaCheckError();

            //Copy data to the GPU side
            cudaSafeCall(cudaMemcpy(x_packed_gpu, x_packed_cpu, packed_matrix_size,
                                    cudaMemcpyHostToDevice));
            cudaSafeCall(cudaMemcpy(z_packed_gpu, z_packed_cpu, packed_matrix_size,
                                    cudaMemcpyHostToDevice));
            cudaSafeCall(cudaMemcpy(r_packed_gpu, r_packed_cpu, packed_r_size, cudaMemcpyHostToDevice));
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
        }

        void pack_to_gpu()
        {
            for(int i = 0; i < rows; i++)
            {
                mask = i % packed_bits;
                r_packed_cpu[i/packed_bits] |= r[i] << (mask);
                for(int j = 0; j < cols; j++)
                {
                    x_packed_cpu[i/packed_bits][j] |= x[i][j] << (mask);
                    x_packed_cpu[i/packed_bits][j] |= x[i][j] << (mask);
                }
            }
            //Copy data to the GPU side
            cudaSafeCall(cudaMemcpy(x_packed_gpu, x_packed_cpu, packed_matrix_size,
                                    cudaMemcpyHostToDevice));
            cudaSafeCall(cudaMemcpy(z_packed_gpu, z_packed_cpu, packed_matrix_size,
                                    cudaMemcpyHostToDevice));
            cudaSafeCall(cudaMemcpy(r_packed_gpu, r_packed_cpu, packed_r_size, cudaMemcpyHostToDevice));
        }

        //Unpacks packed CPU arrays back to bit values
        void unpack_to_cpu()
        {
            //Copy data to the CPU side
            cudaSafeCall(cudaMemcpy(x_packed_cpu, x_packed_gpu, packed_matrix_size,
                                    cudaMemcpyDeviceToHost));
            cudaSafeCall(cudaMemcpy(z_packed_cpu, z_packed_gpu, packed_matrix_size,
                                    cudaMemcpyDeviceToHost));
            cudaSafeCall(cudaMemcpy(r_packed_cpu, r_packed_gpu, packed_r_size, cudaMemcpyDeviceToHost));
            for(int i = 0; i < rows; i++)
            {
                mask = i % packed_bits;
                r[i] = (r_packed_gpu[i/packed_bits] >> mask) & 1;
                for(int j = 0; j < cols; j++)
                {
                    x[i][j] = (x_packed_gpu[i/packed_bits][j] >> mask) & 1;
                    z[i][j] = (z_packed_gpu[i/packed_bits][j] >> mask) & 1;
                }
            }
        }

        // Function to check CUDA errors
        void checkCudaError(cudaError_t error, const char* message) {
            if (error != cudaSuccess) {
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
        //         if (laneid == i ) 
        //             Bval = r0;
        //     }
        //     B[blockIdx.y * height + blockIdx.x * WARP_SIZE + laneid] = Bval;
        // }
        // __global__ void PackTo32Row(const int* input, unsigned* output, int width, int height) 
        // {
        //     int laneId = threadIdx.x; // Warp lane ID
        //     int blockRow = blockIdx.x * 32; // Block row start
        //     int blockCol = blockIdx.y * 32; // Block column start
        //     if (blockRow + laneId < height) {
        //         unsigned packed = 0;
        //         for (int i = 0; i < 32; ++i) {
        //             int colIdx = blockCol + i;
        //             if (colIdx < width) {
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
            /*End initialization*/
        }
        

        void set_seed(IdxType s) override
        {
            rng.seed(s);
        }

        //Prints the tableau including phase
        void print_res_state() override
        {
            for(int i = 0; i < rows; i++)
            {
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

        //Provides a bit/shot measurement output without affecting the original tableau
        IdxType *measure_all(IdxType shots = 2048) override
        {
            //Each index of shotResults is a possible result from a full measurement, ex. 01100101 in an 8 qubit system
            //The integer at that position is the number of times that result occured
            SAFE_ALOC_HOST(totalResults, sizeof(IdxType) * shots);
            memset(totalResults, 0, sizeof(IdxType) * shots);

            int half_rows = rows >> 1;
            
            std::vector<std::vector<int>> temp_x;
            std::vector<std::vector<int>> temp_z;
            std::vector<int> temp_r;
            for(int i = 0; i < shots; i++)
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

                        totalResults[i] |= (temp_r[p] << a);
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
                        totalResults[i] |=  (temp_r[rows-1] << a);
                    } //End if else
                } //End single shot for all qubits
            }//End shots
            return totalResults;
        }//End measure_all

        //Simulate the gates from a circuit in the tableau
        void sim(std::shared_ptr<NWQSim::Circuit> circuit, double &sim_time) override
        {
            std::vector<Gate> gates = circuit->get_gates();
            n_gates = gates.size();
            assert(circuit->num_qubits() == n);

            //Copy gates to the gpu side
            Gate* gpu_gates;
            SAFE_ALOC_GPU(gates, n_gates * sizeof(Gate));
            cudaSafeCall(cudaMemcpy(gpu_gates, gates.data(), n_gates * sizeof(Gate), cudaMemcpyHostToDevice));

            //Start a timer
            gpu_timer sim_timer;
            sim_timer.start_timer();

            //Calculate threads and blocks
            int threadsPerBlock = 256;
            int blocksPerGrid = (rows + threadsPerBlock - 1) / threadsPerBlock;
            
            //Simulate
            stab_simulation_kernel_cuda<<<blocksPerGrid, threadsPerBlock>>>(gpu_gates, n_gates, x_packed_gpu, z_packed_gpu, r_packed_gpu, rows, cols);
            cudaDeviceSynchronize();
            
            //End timer
            sim_timer.stop_timer();
            sim_time = sim_timer.measure();

            //Copy result back to host
            unpack_to_cpu();

            if (Config::PRINT_SIM_TRACE)
            {
                printf("\n============== STAB-Sim ===============\n");
                printf("n_qubits:%lld, n_gates:%lld, ncpus:%d, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms\n",
                       n, n_gates, 1, sim_time, 0.,
                       sim_time);
                printf("=====================================\n");
            }

            SAFE_FREE_GPU(gpu_gates);

            //=========================================
        }

        IdxType measure(IdxType qubit) override
        {
            throw std::logic_error("measure Not implemented (STAB_CPU)");
        }
        void set_initial(std::string fpath, std::string format) override
        {
            throw std::logic_error("set_initial Not implemented (STAB_CPU)");
        }
        void dump_res_state(std::string outfile) override
        {
            throw std::logic_error("dump_res_state Not implemented (STAB_CPU)");
        }
        ValType *get_real() const override
        {
            throw std::logic_error("get_real Not implemented (STAB_CPU)");
        }
        ValType *get_imag() const override
        {
            throw std::logic_error("get_imag Not implemented (STAB_CPU)");
        }
        ValType get_exp_z(const std::vector<size_t> &in_bits) override
        {
            throw std::logic_error("get_exp_z Not implemented (STAB_CPU)");
        }
        ValType get_exp_z() override
        {
            throw std::logic_error("get_exp_z Not implemented (STAB_CPU)");
        }

    protected:
        IdxType n;
        IdxType n_gates;
        int stabCounts;
        int mask;
        IdxType rows;
        IdxType packed_rows;
        IdxType cols;
        int packed_bits;
        IdxType packed_r_size;
        IdxType packed_matrix_size;
        std::vector<std::vector<int>> x;
        std::vector<std::vector<int>> z;
        std::vector<int> r;
        //CPU Arrays
        uint32_t** x_packed_cpu;
        uint32_t** z_packed_cpu;
        uint32_t* r_packed_cpu;
        //GPU Arrays
        uint32_t** x_packed_gpu;
        uint32_t** z_packed_gpu;
        uint32_t* r_packed_gpu;

        IdxType* totalResults = NULL;

        //Random
        std::mt19937 rng;
        std::uniform_int_distribution<int> dist;   
        

        __global__ void stab_simulation_kernel_cuda(Gate* gates, int g, uint32_t** x, uint32_t** z, uint32_t* r, IdxType rows, IdxType cols)
        {
            int row = blockIdx.x * blockDim.x + threadIdx.x;
            if (row >= rows) return;

            for (int k = 0; k < g; k++) 
            {
                Gate gate = gates[k];
                int a = gate.qubit;
                assert(a < cols);

                switch (gate.op_name) 
                {
                    case OP::H:
                        r[row] ^= (x[row][a] & z[row][a]);
                        uint32_t temp = x[row][a];
                        x[row][a] = z[row][a];
                        z[row][a] = temp;
                        break;

                    case OP::S:
                        r[row] ^= (x[row][a] & z[row][a]);
                        z[row][a] ^= x[row][a];
                        break;

                    case OP::CX:
                        int b = gate.qubit;
                        int ctrl = gate.ctrl;
                        r[row] ^= ((x[row][ctrl] & z[row][b]) & (x[row][b] ^ z[row][ctrl] ^ 1));
                        x[row][b] ^= x[row][ctrl];
                        z[row][ctrl] ^= z[row][b];
                        break;
                    
                    case OP::X:
                        r[row] ^= z[row][a];
                        break;

                    case OP::Y: 
                        r[row] ^= x[row][a] ^ z[row][a];
                        break;

                    case OP::Z:
                        r[row] ^= x[row][a];
                        break;
                    
                    default:
                        printf("Non-Clifford or unrecognized gate: %d\n", gate.op_name);
                        assert(false);
                }
            }
        }
    }; //End tableau class
} // namespace NWQSim

// #endif