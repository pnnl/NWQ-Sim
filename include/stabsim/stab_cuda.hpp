#pragma once

// #ifndef tableau_addons
// #define tableau_addons

//Eigen and pauli math
#include "pauli_math.hpp"

#include "../state.hpp"

#include "../nwq_util.hpp"
#include "../gate.hpp"
#include "../circuit.hpp"
#include "../config.hpp"
#include "../private/exp_gate_declarations_host.hpp"

#include "../circuit_pass/fusion.hpp"
#include "../private/macros.hpp"
#include "../private/sim_gate.hpp"

#include <random>
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <set>

namespace NWQSim
{
    class STAB_CUDA : public QuantumState
    {

    public:
        //Default identity constructor
        STAB_CUDA(IdxType _n_qubits) : QuantumState(SimType::STAB)
        {
            n = _n_qubits;
            rows = 2*n+1;
            cols = n;
            stabCounts = n;
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
            //(2n+1 x 2n+1) 
            //2n+1 rows comes from destab+stab+scratch space
            //2n+1 cols comes from each qubit getting an x and z matrix, plus the r column
            cudaMalloc(&host_XZR, (2*n+1)/*rows*/ * (2*n+1)/*cols*/ * sizeof(uint32_t)/32);
            /*we packed the uint32 so we can divide out 32 bits*/

            combine_words_and_interleave();

            //Packs x, z, and r bits into interleaved gpu warps (2D array of uint 32's)
            //Columns in the new format are stabilizers, rows are qubits. Each qubit gets a row of x and a row of z
            PackTo32Col();

            cudaMemcpy(gpu_XZR, host_XZR, (2n) * (2*n+2) * sizeof(uint32_t)/32, cudaMemcpyHostToDevice);

            int threads_per_block = 32;
            int num_blocks = (2*n) * (2*n+2) / (32*32);

            rng.seed(Config::RANDOM_SEED);
            dist = std::uniform_int_distribution<int>(0,1);

            SAFE_ALOC_HOST(totalResults, sizeof(IdxType));
            memset(totalResults, 0, sizeof(IdxType));
        }

        ~STAB_CUDA()
        {
            // Release for CPU side
            SAFE_FREE_HOST(totalResults);
        }

        void combine_words_and_interleave()
        {
            host_XZR = nullptr;

            uint32_t* temp_row;
            
            //Number of warps wide (number of lanes in a row) = (x.size()/32)+1
            for(int row = 0; row < (n*2); row++)
            {
                for(int i = 0; i < (x.size()/32)+1; i++)
                {
                    //Pack x bits into a uint 32 array. This array serves as a row of the total grid.
                    for(int j = 0; j < 32; j++)
                    {
                        host_XZR[row][i] |= x[(i * 32) + j] << j;
                    }

                    row++;

                    //Pack z bits into a uint 32 array. This array serves as a row of the total grid.
                    for(int j = 0; j < 32; j++)
                    {
                        host_XZR[row][i] |= z[(i * 32) + j] << j;
                    }
                }
            }
            
            
            //the last two rows -- buffer and r -- are treated seperately and appended to the grid
        }

        template <typename T>
        __global__ void PackTo32Col(const T* X, unsigned* B, int height, int width)
        {
            unsigned Bval, laneid;//fetch lane−id
            asm("mov.u32 %0, %%laneid;":"=r"(laneid));

            #pragma unroll
            for(int i = 0; i < WARP_SIZE; i++)
            {
                T f0 = X[(blockIdx.x * WARP_SIZE + i) * width + blockIdx.y * WARP_SIZE + laneid];

                //rotate anticlockwise for Col packing
                unsigned r0 = __brev(__ballot(f0> = 0));

                if (laneid == i ) 
                    Bval = r0;
            }
            
            B[blockIdx.y * height + blockIdx.x * WARP_SIZE + laneid] = Bval;
        }

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
            //Tableau format: 2n x n X matrix, 2n x n Z matrix and 2n x 1 phase column. All have an extra 1 x n scratch space (1x1 for phase).
            //n is the number of qubits. First n rows represent destabilizer pauli strings, second n rows represent stabilizers
            //columns represent qubits 
            //Destabilizer/stabilizer rows can have more stabilizers added later ti become 2*n x m
            rows = 2*n+1;
            cols = n;

            stabCounts = n;

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
                        //std::cout << "x at [" << p_index << "][" << a << "] = " << x[p_index][a] << std::endl;
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
            IdxType n_gates = gates.size();

            //Check that the circuit object matches the tableau
            assert(circuit->num_qubits() == n);

            //Start a timer
            cpu_timer sim_timer;
            sim_timer.start_timer();

            //Perform the tableau simulation
            simulation_kernel(gates);

            //End timer
            sim_timer.stop_timer();
            sim_time = sim_timer.measure();

            if (Config::PRINT_SIM_TRACE)
            {
                printf("\n============== STAB-Sim ===============\n");
                printf("n_qubits:%lld, n_gates:%lld, ncpus:%d, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms\n",
                       n, n_gates, 1, sim_time, 0.,
                       sim_time);
                printf("=====================================\n");
            }
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

        // Kernel for parallelized row operations
        __global__ void S_kernel(uint32_t* x, uint32_t* z, uint32_t* r, int rows, int cols, int a) {
            int tid = threadIdx.x + blockIdx.x * blockDim.x; // Thread ID
            int warp_id = tid / 32;  // Warp ID (32 threads per warp)
            int lane = tid % 32;     // Lane ID within the warp

            int rows_per_warp = (rows + gridDim.x * blockDim.x / 32 - 1) / (gridDim.x * blockDim.x / 32);
            int start_row = warp_id * rows_per_warp;
            int enRow = min(start_row + rows_per_warp, rows);

            for (int i = start_row + lane; i < enRow; i += 32) {
                // Load data for the row
                uint32_t x_row = x[i];
                uint32_t z_row = z[i];
                uint32_t r_row = r[i];

                // Perform operations on the row
                uint32_t xa = (x_row >> a) & 1; // Extract bit `a` from x_row
                uint32_t za = (z_row >> a) & 1; // Extract bit `a` from z_row

                // Phase update
                r_row ^= (xa & za);

                // Entry update
                z_row ^= (xa << a);

                // Store updated data back to global memory
                z[i] = z_row;
                r[i] = r_row;
            }
        }


    protected:
        IdxType n;
        int stabCounts;
        IdxType rows;
        IdxType cols;
        std::vector<std::vector<int>> x;
        std::vector<std::vector<int>> z;
        std::vector<int> r;
        IdxType* totalResults = NULL;
        //Device pointers for the gpu
        uint32_t** host_XZR = nullptr;
        uint32_t** gpu_XZR = nullptr;


        std::mt19937 rng;
        std::uniform_int_distribution<int> dist;   

        void simulation_kernel(std::vector<Gate>& gates)
        {
            int g = gates.size();

            //For swapping rows
            int tempVal;
            int half_rows = rows >> 1;
            //Loop over every gate in the circuit and apply them
            for (int k = 0; k < g; k++)
            {
                auto gate = gates[k];
                int a = gate.qubit;
                assert(a < n);

                if (gate.op_name == OP::H)
                {
                    
                    for(int i = 0; i < rows-1; i++)
                    {
                        //Phase
                        r[i] ^= (x[i][a] & z[i][a]);
                        //Entry -- swap x and z bits at qubit a
                        std::swap(x[i][a], z[i][a]);
                    } 
                }
                else if (gate.op_name == OP::S)
                {
                    int a = gate.qubit;
                    Skernel<<<num_blocks, threads_per_block>>>(X, Z, R, rows, cols, a);


                }
                else if (gate.op_name == OP::CX)
                {  
                    int a = gate.ctrl;
                    int b = gate.qubit;
                    for(int i = 0; i < rows-1; i++)
                    {
                        //Phase
                        r[i] ^= ((x[i][a] & z[i][b]) & (x[i][b]^z[i][a]^1));

                        //Entry
                        x[i][b] ^= x[i][a];
                        z[i][a] ^= z[i][b];
                    }
                }
                else if (gate.op_name == OP::M)
                {  
                    

                    int a = gate.qubit;
                    int p = -1;
                    for(int p_index = half_rows; p_index < rows-1; p_index++)
                    {  
                        //std::cout << "x at [" << p_index << "][" << a << "] = " << x[p_index][a] << std::endl;
                        if(x[p_index][a])
                        {
                            p = p_index;
                            break;
                        }
                    }
                    //A p such that x[p][a] = 1 exists
                    if(p > -1)
                    {
                        for(int i = 0; i < rows-1; i++)
                        {
                            if((x[i][a]) && (i != p))
                            {
                                rowsum(i, p);
                            }
                        }
                        x[p-half_rows] = x[p];
                        z[p-half_rows] = z[p];
                        //Change all the columns in row p to be 0
                        for(int i = 0; i < n; i++)
                        {
                            x[p][i] = 0;
                            z[p][i] = 0;                        
                        }

                        //Generate and display a random number
                        int randomBit = dist(rng);
                        
                        if(randomBit)
                        {
                            //std::cout << "Random result of 1" << std::endl;
                            r[p] = 1;
                        }
                        else
                        {
                            //std::cout << "Random result of 0" << std::endl;
                            r[p] = 0;
                        }
                        z[p][a] = 1;

                        totalResults[0] += r[p] << a;
                        // std::cout << "Random measurement at qubit " << a << " value: " << (r[p] << a) << std::endl;
                    }
                    else
                    {

                        //Set the scratch space row to be 0
                        //i is the column indexer in this case
                        for(int i = 0; i < n; i++)
                        {
                            x[rows-1][i] = 0;
                            z[rows-1][i] = 0;
                        }
                        r[rows-1] = 0;

                        //Run rowsum subroutine
                        for(int i = 0; i < half_rows; i++)
                        {
                            if(x[i][a] == 1)
                            {
                                //std::cout << "Perform rowsum at " << i << " + n" << std::endl;
                                rowsum(rows-1, i+half_rows);
                            }
                        }
                        // std::cout << "Deterministc measurement at qubit " << a << " value: " << (r[rows-1] << a) << std::endl;
                        totalResults[0] += r[rows-1] << a;
                    }
                } //End M
                else if (gate.op_name == OP::X)
                {
                    //equiv to H S S H or H Z H
                    for(int i = 0; i < rows-1; i++)
                    {
                        r[i] ^= z[i][a];
                    } 
                }
                else if (gate.op_name == OP::Y)
                {   
                    //equiv to Z X
                    for(int i = 0; i < rows-1; i++)
                    {
                        r[i] = r[i] ^ x[i][a] ^ z[i][a];
                    }
                }
                else if (gate.op_name == OP::Z)
                {
                    //equiv to S S gates
                    for(int i = 0; i < rows-1; i++)
                    {
                        //Phase
                        r[i] ^= x[i][a];
                    }
                }
                else if (gate.op_name == OP::MA)
                {
                    measure_all();
                }
                else    
                {
                    std::cout << "Non-Clifford or unrecognized gate: "
                                << OP_NAMES[gate.op_name] << std::endl;
                    std::logic_error("Invalid gate type");
                    exit(1);
                }
            } //End gates for loop
        } //End simulate
    }; //End tableau class
} // namespace NWQSim

// #endif