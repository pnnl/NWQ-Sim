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
    class STAB_AVX : public QuantumState
    {

    public:
        //Default identity constructor
        STAB_AVX(IdxType _n_qubits) : QuantumState(SimType::STAB)
        {
            n = _n_qubits;
            rows = 2*n+1;
            cols = n;
            stabCounts = n;
            bit_x.resize(rows, std::vector<int>(cols,0)); //first 2n+1 x n block. first n represents destabilizers
                                                       //second n represents stabilizers + 1 extra row
            bit_z.resize(rows, std::vector<int>(cols,0)); //second 2n+1 x n block to form the 2n+1 x 2n sized tableau
            bit_r.resize(rows, 0); //column on the right with 2n+1 rows
            //The 2n+1 th row is scratch space

            //Intialize the identity tableau
            for(int i = 0; i < n; i++)
            {
                bit_x[i][i] = 1;
                bit_z[i+n][i] = 1;
            }

            x = packBitsToAVX(bit_x);
            z = packBitsToAVX(bit_z);
            r = packBitsToAVX(bit_r);
            
            rng.seed(Config::RANDOM_SEED);
            dist = std::uniform_int_distribution<int>(0,1);

            has_destabilizers = true;

            SAFE_ALOC_HOST(totalResults, sizeof(IdxType));
            memset(totalResults, 0, sizeof(IdxType));
        }

        ~STAB_AVX()
        {
            // Release for CPU side
            SAFE_FREE_HOST(totalResults);
        }

        //Pack bits to a 2D AVX array
        std::vector<std::vector<__m256i>> packBitsToAVX(const std::vector<std::vector<int>>& input) 
        {
            size_t avx_cols = (input[0].size()/256) + 1; //round up

            //Resulting AVX 2D vector
            std::vector<std::vector<__m256i>> avxData(rows, std::vector<__m256i>(avx_cols, _mm256_setzero_si256()));

            for (size_t i = 0; i < rows; i++) 
            {
                for (size_t j = 0; j < avx_cols; j++)
                {
                    for (size_t k = 0; k < 256; k++) 
                    {
                        size_t target_col = (j * 256) + k;
                        while(cols > target_col)
                        {
                            if(input[i][target_col] == 1) 
                            {
                                avxData[i][j] = _mm256_or_si256(avxData[i][j], _mm256_set1_epi64x(1ULL << k));
                            }
                        }
                    }
                }
            }
            return avxData;
        }

        void delete_all_rows() override
        {
            rows = 0;
            stabCounts = 0;
            x.resize(rows, std::vector<__m256i>(cols,0));
            z.resize(rows, std::vector<__m256i>(cols,0)); 
            r.resize(rows, 0);   
        }

        void remove_destabilizers() override
        {
            has_destabilizers = false;
        
            x.erase(x.begin(), x.begin()+rows/2);
            z.erase(z.begin(), z.begin()+rows/2);
            r.erase(r.begin(), r.begin()+rows/2);

            rows = (rows/2) + 1;
        }

        int get_num_rows() override
        {
            return rows;
        }

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
            bit_x.resize(rows, std::vector<int>(cols,0)); //first 2n+1 x n block. first n represents destabilizers
                                                       //second n represents stabilizers + 1 extra row
            bit_z.resize(rows, std::vector<int>(cols,0)); //second 2n+1 x n block to form the 2n+1 x 2n sized tableau
            bit_r.resize(rows, 0); //column on the right with 2n+1 rows
            //The 2n+1 th row is scratch space

            //Intialize the identity tableau
            for(int i = 0; i < n; i++)
            {
                bit_x[i][i] = 1;
                bit_z[i+n][i] = 1;
            }

            x = packBitsToAVX(bit_x);
            z = packBitsToAVX(bit_z);
            r = packBitsToAVX(bit_r);
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

        //Get the stabilizers in the tableau, i.e. all of the Pauli strings that stabilize the circuit
        std::vector<std::string> get_destabilizers()
        {
            int x_val;
            int z_val;
            std::string stabilizers;
            std::vector<std::string> pauliStrings;
            for(int i = 0; i < (rows>>1); i++) //rows of stabilizers
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
                            stabilizers += "X";
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

        //Get a pauli string stabilizer and phase bit
        std::pair<std::string,int> get_stabilizer_line(int row) override
        {
            std::string stabilizers;
            for(int j = 0; j < cols; j++) //qubits/cols
            {
                int x_val = x[row][j];
                int z_val = z[row][j];
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
            return std::make_pair(stabilizers, r[row]);
        }

        std::vector<std::pair<std::string,int>> get_all_lines()
        {
            std::vector<std::pair<std::string,int>> full_paulis;
            for(int row = (rows>>1); row < rows-1; row++) //rows of stabilizers
            {
                std::string stabilizers;
                for(int j = 0; j < cols; j++) //qubits/cols
                {
                    int x_val = x[row][j];
                    int z_val = z[row][j];
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
                full_paulis.push_back(std::make_pair(stabilizers, r[row]));
            }
            return full_paulis;
        }

        //Creates a sparse density matrix out of the Pauli strings stabilizers
        std::vector<std::vector<double>> get_density_matrix() override
        {
            std::vector<std::pair<std::string,int>> full_paulis = get_all_lines();
            std::pair<std::string,int> line;
            ComplexSparseMatrix DM = createSparseIdentity(1 << n);
            for(int i = 0; i < full_paulis.size(); i++)
            {
                line = full_paulis[i];
                std::cout << "Line.first " << line.first << std::endl;
                ComplexSparseMatrix I = createSparseIdentity(1 << n);
                ComplexSparseMatrix S = createSparseIdentity(1);
                for(int j = 0; j < line.first.size(); j++)
                {
                    switch(line.first[j])
                    {
                        case 'I':
                            S = kroneckerProduct(S, sparsePauliI());
                            break;
                        case 'X':
                            S= kroneckerProduct(S, sparsePauliX());
                            break;
                        case 'Y':   
                            S = kroneckerProduct(S, sparsePauliY());
                            break;
                        case 'Z':
                            S = kroneckerProduct(S, sparsePauliZ());
                            break;
                        default:
                            std::logic_error("Invalid stabilizer");
                            break;
                    }
                }
                
                if(line.second == 1)
                    S = S * -1.0;

                DM = DM * .5 * (I + S);
            }

            std::vector<std::vector<double>> densityMatrix;
            densityMatrix.resize(DM.rows(), std::vector<double>(DM.cols(),0));

            for(int i = 0; i < DM.rows(); i++)
            {
                for(int j = 0; j < DM.cols(); j++)
                {
                    densityMatrix[i][j] = DM.coeff(i, j).real();
                }
            }

            return densityMatrix;
        }


        //Get the stabilizers in the tableau, i.e. all of the Pauli strings that stabilize a certain circuit
        std::vector<std::string> get_stabilizers() override
        {
            int x_val;
            int z_val;
            std::string stabilizers;
            std::vector<std::string> pauliStrings;
            int start;
            if(has_destabilizers)
                start = (rows/2);
            else
                start = 0;
            for(int i = start; i < rows-1; i++) //rows of stabilizers
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

        //Removes all but one repitions of a stabilizer
        void remove_repetitions(std::string stab, int num_s)
        {
            std::vector<int> temp_x;
            std::vector<int> temp_z;

            for (int j = 0; j < stab.size(); j++) 
            {
                switch(stab[j])
                {
                    case 'I':
                        temp_x.push_back(0);
                        temp_z.push_back(0);
                        break;
                    case 'X':
                        temp_x.push_back(1);
                        temp_z.push_back(0);
                        break;
                    case 'Y':   
                        temp_x.push_back(1);
                        temp_z.push_back(1);
                        break;
                    case 'Z':
                        temp_x.push_back(0);
                        temp_z.push_back(1);
                        break;
                    default:
                        std::logic_error("Invalid stabilizer");
                        break;
                }

                //remove 2 T's from the given T tableau for every S gate
                while(num_s > 0)
                {
                    for(int j = 0; j < rows; j++)
                    {
                        if((x[j] == temp_x) && (z[j] == temp_z))
                        {
                            remove_stabilizer(j);
                            num_s--;
                            break;
                        }
                    }
                }
            }
        }

        //Returns a map of stabilizers and the number of times they ocurr in the tableau
        std::unordered_map<std::string, int> stabilizer_count() override
        {
            std::unordered_map<std::string, int> map;
            std::vector<std::string> stabs = get_stabilizers();
            for(int i = 0; i < stabs.size(); i++)
                map[stabs[i]]++;            
            return map;
        }

        bool check_commutation(std::string pauliString) override
        {
            int new_x;
            int new_z;
            int start = 0;

            //Loop over every stabilizer of the tableau
            if(has_destabilizers)
                start = rows/2;
            for (int i = start; i < x.size()-1; i++) 
            {
                //Compute inner product
                int product = 0;
                
                //loop over every column
                for (int j = 0; j < x[i].size(); j++) 
                {
                    switch(pauliString[j])
                    {
                        case 'I':
                            new_x = 0;
                            new_z = 0;
                            break;
                        case 'X':
                            new_x = 1;
                            new_z = 0;
                            break;
                        case 'Y':   
                            new_x = 1;
                            new_z = 1;
                            break;
                        case 'Z':
                            new_x = 0;
                            new_z = 1;
                            break;
                        default:
                            std::logic_error("Invalid stabilizer");
                            break;
                    }
                    product += (x[i][j] * new_z + z[i][j] * new_x) % 2;
                }
                if (product % 2 != 0) {
                    return false; //Anti-commutation somewhere in the Pauli string
                }
            }
            return true; //Commutes with all stabilizers
        }
        
        //True if pauliString and target row stabilizer commute
        //False if not
        bool check_row_commutation(std::string pauliString, int row) override
        {
            int temp_x;
            int temp_z;
            int product = 0;
            for(int col = 0; col < cols; col++)
            {
                switch(pauliString[col])
                {
                    case 'I':
                        temp_x = 0;
                        temp_z = 0;
                        break;
                    case 'X':
                        temp_x = 1;
                        temp_z = 0;
                        break;
                    case 'Y':   
                        temp_x = 1;
                        temp_z = 1;
                        break;
                    case 'Z':
                        temp_x = 0;
                        temp_z = 1;
                        break;
                    default:
                        std::logic_error("Invalid stabilizer");
                        break;
                }
                product += (x[row][col] * temp_z + z[row][col] * temp_x);
            }
            if(product%2)
                return true;
            return false;
        }

        void add_stabilizer(std::string pauliString, int phase_bit = 0) override
        {
            assert(pauliString.length() == n);
            if(has_destabilizers)
            {
                if(check_commutation(pauliString))
                {
                    //Start by adding a row of destabilizers and stabilizers to T
                    stabCounts++;
                    rows+=2;
                    x.insert(x.begin() + x.size()/2, std::vector<int>(cols,0));
                    x.insert(x.end()-1, std::vector<int>(cols,0));
                    z.insert(z.begin() + z.size()/2, std::vector<int>(cols,0));
                    z.insert(z.end()-1, std::vector<int>(cols,0));
                    r.insert(r.begin() + r.size()/2, phase_bit);
                    r.insert(r.end()-1, phase_bit);
                    
                    //Stabilizer and destabilizer addition
                    for(int i = 0; i < pauliString.length(); i++)
                    {
                        switch(pauliString[i])
                        {
                            case 'I':
                                x[rows-2][i] = 0;
                                z[rows-2][i] = 0;
                                x[(rows>>1)-1][i] = 0;
                                z[(rows>>1)-1][i] = 0;
                                break;
                            case 'X':
                                x[rows-2][i] = 1;
                                z[rows-2][i] = 0;
                                x[(rows>>1)-1][i] = 0;
                                z[(rows>>1)-1][i] = 1;
                                break;
                            case 'Y':   
                                x[rows-2][i] = 1;
                                z[rows-2][i] = 1;
                                //make the destabilizer X to anticommute with Y
                                x[(rows>>1)-1][i] = 1;
                                z[(rows>>1)-1][i] = 0;
                                //add an i to the 2 bit phase representation at the stabilizer row
                                r[rows-1] = (r[rows-1] + 1) % 4;
                                break;
                            case 'Z':
                                x[rows-2][i] = 0;
                                z[rows-2][i] = 1;
                                x[(rows>>1)-1][i] = 1;
                                z[(rows>>1)-1][i] = 0;
                                break;
                            default:
                                std::logic_error("Invalid stabilizer");
                                break;
                        }
                    }
                }
                else
                {
                    std::logic_error("Stabilizer fails commutation check" + pauliString);
                }
            }
            else
            {
                if(check_commutation(pauliString))
                {
                    //Start by adding a row of stabilizers to T
                    stabCounts++;
                    rows++;
                    x.insert(x.end()-1, std::vector<int>(cols,0));
                    z.insert(z.end()-1, std::vector<int>(cols,0));
                    r.insert(r.end()-1, phase_bit);
                    
                    //Stabilizer and destabilizer addition
                    for(int i = 0; i < pauliString.length(); i++)
                    {
                        switch(pauliString[i])
                        {
                            case 'I':
                                x[rows-2][i] = 0;
                                z[rows-2][i] = 0;
                                x[(rows>>1)-1][i] = 0;
                                z[(rows>>1)-1][i] = 0;
                                break;
                            case 'X':
                                x[rows-2][i] = 1;
                                z[rows-2][i] = 0;
                                x[(rows>>1)-1][i] = 0;
                                z[(rows>>1)-1][i] = 1;
                                break;
                            case 'Y':   
                                x[rows-2][i] = 1;
                                z[rows-2][i] = 1;
                                //make the destabilizer X to anticommute with Y
                                x[(rows>>1)-1][i] = 1;
                                z[(rows>>1)-1][i] = 0;
                                //add an i to the 2 bit phase representation at the stabilizer row
                                r[rows-1] = (r[rows-1] + 1) % 4;
                                break;
                            case 'Z':
                                x[rows-2][i] = 0;
                                z[rows-2][i] = 1;
                                x[(rows>>1)-1][i] = 1;
                                z[(rows>>1)-1][i] = 0;
                                break;
                            default:
                                std::logic_error("Invalid stabilizer");
                                break;
                        }
                    }
                }
                else
                {
                    std::logic_error("Stabilizer fails commutation check" + pauliString);
                }
            }
        }

        //Replaces a stabilizer pauli string at some row in the Tableau. Useful for initializing a
        // new Tableau in a for loop without circuit initialization
        void replace_stabilizer(std::string pauliString, int stabPos) override
        {
            assert(pauliString.length() == n);
            //Start by adding a row of destabilizers and stabilizers to T
            
            //Stabilizer and destabilizer addition
            for(int i = 0; i < pauliString.length(); i++)
            {
                //paulistring length is the same as number of qubits, so we use it as the column iterator
                switch(pauliString[i])
                {
                    case 'I':
                        //stab
                        x[(rows>>1) + stabPos][i] = 0;
                        z[(rows>>1) + stabPos][i] = 0;
                        //destab
                        x[stabPos][i] = 0;
                        z[stabPos][i] = 0;
                        break;
                    case 'X':
                        //stab
                        x[(rows>>1) + stabPos][i] = 1;
                        z[(rows>>1) + stabPos][i] = 0;
                        //destab
                        x[stabPos][i] = 0;
                        z[stabPos][i] = 1;
                        break;
                    case 'Y':
                        //stab
                        x[(rows>>1) + stabPos][i] = 1;
                        z[(rows>>1) + stabPos][i] = 1;
                        //make the destabilizer X to anticommute with Y
                        x[stabPos][i] = 1;
                        z[stabPos][i] = 0;
                        //add an i to the 2 bit phase representation at the stabilizer row
                        r[(rows>>1) + stabPos] = (r[(rows>>1) + stabPos] + 1) % 4;
                        break;
                    case 'Z':
                        x[(rows>>1) + stabPos][i] = 0;
                        z[(rows>>1) + stabPos][i] = 1;
                        x[stabPos][i] = 1;
                        z[stabPos][i] = 0;
                        break;
                    default:
                        std::logic_error("Invalid stabilizer");
                        break;
                }
            }
        }

        void remove_stabilizer(int row_index)
        {
            x.erase(x.begin()+row_index);
            z.erase(z.begin()+row_index);
            r.erase(r.begin()+row_index);
        }
        
        void transpose(std::vector<std::vector<int>>& M)
        {
            int rowSize = M.size();
            int colSize = M[0].size();
            std::vector<std::vector<int>> MT = M;


            for(int i = 0; i < rowSize; i++)
            {
                for(int j = 0; j < colSize; j++)
                {
                    MT[i][j] = M[j][i];
                }
            }

            M = MT;
        }

        //Algorithm to reduce the tableau to rref        
        void gaussianElimination() 
        {
            for (int col = 0; col < cols; col++) {
                int pivotRow = -1;
                for (int row = col; row < rows-1; row++) {
                    if (x[row][col] == 1 || z[row][col] == 1) {
                        pivotRow = row;
                        break;
                    }
                }
                if (pivotRow == -1) continue;
                if (pivotRow != col) {
                    swapRows(pivotRow, col);
                }
                for (int row = 0; row < rows-1; row++) {
                    if (row != col && (x[row][col] == 1 || z[row][col] == 1)) {
                        addRows(row, col);
                    }
                }
            }
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
            sum = sum + 2*r[h] + 2*r[i];

            if(sum % 4 == 0)
                r[h] = 0;
            else
                r[h] = 1;
        } //End rowsum

        //For shot based measurement
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
            sum = sum + 2*temp_r[h] + 2*temp_r[i];

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
            SAFE_FREE_HOST(totalResults);
            SAFE_ALOC_HOST(totalResults, sizeof(IdxType) * shots);
            memset(totalResults, 0, sizeof(IdxType) * shots);

            int half_rows = rows >> 1;
            
            for(int i = 0; i < shots; i++)
            {
                //Make a copy of the class being measured so many shots can be performed
                std::vector<std::vector<int>> temp_x = x;
                std::vector<std::vector<int>> temp_z = z;
                std::vector<int> temp_r = r;
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
                    if(p > -1)
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
                    else
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

        std::vector<std::vector<int>> get_graph_matrix() override
        {
            std::cout << "-----Test-----" << std::endl;

            if (x.size() != z.size()) {
                std::cerr << "Error: Matrices must have the same number of rows to append column-wise." << std::endl;
            }

            //Extract just the stabilizers w/o phase
            std::vector<std::vector<int>> graphMatrix;
            graphMatrix.resize((rows>>1), std::vector<int>(2*cols,0));

            for (int i = 0; i < (rows>>1); i++) {
                for(int j = 0; j < cols; j++)
                {
                    graphMatrix[i][j] = x[i+(rows>>1)][j]; //Copy the x row
                    graphMatrix[i][j+cols] = z[i+(rows>>1)][j];
                }
            }
            std::cout << "\n-----Original-----" << std::endl;
            print_table(graphMatrix);
            
            reduction(graphMatrix);

            std::cout << "\n-----Reduction-----" << std::endl;
            print_table(graphMatrix);

            //convert to 'standard form'
            int xRank = nonZeroRows(extractHalf(graphMatrix, true));
            std::cout << "----- X rank after reduction ----- " << xRank << ", stabCounts: " << stabCounts << std::endl;
            if(xRank != stabCounts)
            {
                zReduction(graphMatrix);
                //permute qubits using the qubitIndex

                //Adjust qubit index to maximize the size of the identity block in the Z part
                std::vector<int> qubitIndex(n, 0);
                for (int i = 0; i < xRank; i++) {
                    qubitIndex[i] = i;
                }

                std::vector<std::vector<int>> zMatrix = extractHalf(graphMatrix, false);
            }//end if xRank/standard form



            xRank = nonZeroRows(extractHalf(graphMatrix, true));

            std::cout << "\n-----Standard form-----" << std::endl;
            print_table(graphMatrix);

            /*
                Graph processing after standard form
            */

            if(nonZeroRows(graphMatrix) != n)
                std::cout << "\n\n\nMismatch row/cols\n\n\n";
            
            std::vector<std::vector<int>> xStabs = extractHalf(graphMatrix, true);

            //find the size of the identity
            xRank = nonZeroRows(xStabs);

            //Extract submatrices from graphMatrix
            std::vector<std::vector<int>> I1(xRank, std::vector<int>(xRank));
            std::vector<std::vector<int>> A(xRank, std::vector<int>(n - xRank));
            std::vector<std::vector<int>> B(xRank, std::vector<int>(xRank));
            std::vector<std::vector<int>> sign1(xRank, std::vector<int>(1, 0));

            std::vector<std::vector<int>> D(graphMatrix.size() - xRank, std::vector<int>(xRank));
            std::vector<std::vector<int>> I2(graphMatrix.size() - xRank, std::vector<int>(n - xRank - 1));
            std::vector<std::vector<int>> sign2(graphMatrix.size() - xRank, std::vector<int>(1, 0));

            /*
            If the stabilizer is full rank (stabilizer count = qubit number), the stabilizer tableau can be reduced to
            [[I A | B 0]
             [0 0 | D I]]

            If B has any nonzero diagonal elements, apply S gates to turn Y -> X
            ==> 
            [[I A | B' 0]
             [0 0 | D  I]]
            where B' has zero diagonal elements
            */

            //Apply S on B (0 on diagonal)
            //This is equiv to z = (z + x)mod2 since the x matrix is an I 
            for(int i = 0; i < B.size(); i++)
            {
                B[i][i] = 0;
            }


            //Fill the submatrices with data from graphMatrix
            for (int i = 0; i < xRank; i++) {
                for (int j = 0; j < xRank; j++) {
                    I1[i][j] = graphMatrix[i][j];
                    B[i][j] = graphMatrix[i][n + j];
                }

                //Apply H on second block of qubits (x-z exchange) with I2
                for (int j = 0; j < n - xRank; j++) {
                    A[i][j] = graphMatrix[i][xRank + j];
                }
            }

            for (int i = 0; i < graphMatrix.size() - xRank; i++) {
                for (int j = 0; j < xRank; j++) {
                    D[i][j] = graphMatrix[xRank + i][n + j];
                }
                //Apply H on second block of qubits (x-z exchange) with A
                for (int j = 0; j < n - xRank - 1; j++) {
                    I2[i][j] = graphMatrix[xRank + i][n + xRank + j];
                }
            }

            //Set diagonal elements of B to zero
            for (int i = 0; i < B.size(); i++) {
                if (B[i][i] != 0) {
                    B[i][i] = 0;
                }
            }

            //Reconstruct the new graph matrix
            std::vector<std::vector<int>> newGraph(graphMatrix.size(), std::vector<int>(graphMatrix[0].size(), 0));

            //Fill newGraph with submatrices I1, A, B, I2, D, sign1, and sign2
            for (int i = 0; i < xRank; i++) {
                for (int j = 0; j < xRank; j++) {
                    newGraph[i][j] = I1[i][j];
                }
                for (int j = 0; j < A[0].size(); j++) {
                    newGraph[i][xRank + j] = A[i][j];
                }
                for (int j = 0; j < B[0].size(); j++) {
                    newGraph[i][n + j] = B[i][j];
                }
                newGraph[i][graphMatrix[0].size() - 1] = sign1[i][0];
            }

            for (int i = 0; i < D.size(); i++) {
                for (int j = 0; j < I2[0].size(); j++) {
                    newGraph[xRank + i][n + xRank + j] = I2[i][j];
                }
                for (int j = 0; j < D[0].size(); j++) {
                    newGraph[xRank + i][n + j] = D[i][j];
                }
                newGraph[xRank + i][graphMatrix[0].size() - 1] = sign2[i][0];
            }


            std::cout << "\n-----New graph-----" << std::endl;
            print_table(newGraph);

            return graphMatrix;
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
        bool has_destabilizers;
        const IdxType n;
        int stabCounts;
        IdxType rows;
        IdxType cols;
        std::vector<std::vector<int>> bit_x;
        std::vector<std::vector<int>> bit_z;
        std::vector<int> bit_r;
        std::vector<std::vector<__m256i>> x;
        std::vector<std::vector<__m256i>> z;
        std::vector<__m256i> r;
        IdxType* totalResults = NULL;

        std::mt19937 rng;
        std::uniform_int_distribution<int> dist;

        //Function to swap two rows of the tableau
        void swapRows(int row1, int row2) {
            std::swap(x[row1], x[row2]);
            std::swap(z[row1], z[row2]);
            std::swap(r[row1], r[row2]);
        }

        //Function to add row2 to row1 (mod 2)
        void addRows(int row1, int row2) {
            for (int i = 0; i < n; i++) {
                x[row1][i] ^= x[row2][i];  //XOR for mod 2 addition
                z[row1][i] ^= z[row2][i];  //XOR for mod 2 addition
            }
            r[row1] = (r[row1] + r[row2]) % 4;  //Phase flip
        }

        //Converts a matrix into column echelon and counts the number of pivot rows
        int countPivots(std::vector<std::vector<int>>& matrix) {
            int rows = matrix.size();
            int cols = matrix[0].size();
            int leadCol = 0; //Track the leading column position
            int pivotCount = 0; //To count the number of pivot rows

            for (int col = 0; col < cols && leadCol < rows; col++) {
                //Find a row with a leading 1 in the current column
                int pivotRow = leadCol;
                while (pivotRow < rows && matrix[pivotRow][col] == 0) {
                    pivotRow++;
                }

                //If no leading 1 is found, continue to the next column
                if (pivotRow == rows) {
                    continue;
                }

                //Swap the row with the leading 1 to the top of the current column's range
                std::swap(matrix[leadCol], matrix[pivotRow]);

                //Increment the pivot count
                pivotCount++;

                //Eliminate all other 1s in this column below the leading 1
                for (int row = leadCol + 1; row < rows; row++) {
                    if (matrix[row][col] == 1) {
                        //XOR the current row with the row containing the leading 1 to eliminate the 1
                        for (int j = 0; j < cols; j++) {
                            matrix[row][j] ^= matrix[leadCol][j];
                        }
                    }
                }

                //Move to the next row for the next iteration
                leadCol++;
            }

            return pivotCount;
        }

        int nonZeroRows(const std::vector<std::vector<int>>& matrix) {
            int nonZeroRows = 0;
            for(const auto& row : matrix) 
            {
                bool isNonZero = false;
                for(int value : row) 
                {
                    if (value != 0) 
                    {
                        isNonZero = true;
                        break;
                    }
                }
                if(isNonZero) 
                {
                    nonZeroRows++;
                }
            }
            return nonZeroRows;
        }

        //Function to search for non-zero rows starting from a specified row
        std::vector<int> searchNonZeroRows(const std::vector<std::vector<int>>& matrix, int col, int startRow = 0, int endRow = -1) {
            std::vector<int> nonZeroRows;
            if(endRow == -1) endRow = matrix.size();

            if(col > matrix[0].size())
            {
                std::cout << "Column index out of range";
            }
            
            for(int i = startRow; i < endRow; i++) {
                if (matrix[i][col] != 0) {
                    nonZeroRows.push_back(i);
                }
            }
            return nonZeroRows;
        }

        std::vector<std::vector<int>> extractHalf(const std::vector<std::vector<int>>& matrix, bool left) {
            int numRows = matrix.size();
            int numCols = matrix[0].size();
            int halfCols = numCols / 2;

            std::vector<std::vector<int>> result(numRows, std::vector<int>(halfCols));

            for (int i = 0; i < numRows; i++) {
                for (int j = 0; j < halfCols; j++) {
                    result[i][j] = left ? matrix[i][j] : matrix[i][j + halfCols];
                }
            }
            return result;
        }

        void gaussianElimination(std::vector<std::vector<int>> &matrix) {
            int rows = matrix.size();
            int cols = matrix[0].size();

            //Perform Gaussian elimination
            for (int i = 0; i < rows; i++) {
                //Find the pivot in the current column or further to the right
                int maxCol = i;
                while (maxCol < cols && matrix[i][maxCol] == 0) {
                    //Search for a non-zero element in the current row to use as a pivot
                    for (int k = i + 1; k < rows; k++) {
                        if (matrix[k][maxCol] != 0) {
                            std::swap(matrix[i], matrix[k]);
                            break;
                        }
                    }
                    if (matrix[i][maxCol] == 0) {
                        maxCol++; //Move to the next column if no pivot is found
                    }
                }

                if (maxCol == cols) {
                    //If no pivot is found in the entire row, skip to the next iteration
                    continue;
                }

                //Swap columns if necessary to make the pivot in the current column
                if (maxCol != i) {
                    for (int k = 0; k < rows; k++) {
                        std::swap(matrix[k][i], matrix[k][maxCol]);
                    }
                }

                //Row operations to eliminate entries below the pivot
                for (int j = i + 1; j < rows; j++) {
                    if (matrix[j][i] != 0) {
                        for (int k = i; k < cols; k++) {
                            matrix[j][k] ^= matrix[i][k];
                        }
                    }
                }
            }

            //Eliminate empty rows after reduction
            bool allzero;
            for(int i = 0; i < matrix.size(); i++)
            {
                allzero = true;
                for(int j = 0; j < matrix.size(); j++)
                {
                    if(matrix[i][j])
                    {
                        allzero = false;
                        break;
                    }
                }
                if(allzero)
                {
                    matrix = std::vector<std::vector<int>>(matrix.begin(), matrix.begin() + i);
                    break;
                }
            }
        }

        //Function to perform the row addition and reduction process (reduceUp)
        void reduceUp(std::vector<std::vector<int>>& matrix, const std::vector<int>& columnIndices) {
            int j0 = columnIndices[columnIndices.size()-1];  //The index of the last row in columnIndices for addition
            // std::cout << columnIndices.size()-1<< std::endl;
            
            //Iterate over all rows in columnIndices except the last one
            for (int i = 0; i < columnIndices.size() - 1; i++) {
                int j = columnIndices[i];
                //Perform row addition: matrix[j] = matrix[j] ⊕ matrix[j0]
                for (int k = 0; k < matrix[j].size(); k++) {
                    matrix[j][k] ^= matrix[j0][k];  //XOR operation for binary addition
                }
            }
        }

        //Function to perform the row addition and reduction process (reduceDown)
        void reduceDown(std::vector<std::vector<int>>& matrix, const std::vector<int>& columnIndices) {
            int j0 = columnIndices[0];  //The index of the base row for addition

            for (int i = 1; i < columnIndices.size(); i++) {
                int j = columnIndices[i];

                //Perform row addition: matrix[j] = matrix[j] ⊕ matrix[j0]
                for (int k = 0; k < matrix[j].size(); k++) {
                    matrix[j][k] ^= matrix[j0][k];  //XOR operation for binary addition
                }
            }
        }
        //Function to perform a row insertion (row permutation) operation in the tableau
        void rowInsert(std::vector<std::vector<int>>& matrix, int i, int j) {
            int totalRows = matrix.size();

            //Validate the indices
            if (i >= totalRows) {
                throw std::out_of_range("The inserted row index is out of range.");
            }
            if (i < j) {
                throw std::invalid_argument("The inserted row index is before the insertion place.");
            }
            if (i == j) {
                return;  //No need to do anything if the row is already in place
            }

            //Create the permutation index order
            std::vector<int> indexOrder;
            for (int k = 0; k < j; k++) {
                indexOrder.push_back(k);
            }
            for (int k = j + 1; k <= i; k++) {
                indexOrder.push_back(k);
            }
            indexOrder.push_back(j);
            for (int k = i + 1; k < totalRows; k++) {
                indexOrder.push_back(k);
            }

            //Permute the rows of the matrix according to indexOrder
            std::vector<std::vector<int>> permutedMatrix(totalRows);
            for (int k = 0; k < totalRows; k++) {
                permutedMatrix[k] = matrix[indexOrder[k]];
            }
            matrix = permutedMatrix;  //Update the original matrix with the permuted matrix
        }

        void qubitReorder(std::vector<std::vector<int>>& matrix, std::vector<int> indices)
        {
            //Check if the length of the indices matches the number of qubits
            if (indices.size() != n) {
                throw std::invalid_argument("The permutation list does not match the qubit number.");
            }

            //Create a new vector of indices to map old qubit positions to new positions
            std::vector<int> newIndex(2 * n + 1);
            for (int i = 0; i < n; i++) {
                newIndex[i] = indices[i];
                newIndex[n + i] = indices[i] + n;
            }
            newIndex[2 * n] = 2 * n;

            //Create a new tableau with permuted columns based on the new indices
            std::vector<std::vector<int>> permutedTableau(matrix.size(), std::vector<int>(matrix[0].size(), 0));

            for (int row = 0; row < matrix.size(); ++row) {
                for (int col = 0; col < matrix[row].size(); ++col) {
                    permutedTableau[row][newIndex[col]] = matrix[row][col];
                }
            }

            //Replace the original tableau with the permuted tableau
            matrix = permutedTableau;
        }

        void reduction(std::vector<std::vector<int>>& graphMatrix) {
            int xRank = nonZeroRows(extractHalf(graphMatrix, true));

            int rowPointer = 0;

            for (int col = 0; col < n * 2; col++) {
                std::vector<int> ones = searchNonZeroRows(graphMatrix, col, rowPointer);

                if (ones.empty()) {
                    //No rows have a "1" in this column position, skip it.
                    continue;
                } 
                else 
                {
                    reduceDown(graphMatrix, ones);
                    rowInsert(graphMatrix, ones[0], rowPointer);
                    rowPointer++;

                    std::vector<int> onesUp = searchNonZeroRows(graphMatrix, col, 0, rowPointer);
                    // std::cout << "\nIs empty: " << onesUp.empty() << std::endl;
                    
                    reduceUp(graphMatrix, onesUp);

                    if (rowPointer == graphMatrix.size()) 
                    {
                        break;
                    }
                }
            }
        }

        void zReduction(std::vector<std::vector<int>>& graphMatrix) {
            int xRank = nonZeroRows(extractHalf(graphMatrix, true));

            int rowPointer = xRank;

            for (int col = n + xRank; col < n * 2; col++) {
                std::vector<int> ones = searchNonZeroRows(graphMatrix, col, rowPointer);

                if (ones.empty()) {
                    //No rows have a "1" in this column position, skip it.
                    continue;
                }
                else 
                {
                    reduceDown(graphMatrix, ones);
                    rowInsert(graphMatrix, ones[0], rowPointer);
                    rowPointer++;

                    std::vector<int> onesUp = searchNonZeroRows(graphMatrix, col, 0, rowPointer);
                    if(!onesUp.empty())
                        reduceUp(graphMatrix, onesUp);

                    if (rowPointer == graphMatrix.size()) 
                    {
                        break;
                    }
                }
            }

            //Adjust qubit index to maximize the size of the identity block in the Z part
            std::vector<int> qubitIndex(n, 0);
            for (int i = 0; i < xRank; i++) {
                qubitIndex[i] = i;
            }

            std::vector<std::vector<int>> zMatrix = extractHalf(graphMatrix, false);

            std::set<int> filledSet;
            for(int i = 0; i < xRank; i++)
            {
                filledSet.insert(i);
            }
            for (int j = xRank; j < stabCounts; j++) {
                std::vector<int> col = zMatrix[j];
                if(col.size() == 0)
                {
                    j--;
                    break;
                }
                else
                {
                    qubitIndex[col[0]+xRank] = j;
                    filledSet.insert(col[0]+xRank);
                }
            }

            int j = xRank;
            for (int i = xRank; i < n; i++) {
                if (filledSet.find(i) == filledSet.end()) {
                    j++;
                    qubitIndex[i] = j;
                }
            }
            qubitReorder(graphMatrix, qubitIndex);
        }

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
                        //Entry -- swap x and z bits
                        tempVal = x[i][a];
                        x[i][a] = z[i][a];
                        z[i][a] = tempVal; 
                    } 
                }
                else if (gate.op_name == OP::S)
                {
                    int a = gate.qubit;

                    for(int i = 0; i < rows-1; i++)
                    {
                        //Phase
                        r[i] ^= (x[i][a] & z[i][a]);

                        //Entry
                        z[i][a] ^= x[i][a];
                    }

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