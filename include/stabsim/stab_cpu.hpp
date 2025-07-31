#pragma once

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
    class STAB_CPU : public QuantumState
    {

    public:
        //Default identity constructor
        bool has_destabilizers;

        STAB_CPU(IdxType _n_qubits) : QuantumState(SimType::STAB)
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
            
            std::random_device rd;
            rng.seed(rd());
            dist = std::uniform_int_distribution<int>(0,1);
            random_float = std::uniform_real_distribution<double>(0.0, 1.0);


            has_destabilizers = true;

            SAFE_ALOC_HOST(totalResults, sizeof(IdxType));
            memset(totalResults, 0, sizeof(IdxType));
        }

        ~STAB_CPU()
        {
            // Release for CPU side
            SAFE_FREE_HOST(totalResults);
        }

        IdxType get_qubits() override
        {
            return n;
        }

        //Delete all stabilizers and leave an empty tableau object
        void delete_all_rows() override
        {
            rows = 0;
            stabCounts = 0;
            x.resize(rows, std::vector<int>(cols,0));
            z.resize(rows, std::vector<int>(cols,0)); 
            r.resize(rows, 0);   
        }

        //Remove only the destabilizers for cases where full states aren't necessary
        void remove_destabilizers() override
        {
            has_destabilizers = false;

            if(rows > 0)
            {
                rows = (rows/2);
                x.erase(x.begin(), x.begin()+rows);
                x.pop_back();
                z.erase(z.begin(), z.begin()+rows);
                z.pop_back();
                r.erase(r.begin(), r.begin()+rows);
                r.pop_back();
            }
        }

        int get_num_rows() override
        {
            return rows;
        }

        //resets the tableau to a full identity
        void reset_state() override
        {
            rows = 2*n+1;
            cols = n;
            stabCounts = n;
            for (auto& row : x) {
                std::fill(row.begin(), row.end(), 0);
            }
            for (auto& row : z) {
                std::fill(row.begin(), row.end(), 0);
            }
            std::fill(r.begin(), r.end(), 0);
            //The 2n+1 th row is scratch space

            //Intialize the identity tableau
            for(int i = 0; i < n; i++)
            {
                x[i][i] = 1;
                z[i+n][i] = 1;
            }
            
            std::random_device rd;
            rng.seed(rd());
            dist = std::uniform_int_distribution<int>(0,1);
            random_float = std::uniform_real_distribution<double>(0.0, 1.0);


            has_destabilizers = true;

            memset(totalResults, 0, sizeof(IdxType));
        }

        //Apply gates to a stabilizer-only tableau
        void apply_gate(std::string gate, int a, int b = -1) override
        {
            int tempVal;
            if(gate == "H")
            {
                //std::cout << "APPLYING H(" << a << ")" << std::endl;
                for(int i = 0; i < rows; i++)
                {
                    //Phase
                    r[i] ^= (x[i][a] & z[i][a]);
                    //Entry -- swap x and z bits
                    tempVal = x[i][a];
                    x[i][a] = z[i][a];
                    z[i][a] = tempVal; 
                }
                //std::cout << "After H --------" << std::endl;
               // print_res_state();
            }
            else if(gate == "S")
            {
                for(int i = 0; i < rows; i++)
                {
                    //Phase
                    r[i] ^= (x[i][a] & z[i][a]);

                    //Entry
                    z[i][a] ^= x[i][a];
                }
            }
            else if(gate == "SDG")
            {
                for(int i = 0; i < rows; i++)
                {
                    r[i] ^= x[i][a] ^ (x[i][a] & z[i][a]);

                    //Entry
                    z[i][a] ^= x[i][a];
                }
            }
                            
            else if(gate == "CX")
            {
                //std::cout << "APPLYING CX(" << a << ", " << b << ")" << std::endl;
                for(int i = 0; i < rows; i++)
                {
                    //Phase
                    r[i] ^= ((x[i][a] & z[i][b]) & (x[i][b]^z[i][a]^1));

                    //Entry
                    x[i][b] ^= x[i][a];
                    z[i][a] ^= z[i][b];
                }
                //std::cout << "After CX --------" << std::endl;
                //print_res_state();
            }
            else
            {
                std::cout << "Non-Clifford Gate in apply_gate!" << std::endl;
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
                if(has_destabilizers && ((i == (rows/2)) || (i == rows-1)))
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

        // void initialize_erasure_mask() override
        // {
        //     x_erasure.resize(rows, std::vector<int>(cols,0)); //first 2n+1 x n block. first n represents destabilizers
        //                                                //second n represents stabilizers + 1 extra row
        //     z_erasure.resize(rows, std::vector<int>(cols,0)); //second 2n+1 x n block to form the 2n+1 x 2n sized tableau
        //     r_erasure.resize(rows, 0); //column on the right with 2n+1 rows
        //     //The 2n+1 th row is scratch space

        //     //Intialize the identity tableau
        //     for(int i = 0; i < n; i++)
        //     {
        //         x_erasure[i][i] = 1;
        //         z_erasure[i+n][i] = 1;
        //     }
        // }

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

#ifdef EIGEN
        //Creates a sparse density matrix out of Pauli stabilizers
        ComplexMatrix get_density_matrix() override
        {
            std::vector<std::pair<std::string, int>> full_paulis = get_all_lines();
            int dim = 1 << n;

            ComplexMatrix rho = createSparseIdentity(dim);  //Start with identity

            for (const auto& line : full_paulis) {
                const std::string& pauli_string = line.first;
                int sign = line.second;

                ComplexMatrix S = Matrix::Ones(1, 1);  //Start with scalar 1

                
                //Build the stabilizer operator
                for (char p : pauli_string) {
                    if (p == 'I')      S = kroneckerProduct(S, pauliI()).eval();
                    else if (p == 'X') S = kroneckerProduct(S, pauliX()).eval();
                    else if (p == 'Y') S = kroneckerProduct(S, pauliY()).eval();
                    else if (p == 'Z') S = kroneckerProduct(S, pauliZ()).eval();
                    else throw std::logic_error("Invalid Pauli character in stabilizer");
                }

                //Flip sign if stabilizer has -1 eigenvalue
                if (sign == 1)
                    S *= -1.0;

                //Projector: (I + S) / 2
                ComplexMatrix proj = (createIdentity(dim) + S) * 0.5;

                //Multiply into density matrix
                rho = (rho * proj);
                // std::cout << "Matrix has dimensions: " 
                //     << rho.rows() << " x " << rho.cols() << std::endl;
            }

            // Convert to dense for return
            return ComplexMatrix(rho);
        }   
#endif


        //Get the stabilizers in the tableau, i.e. all of the Pauli strings that stabilize a certain circuit
        std::vector<std::string> get_stabilizers() override
        {
            int x_val;
            int z_val;
            std::vector<std::string> pauliStrings;
            int start;
            int end;
            if(has_destabilizers)
            {
                start = (rows/2);
                end = rows-1;
            }
            else
            {
                start = 0;
                end = rows;
            }
            for(int i = start; i < end; i++) //rows of stabilizers
            {
                std::string stabilizers;
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

        //Returns a map of stabilizers and the number of times they ocurr in the tableau
        void stabilizer_count(std::unordered_map<std::string, std::pair<int, int>>& stab_counts) override
        {
            for(int i = 0; i < rows; i++)
            {
                std::pair<std::string, int> stab = get_stabilizer_line(i);
                std::string stabilizer = stab.first;
                stab_counts[stabilizer].first++;

                //Keep track of how many +/- rotations for later
                if(stab.second)
                    stab_counts[stabilizer].second--;
                else
                    stab_counts[stabilizer].second++;
            }
            //Go back through and remove excess stabilizers, then add any odd number back
            for(const auto& [stabilizer, count] : stab_counts)
            {
                if(count.first > 1)
                {
                    // std::cout << "Stabilizer: " << stabilizer << std::endl;
                    // std::cout << "Count.first: " << count.first << std::endl;
                    // std::cout << "Count.second: " << count.first << std::endl;

                    for(int j = 0; j < rows; j++)
                    {                    
                        //Remove all the stabilizers that repeat
                        if(rows > 0)
                        {
                            if(get_stabilizer_line(j).first == stabilizer)
                            {
                                // std::cout << "Removed: " << stabilizer << " " << get_stabilizer_line(j).second << std::endl;
                                remove_stabilizer(j);
                                j--;
                            }
                        }
                    }
                    //Add back if there were an odd number of rotations left over
                    //T seperation will take care of the rotations
                    if((count.second % 2) == 1)
                    {
                        add_stabilizer(stabilizer, 0);
                        // std::cout << "Added: " << stabilizer << " 0" << std::endl;
                    }
                    else if((count.second % 2) == -1)
                    {
                        add_stabilizer(stabilizer, 1);
                        // std::cout << "Added: " << stabilizer << " 1" << std::endl;
                    }
                }
                // std::cout << std::endl;
            }
            // std::cout << "--- End current state in stab count ---" << std::endl;
        }

        //Check every stabilizer of the tableau
        //If one of the stabilizers doesn't commute with the provide pauli string, return false
        bool check_commutation(std::string& pauliString) override
        {
            int new_x;
            int new_z;
            int start = 0;
            
            // if(has_destabilizers)
            //     start = rows/2;

            for (int i = start; i < rows; i++) 
            {
                //Compute inner product
                int product = 0;
                
                //Loop over every column
                for(int j = 0; j < cols; j++) 
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
                    product ^= ((x[i][j] & new_z) ^ (z[i][j] & new_x));
                }
                if(product > 0) {
                    return false; //Anti-commutation in the Pauli string at row i
                }
            }
            return true; //Commutes with all stabilizers in the tableau
        }
        
        //True if pauliString and target row stabilizer commute
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
                product ^= ((x[row][col] & temp_z) ^ (z[row][col] & temp_x));
            }
            if(product > 0)
                return false;
            return true;
        }

        void add_stabilizer(std::string pauliString, int phase_bit = 0) override
        {
            assert(pauliString.length() == n);

            //Full stabilizer/destabilizer tableau
            if(has_destabilizers)
            {
                std::cout << "Shouldn't be here for T separation case" << std::endl;
                //Start by adding a row of destabilizers and stabilizers to T
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
                            break;
                        case 'X':
                            x[rows-2][i] = 1;
                            z[rows-2][i] = 0;
                            break;
                        case 'Y':   
                            x[rows-2][i] = 1;
                            z[rows-2][i] = 1;
                            break;
                        case 'Z':
                            x[rows-2][i] = 0;
                            z[rows-2][i] = 1;
                            break;
                        default:
                            std::logic_error("Invalid stabilizer");
                            break;
                    }
                }
            }

            //Only stabilizer tableau
            else
            {
                //Start by adding a row of stabilizers to T
                stabCounts++;
                rows++;
                x.push_back(std::vector<int>(cols,0));
                z.push_back(std::vector<int>(cols,0));
                r.push_back(phase_bit);
                
                //Stabilizer only addition (no temp row)
                for(int i = 0; i < cols; i++)
                {
                    switch(pauliString[i])
                    {
                        case 'I':
                            x[rows-1][i] = 0;
                            z[rows-1][i] = 0;
                            break;
                        case 'X':
                            x[rows-1][i] = 1;
                            z[rows-1][i] = 0;
                            break;
                        case 'Y':   
                            x[rows-1][i] = 1;
                            z[rows-1][i] = 1;
                            r[rows-1] = 1;
                            break;
                        case 'Z':
                            x[rows-1][i] = 0;
                            z[rows-1][i] = 1;
                            break;
                        default:
                            std::logic_error("Invalid Pauli");
                            break;
                    }
                }
            }
        }

        void add_stabilizer_bits(std::vector<int> new_x, std::vector<int> new_z, int phase_bit = 0) override
        {
            stabCounts++;
            rows++;
            x.push_back(new_x);
            z.push_back(new_z);
            r.push_back(phase_bit);
        }

        //Replaces a stabilizer pauli string at some row in the Tableau. Useful for initializing a
        //new Tableau in a for loop without circuit initialization
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

        void remove_stabilizer(int row_index) override
        {
            x.erase(x.begin()+row_index);
            z.erase(z.begin()+row_index);
            r.erase(r.begin()+row_index);
            stabCounts--;
            rows--;
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
            sum += 2*r[i] + 2*r[h];

            if(sum % 4 == 0)
                r[h] = 0;
            else
                r[h] = 1;

            // std::cout << "r[" << h << "]" <<  " after if = " <<  r[h] << std::endl;


        } //End rowsum

        void i_rowsum(int h, int i) override
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
            sum += 1 + 2*r[h] + 2*r[i];

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
            
            for(int shot = 0; shot < shots; shot++)
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

                        totalResults[shot] |= (temp_r[p] << a);
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
        void sim(std::shared_ptr<NWQSim::Circuit> circuit, double &sim_time) override
        {
            std::vector<Gate> gates = circuit->get_gates();
            IdxType n_gates = gates.size();

            //Check that the circuit object matches the tableau
            assert((circuit->num_qubits() == n) && "Circuit does not match the qubit number of the state!");

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
        IdxType n;
        int stabCounts;
        IdxType rows;
        IdxType cols;
        std::vector<std::vector<int>> x;
        std::vector<std::vector<int>> z;
        std::vector<int> r;
        std::vector<std::vector<int>> x_erasure;
        std::vector<std::vector<int>> z_erasure;
        std::vector<int> r_erasure;
        IdxType* totalResults = NULL;
        IdxType** totalResultsLong = NULL;

        double p = NULL;

        std::mt19937 rng;
        std::uniform_int_distribution<int> dist;
        std::uniform_real_distribution<double> random_float;

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

        int damping_generator(double p, double gamma)
        {   //Remove negative and renormalize
            //Range [0.0, 1.0)
            double tolerance = 1e-9;
            double monte = random_float(rng);
            double c0 = 0.5 * (1.0-gamma+sqrt(1.0-gamma));
            double c1 = 0.5 * (1.0-gamma-sqrt(1.0-gamma));
            double c2 = gamma;
            // double total = c0 + c1 + c2;
            // std::cout << "\n C0 " << c0;
            // std::cout << " C1 " << c1;
            // std::cout << " C2 " << c2;
            // double

            // if(c0 < 0){
            //     std::cerr << "Assertion failed: c0 = " << c0 << std::endl;                
            //     std::cerr << "p = " << p << ", gamma = " << gamma << std::endl;
            //     assert(false);
            // }
            // if(c1 < 0){
            //     std::cerr << "c1 = " << c1 << std::endl;                
            //     std::cerr << "p = " << p << ", gamma = " << gamma << std::endl;
            //     assert(false);
            // }

            // if (std::abs(total - 1.0) > tolerance) {
            //     std::cerr << "Assertion failed: c0 + c1 + c2 = " << total << " (expected ~1)\n";
            //     std::cerr << "p = " << p << ", gamma = " << gamma << std::endl;
            //     assert(false);
            // }

            double c0_tot = (1 - p) * c0 + p * c1;
            double c1_tot = p * c0 + (1 - p) * c1;
            double c2_tot = c2;

            // if(c1_tot < 0)
            // {
            //     double tot = c0_tot + c2_tot;
            //     c1_tot = 0;
            //     c0_tot = c0_tot/tot;
            //     c2_tot = c2_tot/tot;
            // }
           

            // if (std::abs(tot_total - 1.0) > tolerance) {
            //     std::cerr << "Assertion failed: c0_tot + c1_tot + c2_tot = " << tot_total << " (expected ~1)\n";
            //     std::cerr << "p = " << p << ", gamma = " << gamma << std::endl;
            //     assert(false);
            // }
            // if(c0_tot < 0){
            //     std::cerr << "Assertion failed: c0_tot = " << c0_tot << std::endl;
            //     std::cerr << "c1_tot = " << c1_tot << std::endl;    
            //     std::cerr << "c0 = " << c0 << std::endl;          
            //     std::cerr << "c1 = " << c1 << std::endl;          
            //     std::cerr << "p = " << p << ", gamma = " << gamma << std::endl;
            //     assert(false);
            // }
            // if(c1_tot < 0){
            //     std::cerr << "Assertion failed: c1_tot = " << c1_tot << std::endl;  
            //     std::cerr << "c0_tot = " << c0_tot << std::endl;               
            //     std::cerr << "c0 = " << c0 << std::endl;          
            //     std::cerr << "c1 = " << c1 << std::endl;          
            //     std::cerr << "p = " << p << ", gamma = " << gamma << std::endl;
            //     assert(false);
            // }

            //Remormalization if C1~ < 0

            // if(c0_tot < 0){
            //     c0_tot = 0;
            //     double sum = c1_tot+c2_tot;
            //     c1_tot = c1_tot/sum;
            //     c2_tot = c2_tot/sum;
            // }
            // if(c1_tot < 0){
            //     c2_tot = 0;
            //     double sum = c0_tot+c2_tot;
            //     c0_tot = c0_tot/sum;
            //     c2_tot = c2_tot/sum;
            // }
            // else if(c2_tot < 0){
            //     c2_tot = 0;
            //     double sum = c1_tot+c0_tot;
            //     c1_tot = c1_tot/sum;
            //     c0_tot = c0_tot/sum;
            // }

            // std::cout << " Monte: " << monte << std::endl;
            
            // std::cout << "c: " << c0 << " " << c1 << " " << c2 << std::endl;
            // std::cout << "c~: " << c0_tot << " " << c1_tot << " " << c2_tot << std::endl;

            double temp = monte - c0_tot;
            if(temp < 0) 
            {
                // std::cout << " C0_tot chosen: " << c0_tot;
                return 0;
            }
            
            temp = temp - c1_tot;
            if(temp < 0) 
            {
                // std::cout << " C1_tot chosen: " << c1_tot;
                return 1;
            }
            temp = temp - c2_tot;
            if(temp < 0) 
            {
                // std::cout << " C2_tot chosen: " << c2_tot;
                return 2;
            }
            else 
            {
                std::cerr << "Impossible monte" << std::endl;
                return 3;
            }
        }

        void reset_routine(int a)
        {
            int p = -1;
            int temp_result = 0;
            int half_rows = rows/2;

            for(int p_index = half_rows; p_index < rows-1; p_index++)
            {  
                //std::cout << "x at [" << p_index << "][" << a << "] = " << x[p_index][a] << std::endl;
                if(x[p_index][a])
                {
                    p = p_index;
                    break;
                }
            }
            // std::cout << "p = " << p << std::endl;
            //A p such that x[p][a] = 1 exists
            //Random
            if(p > -1)
            {
                for(int i = 0; i < rows-1; i++)
                {
                    // std::cout << "x = " << x[i][a] << std::endl;
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

                temp_result = r[p];
                // std::cout << "Random measurement at qubit " << a << " value: " << (r[p] << a) << std::endl;
            }
            //Deterministic
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
                        rowsum(rows-1, i+half_rows);
                    }
                }
                // std::cout << "Deterministc measurement at qubit " << a << " value: " << (r[rows-1] << a) << std::endl;
                temp_result = r[rows-1];
            }
            measurement_results.push_back(temp_result);
            if(temp_result == 1) //Apply X to flip back to 0
            {
                for(int i = 0; i < rows-1; i++)
                {
                    // z[i][a] ^= x[i][a]; 
                    r[i] ^= z[i][a];

                    // r[i] ^= z[i][a];
                }
            }
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

                switch (gate.op_name)
                {
                    case OP::CX:
                    {
                        int a = gate.ctrl;
                        int b = gate.qubit;
                        for(int i = 0; i < rows-1; i++)
                        {
                            //Phase
                            r[i] = r[i]^ (x[i][a] & z[i][b]) & (x[i][b]^z[i][a]^1);

                            //Entry
                            x[i][b] ^= x[i][a];
                            z[i][a] ^= z[i][b];
                        }
                        break;
                    }

                    case OP::H: 
                        for(int i = 0; i < rows-1; i++)
                        {
                            //Phase
                            r[i] ^= (x[i][a] & z[i][a]);
                            //Entry -- swap x and z bits
                            tempVal = x[i][a];
                            x[i][a] = z[i][a];
                            z[i][a] = tempVal; 
                        }
                        break;
        
                    case OP::S:
                        for(int i = 0; i < rows-1; i++)
                        {
                            //Phase
                            r[i] ^= (x[i][a] & z[i][a]);

                            //Entry
                            z[i][a] ^= x[i][a];
                        }
                        break;

                    case OP::SDG:
                        for(int i = 0; i < rows-1; i++)
                        {
                            //Phase -- Equal to Z S or x & !z
                            r[i] ^= x[i][a] ^ (x[i][a] & z[i][a]);

                            //Entry
                            z[i][a] ^= x[i][a];
                        }
                        break;

                    case OP::DAMP:
                        // std::cout << "Gamma " << gate.gamma << std::endl;
                        switch(damping_generator(gate.lam, gate.gamma))
                        {
                            case 0: //Do nothing
                                break;
                            case 1: //Apply Z
                                for(int i = 0; i < rows-1; i++)
                                {
                                    //Phase
                                    r[i] ^= x[i][a];
                                }
                                break;
                            case 2: //Reset to |0>
                                reset_routine(a);
                                break;
                            default:
                                std::logic_error("Invalid damping result");
                                exit(1);
                        }
                        break;
                    
                    case OP::RX:
                        //H SDG H
                        if(gate.theta == PI/2)
                        {
                            for(int i = 0; i < rows-1; i++)
                            {
                                //Phase
                                r[i] ^= (x[i][a] & z[i][a]);
                                //Entry -- swap x and z bits
                                tempVal = x[i][a];
                                x[i][a] = z[i][a];
                                z[i][a] = tempVal; 

                                //Phase -- Equal to Z S or x & !z
                                r[i] ^= x[i][a] ^ (x[i][a] & z[i][a]);

                                //Entry
                                z[i][a] ^= x[i][a];

                                r[i] ^= (x[i][a] & z[i][a]);
                                //Entry -- swap x and z bits
                                tempVal = x[i][a];
                                x[i][a] = z[i][a];
                                z[i][a] = tempVal; 
                            }
                        }
                        //H S
                        else if(gate.theta == -PI/2)
                        {
                            for(int i = 0; i < rows-1; i++)
                            {
                                //Phase
                                r[i] ^= (x[i][a] & z[i][a]);
                                //Entry -- swap x and z bits
                                tempVal = x[i][a];
                                x[i][a] = z[i][a];
                                z[i][a] = tempVal; 

                                //S
                                //Phase
                                r[i] ^= (x[i][a] & z[i][a]);

                                //Entry
                                z[i][a] ^= x[i][a];

                                //Phase
                                r[i] ^= (x[i][a] & z[i][a]);
                                //Entry -- swap x and z bits
                                tempVal = x[i][a];
                                x[i][a] = z[i][a];
                                z[i][a] = tempVal; 
                            }
                        }
                        else
                        {
                            std::cout << "Non-Clifford angle in RX gate! "
                                    << OP_NAMES[gate.op_name] << "(" << gate.theta << ")" << std::endl;
                            std::logic_error("Invalid gate type");
                            exit(1);
                        }
                        break;
                    
                    case OP::RY:
                        //H X -- X : r[i] ^= z[i][a]
                        if(gate.theta == PI/2)
                        {
                            for(int i = 0; i < rows-1; i++)
                            {
                                //Phase -- Equal to Z S or x & !z
                                r[i] ^= x[i][a] ^ (x[i][a] & z[i][a]);

                                //Entry
                                z[i][a] ^= x[i][a];

                                //Phase
                                r[i] ^= (x[i][a] & z[i][a]);
                                //Entry -- swap x and z bits
                                tempVal = x[i][a];
                                x[i][a] = z[i][a];
                                z[i][a] = tempVal; 

                                //Phase
                                r[i] ^= (x[i][a] & z[i][a]);

                                //Entry
                                z[i][a] ^= x[i][a];
                            }
                        }
                        //X H
                        else if(gate.theta == -PI/2)
                        {
                            for(int i = 0; i < rows-1; i++)
                            {
                                //Phase -- Equal to Z S or x & !z
                                r[i] ^= (x[i][a] & z[i][a]);

                                //Entry
                                z[i][a] ^= x[i][a];

                                //Phase
                                r[i] ^= (x[i][a] & z[i][a]);
                                //Entry -- swap x and z bits
                                tempVal = x[i][a];
                                x[i][a] = z[i][a];
                                z[i][a] = tempVal; 

                                //Phase
                                r[i] ^= x[i][a] ^ (x[i][a] & z[i][a]);

                                //Entry
                                z[i][a] ^= x[i][a];
                            }
                        }
                        else
                        {
                            std::cout << "Non-Clifford angle in RY gate! "
                                    << OP_NAMES[gate.op_name] << "(" << gate.phi << ")" << std::endl;
                            std::logic_error("Invalid gate type");
                            exit(1);
                        }
                        break;
                    case OP::CY:
                    {
                        int a = gate.ctrl;
                        int b = gate.qubit;
                        for(int i = 0; i < rows-1; i++)
                        {
                            //Phase -- Equal to Z S or x & !z
                            r[i] ^= x[i][b] ^ (x[i][b] & z[i][b]);

                            //Entry
                            z[i][b] ^= x[i][b];
                        }
                        for(int i = 0; i < rows-1; i++)
                        {
                            //Phase
                            r[i] ^= ((x[i][a] & z[i][b]) & (x[i][b]^z[i][a]^1));

                            //Entry
                            x[i][b] ^= x[i][a];
                            z[i][a] ^= z[i][b];
                        }
                        for(int i = 0; i < rows-1; i++)
                        {
                            //Phase
                            r[i] ^= (x[i][b] & z[i][b]);

                            //Entry
                            z[i][b] ^= x[i][b];
                        }
                        break;
                    }
                    case OP::CZ:
                    {
                        int a = gate.ctrl;
                        int b = gate.qubit;

                        for(int i = 0; i < rows-1; i++)
                        {
                            //Phase
                            r[i] ^= (x[i][b] & z[i][b]);
                            //Entry -- swap x and z bits
                            tempVal = x[i][b];
                            x[i][b] = z[i][b];
                            z[i][b] = tempVal; 
                        }
                        for(int i = 0; i < rows-1; i++)
                        {
                            //Phase
                            r[i] ^= ((x[i][a] & z[i][b]) & (x[i][b]^z[i][a]^1));

                            //Entry
                            x[i][b] ^= x[i][a];
                            z[i][a] ^= z[i][b];
                        }
                        for(int i = 0; i < rows-1; i++)
                        {
                            //Phase
                            r[i] ^= (x[i][b] & z[i][b]);
                            //Entry -- swap x and z bits
                            tempVal = x[i][b];
                            x[i][b] = z[i][b];
                            z[i][b] = tempVal; 
                        }
                        break;
                    }

                    case OP::M:
                    {
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
                        // std::cout << "p = " << p << std::endl;
                        //A p such that x[p][a] = 1 exists
                        //Random
                        if(p > -1)
                        {
                            for(int i = 0; i < rows-1; i++)
                            {
                                // std::cout << "x = " << x[i][a] << std::endl;
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
                        //Deterministic
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
                                    rowsum(rows-1, i+half_rows);
                                }
                            }
                            // std::cout << "Deterministc measurement at qubit " << a << " value: " << (r[rows-1] << a) << std::endl;
                            totalResults[0] += r[rows-1] << a;
                        }
                        break;
                    }
                    case OP::RESET:
                    {
                        reset_routine(a);
                        break;
                    }

                    case OP::X:
                        //equiv to H S S H or H Z H
                        for(int i = 0; i < rows-1; i++)
                        {
                            r[i] ^= z[i][a];
                        }
                        break;
        
                    case OP::Y:
                        //equiv to Z X
                        for(int i = 0; i < rows-1; i++)
                        {
                            r[i] = r[i] ^ x[i][a] ^ z[i][a];
                        }
                        break;

                    case OP::Z:
                        //equiv to S S gates
                        for(int i = 0; i < rows-1; i++)
                        {
                            //Phase
                            r[i] ^= x[i][a];
                        }
                        break;
                    

                        
                    case OP::MA:
                        measure_all();
                        break;
                    
                    default:
                        std::cout << "Non-Clifford or unrecognized gate: "
                                    << OP_NAMES[gate.op_name] << std::endl;
                        std::logic_error("Invalid gate type");
                        exit(1);
                }//End switch
            } //End gates for loop
        } //End simulate
    }; //End tableau class
} //namespace NWQSim

//#endif