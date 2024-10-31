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

namespace NWQSim
{
    class STAB_CPU : public QuantumState
    {

    public:

        //Default identity constructor
        STAB_CPU(uint _n_qubits) : QuantumState(SimType::STAB)
        {
            g = 0;
            n = _n_qubits;
            rows = 2*n+1;
            stabCounts = n;
            cols = n;

            outcomes.resize(cols);
            outcomes.assign(outcomes.size(), 0);
            x.resize(rows, std::vector<uint>(cols,0)); //first 2n+1 x n block. first n represents destabilizers
                                                       //second n represents stabilizers + 1 extra row
            z.resize(rows, std::vector<uint>(cols,0)); //second 2n+1 x n block to form the 2n+1 x 2n sized tableau
            r.resize(rows, 0); //column on the right with 2n+1 rows
            //The 2n+1 th row is scratch space

            //Intialize the identity tableau
            for(int i = 0; i < n; i++)
            {
                x[i][i] = 1;
                z[i+n][i] = 1;
            }
            simulate();
        }

        //Constructor with a full set of gates
        STAB_CPU(std::vector<Gate> &_gates, uint _n_qubits) : QuantumState(SimType::STAB)
        {
            gates = _gates;
            g = gates.size();
            n = _n_qubits;
            rows = 2*n+1;
            stabCounts = n;
            cols = n;

            outcomes.resize(cols);
            outcomes.assign(outcomes.size(), 0);
            x.resize(rows, std::vector<uint>(cols,0)); //first 2n+1 x n block. first n represents destabilizers
                                                       //second n represents stabilizers + 1 extra row
            z.resize(rows, std::vector<uint>(cols,0)); //second 2n+1 x n block to form the 2n+1 x 2n sized tableau
            r.resize(rows, 0); //column on the right with 2n+1 rows
            //The 2n+1 th row is scratch space

            //Intialize the identity tableau
            for(int i = 0; i < n; i++)
            {
                x[i][i] = 1;
                z[i+n][i] = 1;
            }
            simulate();
        }

        //Constructor from a circuit containing gates
        STAB_CPU(std::shared_ptr<Circuit>& circuit) : QuantumState(SimType::STAB)
        {
            gates = circuit->get_gates();
            g = gates.size();
            n = circuit->num_qubits();
            rows = 2*n+1;
            stabCounts = n;
            cols = n;

            outcomes.resize(cols);
            outcomes.assign(outcomes.size(), 0);
            x.resize(rows, std::vector<uint>(cols,0)); //first 2n+1 x n block. first n represents destabilizers
                                                     //second n represents stabilizers + 1 extra row
            z.resize(rows, std::vector<uint>(cols,0)); //second 2n+1 x n block to form the 2n+1 x 2n sized tableau
            r.resize(rows, 0); //column on the right with 2n+1 rows
            //The 2n+1 th row is scratch space

            //Intialize the identity tableau
            for(int i = 0; i < n; i++)
            {
                x[i][i] = 1;
                z[i+n][i] = 1;
            }
            simulate();
        }  

        //Prints the tableau including phase
        void print_table()
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
        void print_table(std::vector<std::vector<uint>>& M)
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

        //Add new gates to the tableau and simulate their evolution
        void add_gates(std::shared_ptr<Circuit>& new_circ)
        {
            gates = new_circ->get_gates();
            g = gates.size();
            simulate();
        }

        //Convert a vector of 1's and 0's to an IdxType decimal number
        IdxType *get_results()
        {
            int conversion = 0;

            for (int i = 0; i < outcomes.size(); i++) {
                conversion = conversion | outcomes[i] << i;  //Left shift and add the current bit
            } 
            
            std::cout << "Conversion: " << conversion << std::endl;

            result = new long long(static_cast<long long>(conversion));

            return result;
        }

        //Get the stabilizers in the tableau, i.e. all of the Pauli strings that stabilize the circuit
        std::vector<std::string> get_destabilizers()
        {
            int x_val;
            int z_val;
            std::string stabilizers;
            std::vector<std::string> pauliStrings;
            for(int i = 0; i < rows/2; i++) //rows of stabilizers
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
        std::pair<std::string,int> get_stabilizier_line(int row)
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
            for(int row = rows/2; row < rows-1; row++) //rows of stabilizers
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
        ComplexSparseMatrix get_DM()
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
                switch(line.second)
                {
                    case 0:
                        S = S * std::complex(1.0, 0.0);
                        break;
                    case 1:
                        S = S * std::complex(0.0, 1.0);
                        break;
                    case 2:   
                        S = S * std::complex(-1.0, 0.0);
                        break;
                    case 3:
                        S = S * std::complex(0.0, -1.0);
                        break;
                    default:
                        std::logic_error("Invalid phase");
                        break;
                }
                DM = DM * .5 * (I + S);
            }
            return DM;
        }

        //Get the stabilizers in the tableau, i.e. all of the Pauli strings that stabilize a certain circuit
        std::vector<std::string> get_stabilizers()
        {
            int x_val;
            int z_val;
            std::string stabilizers;
            std::vector<std::string> pauliStrings;
            for(int i = rows/2; i < rows-1; i++) //rows of stabilizers
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

        bool check_commutation(std::string pauliString)
        {
            int new_x;
            int new_z;
            for (int i = 0; i < x.size(); i++) 
            {
                //Compute inner product
                int product = 0;
                
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

        //Takes a default (or any) tableau and sets its stabilizers according to
        //a Pauli string provided
        void add_stabilizer(std::string pauliString)
        {
            assert(pauliString.length() == n);
            if(check_commutation(pauliString))
            {
                //Start by adding a row of destabilizers and stabilizers to T
                stabCounts++;
                rows+=2;
                x.insert(x.begin() + x.size()/2, std::vector<uint>(cols,0));
                x.insert(x.end()-1, std::vector<uint>(cols,0));
                z.insert(z.begin() + z.size()/2, std::vector<uint>(cols,0));
                z.insert(z.end()-1, std::vector<uint>(cols,0));
                r.insert(r.begin() + r.size()/2, 0);
                r.insert(r.end()-1, 0);
                
                //Stabilizer and destabilizer addition
                for(int i = 0; i < pauliString.length(); i++)
                {
                    switch(pauliString[i])
                    {
                        case 'I':
                            x[rows-2][i] = 0;
                            z[rows-2][i] = 0;
                            x[rows/2-1][i] = 0;
                            z[rows/2-1][i] = 0;
                            break;
                        case 'X':
                            x[rows-2][i] = 1;
                            z[rows-2][i] = 0;
                            x[rows/2-1][i] = 0;
                            z[rows/2-1][i] = 1;
                            break;
                        case 'Y':   
                            x[rows-2][i] = 1;
                            z[rows-2][i] = 1;
                            //make the destabilizer X to anticommute with Y
                            x[rows/2-1][i] = 1;
                            z[rows/2-1][i] = 0;
                            //add an i to the 2 bit phase representation at the stabilizer row
                            r[rows-1] = (r[rows-1] + 1) % 4;
                            break;
                        case 'Z':
                            x[rows-2][i] = 0;
                            z[rows-2][i] = 1;
                            x[rows/2-1][i] = 1;
                            z[rows/2-1][i] = 0;
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

        //Replaces a stabilizer pauli string at some row in the Tableau. Useful for initializing a
        // new Tableau in a for loop without circuit initialization
        void replace_stabilizer(std::string pauliString, int stabPos)
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
                        x[rows/2 + stabPos][i] = 0;
                        z[rows/2 + stabPos][i] = 0;
                        //destab
                        x[stabPos][i] = 0;
                        z[stabPos][i] = 0;
                        break;
                    case 'X':
                        //stab
                        x[rows/2 + stabPos][i] = 1;
                        z[rows/2 + stabPos][i] = 0;
                        //destab
                        x[stabPos][i] = 0;
                        z[stabPos][i] = 1;
                        break;
                    case 'Y':
                        //stab
                        x[rows/2 + stabPos][i] = 1;
                        z[rows/2 + stabPos][i] = 1;
                        //make the destabilizer X to anticommute with Y
                        x[stabPos][i] = 1;
                        z[stabPos][i] = 0;
                        //add an i to the 2 bit phase representation at the stabilizer row
                        r[rows/2 + stabPos] = (r[rows/2 + stabPos] + 1) % 4;
                        break;
                    case 'Z':
                        x[rows/2 + stabPos][i] = 0;
                        z[rows/2 + stabPos][i] = 1;
                        x[stabPos][i] = 1;
                        z[stabPos][i] = 0;
                        break;
                    default:
                        std::logic_error("Invalid stabilizer");
                        break;
                }
            }
        }
        
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
        
        void transpose(){};
        void transpose(std::vector<std::vector<uint>>& M)
        {
            int rowSize = M.size();
            int colSize = M[0].size();
            std::vector<std::vector<uint>> MT = M;


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
            //Perform Gaussian elimination on the destabilizers (first half of the rows)
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
                assert((x[i][j] == 0 || x[i][j] == 1) && (z[i][j] == 0 || z[i][j] == 1));
                //Sum every column in the row
                if(x[i][j] == 1)
                {
                    if(z[i][j] == 1)
                        sum += z[h][j] - x[h][j];
                    else
                        sum += z[h][j] * (2*x[h][j]-1);
                }
                else if(z[i][j] == 1)
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

        //Provides a bit/shot measurement output without affecting the original tableau
        std::vector<int> measureShots(int shots = 2048)
        {
            std::vector<int> shotCounts(1<<n); // <single_outcome, # times>
            int singleOutcome = -1;

            for(int shot = 0; shot < shots; shot++)
            {
                //Make a copy of the class being measured so many shots can be performed
                Tableau temp = *this;
                for(int a = 0; a < n; a++)
                {  
                    //Store measurement outcomes of each bit
                    // int single_outcome;

                    int p = -1;
                    for(int p_index = temp.rows/2; p_index < temp.rows-1; p_index++)
                    {  
                        //std::cout << "x at [" << p_index << "][" << a << "] = " << x[p_index][a] << std::endl;
                        if(temp.x[p_index][a] != 0)
                        {
                            p = p_index;
                            break;
                        }
                    }
                    //A p such that x[p][a] = 1 exists
                    if(p > -1)
                    {
                        std::cout << "Random measurement ";

                        for(int i = 0; i < temp.rows-1; i++)
                        {
                            if((i != p) && (temp.x[i][a] == 1))
                            {
                                temp.rowsum(i, p);
                            }
                        }
                        temp.x[p-(temp.rows/2)] = temp.x[p];
                        temp.z[p-(temp.rows/2)] = temp.z[p];
                        //Change all the columns in row p to be 0
                        for(int i = 0; i < temp.n; i++)
                        {
                            temp.x[p][i] = 0;
                            temp.z[p][i] = 0;                        
                        }

                        std::random_device rd;
                        std::mt19937 gen(rd());  //Mersenne Twister engine

                        //Define the range for random numbers
                        std::uniform_int_distribution<> distr(0, 1);

                        //Generate and display a random number
                        int randomBit = distr(gen);
                        
                        if(randomBit)
                        {
                            //std::cout << "Random result of 1" << std::endl;
                            temp.r[p] = 1;
                        }
                        else
                        {
                            //std::cout << "Random result of 0" << std::endl;
                            temp.r[p] = 0;
                        }
                        temp.z[p][a] = 1;

                        temp.outcomes[a] = temp.r[p];

                        std::cout << temp.outcomes[a] << std::endl;                
                    }
                    else
                    {
                        std::cout << "Deterministic measurement ";

                        //Set the scratch space row to be 0
                        //i is the column indexer in this case
                        for(int i = 0; i < temp.n; i++)
                        {
                            temp.x[temp.rows-1][i] = 0;
                            temp.z[temp.rows-1][i] = 0;
                        }

                        //Run rowsum subroutine
                        for(int i = 0; i < temp.rows/2; i++)
                        {
                            if(temp.x[i][a] == 1)
                            {
                                //std::cout << "Perform rowsum at " << i << " + n" << std::endl;
                                temp.rowsum(temp.rows-1, i+(temp.rows/2));
                            }
                        }
                        temp.outcomes[a] = temp.r[rows-1];
                        std::cout << temp.outcomes[a] << std::endl;
                    }
                    //Convert the bit array outcomes[] to an int singleOutcome
                    singleOutcome += temp.outcomes[a] << a;
                //std::cout << "Result at qubit " << a << " = "  << outcomes[a] << std::endl;
                } //End M
                //Increment the count when an outcome is found
                shotCounts[singleOutcome]++;
            }//End shots
            return shotCounts;
        }

        //Simulate the gates from a circuit in the tableau
        void simulate()
        {
            //For swapping rows
            uint tempVal;

            for (int k = 0; k < g; k++)
            {
                auto gate = gates[k];
                //std::cout << "Current gate " << k << std::endl;
                int a = gate.qubit;

                if (gate.op_name == OP::H)
                {
                    
                    for(int i = 0; i < rows-1; i++)
                    {
                        //Phase
                        r[i] = r[i] ^ ((x[i][a] << 1) + z[i][a]);
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
                        r[i] = r[i] ^ ((x[i][a] << 1) + z[i][a]);

                        //Entry
                        z[i][a] = z[i][a] ^ x[i][a];
                    }

                }
                else if (gate.op_name == OP::CX)
                {  
                    int a = gate.ctrl;
                    int b = gate.qubit;
                    for(int i = 0; i < rows-1; i++)
                    {
                        //Phase
                        r[i] = r[i] ^ (((x[i][a] << 1) + z[i][b])*(x[i][b]^z[i][a]^1));

                        //Entry
                        x[i][b] = x[i][b] ^ x[i][a];
                        z[i][a] = z[i][a] ^ z[i][b];
                    }
                }
                else if (gate.op_name == OP::M)
                {  
                    int a = gate.qubit;
                    int p = -1;
                    for(int p_index = rows/2; p_index < rows-1; p_index++)
                    {  
                        //std::cout << "x at [" << p_index << "][" << a << "] = " << x[p_index][a] << std::endl;
                        if(x[p_index][a] != 0)
                        {
                            p = p_index;
                            break;
                        }
                    }
                    //A p such that x[p][a] = 1 exists
                    if(p > -1)
                    {
                        std::cout << "Random measurement ";

                        for(int i = 0; i < rows-1; i++)
                        {
                            if((i != p) && (x[i][a] == 1))
                            {
                                rowsum(i, p);
                            }
                        }
                        x[p-(rows/2)] = x[p];
                        z[p-(rows/2)] = z[p];
                        //Change all the columns in row p to be 0
                        for(int i = 0; i < n; i++)
                        {
                            x[p][i] = 0;
                            z[p][i] = 0;                        
                        }

                        std::random_device rd;
                        std::mt19937 gen(rd());  //Mersenne Twister engine

                        //Define the range for random numbers
                        std::uniform_int_distribution<> distr(0, 1);

                        //Generate and display a random number
                        int randomBit = distr(gen);
                        
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

                        outcomes[a] = r[p];

                        std::cout << outcomes[a] << std::endl;                
                    }
                    else
                    {
                        std::cout << "Deterministic measurement ";

                        //Set the scratch space row to be 0
                        //i is the column indexer in this case
                        for(int i = 0; i < n; i++)
                        {
                            x[rows-1][i] = 0;
                            z[rows-1][i] = 0;
                        }

                        //Run rowsum subroutine
                        for(int i = 0; i < rows/2; i++)
                        {
                            if(x[i][a] == 1)
                            {
                                //std::cout << "Perform rowsum at " << i << " + n" << std::endl;
                                rowsum(rows-1, i+(rows/2));
                            }
                        }
                        //The result is the 
                        outcomes[a] = r[rows-1];
                        std::cout << outcomes[a] << std::endl;
                    }
                //std::cout << "Result at qubit " << a << " = "  << outcomes[a] << std::endl;
                } //End M
                else    
                {
                    std::cout << "Non-Clifford or unrecognized gate: "
                                << OP_NAMES[gate.op_name] << std::endl;
                    std::logic_error("Invalid gate type");
                }
            } //End gates for loop
        } //End simulate

        //Converts a matrix into column echelon and counts the number of pivot rows
        int countPivots(std::vector<std::vector<uint>>& matrix) {
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

        int nonZeroRows(const std::vector<std::vector<uint>>& matrix) {
            int nonZeroRows = 0;
            for(const auto& row : matrix) 
            {
                bool isNonZero = false;
                for(uint value : row) 
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
        std::vector<int> searchNonZeroRows(const std::vector<std::vector<uint>>& matrix, int col, int startRow = 0, int endRow = -1) {
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

        std::vector<std::vector<uint>> extractHalf(const std::vector<std::vector<uint>>& matrix, bool left) {
            int numRows = matrix.size();
            int numCols = matrix[0].size();
            int halfCols = numCols / 2;

            std::vector<std::vector<uint>> result(numRows, std::vector<uint>(halfCols));

            for (int i = 0; i < numRows; i++) {
                for (int j = 0; j < halfCols; j++) {
                    result[i][j] = left ? matrix[i][j] : matrix[i][j + halfCols];
                }
            }
            return result;
        }

        void gaussianElimination(std::vector<std::vector<uint>> &matrix) {
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
                    matrix = std::vector<std::vector<uint>>(matrix.begin(), matrix.begin() + i);
                    break;
                }
            }
        }
// 







// Segmentation fault here
        //Function to perform the row addition and reduction process (reduceUp)
        void reduceUp(std::vector<std::vector<uint>>& matrix, const std::vector<int>& columnIndices) {
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
        void reduceDown(std::vector<std::vector<uint>>& matrix, const std::vector<int>& columnIndices) {
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
        void rowInsert(std::vector<std::vector<uint>>& matrix, int i, int j) {
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
            std::vector<std::vector<uint>> permutedMatrix(totalRows);
            for (int k = 0; k < totalRows; k++) {
                permutedMatrix[k] = matrix[indexOrder[k]];
            }
            matrix = permutedMatrix;  //Update the original matrix with the permuted matrix
        }

        void qubitReorder(std::vector<std::vector<uint>>& matrix, std::vector<int> indices)
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
            std::vector<std::vector<uint>> permutedTableau(matrix.size(), std::vector<uint>(matrix[0].size(), 0));

            for (int row = 0; row < matrix.size(); ++row) {
                for (int col = 0; col < matrix[row].size(); ++col) {
                    permutedTableau[row][newIndex[col]] = matrix[row][col];
                }
            }

            //Replace the original tableau with the permuted tableau
            matrix = permutedTableau;
        }

        void reduction(std::vector<std::vector<uint>>& graphMatrix) {
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

        void zReduction(std::vector<std::vector<uint>>& graphMatrix) {
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

            std::vector<std::vector<uint>> zMatrix = extractHalf(graphMatrix, false);

            std::set<uint> filledSet;
            for(int i = 0; i < xRank; i++)
            {
                filledSet.insert(i);
            }
            for (int j = xRank; j < stabCounts; j++) {
                std::vector<uint> col = zMatrix[j];
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

        std::vector<std::vector<uint>> convert_to_graph()
        {
            std::cout << "-----Test-----" << std::endl;

            if (x.size() != z.size()) {
                std::cerr << "Error: Matrices must have the same number of rows to append column-wise." << std::endl;
            }

            //Extract just the stabilizers w/o phase
            std::vector<std::vector<uint>> graphMatrix;
            graphMatrix.resize(rows/2, std::vector<uint>(2*cols,0));

            for (int i = 0; i < rows/2; i++) {
                for(int j = 0; j < cols; j++)
                {
                    graphMatrix[i][j] = x[i+rows/2][j]; //Copy the x row
                    graphMatrix[i][j+cols] = z[i+rows/2][j];
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

                std::vector<std::vector<uint>> zMatrix = extractHalf(graphMatrix, false);
            }//end if xRank/standard form



            xRank = nonZeroRows(extractHalf(graphMatrix, true));

            std::cout << "\n-----Standard form-----" << std::endl;
            print_table(graphMatrix);

            /*
                Graph processing after standard form
            */

            if(nonZeroRows(graphMatrix) != n)
                std::cout << "\n\n\nMismatch row/cols\n\n\n";
            
            std::vector<std::vector<uint>> xStabs = extractHalf(graphMatrix, true);

            //find the size of the identity
            xRank = nonZeroRows(xStabs);

            //Extract submatrices from graphMatrix
            std::vector<std::vector<uint>> I1(xRank, std::vector<uint>(xRank));
            std::vector<std::vector<uint>> A(xRank, std::vector<uint>(n - xRank));
            std::vector<std::vector<uint>> B(xRank, std::vector<uint>(xRank));
            std::vector<std::vector<uint>> sign1(xRank, std::vector<uint>(1, 0));

            std::vector<std::vector<uint>> D(graphMatrix.size() - xRank, std::vector<uint>(xRank));
            std::vector<std::vector<uint>> I2(graphMatrix.size() - xRank, std::vector<uint>(n - xRank - 1));
            std::vector<std::vector<uint>> sign2(graphMatrix.size() - xRank, std::vector<uint>(1, 0));

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
            std::vector<std::vector<uint>> newGraph(graphMatrix.size(), std::vector<uint>(graphMatrix[0].size(), 0));

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

    protected:
        int g;
        int n;
        int stabCounts;
        int rows;
        int cols;
        std::vector<int> outcomes; //Vector of binary digits
        std::vector<Gate> gates;
        std::vector<std::vector<uint>> x;
        std::vector<std::vector<uint>> z;
        std::vector<uint> r;
    }; //End tableau class

} // namespace NWQSim
