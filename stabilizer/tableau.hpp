#include "../include/state.hpp"

#include "../include/nwq_util.hpp"
#include "../include/gate.hpp"
#include "../include/circuit.hpp"
#include "../include/config.hpp"
#include "../include/private/exp_gate_declarations_host.hpp"

#include "../include/circuit_pass/fusion.hpp"
#include "../include/private/macros.hpp"
#include "../include/private/sim_gate.hpp"

#include <random>
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <vector>

namespace NWQSim
{
    class Tableau
    {
    public:
        //Constructor with a full set of gates
        Tableau(std::vector<Gate> &_gates, uint _numQubits)
        {
            gates = _gates;
            g = gates.size();
            n = _numQubits;
            result = 0;

            outcomes.resize(n);
            outcomes.assign(outcomes.size(), 0);
            x.resize(2*n+1, std::vector<uint>(n,0)); //first 2n+1 x n block. first n represents destabilizers
                                                     //second n represents stabilizers + 1 extra row
            z.resize(2*n+1, std::vector<uint>(n,0)); //second 2n+1 x n block to form the 2n+1 x 2n sized tableau
            r.resize(2*n+1, 0); //column on the right with 2n+1 rows
            //The 2n+1 th row is scratch space

            //Intialize the identity tableau
            for(int i = 0; i < n; i++)
            {
                x[i][i] = 1;
                z[i+n][i] = 1;
            }
            
        }

        //Manual constructor. Add gates later if necessary
        Tableau(std::shared_ptr<Circuit>& circuit, uint _numQubits)
        {
            gates = circuit->get_gates();
            g = gates.size();
            n = _numQubits;
            result = 0;

            outcomes.resize(n);
            outcomes.assign(outcomes.size(), 0);
            x.resize(2*n+1, std::vector<uint>(n,0)); //first 2n+1 x n block. first n represents destabilizers
                                                     //second n represents stabilizers + 1 extra row
            z.resize(2*n+1, std::vector<uint>(n,0)); //second 2n+1 x n block to form the 2n+1 x 2n sized tableau
            r.resize(2*n+1, 0); //column on the right with 2n+1 rows
            //The 2n+1 th row is scratch space
            //Intialize the identity tableau
            for(int i = 0; i < n; i++)
            {
                x[i][i] = 1;
                z[i+n][i] = 1;
            }
        }   

        void tableau_resize()
        {
            outcomes.resize(n);
            outcomes.assign(outcomes.size(), 0);

            //Add rows to x, z, r/
            x.resize(2*n+1); //first 2n+1 x n block. first n represents destabilizers
                                                     //second n represents stabilizers + 1 extra row
            z.resize(2*n+1); //second 2n+1 x n block to form the 2n+1 x 2n sized tableau
            r.resize(2*n+1); //column on the right with 2n+1 rows
            //The 2n+1 th row is scratch space
        }

        void add_gates(std::shared_ptr<Circuit>& new_circ)
        {
            gates = new_circ->get_gates();
            g = gates.size();
        }

        //Convert a vector of 1's and 0's to an IdxType decimal number
        void get_outcomes(IdxType*& result)
        {
            int conversion = 0;

            for (int i = 0; i < outcomes.size(); i++) {
                conversion = (conversion << 1) | outcomes[i];  // Left shift and add the current bit
            } 
            
            std::cout << "Conversion: " << conversion << std::endl;

            result = new long long(static_cast<long long>(conversion));

            //std::cout << "Result: " << result[0] << std::endl;
        }

        //Get the stabilizers in the tableau, i.e. all of the Pauli strings that stabilize the circuit
        std::vector<std::string> get_stabilizers()
        {
            int x_val;
            int z_val;
            std::string stabilizers;
            std::vector<std::string> pauliStrings;
            for(int i = n; i < 2*n; i++) //rows of stabilizers
            {
                stabilizers.clear(); //clear the temporary stabilizers string
                for(int j = 0; j < n; j++) //qubits/cols
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

        //Takes a default (or any) tableau and sets its stabilizers according to
        //a Pauli string provided
        void set_stabilizers(std::string pauliString)
        {
            for(int qubit = 0; qubit < T.n; qubit++)
            {
                for(int i = 0; i < 2*T.n+1; i++)
                {
                    switch(pauliString[qubit])
                    {
                        case 'I':
                            T.x[i][qubit] = 0;
                            T.z[i][qubit] = 0;
                            break;
                        case 'X':
                            T.x[i][qubit] = 0;
                            T.z[i][qubit] = 1;
                            break;
                        case 'Y':
                            T.x[i][qubit] = 1;
                            T.z[i][qubit] = 1;
                            break;
                        case 'Z':
                            T.x[i][qubit] = 1;
                            T.z[i][qubit] = 0;
                            break;
                        default:
                            std::logic_error("Invalid stabilizer");
                            break;
                    }
                }//All rows (2*n rows, destabilizers and stabilizers)
            }//All columns (qubits)
        }

        //Takes a default (or any) tableau and sets its stabilizers according to
        //a Pauli string provided
        void add_stabilizer(std::string pauliString)
        {
            //Start by adding a row to T
            n++;
            tableau_resize();

            assert(pauliString.length() <= n);

            for(int i = 0; i < pauliString.length(); i++)
            {
                switch(pauliString[i])
                {
                    case 'I':
                        x[x.size()][i] = 0;
                        z[z.size()][i] = 0;
                        break;
                    case 'X':
                        x[x.size()][i] = 0;
                        z[z.size()][i] = 1;
                        break;
                    case 'Y':
                        x[x.size()][i] = 1;
                        z[z.size()][i] = 1;
                        break;
                    case 'Z':
                        x[x.size()][i] = 1;
                        z[z.size()][i] = 0;
                        break;
                    default:
                        std::logic_error("Invalid stabilizer");
                        break;
                }
            }//All columns represented in the pauli string (qubits)
        }
        
        //Function to swap two rows of the tableau
        void swapRows(int row1, int row2) {
            std::swap(x[row1], x[row2]);
            std::swap(z[row1], z[row2]);
            std::swap(r[row1], r[row2]);
        }

        //Function to add row2 to row1 (mod 2)
        void addRows(int row1, int row2) {
            for (int i = 0; i < n; ++i) {
                x[row1][i] ^= x[row2][i];  //XOR for mod 2 addition
                z[row1][i] ^= z[row2][i];  //XOR for mod 2 addition
            }
            r[row1] = (r[row1] + r[row2]) % 4;  //Phase flip
        }

        //Algorithm to reduce the tableau to rref        
        void gaussianElimination() 
        {
            //Perform Gaussian elimination on the destabilizers (first half of the rows)
            for (int col = 0; col < n; col++) {
                int pivotRow = -1;
                for (int row = col; row < n; row++) {
                    if (x[row][col] == 1 || z[row][col] == 1) {
                        pivotRow = row;
                        break;
                    }
                }
                if (pivotRow == -1) continue;
                if (pivotRow != col) {
                    swapRows(pivotRow, col);
                }
                for (int row = 0; row < n; row++) {
                    if (row != col && (x[row][col] == 1 || z[row][col] == 1)) {
                        addRows(row, col);
                    }
                }
            }

            //Perform Gaussian elimination on the stabilizers (second half of the rows)
            for (int col = 0; col < n; col++) {
                int pivotRow = -1;
                for (int row = n + col; row < 2 * n; row++) {
                    if (x[row][col] == 1 || z[row][col] == 1) {
                        pivotRow = row;
                        break;
                    }
                }
                if (pivotRow == -1) continue;
                if (pivotRow != (n + col)) {
                    swapRows(pivotRow, n + col);
                }
                for (int row = n; row < 2 * n; row++) {
                    if (row != (n + col) && (x[row][col] == 1 || z[row][col] == 1)) {
                        addRows(row, n + col);
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
                    
                    for(int i = 0; i < 2*n; i++)
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

                    for(int i = 0; i < 2*n; i++)
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
                    for(int i = 0; i < 2*n; i++)
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
                    for(int p_index = n; p_index < 2*n; p_index++)
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

                        for(int i = 0; i < 2 * n; i++)
                        {
                            if((i != p) && (x[i][a] == 1))
                            {
                                rowsum(i, p);
                            }
                        }
                        x[p-n] = x[p];
                        z[p-n] = z[p];
                        //Change all the columns in row p to be 0
                        for(int i = 0; i < n; i++)
                        {
                            x[p][i] = 0;
                            z[p][i] = 0;                        
                        }

                        std::random_device rd;
                        std::mt19937 gen(rd());  //Mersenne Twister engine

                        // Define the range for random numbers
                        std::uniform_int_distribution<> distr(0, 1);

                        // Generate and display a random number
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
                            x[2*n][i] = 0;
                            z[2*n][i] = 0;
                        }

                        //Run rowsum subroutine
                        for(int i = 0; i < n; i++)
                        {
                            if(x[i][a] == 1)
                            {
                                //std::cout << "Perform rowsum at " << i << " + n" << std::endl;
                                rowsum(2*n, i+n);
                            }
                        }
                        //The result is the 
                        outcomes[a] = r[2*n];
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

                gaussianElimination();
            } //End gates for loop

        } //End tableau_simulation
    protected:
        int g;
        int n;
        std::vector<int> outcomes; //Basically a vector of binary numbers
        IdxType* result;
        std::vector<Gate> gates;
        std::vector<std::vector<uint>> x;
        std::vector<std::vector<uint>> z;
        std::vector<uint> r;
    }; //End tableau class
} //End namespace NWQSim