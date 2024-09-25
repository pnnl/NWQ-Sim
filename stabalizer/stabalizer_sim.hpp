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
    class tableau
    {
    public:
        tableau(const std::vector<SVGate> &_gates, uint _numQubits)
        {
            gates = _gates;
            g = gates.size();
            n = _numQubits;
            x.resize(2*n+1, std::vector<uint>(n,0));
            z.resize(2*n+1, std::vector<uint>(n,0));
            r.resize(2*n+1, 0);
        }

        void simulate()
        {
            //For swapping rows
            std::vector<uint> tempRow(n);

            for (int j = 0; j < g; j++)
            {
                auto gate = gates[j];
                int a = gate.qubit;
                if (gate.op_name == OP::X)
                {
                    for(int i = 0; i < 2*n; i++)
                        z[a][i] = z[a][i] ^ 1;
                }
                else if (gate.op_name == OP::Y)
                {
                    for(int i = 0; i < 2*n; i++)
                    {
                        x[a][i] = x[a][i] ^ 1;
                        z[a][i] = z[a][i] ^ 1;
                    }
                }
                else if (gate.op_name == OP::Z)
                {
                    for(int i = 0; i < 2*n; i++)
                        x[a][i] = x[a][i] ^ 1;
                }
                else if (gate.op_name == OP::H)
                {
                    //Phase
                    for(int i = 0; i < 2*n; i++)
                        r[i] = r[i] ^ ((x[i][a] >> 1) + z[i][a]);

                    //Entry
                    tempRow = x[a];
                    x[a] = z[a];
                    z[a] = tempRow;    

                }
                else if (gate.op_name == OP::S)
                {
                    int a = gate.qubit;

                    for(int i = 0; i < 2*n; i++)
                    {
                        //Phase
                        r[i] = r[i] ^ ((x[i][a] >> 1) + z[i][a]);

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
                        r[i] = r[i] ^ (((x[i][a] >> 1) + z[i][b])*(x[i][b]^z[i][a]^1));

                        //Entry
                        x[b][i] = x[b][i] ^ x[a][i];
                        z[a][i] = z[a][i] ^ z[b][i];
                    }
                }
                else if (gate.op_name == OP::M)
                {  
                    int a = gate.ctrl;
                    int p = -1;
                    for(int i = n; i < n; i++)
                    {
                        if(x[i][a])
                        {
                            p = i;
                            break;
                        }
                    }
                    if(p > -1)
                    {
                        for(int i = 0; i < 2 * n; i++)
                        {
                            if((i != p) && (x[i][a] == 1))
                            {
                                rowsum(i, p);
                            }
                        }
                        x[p-n] = x[p];
                        z[p-n] = z[p];
                        std::fill(x[p].begin(), x[p].end(), 0);
                        std::fill(z[p].begin(), z[p].end(), 0);

                        std::srand(std::time(0));
                        if(std::rand() % 2)
                            r[p] = 1;
                        else
                            r[p] = 0;
                        z[p][a] = 1;

                        int outcome = r[p];                
                    }
                    else
                    {
                        x[2*n + 1].assign(x[2*n + 1].size(),0);
                        for(int i = 0; i < n; i++)
                        {
                            if(x[i][a] == 1)
                            {
                                rowsum(2*n+1, i+n);
                            }
                        }
                        int outcome = r[2*n+1];
                    }

                }
                else    
                {
                    std::cout << "../Non-Clifford or unrecognized gate" << std::endl
                                << OP_NAMES[gate.op_name] << std::endl;
                    std::logic_error("../Invalid gate type");
                }
            } //End for g
        } //End tableau_simulation
    protected:
        int g;
        int n;
        std::vector<SVGate> gates;
        std::vector<std::vector<uint>> x;
        std::vector<std::vector<uint>> z;
        std::vector<uint> r;

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
    }; //End tableau class
} //End namespace NWQSim
