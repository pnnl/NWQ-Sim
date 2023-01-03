// ---------------------------------------------------------------------------
// NWQsim: Northwest Quantum Circuit Simulation Environment
// ---------------------------------------------------------------------------
// Ang Li, Senior Computer Scientist
// Pacific Northwest National Laboratory(PNNL), U.S.
// Homepage: http://www.angliphd.com
// GitHub repo: http://www.github.com/pnnl/DM-Sim
// PNNL-IPID: 31919-E, ECCN: EAR99, IR: PNNL-SA-143160
// GitHub repo: http://www.github.com/pnnl/DM-Sim
// PNNL-IPID: 32166, ECCN: EAR99, IR: PNNL-SA-161181
// BSD Lincese.
// ---------------------------------------------------------------------------
// File: fusion.h
// Define gate fusion functions for SV-Sim
// ---------------------------------------------------------------------------
#ifndef FUSION_H_
#define FUSION_H_

#include "gate.h"

namespace NWQSim {

//This is essentially matrix multiplication
void fuse_gate(const Gate& g0, const Gate& g1, Gate& g, const IdxType dim)
{
    for (IdxType m=0; m<dim; m++)
    {
        for (IdxType n=0; n<dim; n++)
        {
            g.gm_real[m*dim+n] = 0;
            g.gm_imag[m*dim+n] = 0;
            for (IdxType k=0; k<dim; k++)
            {
                g.gm_real[m*dim+n] += g0.gm_real[m*dim+k] * g1.gm_real[k*dim+n] - 
                    g0.gm_imag[m*dim+k] * g1.gm_imag[k*dim+n];
                g.gm_imag[m*dim+n] += g0.gm_real[m*dim+k] * g1.gm_imag[k*dim+n] +
                    g0.gm_imag[m*dim+k] * g1.gm_real[k*dim+n];
            }
        }
    }
}
//This is kronecker product
void kron(const Gate& g0, const Gate& g1, Gate& g, const IdxType dim)
{
    for (IdxType r = 0; r < dim; r++) 
    {
        for (IdxType s = 0; s < dim; s++) 
        {
            for (IdxType v = 0; v < dim; v++) 
            {
                for (IdxType w = 0; w < dim; w++) 
                {
                    g.gm_real[(dim*r+v)*dim*dim+(dim*s+w)] = g0.gm_real[r*dim+s] * g1.gm_real[v*dim+w] - 
                        g0.gm_imag[r*dim+s] * g1.gm_imag[v*dim+w];
                    g.gm_imag[(dim*r+v)*dim*dim+(dim*s+w)] = g0.gm_real[r*dim+s] * g1.gm_imag[v*dim+w] + 
                        g0.gm_imag[r*dim+s] * g1.gm_real[v*dim+w];
                }
            }
        }
    }
}
//This function is used for switch ctrl and qubit 
void reverse_ctrl_qubit(const Gate& g0, Gate& g)
{
    Gate SWAP(OP::C2, -1, -1, 0);
    SWAP.gm_real[0] = SWAP.gm_real[6] = SWAP.gm_real[9] = SWAP.gm_real[15] = 1;
    Gate tmp_g(OP::C2, -1, -1, 0);
    fuse_gate(g0, SWAP, tmp_g, 4);
    fuse_gate(SWAP, tmp_g, g, 4);
}
//This is for 2-qubit gate absorbing 1-qubit gate
void expand_ctrl(const Gate& g0, Gate& g, const IdxType ctrl) //This is to krok I for the ctrl qubit position
{
    Gate gI(OP::C1, ctrl, -1, 0);
    gI.gm_real[0] = gI.gm_real[3] = 1; //set I
    kron(gI, g0, g, 2);
}
void expand_qubit(const Gate& g0, Gate& g, const IdxType qubit)
{
    Gate gI(OP::C1, qubit, -1, 0);
    gI.gm_real[0] = gI.gm_real[3] = 1; //set I
    kron(g0, gI, g, 2);
}
//This function is to fuse C2 gates (allows switching ctrl/qubit, e.g., in a SWAP gate)
void gate_fusion_2q(const vector<Gate>& circuit_in, vector<Gate>& circuit_out, const IdxType n_qubits)
{
    //prepare
    IdxType* table = new IdxType[n_qubits*n_qubits];
    bool* canfuse = new bool[n_qubits*n_qubits];
    for (IdxType i=0; i<n_qubits*n_qubits; i++) table[i] = -1;
    for (IdxType i=0; i<n_qubits*n_qubits; i++) canfuse[i] = false;
    circuit_out.clear();

    for (IdxType i=0; i<circuit_in.size(); i++)
    {
        if (circuit_in[i].op_name == C1) //1-qubit gate
        {
            IdxType qubit = circuit_in[i].qubit;
            Gate g(circuit_in[i]);
            for (IdxType j=0; j<n_qubits; j++)
            {
                canfuse[j*n_qubits+qubit] = false;
                canfuse[qubit*n_qubits+j] = false;
            }
            circuit_out.push_back(g);
        }
        if (circuit_in[i].op_name == C2) //2-qubit gate
        {
            IdxType qubit = circuit_in[i].qubit;
            IdxType ctrl = circuit_in[i].ctrl;
            Gate g0(circuit_in[i]);

            //we reverse ctrl-qubit to harvest more fusion opportunities
            // M = SWAP N SWAP
            if (ctrl > qubit) 
            {
                reverse_ctrl_qubit(circuit_in[i], g0);
                g0.ctrl = circuit_in[i].qubit;
                ctrl = circuit_in[i].qubit;
                g0.qubit = circuit_in[i].ctrl;
                qubit = circuit_in[i].ctrl;
            }
            if (canfuse[ctrl*n_qubits+qubit] == false) //cannot fuse
            {
                Gate g(g0);
                circuit_out.push_back(g);
                table[ctrl*n_qubits+qubit] = circuit_out.size()-1; //point for later fusion
                for (IdxType j=0; j<n_qubits; j++)
                {
                    canfuse[qubit*n_qubits+j] = false;
                    canfuse[ctrl*n_qubits+j] = false;
                    canfuse[j*n_qubits+qubit] = false;
                    canfuse[j*n_qubits+ctrl] = false;
                }
                canfuse[ctrl*n_qubits+qubit] = true;
            }
            else //able to fuse
            {
                Gate& g1 = circuit_out[table[ctrl*n_qubits+qubit]];
                Gate final_g(OP::C2, -1, -1, 0);
                fuse_gate(g0, g1, final_g, 4);
                memcpy(g1.gm_real, final_g.gm_real, 16*sizeof(ValType));
                memcpy(g1.gm_imag, final_g.gm_imag, 16*sizeof(ValType));
            }

        }
        if (circuit_in[i].op_name == M) //1-qubit measure gate
        {
            IdxType qubit = circuit_in[i].qubit;
            Gate g(circuit_in[i]);
            for (IdxType j=0; j<n_qubits; j++)
            {
                canfuse[j*n_qubits+qubit] = false;
                canfuse[qubit*n_qubits+j] = false;
            }
            circuit_out.push_back(g);
        }
        if (circuit_in[i].op_name == RESET) //1-qubit reset gate
        {
            IdxType qubit = circuit_in[i].qubit;
            Gate g(circuit_in[i]);
            for (IdxType j=0; j<n_qubits; j++)
            {
                canfuse[j*n_qubits+qubit] = false;
                canfuse[qubit*n_qubits+j] = false;
            }
            circuit_out.push_back(g);
        }
        if (circuit_in[i].op_name == MA) //all-qubit measure gate
        {
            Gate g(circuit_in[i]);
            for (IdxType j=0; j<n_qubits; j++)
                for (IdxType k=0; k<n_qubits; k++)
                    canfuse[j*n_qubits+k] = false;
            circuit_out.push_back(g);
        }
    }
    //clean
    delete[] table;
    delete[] canfuse;
}

//This function is to fuse C1 gates in a circuit
void gate_fusion_1q(const vector<Gate>& circuit_in, vector<Gate>& circuit_out, const IdxType n_qubits)
{
    //prepare
    IdxType* table = new IdxType[n_qubits];
    bool* canfuse = new bool[n_qubits];
    for (IdxType i=0; i<n_qubits; i++) table[i] = -1;
    for (IdxType i=0; i<n_qubits; i++) canfuse[i] = false;
    circuit_out.clear();
    //parse
    for (IdxType i=0; i<circuit_in.size(); i++)
    {
        if (circuit_in[i].op_name == C1) //1-qubit gate
        {
            IdxType qubit = circuit_in[i].qubit;
            if (canfuse[qubit] == false) //cannot fuse
            {
                Gate g(circuit_in[i]);
                circuit_out.push_back(g);
                canfuse[qubit] = true;
                table[qubit] = circuit_out.size()-1; //point to this gate for later fusion
            }
            else //able to fuse
            {
                //start to fuse circuit_in[i] and circuit_out[table[qubit]]
                const Gate& g0 = circuit_in[i];
                Gate& g1 = circuit_out[table[qubit]];

                ValType res_real[4] = {0};
                ValType res_imag[4] = {0};
                for (int m=0; m<2; m++)
                    for (int n=0; n<2; n++)
                        for (int k=0; k<2; k++)
                        {
                            res_real[m*2+n] += g0.gm_real[m*2+k] * g1.gm_real[k*2+n] - 
                                g0.gm_imag[m*2+k] * g1.gm_imag[k*2+n];
                            res_imag[m*2+n] += g0.gm_real[m*2+k] * g1.gm_imag[k*2+n] +
                                g0.gm_imag[m*2+k] * g1.gm_real[k*2+n];
                        }
                memcpy(g1.gm_real, res_real, 4*sizeof(ValType));
                memcpy(g1.gm_imag, res_imag, 4*sizeof(ValType));
            }
        }
        if (circuit_in[i].op_name == C2) //2-qubit gate
        {
            canfuse[circuit_in[i].qubit] = false;
            canfuse[circuit_in[i].ctrl] = false;
            Gate g(circuit_in[i]);
            circuit_out.push_back(g);
        }
        if (circuit_in[i].op_name == M) //1-qubit measure gate
        {
            canfuse[circuit_in[i].qubit] = false;
            Gate g(circuit_in[i]);
            circuit_out.push_back(g);
        }
        if (circuit_in[i].op_name == RESET) //1-qubit reset gate
        {
            canfuse[circuit_in[i].qubit] = false;
            Gate g(circuit_in[i]);
            circuit_out.push_back(g);
        }
        if (circuit_in[i].op_name == MA) //all-qubit measure gate
        {
            Gate g(circuit_in[i]);
            circuit_out.push_back(g);
            for (IdxType q=0; q<n_qubits; q++) canfuse[q] = false;
        }
    }
    //clean
    delete[] table;
    delete[] canfuse;
}

void gate_fusion_2q_absorb_1q_forward(const vector<Gate>& circuit_in, vector<Gate>& circuit_out, const IdxType n_qubits)
{
    //prepare
    IdxType* table = new IdxType[n_qubits*n_qubits];
    bool* canfuse = new bool[n_qubits*n_qubits];
    for (IdxType i=0; i<n_qubits*n_qubits; i++) table[i] = -1;
    for (IdxType i=0; i<n_qubits*n_qubits; i++) canfuse[i] = false;
    circuit_out.clear();
    //forward parse
    for (IdxType i=0; i<circuit_in.size(); i++)
    {
        if (circuit_in[i].op_name == C1) //1-qubit gate
        {
            IdxType qubit = circuit_in[i].qubit;
            bool fused = false;
            for (IdxType j=0; j<n_qubits; j++)
            {
                //target qubit can merge
                if (canfuse[j*n_qubits+qubit] == true)
                {
                    const Gate& g0 = circuit_in[i];
                    Gate& g1 = circuit_out[table[j*n_qubits+qubit]];
                    Gate expand_g(OP::C2,-1,-1,0);
                    Gate final_g(OP::C2,-1,-1,0);
                    expand_ctrl(g0, expand_g, j);
                    //now ready to fuse
                    fuse_gate(expand_g, g1, final_g, 4); 
                    memcpy(g1.gm_real, final_g.gm_real, 16*sizeof(ValType));
                    memcpy(g1.gm_imag, final_g.gm_imag, 16*sizeof(ValType));
                    fused = true;
                    break;
                }
                //ctrl qubit can merge
                if (canfuse[qubit*n_qubits+j] == true)
                {
                    const Gate& g0 = circuit_in[i];
                    Gate& g1 = circuit_out[table[qubit*n_qubits+j]];
                    Gate expand_g(OP::C2,-1,-1,0);
                    expand_qubit(g0, expand_g, j);
                    //now ready to fuse
                    Gate final_g(OP::C2,-1,-1,0);
                    fuse_gate(expand_g, g1, final_g, 4); 
                    memcpy(g1.gm_real, final_g.gm_real, 16*sizeof(ValType));
                    memcpy(g1.gm_imag, final_g.gm_imag, 16*sizeof(ValType));
                    fused = true;
                    break;
                }
            }
            if (fused == false) //did not fuse with any C4 gate, add it
            {
                for (IdxType j=0; j<n_qubits; j++)
                {
                    canfuse[j*n_qubits+qubit] = false;
                    canfuse[qubit*n_qubits+j] = false;
                }
                Gate g(circuit_in[i]);
                circuit_out.push_back(g);
            }
        }
        if (circuit_in[i].op_name == C2) //2-qubit gate
        {
            IdxType qubit = circuit_in[i].qubit;
            IdxType ctrl = circuit_in[i].ctrl;
            for (IdxType j=0; j<n_qubits; j++)
            {
                canfuse[j*n_qubits+qubit] = false;
                canfuse[qubit*n_qubits+j] = false;
                canfuse[j*n_qubits+ctrl] = false;
                canfuse[ctrl*n_qubits+j] = false;
            }
            Gate g(circuit_in[i]);
            circuit_out.push_back(g);
            canfuse[ctrl*n_qubits+qubit] = true;
            table[ctrl*n_qubits+qubit] = circuit_out.size()-1; //this gate can fuse
        }
        if (circuit_in[i].op_name == M) //1-qubit measure gate
        {
            IdxType qubit = circuit_in[i].qubit;
            for (IdxType j=0; j<n_qubits; j++)
            {
                canfuse[j*n_qubits+qubit] = false;
                canfuse[qubit*n_qubits+j] = false;
            }
            Gate g(circuit_in[i]);
            circuit_out.push_back(g);
        }
        if (circuit_in[i].op_name == RESET) //1-qubit reset gate
        {
            IdxType qubit = circuit_in[i].qubit;
            for (IdxType j=0; j<n_qubits; j++)
            {
                canfuse[j*n_qubits+qubit] = false;
                canfuse[qubit*n_qubits+j] = false;
            }
            Gate g(circuit_in[i]);
            circuit_out.push_back(g);
        }
        if (circuit_in[i].op_name == MA) //all-qubit measure gate
        {
            for (IdxType j=0; j<n_qubits; j++)
                for (IdxType k=0; k<n_qubits; k++)
                    canfuse[j*n_qubits+k] = false;
            Gate g(circuit_in[i]);
            circuit_out.push_back(g);
        }
    }
    //clean
    delete[] table;
    delete[] canfuse;
}

void gate_fusion_2q_absorb_1q_backward(const vector<Gate>& circuit_in, vector<Gate>& circuit_out, const IdxType n_qubits)
{
    //prepare
    IdxType* table = new IdxType[n_qubits*n_qubits];
    bool* canfuse = new bool[n_qubits*n_qubits];
    for (IdxType i=0; i<n_qubits*n_qubits; i++) table[i] = -1;
    for (IdxType i=0; i<n_qubits*n_qubits; i++) canfuse[i] = false;
    circuit_out.clear();
    //forward parse
    for (IdxType i=circuit_in.size()-1; i>=0; i--)
    {
        if (circuit_in[i].op_name == C1) //1-qubit gate
        {
            IdxType qubit = circuit_in[i].qubit;
            bool fused = false;
            for (IdxType j=0; j<n_qubits; j++)
            {
                //target qubit can merge
                if (canfuse[j*n_qubits+qubit] == true)
                {
                    const Gate& g0 = circuit_in[i];
                    Gate& g1 = circuit_out[table[j*n_qubits+qubit]];
                    Gate expand_g(OP::C2,-1,-1,0);
                    Gate final_g(OP::C2,-1,-1,0);
                    expand_ctrl(g0, expand_g, j);
                    //now ready to fuse
                    fuse_gate(g1, expand_g, final_g, 4); 
                    memcpy(g1.gm_real, final_g.gm_real, 16*sizeof(ValType));
                    memcpy(g1.gm_imag, final_g.gm_imag, 16*sizeof(ValType));
                    fused = true;
                    break;
                }
                //ctrl qubit can merge
                if (canfuse[qubit*n_qubits+j] == true)
                {
                    const Gate& g0 = circuit_in[i];
                    Gate& g1 = circuit_out[table[qubit*n_qubits+j]];
                    Gate expand_g(OP::C2,-1,-1,0);
                    expand_qubit(g0, expand_g, j);
                    //now ready to fuse
                    Gate final_g(OP::C2,-1,-1,0);
                    fuse_gate(g1, expand_g, final_g, 4); 
                    memcpy(g1.gm_real, final_g.gm_real, 16*sizeof(ValType));
                    memcpy(g1.gm_imag, final_g.gm_imag, 16*sizeof(ValType));
                    fused = true;
                    break;
                }
            }
            if (fused == false) //did not fuse with any C4 gate, add it
            {
                for (IdxType j=0; j<n_qubits; j++)
                {
                    canfuse[j*n_qubits+qubit] = false;
                    canfuse[qubit*n_qubits+j] = false;
                }
                Gate g(circuit_in[i]);
                circuit_out.push_back(g);
            }
        }
        if (circuit_in[i].op_name == C2) //2-qubit gate
        {
            IdxType qubit = circuit_in[i].qubit;
            IdxType ctrl = circuit_in[i].ctrl;
            for (IdxType j=0; j<n_qubits; j++)
            {
                canfuse[j*n_qubits+qubit] = false;
                canfuse[qubit*n_qubits+j] = false;
                canfuse[j*n_qubits+ctrl] = false;
                canfuse[ctrl*n_qubits+j] = false;
            }
            Gate g(circuit_in[i]);
            circuit_out.push_back(g);

            canfuse[ctrl*n_qubits+qubit] = true;
            table[ctrl*n_qubits+qubit] = circuit_out.size()-1; //this gate can fuse
        }
        if (circuit_in[i].op_name == M) //1-qubit measure gate
        {
            IdxType qubit = circuit_in[i].qubit;
            for (IdxType j=0; j<n_qubits; j++)
            {
                canfuse[j*n_qubits+qubit] = false;
                canfuse[qubit*n_qubits+j] = false;
            }
            Gate g(circuit_in[i]);
            circuit_out.push_back(g);
        }
        if (circuit_in[i].op_name == RESET) //1-qubit reset gate
        {
            IdxType qubit = circuit_in[i].qubit;
            for (IdxType j=0; j<n_qubits; j++)
            {
                canfuse[j*n_qubits+qubit] = false;
                canfuse[qubit*n_qubits+j] = false;
            }
            Gate g(circuit_in[i]);
            circuit_out.push_back(g);
        }
        if (circuit_in[i].op_name == MA) //all-qubit measure gate
        {
            for (IdxType j=0; j<n_qubits; j++)
                for (IdxType k=0; k<n_qubits; k++)
                    canfuse[j*n_qubits+k] = false;
            Gate g(circuit_in[i]);
            circuit_out.push_back(g);
        }
    }
    //clean
    delete[] table;
    delete[] canfuse;
    reverse(circuit_out.begin(), circuit_out.end());
}


}//end of NWQSim
#endif
