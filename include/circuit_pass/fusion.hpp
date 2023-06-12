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
#pragma once

#include <vector>
#include <cstring>
#include <algorithm>

#include "../private/config.hpp"
#include "../public/circuit.hpp"
#include "../public/util.hpp"

#include "../private/sim_gate.hpp"
#include "../private/gate_factory/sv_gates.hpp"
#include "../private/gate_factory/dm_gates.hpp"

#define vec_dim_1q(sim_type) (sim_type == SV ? 2 : 4)
#define vec_dim_2q(sim_type) (sim_type == SV ? 4 : 16)
#define vec_length_1q(sim_type) (sim_type == SV ? 4 : 16)
#define vec_length_2q(sim_type) (sim_type == SV ? 16 : 256)

namespace NWQSim
{
    enum SimType
    {
        SV,
        DM
    };

    enum OP get_op(IdxType n_qubits, SimType sim_type)
    {
        if (sim_type == SV)
        {
            return n_qubits == 1 ? OP::C1 : OP::C2;
        }
        else
        {
            return n_qubits == 1 ? OP::C2 : OP::C4;
        }
    }

    // This is essentially matrix multiplication
    template <typename GateType>
    void fuse_gate(const GateType &g0, const GateType &g1, GateType &g, const IdxType dim)
    {
        for (IdxType m = 0; m < dim; m++)
        {
            for (IdxType n = 0; n < dim; n++)
            {
                g.gm_real[m * dim + n] = 0;
                g.gm_imag[m * dim + n] = 0;
                for (IdxType k = 0; k < dim; k++)
                {
                    g.gm_real[m * dim + n] += g0.gm_real[m * dim + k] * g1.gm_real[k * dim + n] -
                                              g0.gm_imag[m * dim + k] * g1.gm_imag[k * dim + n];
                    g.gm_imag[m * dim + n] += g0.gm_real[m * dim + k] * g1.gm_imag[k * dim + n] +
                                              g0.gm_imag[m * dim + k] * g1.gm_real[k * dim + n];
                }
            }
        }
    }
    // This is kronecker product
    template <typename GateType>
    void kron(const GateType &g0, const GateType &g1, GateType &g, const IdxType dim)
    {
        for (IdxType r = 0; r < dim; r++)
        {
            for (IdxType s = 0; s < dim; s++)
            {
                for (IdxType v = 0; v < dim; v++)
                {
                    for (IdxType w = 0; w < dim; w++)
                    {
                        g.gm_real[(dim * r + v) * dim * dim + (dim * s + w)] = g0.gm_real[r * dim + s] * g1.gm_real[v * dim + w] -
                                                                               g0.gm_imag[r * dim + s] * g1.gm_imag[v * dim + w];
                        g.gm_imag[(dim * r + v) * dim * dim + (dim * s + w)] = g0.gm_real[r * dim + s] * g1.gm_imag[v * dim + w] +
                                                                               g0.gm_imag[r * dim + s] * g1.gm_real[v * dim + w];
                    }
                }
            }
        }
    }

    // This function is used for switch ctrl and qubit
    template <typename GateType>
    void reverse_ctrl_qubit(const GateType &g0, GateType &g, const SimType sim_type)
    {

        GateType SWAP_SV(get_op(2, sim_type), -1, -1);
        SWAP_SV.gm_real[0] = SWAP_SV.gm_real[6] = SWAP_SV.gm_real[9] = SWAP_SV.gm_real[15] = 1;

        GateType SWAP(SWAP_SV);

        if (sim_type == DM)
            kron(SWAP_SV, SWAP_SV, SWAP, 4); // create SWAP superop for DM

        GateType tmp_g(get_op(2, sim_type), -1, -1);
        fuse_gate(g0, SWAP, tmp_g, vec_dim_2q(sim_type));
        fuse_gate(SWAP, tmp_g, g, vec_dim_2q(sim_type));
    }

    // This is for 2-qubit gate absorbing 1-qubit gate
    template <typename GateType>
    void expand_ctrl(const GateType &g0, GateType &g, const IdxType ctrl, const SimType sim_type) // This is to krok I for the ctrl qubit position
    {
        GateType gI(get_op(1, sim_type), ctrl, -1);

        if (sim_type == SV)
            gI.gm_real[0] = gI.gm_real[3] = 1; // set I for SV
        else
            gI.gm_real[0] = gI.gm_real[5] = gI.gm_real[10] = gI.gm_real[15] = 1; // set I for DM

        kron(gI, g0, g, vec_dim_1q(sim_type));
    }

    template <typename GateType>
    void expand_qubit(const GateType &g0, GateType &g, const IdxType qubit, const SimType sim_type)
    {
        GateType gI(get_op(1, sim_type), qubit, -1);

        if (sim_type == SV)
            gI.gm_real[0] = gI.gm_real[3] = 1; // set I for SV
        else
            gI.gm_real[0] = gI.gm_real[5] = gI.gm_real[10] = gI.gm_real[15] = 1; // set I for DM

        kron(g0, gI, g, vec_dim_1q(sim_type));
    }

    // This function is to fuse C2 gates (allows switching ctrl/qubit, e.g., in a SWAP gate)
    template <typename GateType>
    void gate_fusion_2q(const std::vector<GateType> &circuit_in, std::vector<GateType> &circuit_out, const IdxType n_qubits, const SimType sim_type)
    {
        // prepare
        IdxType *table = new IdxType[n_qubits * n_qubits];
        bool *canfuse = new bool[n_qubits * n_qubits];
        for (IdxType i = 0; i < n_qubits * n_qubits; i++)
            table[i] = -1;
        for (IdxType i = 0; i < n_qubits * n_qubits; i++)
            canfuse[i] = false;
        circuit_out.clear();

        for (size_t i = 0; i < circuit_in.size(); i++)
        {

            if (circuit_in[i].op_name == M) // 1-qubit measure gate
            {
                IdxType qubit = circuit_in[i].qubit;
                GateType g(circuit_in[i]);
                for (IdxType j = 0; j < n_qubits; j++)
                {
                    canfuse[j * n_qubits + qubit] = false;
                    canfuse[qubit * n_qubits + j] = false;
                }
                circuit_out.push_back(g);
            }
            else if (circuit_in[i].op_name == MA) // all-qubit measure gate
            {
                GateType g(circuit_in[i]);
                for (IdxType j = 0; j < n_qubits; j++)
                    for (IdxType k = 0; k < n_qubits; k++)
                        canfuse[j * n_qubits + k] = false;
                circuit_out.push_back(g);
            }
            else if (circuit_in[i].op_name == RESET) // 1-qubit reset gate
            {
                IdxType qubit = circuit_in[i].qubit;
                GateType g(circuit_in[i]);
                for (IdxType j = 0; j < n_qubits; j++)
                {
                    canfuse[j * n_qubits + qubit] = false;
                    canfuse[qubit * n_qubits + j] = false;
                }
                circuit_out.push_back(g);
            }
            else if (circuit_in[i].op_name == get_op(1, sim_type)) // 1-qubit gate
            {
                IdxType qubit = circuit_in[i].qubit;
                GateType g(circuit_in[i]);
                for (IdxType j = 0; j < n_qubits; j++)
                {
                    canfuse[j * n_qubits + qubit] = false;
                    canfuse[qubit * n_qubits + j] = false;
                }
                circuit_out.push_back(g);
            }
            else if (circuit_in[i].op_name == get_op(2, sim_type)) // 2-qubit gate
            {
                IdxType qubit = circuit_in[i].qubit;
                IdxType ctrl = circuit_in[i].ctrl;
                GateType g0(circuit_in[i]);

                // we reverse ctrl-qubit to harvest more fusion opportunities
                //  M = SWAP N SWAP
                if (ctrl > qubit)
                {
                    reverse_ctrl_qubit(circuit_in[i], g0, sim_type);
                    g0.ctrl = circuit_in[i].qubit;
                    ctrl = circuit_in[i].qubit;
                    g0.qubit = circuit_in[i].ctrl;
                    qubit = circuit_in[i].ctrl;
                }
                if (canfuse[ctrl * n_qubits + qubit] == false) // cannot fuse
                {
                    GateType g(g0);
                    circuit_out.push_back(g);
                    table[ctrl * n_qubits + qubit] = circuit_out.size() - 1; // point for later fusion
                    for (IdxType j = 0; j < n_qubits; j++)
                    {
                        canfuse[qubit * n_qubits + j] = false;
                        canfuse[ctrl * n_qubits + j] = false;
                        canfuse[j * n_qubits + qubit] = false;
                        canfuse[j * n_qubits + ctrl] = false;
                    }
                    canfuse[ctrl * n_qubits + qubit] = true;
                }
                else // able to fuse
                {
                    GateType &g1 = circuit_out[table[ctrl * n_qubits + qubit]];
                    GateType final_g(get_op(2, sim_type), -1, -1);
                    fuse_gate(g0, g1, final_g, vec_dim_2q(sim_type));
                    memcpy(g1.gm_real, final_g.gm_real, vec_length_2q(sim_type) * sizeof(ValType));
                    memcpy(g1.gm_imag, final_g.gm_imag, vec_length_2q(sim_type) * sizeof(ValType));
                }
            }
        }
        // clean
        delete[] table;
        delete[] canfuse;
    }

    // This function is to fuse C1 gates in a circuit
    template <typename GateType>
    void gate_fusion_1q(const std::vector<GateType> &circuit_in, std::vector<GateType> &circuit_out, const IdxType n_qubits, const SimType sim_type)
    {
        // prepare
        IdxType *table = new IdxType[n_qubits];
        bool *canfuse = new bool[n_qubits];
        for (IdxType i = 0; i < n_qubits; i++)
            table[i] = -1;
        for (IdxType i = 0; i < n_qubits; i++)
            canfuse[i] = false;
        circuit_out.clear();
        // parse
        for (size_t i = 0; i < circuit_in.size(); i++)
        {
            if (circuit_in[i].op_name == M) // 1-qubit measure gate
            {
                canfuse[circuit_in[i].qubit] = false;
                GateType g(circuit_in[i]);
                circuit_out.push_back(g);
            }
            else if (circuit_in[i].op_name == RESET) // 1-qubit reset gate
            {
                canfuse[circuit_in[i].qubit] = false;
                GateType g(circuit_in[i]);
                circuit_out.push_back(g);
            }
            else if (circuit_in[i].op_name == MA) // all-qubit measure gate
            {
                GateType g(circuit_in[i]);
                circuit_out.push_back(g);
                for (IdxType q = 0; q < n_qubits; q++)
                    canfuse[q] = false;
            }
            else if (circuit_in[i].op_name == get_op(1, sim_type)) // 1-qubit gate
            {
                IdxType qubit = circuit_in[i].qubit;
                if (canfuse[qubit] == false) // cannot fuse
                {
                    GateType g(circuit_in[i]);
                    circuit_out.push_back(g);
                    canfuse[qubit] = true;
                    table[qubit] = circuit_out.size() - 1; // point to this gate for later fusion
                }
                else // able to fuse
                {
                    // start to fuse circuit_in[i] and circuit_out[table[qubit]]
                    const GateType &g0 = circuit_in[i];
                    GateType &g1 = circuit_out[table[qubit]];

                    ValType res_real[16] = {0};
                    ValType res_imag[16] = {0};

                    IdxType dim = vec_dim_1q(sim_type);
                    for (int m = 0; m < dim; m++)
                        for (int n = 0; n < dim; n++)
                            for (int k = 0; k < dim; k++)
                            {
                                res_real[m * dim + n] += g0.gm_real[m * dim + k] * g1.gm_real[k * dim + n] -
                                                         g0.gm_imag[m * dim + k] * g1.gm_imag[k * dim + n];
                                res_imag[m * dim + n] += g0.gm_real[m * dim + k] * g1.gm_imag[k * dim + n] +
                                                         g0.gm_imag[m * dim + k] * g1.gm_real[k * dim + n];
                            }

                    memcpy(g1.gm_real, res_real, vec_length_1q(sim_type) * sizeof(ValType));
                    memcpy(g1.gm_imag, res_imag, vec_length_1q(sim_type) * sizeof(ValType));
                }
            }
            else if (circuit_in[i].op_name == get_op(2, sim_type)) // 2-qubit gate
            {
                canfuse[circuit_in[i].qubit] = false;
                canfuse[circuit_in[i].ctrl] = false;
                GateType g(circuit_in[i]);
                circuit_out.push_back(g);
            }
        }
        // clean
        delete[] table;
        delete[] canfuse;
    }

    template <typename GateType>
    void gate_fusion_2q_absorb_1q_forward(const std::vector<GateType> &circuit_in, std::vector<GateType> &circuit_out, const IdxType n_qubits, const SimType sim_type)
    {
        // prepare
        IdxType *table = new IdxType[n_qubits * n_qubits];
        bool *canfuse = new bool[n_qubits * n_qubits];
        for (IdxType i = 0; i < n_qubits * n_qubits; i++)
            table[i] = -1;
        for (IdxType i = 0; i < n_qubits * n_qubits; i++)
            canfuse[i] = false;
        circuit_out.clear();
        // forward parse
        for (size_t i = 0; i < circuit_in.size(); i++)
        {

            if (circuit_in[i].op_name == M) // 1-qubit measure gate
            {
                IdxType qubit = circuit_in[i].qubit;
                for (IdxType j = 0; j < n_qubits; j++)
                {
                    canfuse[j * n_qubits + qubit] = false;
                    canfuse[qubit * n_qubits + j] = false;
                }
                GateType g(circuit_in[i]);
                circuit_out.push_back(g);
            }
            else if (circuit_in[i].op_name == RESET) // 1-qubit reset gate
            {
                IdxType qubit = circuit_in[i].qubit;
                for (IdxType j = 0; j < n_qubits; j++)
                {
                    canfuse[j * n_qubits + qubit] = false;
                    canfuse[qubit * n_qubits + j] = false;
                }
                GateType g(circuit_in[i]);
                circuit_out.push_back(g);
            }
            else if (circuit_in[i].op_name == MA) // all-qubit measure gate
            {
                for (IdxType j = 0; j < n_qubits; j++)
                    for (IdxType k = 0; k < n_qubits; k++)
                        canfuse[j * n_qubits + k] = false;
                GateType g(circuit_in[i]);
                circuit_out.push_back(g);
            }
            else if (circuit_in[i].op_name == get_op(1, sim_type)) // 1-qubit gate
            {
                IdxType qubit = circuit_in[i].qubit;
                bool fused = false;
                for (IdxType j = 0; j < n_qubits; j++)
                {
                    // target qubit can merge
                    if (canfuse[j * n_qubits + qubit] == true)
                    {
                        const GateType &g0 = circuit_in[i];
                        GateType &g1 = circuit_out[table[j * n_qubits + qubit]];

                        GateType expand_g(get_op(2, sim_type), -1, -1);
                        GateType final_g(get_op(2, sim_type), -1, -1);

                        expand_ctrl(g0, expand_g, j, sim_type);
                        // now ready to fuse
                        fuse_gate(expand_g, g1, final_g, vec_dim_2q(sim_type));
                        memcpy(g1.gm_real, final_g.gm_real, vec_length_2q(sim_type) * sizeof(ValType));
                        memcpy(g1.gm_imag, final_g.gm_imag, vec_length_2q(sim_type) * sizeof(ValType));
                        fused = true;
                        break;
                    }
                    // ctrl qubit can merge
                    if (canfuse[qubit * n_qubits + j] == true)
                    {
                        const GateType &g0 = circuit_in[i];
                        GateType &g1 = circuit_out[table[qubit * n_qubits + j]];
                        GateType expand_g(get_op(2, sim_type), -1, -1);
                        expand_qubit(g0, expand_g, j, sim_type);
                        // now ready to fuse
                        GateType final_g(get_op(2, sim_type), -1, -1);
                        fuse_gate(expand_g, g1, final_g, vec_dim_2q(sim_type));
                        memcpy(g1.gm_real, final_g.gm_real, vec_length_2q(sim_type) * sizeof(ValType));
                        memcpy(g1.gm_imag, final_g.gm_imag, vec_length_2q(sim_type) * sizeof(ValType));
                        fused = true;
                        break;
                    }
                }
                if (fused == false) // did not fuse with any C4 gate, add it
                {
                    for (IdxType j = 0; j < n_qubits; j++)
                    {
                        canfuse[j * n_qubits + qubit] = false;
                        canfuse[qubit * n_qubits + j] = false;
                    }
                    GateType g(circuit_in[i]);
                    circuit_out.push_back(g);
                }
            }
            else if (circuit_in[i].op_name == get_op(2, sim_type)) // 2-qubit gate
            {
                IdxType qubit = circuit_in[i].qubit;
                IdxType ctrl = circuit_in[i].ctrl;
                for (IdxType j = 0; j < n_qubits; j++)
                {
                    canfuse[j * n_qubits + qubit] = false;
                    canfuse[qubit * n_qubits + j] = false;
                    canfuse[j * n_qubits + ctrl] = false;
                    canfuse[ctrl * n_qubits + j] = false;
                }
                GateType g(circuit_in[i]);
                circuit_out.push_back(g);
                canfuse[ctrl * n_qubits + qubit] = true;
                table[ctrl * n_qubits + qubit] = circuit_out.size() - 1; // this gate can fuse
            }
        }
        // clean
        delete[] table;
        delete[] canfuse;
    }

    template <typename GateType>
    void gate_fusion_2q_absorb_1q_backward(const std::vector<GateType> &circuit_in, std::vector<GateType> &circuit_out, const IdxType n_qubits, const SimType sim_type)
    {
        // prepare
        IdxType *table = new IdxType[n_qubits * n_qubits];
        bool *canfuse = new bool[n_qubits * n_qubits];
        for (IdxType i = 0; i < n_qubits * n_qubits; i++)
            table[i] = -1;
        for (IdxType i = 0; i < n_qubits * n_qubits; i++)
            canfuse[i] = false;
        circuit_out.clear();
        // forward parse
        for (IdxType i = circuit_in.size() - 1; i >= 0; i--)
        {

            if (circuit_in[i].op_name == M) // 1-qubit measure gate
            {
                IdxType qubit = circuit_in[i].qubit;
                for (IdxType j = 0; j < n_qubits; j++)
                {
                    canfuse[j * n_qubits + qubit] = false;
                    canfuse[qubit * n_qubits + j] = false;
                }
                GateType g(circuit_in[i]);
                circuit_out.push_back(g);
            }
            else if (circuit_in[i].op_name == RESET) // 1-qubit reset gate
            {
                IdxType qubit = circuit_in[i].qubit;
                for (IdxType j = 0; j < n_qubits; j++)
                {
                    canfuse[j * n_qubits + qubit] = false;
                    canfuse[qubit * n_qubits + j] = false;
                }
                GateType g(circuit_in[i]);
                circuit_out.push_back(g);
            }
            else if (circuit_in[i].op_name == MA) // all-qubit measure gate
            {
                for (IdxType j = 0; j < n_qubits; j++)
                    for (IdxType k = 0; k < n_qubits; k++)
                        canfuse[j * n_qubits + k] = false;
                GateType g(circuit_in[i]);
                circuit_out.push_back(g);
            }
            else if (circuit_in[i].op_name == get_op(1, sim_type)) // 1-qubit gate
            {
                IdxType qubit = circuit_in[i].qubit;
                bool fused = false;
                for (IdxType j = 0; j < n_qubits; j++)
                {
                    // target qubit can merge
                    if (canfuse[j * n_qubits + qubit] == true)
                    {
                        const GateType &g0 = circuit_in[i];
                        GateType &g1 = circuit_out[table[j * n_qubits + qubit]];
                        GateType expand_g(get_op(2, sim_type), -1, -1);
                        GateType final_g(get_op(2, sim_type), -1, -1);

                        expand_ctrl(g0, expand_g, j, sim_type);
                        // now ready to fuse
                        fuse_gate(g1, expand_g, final_g, vec_dim_2q(sim_type));
                        memcpy(g1.gm_real, final_g.gm_real, vec_length_2q(sim_type) * sizeof(ValType));
                        memcpy(g1.gm_imag, final_g.gm_imag, vec_length_2q(sim_type) * sizeof(ValType));

                        fused = true;
                        break;
                    }
                    // ctrl qubit can merge
                    if (canfuse[qubit * n_qubits + j] == true)
                    {
                        const GateType &g0 = circuit_in[i];
                        GateType &g1 = circuit_out[table[qubit * n_qubits + j]];
                        GateType expand_g(get_op(2, sim_type), -1, -1);
                        expand_qubit(g0, expand_g, j, sim_type);
                        // now ready to fuse
                        GateType final_g(get_op(2, sim_type), -1, -1);
                        fuse_gate(g1, expand_g, final_g, vec_dim_2q(sim_type));
                        memcpy(g1.gm_real, final_g.gm_real, vec_length_2q(sim_type) * sizeof(ValType));
                        memcpy(g1.gm_imag, final_g.gm_imag, vec_length_2q(sim_type) * sizeof(ValType));
                        fused = true;
                        break;
                    }
                }
                if (fused == false) // did not fuse with any C4 gate, add it
                {
                    for (IdxType j = 0; j < n_qubits; j++)
                    {
                        canfuse[j * n_qubits + qubit] = false;
                        canfuse[qubit * n_qubits + j] = false;
                    }
                    GateType g(circuit_in[i]);
                    circuit_out.push_back(g);
                }
            }
            if (circuit_in[i].op_name == get_op(2, sim_type)) // 2-qubit gate
            {
                IdxType qubit = circuit_in[i].qubit;
                IdxType ctrl = circuit_in[i].ctrl;
                for (IdxType j = 0; j < n_qubits; j++)
                {
                    canfuse[j * n_qubits + qubit] = false;
                    canfuse[qubit * n_qubits + j] = false;
                    canfuse[j * n_qubits + ctrl] = false;
                    canfuse[ctrl * n_qubits + j] = false;
                }
                GateType g(circuit_in[i]);
                circuit_out.push_back(g);

                canfuse[ctrl * n_qubits + qubit] = true;
                table[ctrl * n_qubits + qubit] = circuit_out.size() - 1; // this gate can fuse
            }
        }
        // clean
        delete[] table;
        delete[] canfuse;
        reverse(circuit_out.begin(), circuit_out.end());
    }

    template <typename GateType>
    std::vector<GateType> fuse_circuit_gates(std::vector<GateType> gates, IdxType n_qubits, SimType sim_type)
    {
        //====================== Fuse ========================
        std::vector<GateType> tmp1_circuit;
        std::vector<GateType> tmp2_circuit;
        std::vector<GateType> tmp3_circuit;
        std::vector<GateType> fused_circuit;

        gate_fusion_1q(gates, tmp1_circuit, n_qubits, sim_type);
        gate_fusion_2q_absorb_1q_forward(tmp1_circuit, tmp2_circuit, n_qubits, sim_type);
        gate_fusion_2q_absorb_1q_backward(tmp2_circuit, tmp3_circuit, n_qubits, sim_type);
        gate_fusion_2q(tmp3_circuit, fused_circuit, n_qubits, sim_type);

        return fused_circuit;
    }

    std::vector<SVGate> fuse_circuit_sv(Circuit *circuit)
    {
        std::vector<SVGate> gates = GateFactory::getInstance().getSVGates(circuit->get_gates());

        if (!Config::ENABLE_FUSION)
        {
            return gates;
        }
        IdxType n_qubits = circuit->num_qubits();

        return fuse_circuit_gates(gates, n_qubits, SV);
    }

    std::vector<DMGate> fuse_circuit_dm(Circuit *circuit)
    {
        std::vector<DMGate> gates = getDMGates(circuit->get_gates(), circuit->num_qubits());

        if (!Config::ENABLE_FUSION)
        {
            return gates;
        }
        IdxType n_qubits = circuit->num_qubits();

        std::vector<DMGate> tmp1_circuit;
        std::vector<DMGate> fused_circuit;

        gate_fusion_1q(gates, tmp1_circuit, n_qubits, DM);
        gate_fusion_2q(tmp1_circuit, fused_circuit, n_qubits, DM);

        return fused_circuit;
    }

} // end of NWQSim
