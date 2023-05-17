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
// File: circuit.hpp
// Define the quantum circuit class.
// ---------------------------------------------------------------------------

#ifndef CIRCUIT_HPP
#define CIRCUIT_HPP

#include "util.h"
#include "gate.h"
#include "config.h"
#include "fusion.h"

namespace NWQSim
{
    using namespace std;

    class Circuit
    {
    public:
        Circuit(IdxType _n_qubits = 0) : n_qubits(_n_qubits), n_gates(0), circuit_cpu(NULL)
        {
        }
        ~Circuit() { clear(); }
        void append(Gate &g)
        {
#ifdef PRINT_GATE_TRACE
            std::cout << g.gateToString() << std::flush;
#endif
            if (g.qubit >= n_qubits || ((g.op_name != MA) && (g.ctrl >= n_qubits)))
            {
                std::cerr << g.gateToString() << std::flush;
                string msg = "Gate uses qubit out of range of " + to_string(n_qubits) + " !\n";
                throw std::logic_error(msg.c_str());
            }
            circuit.push_back(g);
            delete (&g);
            n_gates++;
        }
        void AllocateQubit()
        {
            n_qubits++;
#ifdef PRINT_QUBIT_ALC_AND_RLS
            std::cout << "allocate 1 qubit, now in total:" << n_qubits << std::endl;
#endif
        }
        void ReleaseQubit()
        {
            n_qubits--;
#ifdef PRINT_QUBIT_ALC_AND_RLS
            std::cout << "release 1 qubit, now in total:" << n_qubits << std::endl;
#endif
        }
        void clear()
        {
            circuit.clear();
            n_gates = 0;
            SAFE_FREE_HOST(circuit_cpu);
        }
        void reset()
        {
            clear();
        }

        Gate *upload()
        {
#ifdef PRINT_CIRCUIT_METRICS
            circuit_metrics(circuit, n_qubits);
#endif
#ifdef DISABLE_GATE_FUSION
            //====================== No Fuse =====================
            SAFE_FREE_HOST(circuit_cpu);
            SAFE_ALOC_HOST(circuit_cpu, n_gates * sizeof(Gate));
            memcpy(circuit_cpu, circuit.data(), n_gates * sizeof(Gate));
            //====================================================
#else
            //====================== Fuse ========================
            vector<Gate> tmp1_circuit;
            vector<Gate> tmp2_circuit;
            vector<Gate> tmp3_circuit;
            tmp1_circuit.clear();
            tmp2_circuit.clear();

            gate_fusion_1q(circuit, tmp1_circuit, n_qubits);
            gate_fusion_2q_absorb_1q_forward(tmp1_circuit, tmp2_circuit, n_qubits);
            gate_fusion_2q_absorb_1q_backward(tmp2_circuit, tmp3_circuit, n_qubits);
            gate_fusion_2q(tmp3_circuit, fused_circuit, n_qubits);

            this->n_gates = fused_circuit.size();
            SAFE_FREE_HOST(circuit_cpu);
            SAFE_ALOC_HOST(circuit_cpu, n_gates * sizeof(Gate));
            memcpy(circuit_cpu, fused_circuit.data(), n_gates * sizeof(Gate));
            //====================================================
#endif
            return circuit_cpu;
        }
        std::string circuitToString()
        {
            stringstream ss;
            for (IdxType t = 0; t < n_gates; t++)
                ss << circuit[t].gateToString();
            return ss.str();
        }

    public:
        // number of qubits
        IdxType n_qubits;
        // number of gates
        IdxType n_gates;
        // user input gate sequence
        vector<Gate> circuit;
        // fused gate sequence
        vector<Gate> fused_circuit;
        Gate *circuit_cpu;
    };

}

#endif // CIRCUIT_HPP