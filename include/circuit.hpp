#pragma once

#include "nwq_util.hpp"
#include "gate.hpp"

#include <string>
#include <vector>
#include <sstream>
#include <memory>
#include <cmath>

namespace NWQSim
{
    struct CircuitMetrics {
        IdxType depth;
        IdxType one_q_gates;
        IdxType two_q_gates;
        ValType gate_density;
        ValType retention_lifespan; 
        ValType measurement_density;
        ValType entanglement_variance;
    };
    
    class Circuit
    {

    private:
        // number of qubits
        IdxType n_qubits;
        IdxType n_expect;

    public:
        // user input gate sequence
        std::shared_ptr<std::vector<Gate>> gates;

        Circuit(IdxType _n_qubits) : n_qubits(_n_qubits), n_expect(0)
        {
            // Implementation of constructor
            gates = std::make_shared<std::vector<Gate>>();
        }
        ~Circuit() {};

        IdxType num_qubits() { return n_qubits; };
        IdxType num_gates() { return gates->size(); };

        bool is_empty() { return gates->empty(); };

        std::vector<Gate> get_gates()
        {
            return *gates;
        }
        IdxType num_expect()
        {
            return n_expect;
        }
        void set_gates(std::vector<Gate> new_gates)
        {
            gates = std::make_shared<std::vector<Gate>>(new_gates);
        }
        void set_num_qubits(IdxType _n_qubits)
        {
            n_qubits = _n_qubits;
        }
        void clear()
        {
            // Implementation of clear function
            gates->clear();
            // n_qubits = 0;
        }
        void reset()
        {
            // Implementation of reset function
            clear();
        }
        std::string to_string()
        {
            // Implementation of to_string function
            std::stringstream ss;
            for (auto gate : *gates)
                ss << gate.gateToString() << std::endl;
            return ss.str();
        }

        CircuitMetrics circuit_metrics()
        {
            // Initialize a vector to track the current depth of each qubit
            std::vector<IdxType> qubit_depth(n_qubits, 0);
            // Initialize a vector to track the number of two-qubit gates applied to each qubit
            std::vector<IdxType> qubit_g2_gates(n_qubits, 0);

            // Track the total number of single-qubit gates and two-qubit gates in the circuit
            IdxType g1_gates = 0;
            IdxType g2_gates = 0;
            IdxType n_measure = 0;

            // Loop over each gate in the circuit and update the depth of the corresponding qubits
            IdxType max_depth = 0;
            IdxType one_q_gates = 0;
            IdxType two_q_gates = 0;
            for (IdxType i = 0; i < gates->size(); i++)
            {
                if (gates->at(i).op_name == OP::MA)
                {
                    n_measure += n_qubits;
                    continue;
                }
                if (gates->at(i).op_name == OP::M)
                {
                    n_measure++;
                    continue;
                }
                if (gates->at(i).op_name == OP::MOD_NOISE)
                {
                    continue; // Skip noise update gate
                }
                IdxType ctrl = gates->at(i).ctrl;
                IdxType target = gates->at(i).qubit;

                // Calculate the depth of this gate based on the depths of the control and target qubits
                IdxType depth;
                if (ctrl == -1)
                {
                    // Single-qubit gate
                    depth = qubit_depth[target] + 1;
                    g1_gates++;
                }
                else
                {
                    // Two-qubit gate
                    depth = std::max(qubit_depth[ctrl], qubit_depth[target]) + 1;
                    g2_gates++;

                    // Increment the number of two-qubit gates applied to each qubit
                    qubit_g2_gates[ctrl]++;
                    qubit_g2_gates[target]++;

                    if (ctrl == target)
                        printf("Exception: target==ctrl\n");
                }
                // Update the depth of the control and target qubits
                if (ctrl == -1)
                {
                    one_q_gates++;
                }
                else
                {
                    qubit_depth[ctrl] = depth;
                    two_q_gates++;
                }
                qubit_depth[target] = depth;
                // Update the maximum depth if the current depth is greater than the previous maximum
                if (depth > max_depth)
                {
                    max_depth = depth;
                }
            }
            // Calculate the gate density, retention lifespan, and entanglement variance of the circuit
            ValType gate_density = (g1_gates + 2 * g2_gates) / (ValType)(max_depth * n_qubits);
            ValType retention_lifespan = log(max_depth);
            ValType measurement_density = log((ValType)max_depth * n_qubits) / (ValType)n_measure;

            IdxType sum_g2_gates = 0;
            for (auto val : qubit_g2_gates)
            {
                sum_g2_gates += val;
            }
            ValType average_g2_gates = sum_g2_gates / (ValType)n_qubits;

            ValType entanglement_var = 0;
            for (auto val : qubit_g2_gates)
            {
                entanglement_var += log(pow(val - average_g2_gates, 2) + 1);
            }
            entanglement_var /= n_qubits;

            // Write to the metrics struct
            CircuitMetrics metrics;
            metrics.depth = max_depth;
            metrics.one_q_gates = one_q_gates;
            metrics.two_q_gates = two_q_gates;
            metrics.gate_density = gate_density;
            metrics.retention_lifespan = retention_lifespan;
            metrics.measurement_density = measurement_density;
            metrics.entanglement_variance = entanglement_var;
            return metrics;
        }

        void print_metrics()
        {
            CircuitMetrics metrics = circuit_metrics();
            IdxType max_depth = metrics.depth;
            IdxType one_q_gates = metrics.one_q_gates;
            IdxType two_q_gates = metrics.two_q_gates;
            ValType gate_density = metrics.gate_density;
            ValType retention_lifespan = metrics.retention_lifespan;
            ValType measurement_density = metrics.measurement_density;
            ValType entanglement_var = metrics.entanglement_variance;
            // Print the results to the console
            printf("Circuit Depth: %lld; One-qubit Gates: %lld; Two-qubit Gates: %lld; Gate Density: %.4lf; Retention Lifespan: %.4lf; Measurement Density: %.4lf; Entanglement Variance: %.4lf\n\n",
                   max_depth, one_q_gates, two_q_gates, gate_density, retention_lifespan, measurement_density, entanglement_var);
        }

        // ===================== Standard Gates =========================
        void MOD_NOISE(std::string mod_op, std::string mod_noise, ValType value, const std::vector<IdxType> &qubit_list)
        {
            Gate G(OP::MOD_NOISE, 0);

            G.mod_op = mod_op;
            G.mod_noise = mod_noise;
            G.mod_value = value;

            for (auto q : qubit_list)
            {
                G.mod_qubits.push_back(q);
            }

            gates->push_back(G);
        }

        void X(IdxType qubit)
        {
            // Pauli X-gate: bit flip
            /** X = [0 1]
                    [1 0]
             */
            Gate G(OP::X, qubit);
            gates->push_back(G);
        }
        void Y(IdxType qubit)
        {
            // Pauli-Y gate: bit and phase flip
            /** Y = [0 -i]
                    [i  0]
             */
            Gate G(OP::Y, qubit);
            gates->push_back(G);
        }
        void Z(IdxType qubit)
        {
            // Pauli-Z gate: phase flip
            /** Z = [1  0]
                    [0 -1]
             */
            Gate G(OP::Z, qubit);
            gates->push_back(G);
        }
        void H(IdxType qubit)
        {
            // Clifford gate: Hadamard
            /** H = 1/sqrt(2) * [1  1]
                                [1 -1]
             */
            Gate G(OP::H, qubit);
            gates->push_back(G);
        }
        void S(IdxType qubit)
        {
            // Clifford gate: sqrt(Z) phase gate
            /** S = [1 0]
                    [0 i]
            */
            Gate G(OP::S, qubit);
            gates->push_back(G);
        }
        void EXPECT(void *obsptr)
        {
            // Compute the expectation value of a list of observables (data stored in the ObservableList struct)
            Gate G(OP::EXPECT, 0, -1, 1, 0, 0, 0, 0, obsptr);
            gates->push_back(G);
            n_expect++;
        }
        void SDG(IdxType qubit)
        {
            // Clifford gate: conjugate of sqrt(Z) phase gate
            /** SDG = [1  0]
                      [0 -i]
            */
            Gate G(OP::SDG, qubit);
            gates->push_back(G);
        }
        void T(IdxType qubit)
        {
            // C3 gate: sqrt(S) phase gate
            /** T = [1 0]
                    [0 s2i+s2i*i]
            */
            Gate G(OP::T, qubit);
            gates->push_back(G);
        }
        void TDG(IdxType qubit)
        {
            // C3 gate: conjugate of sqrt(S) phase gate
            /** TDG = [1 0]
                      [0 s2i-s2i*i]
            */
            Gate G(OP::TDG, qubit);
            gates->push_back(G);
        }
        void RI(ValType theta, IdxType qubit)
        {
            // Global phase gate
            /** RI = [e^(ia) 0] = [cos(a)+i*sin(a) 0]
                     [0 e^(ia)]   [0 cos(a)+i*sin(a)]
            */
            Gate G(OP::RI, qubit, -1, 1, theta);
            gates->push_back(G);
        }
        void RX(ValType theta, IdxType qubit)
        {
            // Rotation around X axis
            /** RX = [cos(a/2) -i*sin(a/2)]
                     [-i*sin(a/2) cos(a/2)]
            */
            Gate G(OP::RX, qubit, -1, 1, theta);
            gates->push_back(G);
        }
        void RY(ValType theta, IdxType qubit)
        {
            // Rotation around Y axis
            /** RY = [cos(a/2) -sin(a/2)]
                     [sin(a/2)  cos(a/2)]
            */
            Gate G(OP::RY, qubit, -1, 1, theta);
            gates->push_back(G);
        }
        void RZ(ValType theta, IdxType qubit)
        {
            // Rotation around Z axis
            /** RZ = [cos(a/2)-i*sin(a/2)  0]
                     [0  cos(a/2)+i*sin(a/2)]
            */
            Gate G(OP::RZ, qubit, -1, 1, theta);
            gates->push_back(G);
        }

        void P(ValType theta, IdxType qubit)
        {
            // Phase gate defined by Qiskit
            /** P = [1, 0     ]  = [1,0]
                    [0, e^(ia)]    [0,cos(a)+i*sin(a)]
            */
            Gate G(OP::P, qubit, -1, 1, theta);
            gates->push_back(G);
        }
        void U(ValType theta, ValType phi, ValType lam, IdxType qubit)
        {
            // Generic single-qubit rotation gate with 3 Euler angles
            /** U = [cos(theta/2), -e^(i*lam)sin(theta/2)]
                    [e^(i*phi)sin(theta/2), e^(i*(phi+lam))cos(theta/2)]
            */
            Gate G(OP::U, qubit, -1, 1, theta, phi, lam);
            gates->push_back(G);
        }
        void CX(IdxType ctrl, IdxType qubit)
        {
            // Controlled-NOT or CNOT
            /**  CX   = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 0 1]
                        [0 0 1 0]
            */
            Gate G(OP::CX, qubit, ctrl, 2);
            gates->push_back(G);
        }
        void CY(IdxType ctrl, IdxType qubit)
        {
            // Controlled-Y
            /**  CY   = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 0 -i]
                        [0 0 i 0]
            */
            Gate G(OP::CY, qubit, ctrl, 2);
            gates->push_back(G);
        }
        void CZ(IdxType ctrl, IdxType qubit)
        {
            // Controlled-Z
            /**  CZ   = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 1 0]
                        [0 0 0 -1]
            */
            Gate G(OP::CZ, qubit, ctrl, 2);
            gates->push_back(G);
        }
        void CH(IdxType ctrl, IdxType qubit)
        {
            // Controlled-H
            /**  CH   = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 s2i s2i]
                        [0 0 s2i -s2i]
            */
            Gate G(OP::CH, qubit, ctrl, 2);
            gates->push_back(G);
        }
        void CS(IdxType ctrl, IdxType qubit)
        {
            // Controlled-S
            /**  CS   = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 1 0]
                        [0 0 0 i]
            */
            Gate G(OP::CS, qubit, ctrl, 2);
        }
        void CSDG(IdxType ctrl, IdxType qubit)
        {
            // Controlled-SDG
            /**  CSDG = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 1 0]
                        [0 0 0 -i]
            */
            Gate G(OP::CSDG, qubit, ctrl, 2);
            gates->push_back(G);
        }
        void CT(IdxType ctrl, IdxType qubit)
        {
            // Controlled-T
            /**  CT   = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 1 0]
                        [0 0 0 s2i+si2*i]
            */
            Gate G(OP::CT, qubit, ctrl, 2);
            gates->push_back(G);
        }
        void CTDG(IdxType ctrl, IdxType qubit)
        {
            // Controlled-TDG
            /**  CTDG = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 1 0]
                        [0 0 0 s2i-si2*i]
            */
            Gate G(OP::CTDG, qubit, ctrl, 2);
            gates->push_back(G);
        }
        void CRX(ValType theta, IdxType ctrl, IdxType qubit)
        {
            // Controlled-RX
            /**  CRX  = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 cos(a/2) -i*sin(a/2)]
                        [0 0 -i*sin(a/2) cos(a/2)]
            */
            Gate G(OP::CRX, qubit, ctrl, 2, theta);
            gates->push_back(G);
        }
        void CRY(ValType theta, IdxType ctrl, IdxType qubit)
        {
            // Controlled-RY
            /**  CRY  = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 cos(a/2) -sin(a/2)]
                        [0 0 sin(a/2)  cos(a/2)]
            */
            Gate G(OP::CRY, qubit, ctrl, 2, theta);
            gates->push_back(G);
        }
        void CRZ(ValType theta, IdxType ctrl, IdxType qubit)
        {
            // Controlled-RZ
            /**  CRZ  = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 cos(a/2)-i*sin(a/2)  0]
                        [0 0 0  cos(a/2)+i*sin(a/2)]
            */
            Gate G(OP::CRZ, qubit, ctrl, 2, theta);
            gates->push_back(G);
        }
        void CSX(IdxType ctrl, IdxType qubit)
        {
            // Controlled-SX
            /**  CSX  = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 (1+i)/2 (1-i)/2]
                        [0 0 (1-i)/2 (1+i)/2]
            */
            Gate G(OP::CSX, qubit, ctrl, 2);
            gates->push_back(G);
        }
        void CP(ValType theta, IdxType ctrl, IdxType qubit)
        {
            // Controlled-P
            /**  CP   = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 1 0]
                        [0 0 0 cos(a)+i*sin(a)]
            */
            Gate G(OP::CP, qubit, ctrl, 2, theta);
            gates->push_back(G);
        }
        void CU(ValType theta, ValType phi, ValType lam, ValType gamma,
                IdxType ctrl, IdxType qubit)
        {
            // Controlled-U, w.s.p. to Qiksik: https://qiskit.org/documentation/stubs/qiskit.circuit.library.CUGate.html
            /**  CU   = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 e^(i*gamma)cos(theta/2), -e^(i*(gamma+lam))sin(theta/2)]
                        [0 0 e^(i*(gamma+phi))sin(theta/2), e^(i*(gamma+phi+lam))cos(theta/2)]
            */
            Gate G(OP::CU, qubit, ctrl, 2, theta, phi, lam);
            gates->push_back(G);
        }
        void ECR(IdxType ctrl, IdxType qubit)
        {
            /******************************************
             * Echoed Cross-Resonance Gate
             * Implements 1/sqrt(2) (IX - XY)
             * ECR = 1/$\sqrt{2}$ * [ 0  1  0  i]
             *                      [ 1  0 -i  0]
             *                      [ 0  i  0  1]
             *                      [-i  0  1  0]
             ******************************************/
            Gate G(OP::ECR, qubit, ctrl, 2);
            gates->push_back(G);
        }
        void RXX(ValType theta, IdxType qubit0, IdxType qubit1)
        {
            // RXX, w.s.p. to Qiskit: https://qiskit.org/documentation/stubs/qiskit.circuit.library.RXXGate.html
            /**  CU   = [cos(theta/2)    0               0               -i*sin(theta/2)]
                        [0               cos(theta/2)    -i*sin(theta/2) 0]
                        [0               -i*sin(theta/2) cos(theta/2)    0]
                        [-i*sin(theta/2) 0               0               cos(theta/2)]
            */
            Gate G(OP::RXX, qubit0, qubit1, 2, theta);
            gates->push_back(G);
        }
        void RYY(ValType theta, IdxType qubit0, IdxType qubit1)
        {
            // RYY, w.s.p. to Qiksik: https://qiskit.org/documentation/stubs/qiskit.circuit.library.RYYGate.html#qiskit.circuit.library.RYYGate
            /**  CU   = [cos(theta/2)    0               0               i*sin(theta/2)]
                        [0               cos(theta/2)    -i*sin(theta/2) 0]
                        [0               -i*sin(theta/2) cos(theta/2)    0]
                        [i*sin(theta/2) 0               0               cos(theta/2)]
            */
            Gate G(OP::RYY, qubit0, qubit1, 2, theta);
            gates->push_back(G);
        }
        void RZZ(ValType theta, IdxType qubit0, IdxType qubit1)
        {
            // RZZ, w.s.p. to Qiksik: https://qiskit.org/documentation/stubs/qiskit.circuit.library.RZZGate.html#qiskit.circuit.library.RZZGate
            /**  CU   = [e^(-i theta/2)    0               0               0]
             *          [0                 e^(i theta/2)  0               0]
             *          [0                 0               e^(i theta/2)  0]
             *          [0                 0               0               e^(-i theta/2)]
             */
            Gate G(OP::RZZ, qubit0, qubit1, 2, theta);
            gates->push_back(G);
        }
        void SX(IdxType qubit)
        {
            // sqrt(X) gate, basis gate for IBMQ
            /** SX = 1/2 [1+i 1-i]
                         [1-i 1+i]
            */
            Gate G(OP::SX, qubit);
            gates->push_back(G);
        }
        void ID(IdxType qubit)
        {
            // Identity gate
            /** ID  = [1 0]
                      [0 1]
            */
            Gate G(OP::ID, qubit);
            gates->push_back(G);
        }
        void SWAP(IdxType ctrl, IdxType qubit)
        {
            // SWAP gate
            /**  SWAP = [1 0 0 0]
                        [0 0 1 0]
                        [0 1 0 0]
                        [0 0 0 1]
            */
            Gate G(OP::SWAP, qubit, ctrl, 2);
            gates->push_back(G);
        }
        void DELAY(ValType gate_len, IdxType qubit)
        {
            // Identity with user-defined delay and relaxation (stored in theta)
            /** Delay  = [1 0]
                         [0 1]
            */
            Gate G(OP::DELAY, qubit, -1, 1, gate_len);
            gates->push_back(G);
        }
        void M(IdxType qubit) // default is pauli-Z
        {
            Gate G(OP::M, qubit);
            gates->push_back(G);
        }
        void MA(IdxType repetition) // default is pauli-Z
        {
            Gate G(OP::MA, -1, -1, 1, 0, 0, 0, repetition);
            gates->push_back(G);
        }
        void RESET(IdxType qubit)
        {
            Gate G(OP::RESET, qubit);
            gates->push_back(G);
        }

        // ============================== Other Gate Definition ================================
        void U3(ValType theta, ValType phi, ValType lam, IdxType qubit)
        {
            U(theta, phi, lam, qubit);
        }
        void U2(ValType phi, ValType lam, IdxType qubit)
        {
            U(PI / 2, phi, lam, qubit);
        }
        void U1(ValType lam, IdxType qubit)
        {
            U(0, 0, lam, qubit);
        }

        void CCX(IdxType qubit0, IdxType qubit1, IdxType qubit2)
        {
            H(qubit2);
            CX(qubit1, qubit2);
            TDG(qubit2);
            CX(qubit0, qubit2);
            T(qubit2);
            CX(qubit1, qubit2);
            T(qubit1);
            TDG(qubit2);
            CX(qubit0, qubit2);
            T(qubit2);
            CX(qubit0, qubit1);
            T(qubit0);
            TDG(qubit1);
            H(qubit2);
            CX(qubit0, qubit1);
        }
        void CSWAP(IdxType qubit0, IdxType qubit1, IdxType qubit2)
        {
            CX(qubit2, qubit1);
            CCX(qubit0, qubit1, qubit2);
            CX(qubit2, qubit1);
        }
    };

} // namespace NWQSim
