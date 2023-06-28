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
    class Circuit
    {

    private:
        // number of qubits
        IdxType n_qubits;

    public:
        // user input gate sequence
        std::shared_ptr<std::vector<Gate>> gates;
        map<string, creg> list_cregs;
        std::vector<int> initial_mapping;
        Circuit(IdxType _n_qubits) : n_qubits(_n_qubits)
        {
            // Implementation of constructor
            gates = std::make_shared<std::vector<Gate>>();
        }
        ~Circuit(){};

        IdxType num_qubits() { return n_qubits; };
        IdxType num_gates() { return gates->size(); };

        bool is_empty() { return gates->empty(); };

        std::vector<Gate> get_gates()
        {
            return *gates;
        }
        void set_gates(std::vector<Gate> new_gates)
        {
            gates = std::make_shared<std::vector<Gate>>(new_gates);
        }
        void set_creg(map<string, creg> list_cregs)
        {
            this->list_cregs = list_cregs;
        }
        void set_mapping(std::vector<int> initial_mapping)
        {
            for (auto ini : initial_mapping)
            {
                this->initial_mapping.push_back(ini);
            }
        }
        std::vector<int> get_mapping()
        {
            return this->initial_mapping;
        }
        map<string, creg> get_cregs()
        {
            return this->list_cregs;
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

        // ===================== Standard Gates =========================

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
        void RXX(ValType theta, IdxType qubit0, IdxType qubit1)
        {
            // RXX, w.s.p. to Qiksik: https://qiskit.org/documentation/stubs/qiskit.circuit.library.RXXGate.html
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