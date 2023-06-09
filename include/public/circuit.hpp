#pragma once

#include "util.hpp"
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
        std::shared_ptr<std::vector<Gate>>
            gates;

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
            Gate *G = new Gate(OP::X, qubit);
            gates->push_back(*G);
        }
        void Y(IdxType qubit)
        {
            Gate *G = new Gate(OP::Y, qubit);
            gates->push_back(*G);
        }
        void Z(IdxType qubit)
        {
            Gate *G = new Gate(OP::Z, qubit);
            gates->push_back(*G);
        }
        void H(IdxType qubit)
        {
            Gate *G = new Gate(OP::H, qubit);
            gates->push_back(*G);
        }
        void S(IdxType qubit)
        {
            Gate *G = new Gate(OP::S, qubit);
            gates->push_back(*G);
        }
        void SDG(IdxType qubit)
        {
            Gate *G = new Gate(OP::SDG, qubit);
            gates->push_back(*G);
        }
        void T(IdxType qubit)
        {
            Gate *G = new Gate(OP::T, qubit);
            gates->push_back(*G);
        }
        void TDG(IdxType qubit)
        {
            Gate *G = new Gate(OP::TDG, qubit);
            gates->push_back(*G);
        }
        void RI(ValType theta, IdxType qubit)
        {
            Gate *G = new Gate(OP::RI, qubit, -1, 1, theta);
            gates->push_back(*G);
        }
        void RX(ValType theta, IdxType qubit)
        {
            Gate *G = new Gate(OP::RX, qubit, -1, 1, theta);
            gates->push_back(*G);
        }
        void RY(ValType theta, IdxType qubit)
        {
            Gate *G = new Gate(OP::RY, qubit, -1, 1, theta);
            gates->push_back(*G);
        }
        void RZ(ValType theta, IdxType qubit)
        {
            Gate *G = new Gate(OP::RZ, qubit, -1, 1, theta);
            gates->push_back(*G);
        }

        void P(ValType theta, IdxType qubit)
        {
            Gate *G = new Gate(OP::P, qubit, -1, 1, theta);
            gates->push_back(*G);
        }
        void U(ValType theta, ValType phi, ValType lam, IdxType qubit)
        {
            Gate *G = new Gate(OP::U, qubit, -1, 1, theta, phi, lam);
            gates->push_back(*G);
        }
        void CX(IdxType ctrl, IdxType qubit)
        {
            Gate *G = new Gate(OP::CX, qubit, ctrl, 2);
            gates->push_back(*G);
        }
        void CY(IdxType ctrl, IdxType qubit)
        {
            Gate *G = new Gate(OP::CY, qubit, ctrl, 2);
            gates->push_back(*G);
        }
        void CZ(IdxType ctrl, IdxType qubit)
        {
            Gate *G = new Gate(OP::CZ, qubit, ctrl, 2);
            gates->push_back(*G);
        }
        void CH(IdxType ctrl, IdxType qubit)
        {
            Gate *G = new Gate(OP::CH, qubit, ctrl, 2);
            gates->push_back(*G);
        }
        void CS(IdxType ctrl, IdxType qubit)
        {
            Gate *G = new Gate(OP::CS, qubit, ctrl, 2);
            gates->push_back(*G);
        }
        void CSDG(IdxType ctrl, IdxType qubit)
        {
            Gate *G = new Gate(OP::CSDG, qubit, ctrl, 2);
            gates->push_back(*G);
        }
        void CT(IdxType ctrl, IdxType qubit)
        {
            Gate *G = new Gate(OP::CT, qubit, ctrl, 2);
            gates->push_back(*G);
        }
        void CTDG(IdxType ctrl, IdxType qubit)
        {
            Gate *G = new Gate(OP::CTDG, qubit, ctrl, 2);
            gates->push_back(*G);
        }
        void CRX(ValType theta, IdxType ctrl, IdxType qubit)
        {
            Gate *G = new Gate(OP::CRX, qubit, ctrl, 2, theta);
            gates->push_back(*G);
        }
        void CRY(ValType theta, IdxType ctrl, IdxType qubit)
        {
            Gate *G = new Gate(OP::CRY, qubit, ctrl, 2, theta);
            gates->push_back(*G);
        }
        void CRZ(ValType theta, IdxType ctrl, IdxType qubit)
        {
            Gate *G = new Gate(OP::CRZ, qubit, ctrl, 2, theta);
            gates->push_back(*G);
        }
        void CSX(IdxType ctrl, IdxType qubit)
        {
            Gate *G = new Gate(OP::CSX, qubit, ctrl, 2);
            gates->push_back(*G);
        }
        void CP(ValType theta, IdxType ctrl, IdxType qubit)
        {
            Gate *G = new Gate(OP::CP, qubit, ctrl, 2, theta);
            gates->push_back(*G);
        }
        void CU(ValType theta, ValType phi, ValType lam, ValType gamma,
                IdxType ctrl, IdxType qubit)
        {
            Gate *G = new Gate(OP::CU, qubit, ctrl, 2, theta, phi, lam, gamma);
            gates->push_back(*G);
        }
        void RXX(ValType theta, IdxType qubit0, IdxType qubit1)
        {
            Gate *G = new Gate(OP::RXX, qubit0, qubit1, 2, theta);
            gates->push_back(*G);
        }
        void RYY(ValType theta, IdxType qubit0, IdxType qubit1)
        {
            Gate *G = new Gate(OP::RYY, qubit0, qubit1, 2, theta);
            gates->push_back(*G);
        }
        void RZZ(ValType theta, IdxType qubit0, IdxType qubit1)
        {
            Gate *G = new Gate(OP::RZZ, qubit0, qubit1, 2, theta);
            gates->push_back(*G);
        }

        void SX(IdxType qubit)
        {
            Gate *G = new Gate(OP::SX, qubit);
            gates->push_back(*G);
        }

        void ID(IdxType qubit)
        {
            Gate *G = new Gate(OP::ID, qubit);
            gates->push_back(*G);
        }

        void SWAP(IdxType ctrl, IdxType qubit)
        {
            Gate *G = new Gate(OP::SWAP, qubit, ctrl, 2);
            gates->push_back(*G);
        }

        void M(IdxType qubit) // default is pauli-Z
        {
            Gate *G = new Gate(OP::M, qubit);
            gates->push_back(*G);
        }
        void MA(IdxType repetition) // default is pauli-Z
        {
            Gate *G = new Gate(OP::MA, -1, -1, 1, 0, 0, 0, 0, repetition);
            gates->push_back(*G);
        }
        void RESET(IdxType qubit)
        {
            Gate *G = new Gate(OP::RESET, qubit);
            gates->push_back(*G);
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
