#include "../include/circuit.hpp"
#include "../include/gate.hpp"
#include "../include/util.hpp"

#include <algorithm>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

namespace NWQSim
{
    Circuit::Circuit()
    {
        // Implementation of constructor
        Circuit::gates = std::make_shared<std::vector<Gate>>();
    }

    Circuit::~Circuit()
    {
        // Implementation of destructor
    }

    IdxType Circuit::num_qubits()
    {
        return Circuit::n_qubits + 1;
    }

    IdxType Circuit::num_gates()
    {
        return Circuit::gates->size();
    }

    std::vector<Gate> Circuit::get_gates()
    {
        return *Circuit::gates;
    }

    void Circuit::set_gates(std::vector<Gate> new_gates)
    {
        Circuit::gates = std::make_shared<std::vector<Gate>>(new_gates);
    }

    void Circuit::clear()
    {
        // Implementation of clear function
        Circuit::gates->clear();
        Circuit::n_qubits = 0;
    }

    void Circuit::reset()
    {
        // Implementation of reset function
        Circuit::clear();
    }

    std::string Circuit::circuitToString()
    {
        // Implementation of circuitToString function
        std::stringstream ss;
        for (auto gate : *Circuit::gates)
            ss << gate.gateToString() << std::endl;
        return ss.str();
    }

    // =============================== Standard Gates ===================================
    void Circuit::X(IdxType qubit)
    {
        // Pauli X-gate: bit flip
        /** X = [0 1]
                [1 0]
         */
        Gate *G = new Gate(OP::X, qubit);
        ValType gm_real[4] = {0, 1, 1, 0};
        ValType gm_imag[4] = {0, 0, 0, 0};
        G->set_gm(gm_real, gm_imag, 2);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(qubit, Circuit::n_qubits);
    }
    void Circuit::Y(IdxType qubit)
    {
        // Pauli-Y gate: bit and phase flip
        /** Y = [0 -i]
                [i  0]
         */
        Gate *G = new Gate(OP::Y, qubit);
        ValType gm_real[4] = {0, 0, 0, 0};
        ValType gm_imag[4] = {0, -1, 1, 0};
        G->set_gm(gm_real, gm_imag, 2);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(qubit, Circuit::n_qubits);
    }
    void Circuit::Z(IdxType qubit)
    {
        // Pauli-Z gate: phase flip
        /** Z = [1  0]
                [0 -1]
         */
        Gate *G = new Gate(OP::Z, qubit);
        ValType gm_real[4] = {1, 0, 0, -1};
        ValType gm_imag[4] = {0, 0, 0, 0};
        G->set_gm(gm_real, gm_imag, 2);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(qubit, Circuit::n_qubits);
    }
    void Circuit::H(IdxType qubit)
    {
        // Clifford gate: Hadamard
        /** H = 1/sqrt(2) * [1  1]
                            [1 -1]
         */
        Gate *G = new Gate(OP::H, qubit);
        ValType gm_real[4] = {S2I, S2I, S2I, -S2I};
        ValType gm_imag[4] = {0, 0, 0, 0};
        G->set_gm(gm_real, gm_imag, 2);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(qubit, Circuit::n_qubits);
    }
    void Circuit::S(IdxType qubit)
    {
        // Clifford gate: sqrt(Z) phase gate
        /** S = [1 0]
                [0 i]
        */
        Gate *G = new Gate(OP::S, qubit);
        ValType gm_real[4] = {1, 0, 0, 0};
        ValType gm_imag[4] = {0, 0, 0, 1};
        G->set_gm(gm_real, gm_imag, 2);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(qubit, Circuit::n_qubits);
    }
    void Circuit::SDG(IdxType qubit)
    {
        // Clifford gate: conjugate of sqrt(Z) phase gate
        /** SDG = [1  0]
                  [0 -i]
        */
        Gate *G = new Gate(OP::SDG, qubit);
        ValType gm_real[4] = {1, 0, 0, 0};
        ValType gm_imag[4] = {0, 0, 0, -1};
        G->set_gm(gm_real, gm_imag, 2);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(qubit, Circuit::n_qubits);
    }
    void Circuit::T(IdxType qubit)
    {
        // C3 gate: sqrt(S) phase gate
        /** T = [1 0]
                [0 s2i+s2i*i]
        */
        Gate *G = new Gate(OP::T, qubit);
        ValType gm_real[4] = {1, 0, 0, S2I};
        ValType gm_imag[4] = {0, 0, 0, S2I};
        G->set_gm(gm_real, gm_imag, 2);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(qubit, Circuit::n_qubits);
    }
    void Circuit::TDG(IdxType qubit)
    {
        // C3 gate: conjugate of sqrt(S) phase gate
        /** TDG = [1 0]
                  [0 s2i-s2i*i]
        */
        Gate *G = new Gate(OP::TDG, qubit);
        ValType gm_real[4] = {1, 0, 0, S2I};
        ValType gm_imag[4] = {0, 0, 0, -S2I};
        G->set_gm(gm_real, gm_imag, 2);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(qubit, Circuit::n_qubits);
    }
    void Circuit::RI(ValType theta, IdxType qubit)
    {
        // Global phase gate
        /** RI = [e^(ia) 0] = [cos(a)+i*sin(a) 0]
                 [0 e^(ia)]   [0 cos(a)+i*sin(a)]
        */
        Gate *G = new Gate(OP::RI, qubit, -1, 1, theta);
        ValType gm_real[4] = {cos(theta), 0, 0, cos(theta)};
        ValType gm_imag[4] = {sin(theta), 0, 0, sin(theta)};
        G->set_gm(gm_real, gm_imag, 2);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(qubit, Circuit::n_qubits);
    }
    void Circuit::RX(ValType theta, IdxType qubit)
    {
        // Rotation around X axis
        /** RX = [cos(a/2) -i*sin(a/2)]
                 [-i*sin(a/2) cos(a/2)]
        */
        Gate *G = new Gate(OP::RX, qubit, -1, 1, theta);
        ValType gm_real[4] = {cos(HALF * theta), 0, 0, cos(HALF * theta)};
        ValType gm_imag[4] = {0, -sin(HALF * theta), -sin(HALF * theta), 0};
        G->set_gm(gm_real, gm_imag, 2);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(qubit, Circuit::n_qubits);
    }
    void Circuit::RY(ValType theta, IdxType qubit)
    {
        // Rotation around Y axis
        /** RY = [cos(a/2) -sin(a/2)]
                 [sin(a/2)  cos(a/2)]
        */
        Gate *G = new Gate(OP::RY, qubit, -1, 1, theta);
        ValType gm_real[4] = {cos(HALF * theta), -sin(HALF * theta), sin(HALF * theta), cos(HALF * theta)};
        ValType gm_imag[4] = {0};
        G->set_gm(gm_real, gm_imag, 2);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(qubit, Circuit::n_qubits);
    }
    void Circuit::RZ(ValType theta, IdxType qubit)
    {
        // Rotation around Z axis
        /** RZ = [cos(a/2)-i*sin(a/2)  0]
                 [0  cos(a/2)+i*sin(a/2)]
        */
        Gate *G = new Gate(OP::RZ, qubit, -1, 1, theta);
        ValType gm_real[4] = {cos(HALF * theta), 0, 0, cos(HALF * theta)};
        ValType gm_imag[4] = {-sin(HALF * theta), 0, 0, sin(HALF * theta)};
        G->set_gm(gm_real, gm_imag, 2);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(qubit, Circuit::n_qubits);
    }

    void Circuit::P(ValType theta, IdxType qubit)
    {
        // Phase gate defined by Qiskit
        /** P = [1, 0     ]  = [1,0]
                [0, e^(ia)]    [0,cos(a)+i*sin(a)]
        */
        Gate *G = new Gate(OP::P, qubit, -1, 1, theta);
        ValType gm_real[4] = {1, 0, 0, cos(theta)};
        ValType gm_imag[4] = {0, 0, 0, sin(theta)};
        G->set_gm(gm_real, gm_imag, 2);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(qubit, Circuit::n_qubits);
    }
    void Circuit::U(ValType theta, ValType phi, ValType lam, IdxType qubit)
    {
        // Generic single-qubit rotation gate with 3 Euler angles
        /** U = [cos(theta/2), -e^(i*lam)sin(theta/2)]
                [e^(i*phi)sin(theta/2), e^(i*(phi+lam))cos(theta/2)]
        */
        Gate *G = new Gate(OP::U, qubit, -1, 1, theta, phi, lam);
        ValType gm_real[4] = {cos(HALF * theta),
                              -cos(lam) * sin(HALF * theta),
                              cos(phi) * sin(HALF * theta),
                              cos(phi + lam) * cos(HALF * theta)};
        ValType gm_imag[4] = {0,
                              -sin(lam) * sin(HALF * theta),
                              sin(phi) * sin(HALF * theta),
                              sin(lam + phi) * cos(HALF * theta)};
        G->set_gm(gm_real, gm_imag, 2);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(qubit, Circuit::n_qubits);
    }
    void Circuit::CX(IdxType ctrl, IdxType qubit)
    {
        // Controlled-NOT or CNOT
        /**  CX   = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 0 1]
                    [0 0 1 0]
        */
        Gate *G = new Gate(OP::CX, qubit, ctrl, 2);
        ValType gm_real[16] = {1, 0, 0, 0,
                               0, 1, 0, 0,
                               0, 0, 0, 1,
                               0, 0, 1, 0};
        ValType gm_imag[16] = {0};
        G->set_gm(gm_real, gm_imag, 4);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(std::max(qubit, ctrl), Circuit::n_qubits);
    }
    void Circuit::CY(IdxType ctrl, IdxType qubit)
    {
        // Controlled-Y
        /**  CY   = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 0 -i]
                    [0 0 i 0]
        */
        Gate *G = new Gate(OP::CY, qubit, ctrl, 2);
        ValType gm_real[16] = {1, 0, 0, 0,
                               0, 1, 0, 0,
                               0, 0, 0, 0,
                               0, 0, 0, 0};
        ValType gm_imag[16] = {0, 0, 0, 0,
                               0, 0, 0, 0,
                               0, 0, 0, -1,
                               0, 0, 1, 0};
        G->set_gm(gm_real, gm_imag, 4);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(std::max(qubit, ctrl), Circuit::n_qubits);
    }
    void Circuit::CZ(IdxType ctrl, IdxType qubit)
    {
        // Controlled-Z
        /**  CZ   = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 1 0]
                    [0 0 0 -1]
        */
        Gate *G = new Gate(OP::CZ, qubit, ctrl, 2);
        ValType gm_real[16] = {1, 0, 0, 0,
                               0, 1, 0, 0,
                               0, 0, 1, 0,
                               0, 0, 0, -1};
        ValType gm_imag[16] = {0};
        G->set_gm(gm_real, gm_imag, 4);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(std::max(qubit, ctrl), Circuit::n_qubits);
    }
    void Circuit::CH(IdxType ctrl, IdxType qubit)
    {
        // Controlled-H
        /**  CH   = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 s2i s2i]
                    [0 0 s2i -s2i]
        */
        Gate *G = new Gate(OP::CH, qubit, ctrl, 2);
        ValType gm_real[16] = {1, 0, 0, 0,
                               0, 1, 0, 0,
                               0, 0, S2I, S2I,
                               0, 0, S2I, -S2I};
        ValType gm_imag[16] = {0};
        G->set_gm(gm_real, gm_imag, 4);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(std::max(qubit, ctrl), Circuit::n_qubits);
    }
    void Circuit::CS(IdxType ctrl, IdxType qubit)
    {
        // Controlled-S
        /**  CS   = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 1 0]
                    [0 0 0 i]
        */
        Gate *G = new Gate(OP::CS, qubit, ctrl, 2);
        ValType gm_real[16] = {1, 0, 0, 0,
                               0, 1, 0, 0,
                               0, 0, 1, 0,
                               0, 0, 0, 0};
        ValType gm_imag[16] = {0, 0, 0, 0,
                               0, 0, 0, 0,
                               0, 0, 0, 0,
                               0, 0, 0, 1};
        G->set_gm(gm_real, gm_imag, 4);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(std::max(qubit, ctrl), Circuit::n_qubits);
    }
    void Circuit::CSDG(IdxType ctrl, IdxType qubit)
    {
        // Controlled-SDG
        /**  CSDG = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 1 0]
                    [0 0 0 -i]
        */
        Gate *G = new Gate(OP::CSDG, qubit, ctrl, 2);
        ValType gm_real[16] = {1, 0, 0, 0,
                               0, 1, 0, 0,
                               0, 0, 1, 0,
                               0, 0, 0, 0};
        ValType gm_imag[16] = {0, 0, 0, 0,
                               0, 0, 0, 0,
                               0, 0, 0, 0,
                               0, 0, 0, -1};
        G->set_gm(gm_real, gm_imag, 4);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(std::max(qubit, ctrl), Circuit::n_qubits);
    }
    void Circuit::CT(IdxType ctrl, IdxType qubit)
    {
        // Controlled-T
        /**  CT   = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 1 0]
                    [0 0 0 s2i+si2*i]
        */
        Gate *G = new Gate(OP::CT, qubit, ctrl, 2);
        ValType gm_real[16] = {1, 0, 0, 0,
                               0, 1, 0, 0,
                               0, 0, 1, 0,
                               0, 0, 0, S2I};
        ValType gm_imag[16] = {0, 0, 0, 0,
                               0, 0, 0, 0,
                               0, 0, 0, 0,
                               0, 0, 0, S2I};
        G->set_gm(gm_real, gm_imag, 4);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(std::max(qubit, ctrl), Circuit::n_qubits);
    }
    void Circuit::CTDG(IdxType ctrl, IdxType qubit)
    {
        // Controlled-TDG
        /**  CTDG = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 1 0]
                    [0 0 0 s2i-si2*i]
        */
        Gate *G = new Gate(OP::CTDG, qubit, ctrl, 2);
        ValType gm_real[16] = {1, 0, 0, 0,
                               0, 1, 0, 0,
                               0, 0, 1, 0,
                               0, 0, 0, S2I};
        ValType gm_imag[16] = {0, 0, 0, 0,
                               0, 0, 0, 0,
                               0, 0, 0, 0,
                               0, 0, 0, -S2I};
        G->set_gm(gm_real, gm_imag, 4);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(std::max(qubit, ctrl), Circuit::n_qubits);
    }
    void Circuit::CRX(ValType theta, IdxType ctrl, IdxType qubit)
    {
        // Controlled-RX
        /**  CRX  = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 cos(a/2) -i*sin(a/2)]
                    [0 0 -i*sin(a/2) cos(a/2)]
        */
        Gate *G = new Gate(OP::CRX, qubit, ctrl, 2, theta);
        ValType gm_real[16] = {1, 0, 0, 0,
                               0, 1, 0, 0,
                               0, 0, cos(HALF * theta), 0,
                               0, 0, 0, cos(HALF * theta)};
        ValType gm_imag[16] = {0, 0, 0, 0,
                               0, 0, 0, 0,
                               0, 0, 0, -sin(HALF * theta),
                               0, 0, -sin(HALF * theta), 0};
        G->set_gm(gm_real, gm_imag, 4);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(std::max(qubit, ctrl), Circuit::n_qubits);
    }
    void Circuit::CRY(ValType theta, IdxType ctrl, IdxType qubit)
    {
        // Controlled-RY
        /**  CRY  = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 cos(a/2) -sin(a/2)]
                    [0 0 sin(a/2)  cos(a/2)]
        */
        Gate *G = new Gate(OP::CRY, qubit, ctrl, 2, theta);
        ValType gm_real[16] = {1, 0, 0, 0,
                               0, 1, 0, 0,
                               0, 0, cos(HALF * theta), -sin(HALF * theta),
                               0, 0, sin(HALF * theta), cos(HALF * theta)};
        ValType gm_imag[16] = {0};
        G->set_gm(gm_real, gm_imag, 4);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(std::max(qubit, ctrl), Circuit::n_qubits);
    }
    void Circuit::CRZ(ValType theta, IdxType ctrl, IdxType qubit)
    {
        // Controlled-RZ
        /**  CRZ  = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 cos(a/2)-i*sin(a/2)  0]
                    [0 0 0  cos(a/2)+i*sin(a/2)]
        */
        Gate *G = new Gate(OP::CRZ, qubit, ctrl, 2, theta);
        ValType gm_real[16] = {1, 0, 0, 0,
                               0, 1, 0, 0,
                               0, 0, cos(HALF * theta), 0,
                               0, 0, 0, cos(HALF * theta)};
        ValType gm_imag[16] = {0, 0, 0, 0,
                               0, 0, 0, 0,
                               0, 0, -sin(HALF * theta), 0,
                               0, 0, 0, sin(HALF * theta)};
        G->set_gm(gm_real, gm_imag, 4);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(std::max(qubit, ctrl), Circuit::n_qubits);
    }
    void Circuit::CSX(IdxType ctrl, IdxType qubit)
    {
        // Controlled-SX
        /**  CSX  = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 (1+i)/2 (1-i)/2]
                    [0 0 (1-i)/2 (1+i)/2]
        */
        Gate *G = new Gate(OP::CSX, qubit, ctrl, 2);
        ValType gm_real[16] = {1, 0, 0, 0,
                               0, 1, 0, 0,
                               0, 0, HALF, HALF,
                               0, 0, HALF, HALF};
        ValType gm_imag[16] = {0, 0, 0, 0,
                               0, 0, 0, 0,
                               0, 0, HALF, -HALF,
                               0, 0, -HALF, HALF};
        G->set_gm(gm_real, gm_imag, 4);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(std::max(qubit, ctrl), Circuit::n_qubits);
    }
    void Circuit::CP(ValType theta, IdxType ctrl, IdxType qubit)
    {
        // Controlled-P
        /**  CP   = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 1 0]
                    [0 0 0 cos(a)+i*sin(a)]
        */
        Gate *G = new Gate(OP::CP, qubit, ctrl, 2, theta);
        ValType gm_real[16] = {1, 0, 0, 0,
                               0, 1, 0, 0,
                               0, 0, 1, 0,
                               0, 0, 0, cos(theta)};
        ValType gm_imag[16] = {0, 0, 0, 0,
                               0, 0, 0, 0,
                               0, 0, 0, 0,
                               0, 0, 0, sin(theta)};
        G->set_gm(gm_real, gm_imag, 4);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(std::max(qubit, ctrl), Circuit::n_qubits);
    }
    void Circuit::CU(ValType theta, ValType phi, ValType lam, ValType gamma,
                     IdxType ctrl, IdxType qubit)
    {
        // Controlled-U, w.s.p. to Qiksik: https://qiskit.org/documentation/stubs/qiskit.circuit.library.CUGate.html
        /**  CU   = [1 0 0 0]
                    [0 1 0 0]
                    [0 0 e^(i*gamma)cos(theta/2), -e^(i*(gamma+lam))sin(theta/2)]
                    [0 0 e^(i*(gamma+phi))sin(theta/2), e^(i*(gamma+phi+lam))cos(theta/2)]
        */
        Gate *G = new Gate(OP::CU, qubit, ctrl, 2, theta, phi, lam);
        ValType gm_real[16] = {1, 0, 0, 0,
                               0, 1, 0, 0,
                               0, 0, cos(gamma) * cos(HALF * theta), -cos(gamma + lam) * sin(HALF * theta),
                               0, 0, cos(gamma + phi) * sin(HALF * theta), cos(gamma + phi + lam) * cos(HALF * theta)};
        ValType gm_imag[16] = {0, 0, 0, 0,
                               0, 0, 0, 0,
                               0, 0, sin(gamma) * cos(HALF * theta), -sin(gamma + lam) * sin(HALF * theta),
                               0, 0, sin(gamma + phi) * sin(HALF * theta), sin(gamma + phi + lam) * cos(HALF * theta)};
        G->set_gm(gm_real, gm_imag, 4);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(std::max(qubit, ctrl), Circuit::n_qubits);
    }
    void Circuit::RXX(ValType theta, IdxType qubit0, IdxType qubit1)
    {
        // RXX, w.s.p. to Qiksik: https://qiskit.org/documentation/stubs/qiskit.circuit.library.RXXGate.html
        /**  CU   = [cos(theta/2)    0               0               -i*sin(theta/2)]
                    [0               cos(theta/2)    -i*sin(theta/2) 0]
                    [0               -i*sin(theta/2) cos(theta/2)    0]
                    [-i*sin(theta/2) 0               0               cos(theta/2)]
        */
        Gate *G = new Gate(OP::RXX, qubit0, qubit1, 2, theta);
        ValType gm_real[16] = {cos(HALF * theta), 0, 0, 0,
                               0, cos(HALF * theta), 0, 0,
                               0, 0, cos(HALF * theta), 0,
                               0, 0, 0, cos(HALF * theta)};

        ValType gm_imag[16] = {0, 0, 0, -sin(HALF * theta),
                               0, 0, -sin(HALF * theta), 0,
                               0, -sin(HALF * theta), 0, 0,
                               -sin(HALF * theta), 0, 0, 0};

        G->set_gm(gm_real, gm_imag, 4);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(std::max(qubit0, qubit1), Circuit::n_qubits);
    }
    void Circuit::RYY(ValType theta, IdxType qubit0, IdxType qubit1)
    {
        // RYY, w.s.p. to Qiksik: https://qiskit.org/documentation/stubs/qiskit.circuit.library.RYYGate.html#qiskit.circuit.library.RYYGate
        /**  CU   = [cos(theta/2)    0               0               i*sin(theta/2)]
                    [0               cos(theta/2)    -i*sin(theta/2) 0]
                    [0               -i*sin(theta/2) cos(theta/2)    0]
                    [i*sin(theta/2) 0               0               cos(theta/2)]
        */
        Gate *G = new Gate(OP::RYY, qubit0, qubit1, 2, theta);
        ValType gm_real[16] = {cos(HALF * theta), 0, 0, 0,
                               0, cos(HALF * theta), 0, 0,
                               0, 0, cos(HALF * theta), 0,
                               0, 0, 0, cos(HALF * theta)};

        ValType gm_imag[16] = {0, 0, 0, sin(HALF * theta),
                               0, 0, -sin(HALF * theta), 0,
                               0, -sin(HALF * theta), 0, 0,
                               sin(HALF * theta), 0, 0, 0};

        G->set_gm(gm_real, gm_imag, 4);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(std::max(qubit0, qubit1), Circuit::n_qubits);
    }
    void Circuit::RZZ(ValType theta, IdxType qubit0, IdxType qubit1)
    {
        // RZZ, w.s.p. to Qiksik: https://qiskit.org/documentation/stubs/qiskit.circuit.library.RZZGate.html#qiskit.circuit.library.RZZGate
        /**  CU   = [e^(-i theta/2)    0               0               0]
         *          [0                 e^(i theta/2)  0               0]
         *          [0                 0               e^(i theta/2)  0]
         *          [0                 0               0               e^(-i theta/2)]
         */
        Gate *G = new Gate(OP::RZZ, qubit0, qubit1, 2, theta);

        ValType gm_real[16] = {cos(-HALF * theta), 0, 0, 0,
                               0, cos(HALF * theta), 0, 0,
                               0, 0, cos(HALF * theta), 0,
                               0, 0, 0, cos(-HALF * theta)};

        ValType gm_imag[16] = {sin(-HALF * theta), 0, 0, 0,
                               0, sin(HALF * theta), 0, 0,
                               0, 0, sin(HALF * theta), 0,
                               0, 0, 0, sin(-HALF * theta)};

        G->set_gm(gm_real, gm_imag, 4);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(std::max(qubit0, qubit1), Circuit::n_qubits);
    }

    void Circuit::SX(IdxType qubit)
    {
        // sqrt(X) gate, basis gate for IBMQ
        /** SX = 1/2 [1+i 1-i]
                     [1-i 1+i]
        */
        Gate *G = new Gate(OP::SX, qubit);
        ValType gm_real[4] = {HALF, HALF, HALF, HALF};
        ValType gm_imag[4] = {HALF, -HALF, -HALF, HALF};
        G->set_gm(gm_real, gm_imag, 2);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(qubit, Circuit::n_qubits);
    }

    void Circuit::ID(IdxType qubit)
    {
        // Identity gate
        /** ID  = [1 0]
                  [0 1]
        */
        Gate *G = new Gate(OP::ID, qubit);
        ValType gm_real[4] = {1, 0, 0, 1};
        ValType gm_imag[4] = {0, 0, 0, 0};
        G->set_gm(gm_real, gm_imag, 2);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(qubit, Circuit::n_qubits);
    }

    void Circuit::SWAP(IdxType ctrl, IdxType qubit)
    {
        // SWAP gate
        /**  SWAP = [1 0 0 0]
                    [0 0 1 0]
                    [0 1 0 0]
                    [0 0 0 1]
        */
        Gate *G = new Gate(OP::SWAP, qubit, ctrl, 2);
        ValType gm_real[16] = {1, 0, 0, 0,
                               0, 0, 1, 0,
                               0, 1, 0, 0,
                               0, 0, 0, 1};
        ValType gm_imag[16] = {0, 0, 0, 0,
                               0, 0, 0, 0,
                               0, 0, 0, 0,
                               0, 0, 0, 0};
        G->set_gm(gm_real, gm_imag, 4);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(std::max(ctrl, qubit), Circuit::n_qubits);
    }

    void Circuit::M(IdxType qubit) // default is pauli-Z
    {
        Gate *G = new Gate(OP::M, qubit);
        Circuit::gates->push_back(*G);
    }
    void Circuit::MA(IdxType repetition) // default is pauli-Z
    {
        Gate *G = new Gate(OP::MA, repetition);
        Circuit::gates->push_back(*G);
    }
    void Circuit::RESET(IdxType qubit)
    {
        Gate *G = new Gate(OP::RESET, qubit);
        Circuit::gates->push_back(*G);

        Circuit::n_qubits = std::max(qubit, Circuit::n_qubits);
    }

    // ============================== Other Gate Definition ================================
    void Circuit::U3(ValType theta, ValType phi, ValType lam, IdxType qubit)
    {
        Circuit::U(theta, phi, lam, qubit);
    }
    void Circuit::U2(ValType phi, ValType lam, IdxType qubit)
    {
        Circuit::U(PI / 2, phi, lam, qubit);
    }
    void Circuit::U1(ValType lam, IdxType qubit)
    {
        Circuit::U(0, 0, lam, qubit);
    }

    void Circuit::CCX(IdxType qubit0, IdxType qubit1, IdxType qubit2)
    {
        Circuit::H(qubit2);
        Circuit::CX(qubit1, qubit2);
        Circuit::TDG(qubit2);
        Circuit::CX(qubit0, qubit2);
        Circuit::T(qubit2);
        Circuit::CX(qubit1, qubit2);
        Circuit::T(qubit1);
        Circuit::TDG(qubit2);
        Circuit::CX(qubit0, qubit2);
        Circuit::T(qubit2);
        Circuit::CX(qubit0, qubit1);
        Circuit::T(qubit0);
        Circuit::TDG(qubit1);
        Circuit::H(qubit2);
        Circuit::CX(qubit0, qubit1);
    }
    void Circuit::CSWAP(IdxType qubit0, IdxType qubit1, IdxType qubit2)
    {
        Circuit::CX(qubit2, qubit1);
        Circuit::CCX(qubit0, qubit1, qubit2);
        Circuit::CX(qubit2, qubit1);
    }
}