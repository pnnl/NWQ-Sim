// ---------------------------------------------------------------------------
// NWQsim: Northwest Quantum Circuit Simulation Environment
// ---------------------------------------------------------------------------
// Ang Li, Senior Computer Scientist
// Pacific Northwest National Laboratory(PNNL), U.S.
// Homepage: http://www.angliphd.com
// GitHub repo: http://www.github.com/pnnl/NWQ-Sim
// PNNL-IPID: 32166, ECCN: EAR99, IR: PNNL-SA-161181
// BSD Lincese.
// ---------------------------------------------------------------------------
// File:

// ---------------------------------------------------------------------------

#ifndef NWQSIM_HPP
#define NWQSIM_HPP

#include "util.hpp"
#include "circuit.hpp"
#include "gate.hpp"

#include <cmath>
#include <string>
#define SVSIM

namespace NWQSim
{

    class Simulation
    {
    public:
        Simulation(IdxType _n_qubits);

        ~Simulation();
        void AllocateQubit();
        void ReleaseQubit();

        void reset_sim();
        void reset_circuit();

        IdxType get_n_qubits();
        IdxType get_n_gates();
        void set_seed(IdxType seed);
        void clear_circuit();
        void update(const IdxType _n_qubits, const IdxType _n_gates);
        std::string circuitToString();
        void sim();
        IdxType measure(IdxType qubit);
        IdxType *measure_all(IdxType repetition = DEFAULT_REPETITIONS);

        // cricuit
        Circuit *circuit_handle;

        // =============================== Standard Gates ===================================
        void X(IdxType qubit)
        {
            // Pauli X-gate: bit flip
            /** X = [0 1]
                    [1 0]
             */
            Gate *G = new Gate(OP::C1, qubit);
            ValType gm_real[4] = {0, 1, 1, 0};
            ValType gm_imag[4] = {0, 0, 0, 0};
            G->set_gm(gm_real, gm_imag, 2);
            circuit_handle->append(*G);
        }
        void SX(IdxType qubit)
        {
            // sqrt(X) gate, basis gate for IBMQ
            /** SX = 1/2 [1+i 1-i]
                         [1-i 1+i]
            */
            Gate *G = new Gate(OP::C1, qubit);
            ValType gm_real[4] = {HALF, HALF, HALF, HALF};
            ValType gm_imag[4] = {HALF, -HALF, -HALF, HALF};
            G->set_gm(gm_real, gm_imag, 2);
            circuit_handle->append(*G);
        }
        void RZ(ValType theta, IdxType qubit)
        {
            // Rotation around Z axis
            /** RZ = [cos(a/2)-i*sin(a/2)  0]
                     [0  cos(a/2)+i*sin(a/2)]
            */
            Gate *G = new Gate(OP::C1, qubit, -1, theta);
            ValType gm_real[4] = {cos(HALF * theta), 0, 0, cos(HALF * theta)};
            ValType gm_imag[4] = {-sin(HALF * theta), 0, 0, sin(HALF * theta)};
            G->set_gm(gm_real, gm_imag, 2);
            circuit_handle->append(*G);
        }
        void CX(IdxType ctrl, IdxType qubit)
        {
            // Controlled-NOT or CNOT
            /**  CX   = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 0 1]
                        [0 0 1 0]
            */
            Gate *G = new Gate(OP::C2, qubit, ctrl);
            ValType gm_real[16] = {1, 0, 0, 0,
                                   0, 1, 0, 0,
                                   0, 0, 0, 1,
                                   0, 0, 1, 0};
            ValType gm_imag[16] = {0};
            G->set_gm(gm_real, gm_imag, 4);
            circuit_handle->append(*G);
        }
        void ID(IdxType qubit)
        {
            // Identity gate
            /** ID  = [1 0]
                      [0 1]
            */
            Gate *G = new Gate(OP::C1, qubit);
            ValType gm_real[4] = {1, 0, 0, 1};
            ValType gm_imag[4] = {0, 0, 0, 0};
            G->set_gm(gm_real, gm_imag, 2);
            circuit_handle->append(*G);
        }

        void M(IdxType qubit) // default is pauli-Z
        {
            Gate *G = new Gate(OP::M, qubit, -1, 0);
            circuit_handle->append(*G);
        }
        void MA(IdxType repetition) // default is pauli-Z
        {
            Gate *G = new Gate(OP::MA, 0, repetition, 0);
            circuit_handle->append(*G);
        }
        void RESET(IdxType qubit)
        {
            Gate *G = new Gate(OP::RESET, qubit);
            circuit_handle->append(*G);
        }

#ifdef SVSIM
        // ===================== Additional Standard Gates for SVSIM =========================
        void Y(IdxType qubit)
        {
            // Pauli-Y gate: bit and phase flip
            /** Y = [0 -i]
                    [i  0]
             */
            Gate *G = new Gate(OP::C1, qubit);
            ValType gm_real[4] = {0, 0, 0, 0};
            ValType gm_imag[4] = {0, -1, 1, 0};
            G->set_gm(gm_real, gm_imag, 2);
            circuit_handle->append(*G);
        }
        void Z(IdxType qubit)
        {
            // Pauli-Z gate: phase flip
            /** Z = [1  0]
                    [0 -1]
             */
            Gate *G = new Gate(OP::C1, qubit);
            ValType gm_real[4] = {1, 0, 0, -1};
            ValType gm_imag[4] = {0, 0, 0, 0};
            G->set_gm(gm_real, gm_imag, 2);
            circuit_handle->append(*G);
        }
        void H(IdxType qubit)
        {
            // Clifford gate: Hadamard
            /** H = 1/sqrt(2) * [1  1]
                                [1 -1]
             */
            Gate *G = new Gate(OP::C1, qubit);
            ValType gm_real[4] = {S2I, S2I, S2I, -S2I};
            ValType gm_imag[4] = {0, 0, 0, 0};
            G->set_gm(gm_real, gm_imag, 2);
            circuit_handle->append(*G);
        }
        void S(IdxType qubit)
        {
            // Clifford gate: sqrt(Z) phase gate
            /** S = [1 0]
                    [0 i]
            */
            Gate *G = new Gate(OP::C1, qubit);
            ValType gm_real[4] = {1, 0, 0, 0};
            ValType gm_imag[4] = {0, 0, 0, 1};
            G->set_gm(gm_real, gm_imag, 2);
            circuit_handle->append(*G);
        }
        void SDG(IdxType qubit)
        {
            // Clifford gate: conjugate of sqrt(Z) phase gate
            /** SDG = [1  0]
                      [0 -i]
            */
            Gate *G = new Gate(OP::C1, qubit);
            ValType gm_real[4] = {1, 0, 0, 0};
            ValType gm_imag[4] = {0, 0, 0, -1};
            G->set_gm(gm_real, gm_imag, 2);
            circuit_handle->append(*G);
        }
        void T(IdxType qubit)
        {
            // C3 gate: sqrt(S) phase gate
            /** T = [1 0]
                    [0 s2i+s2i*i]
            */
            Gate *G = new Gate(OP::C1, qubit);
            ValType gm_real[4] = {1, 0, 0, S2I};
            ValType gm_imag[4] = {0, 0, 0, S2I};
            G->set_gm(gm_real, gm_imag, 2);
            circuit_handle->append(*G);
        }
        void TDG(IdxType qubit)
        {
            // C3 gate: conjugate of sqrt(S) phase gate
            /** TDG = [1 0]
                      [0 s2i-s2i*i]
            */
            Gate *G = new Gate(OP::C1, qubit);
            ValType gm_real[4] = {1, 0, 0, S2I};
            ValType gm_imag[4] = {0, 0, 0, -S2I};
            G->set_gm(gm_real, gm_imag, 2);
            circuit_handle->append(*G);
        }
        void RI(ValType theta, IdxType qubit)
        {
            // Global phase gate
            /** RI = [e^(ia) 0] = [cos(a)+i*sin(a) 0]
                     [0 e^(ia)]   [0 cos(a)+i*sin(a)]
            */
            Gate *G = new Gate(OP::C1, qubit, -1, theta);
            ValType gm_real[4] = {cos(theta), 0, 0, cos(theta)};
            ValType gm_imag[4] = {sin(theta), 0, 0, sin(theta)};
            G->set_gm(gm_real, gm_imag, 2);
            circuit_handle->append(*G);
        }
        void RX(ValType theta, IdxType qubit)
        {
            // Rotation around X axis
            /** RX = [cos(a/2) -i*sin(a/2)]
                     [-i*sin(a/2) cos(a/2)]
            */
            Gate *G = new Gate(OP::C1, qubit, -1, theta);
            ValType gm_real[4] = {cos(HALF * theta), 0, 0, cos(HALF * theta)};
            ValType gm_imag[4] = {0, -sin(HALF * theta), -sin(HALF * theta), 0};
            G->set_gm(gm_real, gm_imag, 2);
            circuit_handle->append(*G);
        }
        void RY(ValType theta, IdxType qubit)
        {
            // Rotation around Y axis
            /** RY = [cos(a/2) -sin(a/2)]
                     [sin(a/2)  cos(a/2)]
            */
            Gate *G = new Gate(OP::C1, qubit, -1, theta);
            ValType gm_real[4] = {cos(HALF * theta), -sin(HALF * theta), sin(HALF * theta), cos(HALF * theta)};
            ValType gm_imag[4] = {0};
            G->set_gm(gm_real, gm_imag, 2);
            circuit_handle->append(*G);
        }

        void P(ValType theta, IdxType qubit)
        {
            // Phase gate defined by Qiskit
            /** P = [1, 0     ]  = [1,0]
                    [0, e^(ia)]    [0,cos(a)+i*sin(a)]
            */
            Gate *G = new Gate(OP::C1, qubit, -1, theta);
            ValType gm_real[4] = {1, 0, 0, cos(theta)};
            ValType gm_imag[4] = {0, 0, 0, sin(theta)};
            G->set_gm(gm_real, gm_imag, 2);
            circuit_handle->append(*G);
        }
        void U(ValType theta, ValType phi, ValType lam, IdxType qubit)
        {
            // Generic single-qubit rotation gate with 3 Euler angles
            /** U = [cos(theta/2), -e^(i*lam)sin(theta/2)]
                    [e^(i*phi)sin(theta/2), e^(i*(phi+lam))cos(theta/2)]
            */
            Gate *G = new Gate(OP::C1, qubit, -1, theta);
            ValType gm_real[4] = {cos(HALF * theta),
                                  -cos(lam) * sin(HALF * theta),
                                  cos(phi) * sin(HALF * theta),
                                  cos(phi + lam) * cos(HALF * theta)};
            ValType gm_imag[4] = {0,
                                  -sin(lam) * sin(HALF * theta),
                                  sin(phi) * sin(HALF * theta),
                                  sin(lam + phi) * cos(HALF * theta)};
            G->set_gm(gm_real, gm_imag, 2);
            circuit_handle->append(*G);
        }
        void CY(IdxType ctrl, IdxType qubit)
        {
            // Controlled-Y
            /**  CY   = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 0 -i]
                        [0 0 i 0]
            */
            Gate *G = new Gate(OP::C2, qubit, ctrl);
            ValType gm_real[16] = {1, 0, 0, 0,
                                   0, 1, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, 0, 0};
            ValType gm_imag[16] = {0, 0, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, 0, -1,
                                   0, 0, 1, 0};
            G->set_gm(gm_real, gm_imag, 4);
            circuit_handle->append(*G);
        }
        void CZ(IdxType ctrl, IdxType qubit)
        {
            // Controlled-Z
            /**  CZ   = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 1 0]
                        [0 0 0 -1]
            */
            Gate *G = new Gate(OP::C2, qubit, ctrl);
            ValType gm_real[16] = {1, 0, 0, 0,
                                   0, 1, 0, 0,
                                   0, 0, 1, 0,
                                   0, 0, 0, -1};
            ValType gm_imag[16] = {0};
            G->set_gm(gm_real, gm_imag, 4);
            circuit_handle->append(*G);
        }
        void CH(IdxType ctrl, IdxType qubit)
        {
            // Controlled-H
            /**  CH   = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 s2i s2i]
                        [0 0 s2i -s2i]
            */
            Gate *G = new Gate(OP::C2, qubit, ctrl);
            ValType gm_real[16] = {1, 0, 0, 0,
                                   0, 1, 0, 0,
                                   0, 0, S2I, S2I,
                                   0, 0, S2I, -S2I};
            ValType gm_imag[16] = {0};
            G->set_gm(gm_real, gm_imag, 4);
            circuit_handle->append(*G);
        }
        void CS(IdxType ctrl, IdxType qubit)
        {
            // Controlled-S
            /**  CS   = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 1 0]
                        [0 0 0 i]
            */
            Gate *G = new Gate(OP::C2, qubit, ctrl);
            ValType gm_real[16] = {1, 0, 0, 0,
                                   0, 1, 0, 0,
                                   0, 0, 1, 0,
                                   0, 0, 0, 0};
            ValType gm_imag[16] = {0, 0, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, 0, 1};
            G->set_gm(gm_real, gm_imag, 4);
            circuit_handle->append(*G);
        }
        void CSDG(IdxType ctrl, IdxType qubit)
        {
            // Controlled-SDG
            /**  CSDG = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 1 0]
                        [0 0 0 -i]
            */
            Gate *G = new Gate(OP::C2, qubit, ctrl);
            ValType gm_real[16] = {1, 0, 0, 0,
                                   0, 1, 0, 0,
                                   0, 0, 1, 0,
                                   0, 0, 0, 0};
            ValType gm_imag[16] = {0, 0, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, 0, -1};
            G->set_gm(gm_real, gm_imag, 4);
            circuit_handle->append(*G);
        }
        void CT(IdxType ctrl, IdxType qubit)
        {
            // Controlled-T
            /**  CT   = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 1 0]
                        [0 0 0 s2i+si2*i]
            */
            Gate *G = new Gate(OP::C2, qubit, ctrl);
            ValType gm_real[16] = {1, 0, 0, 0,
                                   0, 1, 0, 0,
                                   0, 0, 1, 0,
                                   0, 0, 0, S2I};
            ValType gm_imag[16] = {0, 0, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, 0, S2I};
            G->set_gm(gm_real, gm_imag, 4);
            circuit_handle->append(*G);
        }
        void CTDG(IdxType ctrl, IdxType qubit)
        {
            // Controlled-TDG
            /**  CTDG = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 1 0]
                        [0 0 0 s2i-si2*i]
            */
            Gate *G = new Gate(OP::C2, qubit, ctrl);
            ValType gm_real[16] = {1, 0, 0, 0,
                                   0, 1, 0, 0,
                                   0, 0, 1, 0,
                                   0, 0, 0, S2I};
            ValType gm_imag[16] = {0, 0, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, 0, -S2I};
            G->set_gm(gm_real, gm_imag, 4);
            circuit_handle->append(*G);
        }
        void CRX(ValType theta, IdxType ctrl, IdxType qubit)
        {
            // Controlled-RX
            /**  CRX  = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 cos(a/2) -i*sin(a/2)]
                        [0 0 -i*sin(a/2) cos(a/2)]
            */
            Gate *G = new Gate(OP::C2, qubit, ctrl, theta);
            ValType gm_real[16] = {1, 0, 0, 0,
                                   0, 1, 0, 0,
                                   0, 0, cos(HALF * theta), 0,
                                   0, 0, 0, cos(HALF * theta)};
            ValType gm_imag[16] = {0, 0, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, 0, -sin(HALF * theta),
                                   0, 0, -sin(HALF * theta), 0};
            G->set_gm(gm_real, gm_imag, 4);
            circuit_handle->append(*G);
        }
        void CRY(ValType theta, IdxType ctrl, IdxType qubit)
        {
            // Controlled-RY
            /**  CRY  = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 cos(a/2) -sin(a/2)]
                        [0 0 sin(a/2)  cos(a/2)]
            */
            Gate *G = new Gate(OP::C2, qubit, ctrl, theta);
            ValType gm_real[16] = {1, 0, 0, 0,
                                   0, 1, 0, 0,
                                   0, 0, cos(HALF * theta), -sin(HALF * theta),
                                   0, 0, sin(HALF * theta), cos(HALF * theta)};
            ValType gm_imag[16] = {0};
            G->set_gm(gm_real, gm_imag, 4);
            circuit_handle->append(*G);
        }
        void CRZ(ValType theta, IdxType ctrl, IdxType qubit)
        {
            // Controlled-RZ
            /**  CRZ  = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 cos(a/2)-i*sin(a/2)  0]
                        [0 0 0  cos(a/2)+i*sin(a/2)]
            */
            Gate *G = new Gate(OP::C2, qubit, ctrl, theta);
            ValType gm_real[16] = {1, 0, 0, 0,
                                   0, 1, 0, 0,
                                   0, 0, cos(HALF * theta), 0,
                                   0, 0, 0, cos(HALF * theta)};
            ValType gm_imag[16] = {0, 0, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, -sin(HALF * theta), 0,
                                   0, 0, 0, sin(HALF * theta)};
            G->set_gm(gm_real, gm_imag, 4);
            circuit_handle->append(*G);
        }
        void CSX(IdxType ctrl, IdxType qubit)
        {
            // Controlled-SX
            /**  CSX  = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 (1+i)/2 (1-i)/2]
                        [0 0 (1-i)/2 (1+i)/2]
            */
            Gate *G = new Gate(OP::C2, qubit, ctrl);
            ValType gm_real[16] = {1, 0, 0, 0,
                                   0, 1, 0, 0,
                                   0, 0, HALF, HALF,
                                   0, 0, HALF, HALF};
            ValType gm_imag[16] = {0, 0, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, HALF, -HALF,
                                   0, 0, -HALF, HALF};
            G->set_gm(gm_real, gm_imag, 4);
            circuit_handle->append(*G);
        }
        void CP(ValType theta, IdxType ctrl, IdxType qubit)
        {
            // Controlled-P
            /**  CP   = [1 0 0 0]
                        [0 1 0 0]
                        [0 0 1 0]
                        [0 0 0 cos(a)+i*sin(a)]
            */
            Gate *G = new Gate(OP::C2, qubit, ctrl, theta);
            ValType gm_real[16] = {1, 0, 0, 0,
                                   0, 1, 0, 0,
                                   0, 0, 1, 0,
                                   0, 0, 0, cos(theta)};
            ValType gm_imag[16] = {0, 0, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, 0, sin(theta)};
            G->set_gm(gm_real, gm_imag, 4);
            circuit_handle->append(*G);
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
            Gate *G = new Gate(OP::C2, qubit, ctrl, theta);
            ValType gm_real[16] = {1, 0, 0, 0,
                                   0, 1, 0, 0,
                                   0, 0, cos(gamma) * cos(HALF * theta), -cos(gamma + lam) * sin(HALF * theta),
                                   0, 0, cos(gamma + phi) * sin(HALF * theta), cos(gamma + phi + lam) * cos(HALF * theta)};
            ValType gm_imag[16] = {0, 0, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, sin(gamma) * cos(HALF * theta), -sin(gamma + lam) * sin(HALF * theta),
                                   0, 0, sin(gamma + phi) * sin(HALF * theta), sin(gamma + phi + lam) * cos(HALF * theta)};
            G->set_gm(gm_real, gm_imag, 4);
            circuit_handle->append(*G);
        }
        void RXX(ValType theta, IdxType qubit0, IdxType qubit1)
        {
            // RXX, w.s.p. to Qiksik: https://qiskit.org/documentation/stubs/qiskit.circuit.library.RXXGate.html
            /**  CU   = [cos(theta/2)    0               0               -i*sin(theta/2)]
                        [0               cos(theta/2)    -i*sin(theta/2) 0]
                        [0               -i*sin(theta/2) cos(theta/2)    0]
                        [-i*sin(theta/2) 0               0               cos(theta/2)]
            */
            Gate *G = new Gate(OP::C2, qubit0, qubit1, theta);
            ValType gm_real[16] = {cos(HALF * theta), 0, 0, 0,
                                   0, cos(HALF * theta), 0, 0,
                                   0, 0, cos(HALF * theta), 0,
                                   0, 0, 0, cos(HALF * theta)};

            ValType gm_imag[16] = {0, 0, 0, -sin(HALF * theta),
                                   0, 0, -sin(HALF * theta), 0,
                                   0, -sin(HALF * theta), 0, 0,
                                   -sin(HALF * theta), 0, 0, 0};

            G->set_gm(gm_real, gm_imag, 4);
            circuit_handle->append(*G);
        }
        void RYY(ValType theta, IdxType qubit0, IdxType qubit1)
        {
            // RYY, w.s.p. to Qiksik: https://qiskit.org/documentation/stubs/qiskit.circuit.library.RYYGate.html#qiskit.circuit.library.RYYGate
            /**  CU   = [cos(theta/2)    0               0               i*sin(theta/2)]
                        [0               cos(theta/2)    -i*sin(theta/2) 0]
                        [0               -i*sin(theta/2) cos(theta/2)    0]
                        [i*sin(theta/2) 0               0               cos(theta/2)]
            */
            Gate *G = new Gate(OP::C2, qubit0, qubit1, theta);
            ValType gm_real[16] = {cos(HALF * theta), 0, 0, 0,
                                   0, cos(HALF * theta), 0, 0,
                                   0, 0, cos(HALF * theta), 0,
                                   0, 0, 0, cos(HALF * theta)};

            ValType gm_imag[16] = {0, 0, 0, sin(HALF * theta),
                                   0, 0, -sin(HALF * theta), 0,
                                   0, -sin(HALF * theta), 0, 0,
                                   sin(HALF * theta), 0, 0, 0};

            G->set_gm(gm_real, gm_imag, 4);
            circuit_handle->append(*G);
        }
        void RZZ(ValType theta, IdxType qubit0, IdxType qubit1)
        {
            // RZZ, w.s.p. to Qiksik: https://qiskit.org/documentation/stubs/qiskit.circuit.library.RZZGate.html#qiskit.circuit.library.RZZGate
            /**  CU   = [e^(-i theta/2)    0               0               0]
             *          [0                 e^(i theta/2)  0               0]
             *          [0                 0               e^(i theta/2)  0]
             *          [0                 0               0               e^(-i theta/2)]
             */
            Gate *G = new Gate(OP::C2, qubit0, qubit1, theta);

            ValType gm_real[16] = {cos(-HALF * theta), 0, 0, 0,
                                   0, cos(HALF * theta), 0, 0,
                                   0, 0, cos(HALF * theta), 0,
                                   0, 0, 0, cos(-HALF * theta)};

            ValType gm_imag[16] = {sin(-HALF * theta), 0, 0, 0,
                                   0, sin(HALF * theta), 0, 0,
                                   0, 0, sin(HALF * theta), 0,
                                   0, 0, 0, sin(-HALF * theta)};

            G->set_gm(gm_real, gm_imag, 4);
            circuit_handle->append(*G);
        }
        void SWAP(IdxType ctrl, IdxType qubit)
        {
            // SWAP gate
            /**  SWAP = [1 0 0 0]
                        [0 0 1 0]
                        [0 1 0 0]
                        [0 0 0 1]
            */
            Gate *G = new Gate(OP::C2, ctrl, qubit);
            ValType gm_real[16] = {1, 0, 0, 0,
                                   0, 0, 1, 0,
                                   0, 1, 0, 0,
                                   0, 0, 0, 1};
            ValType gm_imag[16] = {0, 0, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, 0, 0,
                                   0, 0, 0, 0};
            G->set_gm(gm_real, gm_imag, 4);
            circuit_handle->append(*G);
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
        // =============================== End of Gate Define ===================================
#endif // SVSIM
    };

} // namespace NWQSim

#endif // _NWQSIM_HPP_