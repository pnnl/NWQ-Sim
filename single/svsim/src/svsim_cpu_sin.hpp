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
// File: svsim_cpu_sin.hpp
// Single thread CPU based state-vector simulation
// ---------------------------------------------------------------------------

#ifndef SVSIM_CPU_SIN_HPP_
#define SVSIM_CPU_SIN_HPP_
#include <assert.h>
#include <random>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>

#include "config.h"
#include "gate.h"

#ifndef DISABLE_GATE_FUSION
#include "fusion.h"
#endif

namespace NWQSim
{
    using namespace std;

    class Gate;
    class Simulation;
    void Purity_Check(const Simulation *sim, const IdxType t, ValType *sv_real, ValType *sv_imag);

    // Simulation runtime
    void simulation_kernel(Simulation *);

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

    class Simulation
    {
    public:
        Simulation(IdxType _n_qubits = N_QUBIT_SLOT) : n_qubits(_n_qubits),
                                                       dim((IdxType)1 << (n_qubits)),
                                                       half_dim((IdxType)1 << (n_qubits - 1)),
                                                       n_gates(0),
                                                       cpu_mem(0),
                                                       sv_real(NULL),
                                                       sv_imag(NULL),
                                                       m_real(NULL),
                                                       randoms(NULL),
                                                       results(NULL)
        {
            sv_size = dim * (IdxType)sizeof(ValType);
            // CPU side initialization
            SAFE_ALOC_HOST(sv_real, sv_size);
            SAFE_ALOC_HOST(sv_imag, sv_size);
            memset(sv_real, 0, sv_size);
            memset(sv_imag, 0, sv_size);
            // State-vector initial state [0..0] = 1
            sv_real[0] = 1.;
            cpu_mem += sv_size * 4;
            // Initialize Circuit
            circuit_handle = new Circuit(n_qubits);
            circuit_handle_cpu = NULL;
            SAFE_ALOC_HOST(m_real, sv_size + sizeof(ValType));
            memset(m_real, 0, sv_size + sizeof(ValType));
            rng.seed(time(0));
#ifdef PRINT_SIM_TRACE
            printf("SVSim_cpu is initialized!\n");
#endif
        }

        ~Simulation()
        {
            // Release circuit
            if (circuit_handle != NULL)
                delete circuit_handle;
            // Release for CPU side
            SAFE_FREE_HOST(sv_real);
            SAFE_FREE_HOST(sv_imag);
            SAFE_FREE_HOST(randoms);
            SAFE_FREE_HOST(results);
            SAFE_FREE_HOST(m_real);
#ifdef PRINT_SIM_TRACE
            printf("SVSim_cpu is finalized!\n\n");
#endif
        }
        void AllocateQubit()
        {
            circuit_handle->AllocateQubit();
        }
        void ReleaseQubit()
        {
            circuit_handle->ReleaseQubit();
        }

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
        void M(IdxType qubit) // default is pauli-Z
        {
            SAFE_FREE_HOST(results);
            SAFE_ALOC_HOST(results, sizeof(IdxType));
            memset(results, 0, sizeof(IdxType));
            ValType rand = uni_dist(rng);
            Gate *G = new Gate(OP::M, qubit, -1, rand);
            circuit_handle->append(*G);
        }
        void MA(IdxType repetition) // default is pauli-Z
        {
            SAFE_FREE_HOST(results);
            SAFE_ALOC_HOST(results, sizeof(IdxType) * repetition);
            memset(results, 0, sizeof(IdxType) * repetition);
            SAFE_FREE_HOST(randoms);
            SAFE_ALOC_HOST(randoms, sizeof(ValType) * repetition);
            for (IdxType i = 0; i < repetition; i++)
                randoms[i] = uni_dist(rng);
            Gate *G = new Gate(OP::MA, 0, repetition, 0);
            circuit_handle->append(*G);
        }
        void RESET(IdxType qubit)
        {
            Gate *G = new Gate(OP::RESET, qubit);
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
        void RXX(ValType theta, IdxType qubit0, IdxType qubit1)
        {
            H(qubit0);
            H(qubit1);
            CX(qubit0, qubit1);
            RZ(theta, qubit1);
            CX(qubit0, qubit1);
            H(qubit0);
            H(qubit1);
        }
        void RYY(ValType theta, IdxType qubit0, IdxType qubit1)
        {
            RX(PI / 2, qubit0);
            RX(PI / 2, qubit1);
            CX(qubit0, qubit1);
            RZ(theta, qubit1);
            CX(qubit0, qubit1);
            RX(PI / 2, qubit0);
            RX(PI / 2, qubit1);
        }
        void RZZ(ValType theta, IdxType qubit0, IdxType qubit1)
        {
            CX(qubit0, qubit1);
            RZ(theta, qubit1);
            CX(qubit0, qubit1);
        }
        // =============================== End of Gate Define ===================================
        void reset_sim()
        {
            // Reset CPU input & output
            memset(sv_real, 0, sv_size);
            memset(sv_imag, 0, sv_size);
            memset(m_real, 0, sv_size + sizeof(ValType));
            // State Vector initial state [0..0] = 1
            sv_real[0] = 1.;
            reset_circuit();
        }
        void reset_circuit()
        {
            circuit_handle->reset();
#ifdef PRINT_CIRCUIT_TRACE
            printf("Circuit is reset!\n");
#endif
        }
        IdxType get_n_qubits()
        {
            return circuit_handle->n_qubits;
        }
        IdxType get_n_gates()
        {
            return circuit_handle->n_gates;
        }
        void set_seed(IdxType seed)
        {
            rng.seed(seed);
        }
        void clear_circuit()
        {
            circuit_handle->clear();
        }
        void update(const IdxType _n_qubits, const IdxType _n_gates)
        {
            this->n_gates = _n_gates;
        }
        std::string circuitToString()
        {
            return circuit_handle->circuitToString();
        }
        void sim()
        {
            IdxType input_gates = circuit_handle->n_gates;
            circuit_handle_cpu = circuit_handle->upload();
            // update should be put after upload where gate fusion is applied
            // which may change the number of gates in the circuit
            update(circuit_handle->n_qubits, circuit_handle->n_gates);
#ifdef PRINT_SIM_TRACE
            printf("SVSim_cpu is running! Requesting %lld qubits.\n", circuit_handle->n_qubits);
#endif
#ifdef PRINT_KERNEL_TRACE
            double sim_time;
            cpu_timer sim_timer;
            sim_timer.start_timer();
#endif
            //=========================================
            simulation_kernel(this);
            //=========================================
#ifdef PRINT_KERNEL_TRACE
            sim_timer.stop_timer();
            sim_time = sim_timer.measure();
            printf("\n============== SV-Sim ===============\n");
#ifdef USE_AVX512
            printf("nqubits:%lld, ngates:%lld, sim_gates:%lld, ncpus(avx512):%d, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms, mem:%.3lf MB, mem_per_cpu:%.3lf MB\n",
                   n_qubits, input_gates, n_gates, 1, sim_time, 0.,
                   sim_time, cpu_mem / 1024 / 1024, cpu_mem / 1024 / 1024);
#else
            printf("nqubits:%lld, ngates:%lld, sim_gates:%lld, ncpus:%d, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms, mem:%.3lf MB, mem_per_cpu:%.3lf MB\n",
                   n_qubits, input_gates, n_gates, 1, sim_time, 0.,
                   sim_time, cpu_mem / 1024 / 1024, cpu_mem / 1024 / 1024);
#endif
            printf("=====================================\n");
            fflush(stdout);
#endif
            clear_circuit();
        }
        IdxType measure(IdxType qubit)
        {
            this->M(qubit);
            this->sim();
            return this->results[0];
        }
        IdxType *measure_all(IdxType repetition = DEFAULT_REPETITIONS)
        {
            this->MA(repetition);
            this->sim();
            return this->results;
        }
        void print_res_sv()
        {
            IdxType num = ((IdxType)1 << n_qubits);
            printf("----- SVSim ------\n");
            for (IdxType i = 0; i < num; i++)
            {
                printf("(%.3lf,%.3lfj) ", sv_real[i], sv_imag[i]);
                if ((i + 1) % 8 == 0)
                    printf("\n");
            }
            printf("\n");
        }

    public:
        // n_qubits is the number of qubits
        IdxType n_qubits;
        IdxType sv_size;
        IdxType dim;
        IdxType half_dim;
        IdxType n_gates;
        // CPU arrays
        ValType *sv_real;
        ValType *sv_imag;
        ValType *m_real;
        // For measurement randoms
        ValType *randoms;
        // For measurement result
        IdxType *results;
        // Random
        std::mt19937 rng;
        std::uniform_real_distribution<ValType> uni_dist;
        // CPU memory usage
        ValType cpu_mem;
        // cricuit
        Circuit *circuit_handle;
        Gate *circuit_handle_cpu;
    };

    void CHECK_TRACE(const Simulation *sim, ValType *sv_real, ValType *sv_imag, IdxType t);

    /***********************************************
     * Simulation Main Kernel
     ***********************************************/
    void simulation_kernel(Simulation *sim)
    {
        for (IdxType t = 0; t < (sim->n_gates); t++)
        {
            // CHECK_TRACE(sim, sim->sv_real, sim->sv_imag, t);
            // printf("==%lld== %s(qubit:%lld, ctrl:%lld, theta:%lf)\n",sim->circuit_handle_cpu[t].op_name, OP_NAMES[sim->circuit_handle_cpu[t].op_name], sim->circuit_handle_cpu[t].qubit, sim->circuit_handle_cpu[t].ctrl, sim->circuit_handle_cpu[t].theta);

            // IdxType t0;
            // if (blockIdx.x == 0 && threadIdx.x == 0) t0 = clock64();

            ((sim->circuit_handle_cpu)[t]).exe_op(sim, sim->sv_real, sim->sv_imag);

            /*
                IdxType time = clock64() - t0;
                printf("G%lld:%s(ctrl:%lld,qubit:%lld) ticks: %lld\n",t,
                        OP_NAMES[sim->circuit_handle_cpu[t].op_name],
                        sim->circuit_handle_cpu[t].ctrl,
                        sim->circuit_handle_cpu[t].qubit,
                        time);
             */
            // CHECK_TRACE(sim, sim->sv_real, sim->sv_imag, t);
#ifdef PURITY_CHECK
            Purity_Check(sim, t, sim->sv_real, sim->sv_imag);
#endif
        }
    }

    //================================= Gate Definition ========================================

#ifdef USE_AVX512
#include "svsim_cpu_avx512.hpp"
#else

/***********************************************
 * Key Macros
 ***********************************************/
#define PUT(arr, i, val) (arr[(i)] = (val))
#define GET(arr, i) (arr[(i)])
#define BARR  \
    while (0) \
    {         \
    };

    //============== Check Trace (debug purpose) ================
    void CHECK_TRACE(const Simulation *sim, ValType *sv_real, ValType *sv_imag, IdxType t)
    {
        ValType trace = 0;
        for (IdxType i = 0; i < ((IdxType)1 << (sim->n_qubits)); i++)
        {
            const ValType val = GET(sv_real, i);
            trace += fabs(val);
        }
        printf("%s: Trace is: %lf\n", OP_NAMES[sim->circuit_handle_cpu[t].op_name], trace);
    }

// For C2 and C4 gates
#define DIV2E(x, y) ((x) >> (y))
#define MOD2E(x, y) ((x) & (((IdxType)1 << (y)) - (IdxType)1))
#define EXP2E(x) ((IdxType)1 << (x))
#define SV16IDX(x) (((x >> 3) & 1) * EXP2E(qubit0) + ((x >> 2) & 1) * EXP2E(qubit1) + ((x >> 1) & 1) * EXP2E(qubit2) + ((x & 1) * EXP2E(qubit3)))

    //============== C1 Gate ================
    // Arbitrary 1-qubit gate
    void C1_GATE(const Simulation *sim, ValType *sv_real, ValType *sv_imag,
                 const ValType *gm_real, const ValType *gm_imag, const IdxType qubit)
    {
        for (IdxType i = 0; i < sim->half_dim; i++)
        {
            IdxType outer = (i >> qubit);
            IdxType inner = (i & (((IdxType)1 << qubit) - 1));
            IdxType offset = (outer << (qubit + 1));
            IdxType pos0 = (offset + inner);
            IdxType pos1 = (offset + inner + ((IdxType)1 << qubit));

            const ValType el0_real = sv_real[pos0];
            const ValType el0_imag = sv_imag[pos0];
            const ValType el1_real = sv_real[pos1];
            const ValType el1_imag = sv_imag[pos1];

            ValType sv_real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag);
            ValType sv_imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real);
            ValType sv_real_pos1 = (gm_real[2] * el0_real) - (gm_imag[2] * el0_imag) + (gm_real[3] * el1_real) - (gm_imag[3] * el1_imag);
            ValType sv_imag_pos1 = (gm_real[2] * el0_imag) + (gm_imag[2] * el0_real) + (gm_real[3] * el1_imag) + (gm_imag[3] * el1_real);

            sv_real[pos0] = sv_real_pos0;
            sv_imag[pos0] = sv_imag_pos0;
            sv_real[pos1] = sv_real_pos1;
            sv_imag[pos1] = sv_imag_pos1;
        }
    }

    //============== C2 Gate ================
    // Arbitrary 2-qubit gate
    inline void C2_GATE(const Simulation *sim, ValType *sv_real, ValType *sv_imag,
                        const ValType *gm_real, const ValType *gm_imag,
                        const IdxType qubit0, const IdxType qubit1)
    {

        const IdxType per_pe_work = ((sim->dim) >> 2);
        assert(qubit0 != qubit1); // Non-cloning

        const IdxType q0dim = ((IdxType)1 << max(qubit0, qubit1));
        const IdxType q1dim = ((IdxType)1 << min(qubit0, qubit1));
        const IdxType outer_factor = ((sim->dim) + q0dim + q0dim - 1) >> (max(qubit0, qubit1) + 1);
        const IdxType mider_factor = (q0dim + q1dim + q1dim - 1) >> (min(qubit0, qubit1) + 1);
        const IdxType inner_factor = q1dim;
        const IdxType qubit0_dim = ((IdxType)1 << qubit0);
        const IdxType qubit1_dim = ((IdxType)1 << qubit1);

        for (IdxType i = 0; i < per_pe_work; i++)
        {
            IdxType outer = ((i / inner_factor) / (mider_factor)) * (q0dim + q0dim);
            IdxType mider = ((i / inner_factor) % (mider_factor)) * (q1dim + q1dim);
            IdxType inner = i % inner_factor;
            IdxType pos0 = outer + mider + inner;
            IdxType pos1 = outer + mider + inner + qubit1_dim;
            IdxType pos2 = outer + mider + inner + qubit0_dim;
            IdxType pos3 = outer + mider + inner + q0dim + q1dim;

            const ValType el0_real = GET(sv_real, pos0);
            const ValType el0_imag = GET(sv_imag, pos0);
            const ValType el1_real = GET(sv_real, pos1);
            const ValType el1_imag = GET(sv_imag, pos1);
            const ValType el2_real = GET(sv_real, pos2);
            const ValType el2_imag = GET(sv_imag, pos2);
            const ValType el3_real = GET(sv_real, pos3);
            const ValType el3_imag = GET(sv_imag, pos3);
            // Real part
            ValType sv_real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag) + (gm_real[2] * el2_real) - (gm_imag[2] * el2_imag) + (gm_real[3] * el3_real) - (gm_imag[3] * el3_imag);
            ValType sv_real_pos1 = (gm_real[4] * el0_real) - (gm_imag[4] * el0_imag) + (gm_real[5] * el1_real) - (gm_imag[5] * el1_imag) + (gm_real[6] * el2_real) - (gm_imag[6] * el2_imag) + (gm_real[7] * el3_real) - (gm_imag[7] * el3_imag);
            ValType sv_real_pos2 = (gm_real[8] * el0_real) - (gm_imag[8] * el0_imag) + (gm_real[9] * el1_real) - (gm_imag[9] * el1_imag) + (gm_real[10] * el2_real) - (gm_imag[10] * el2_imag) + (gm_real[11] * el3_real) - (gm_imag[11] * el3_imag);
            ValType sv_real_pos3 = (gm_real[12] * el0_real) - (gm_imag[12] * el0_imag) + (gm_real[13] * el1_real) - (gm_imag[13] * el1_imag) + (gm_real[14] * el2_real) - (gm_imag[14] * el2_imag) + (gm_real[15] * el3_real) - (gm_imag[15] * el3_imag);
            // Imag part
            ValType sv_imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real) + (gm_real[2] * el2_imag) + (gm_imag[2] * el2_real) + (gm_real[3] * el3_imag) + (gm_imag[3] * el3_real);
            ValType sv_imag_pos1 = (gm_real[4] * el0_imag) + (gm_imag[4] * el0_real) + (gm_real[5] * el1_imag) + (gm_imag[5] * el1_real) + (gm_real[6] * el2_imag) + (gm_imag[6] * el2_real) + (gm_real[7] * el3_imag) + (gm_imag[7] * el3_real);
            ValType sv_imag_pos2 = (gm_real[8] * el0_imag) + (gm_imag[8] * el0_real) + (gm_real[9] * el1_imag) + (gm_imag[9] * el1_real) + (gm_real[10] * el2_imag) + (gm_imag[10] * el2_real) + (gm_real[11] * el3_imag) + (gm_imag[11] * el3_real);
            ValType sv_imag_pos3 = (gm_real[12] * el0_imag) + (gm_imag[12] * el0_real) + (gm_real[13] * el1_imag) + (gm_imag[13] * el1_real) + (gm_real[14] * el2_imag) + (gm_imag[14] * el2_real) + (gm_real[15] * el3_imag) + (gm_imag[15] * el3_real);

            PUT(sv_real, pos0, sv_real_pos0);
            PUT(sv_real, pos1, sv_real_pos1);
            PUT(sv_real, pos2, sv_real_pos2);
            PUT(sv_real, pos3, sv_real_pos3);

            PUT(sv_imag, pos0, sv_imag_pos0);
            PUT(sv_imag, pos1, sv_imag_pos1);
            PUT(sv_imag, pos2, sv_imag_pos2);
            PUT(sv_imag, pos3, sv_imag_pos3);
        }
    }

    //============== C4 Gate ================
    // Arbitrary 4-qubit gate
    inline void C4_GATE(const Simulation *sim, ValType *sv_real, ValType *sv_imag,
                        const ValType *gm_real, const ValType *gm_imag,
                        const IdxType qubit0, const IdxType qubit1,
                        const IdxType qubit2, const IdxType qubit3)
    {
        assert(qubit0 != qubit1); // Non-cloning
        assert(qubit0 != qubit2); // Non-cloning
        assert(qubit0 != qubit3); // Non-cloning
        assert(qubit1 != qubit2); // Non-cloning
        assert(qubit1 != qubit3); // Non-cloning
        assert(qubit2 != qubit3); // Non-cloning
        // need to sort qubits: min->max: p, q, r, s
        const IdxType v0 = min(qubit0, qubit1);
        const IdxType v1 = min(qubit2, qubit3);
        const IdxType v2 = max(qubit0, qubit1);
        const IdxType v3 = max(qubit2, qubit3);
        const IdxType p = min(v0, v1);
        const IdxType q = min(min(v2, v3), max(v0, v1));
        const IdxType r = max(min(v2, v3), max(v0, v1));
        const IdxType s = max(v2, v3);
        for (IdxType i = 0; i < ((sim->dim) >> 4); i++)
        {
            const IdxType term0 = MOD2E(i, p);
            const IdxType term1 = MOD2E(DIV2E(i, p), q - p - 1) * EXP2E(p + 1);
            const IdxType term2 = MOD2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1) * EXP2E(q + 1);
            const IdxType term3 = MOD2E(DIV2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(r + 1);
            const IdxType term4 = DIV2E(DIV2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(s + 1);
            const IdxType term = term4 + term3 + term2 + term1 + term0;
            const ValType el_real[16] = {
                GET(sv_real, term + SV16IDX(0)), GET(sv_real, term + SV16IDX(1)),
                GET(sv_real, term + SV16IDX(2)), GET(sv_real, term + SV16IDX(3)),
                GET(sv_real, term + SV16IDX(4)), GET(sv_real, term + SV16IDX(5)),
                GET(sv_real, term + SV16IDX(6)), GET(sv_real, term + SV16IDX(7)),
                GET(sv_real, term + SV16IDX(8)), GET(sv_real, term + SV16IDX(9)),
                GET(sv_real, term + SV16IDX(10)), GET(sv_real, term + SV16IDX(11)),
                GET(sv_real, term + SV16IDX(12)), GET(sv_real, term + SV16IDX(13)),
                GET(sv_real, term + SV16IDX(14)), GET(sv_real, term + SV16IDX(15))};
            const ValType el_imag[16] = {
                GET(sv_imag, term + SV16IDX(0)), GET(sv_imag, term + SV16IDX(1)),
                GET(sv_imag, term + SV16IDX(2)), GET(sv_imag, term + SV16IDX(3)),
                GET(sv_imag, term + SV16IDX(4)), GET(sv_imag, term + SV16IDX(5)),
                GET(sv_imag, term + SV16IDX(6)), GET(sv_imag, term + SV16IDX(7)),
                GET(sv_imag, term + SV16IDX(8)), GET(sv_imag, term + SV16IDX(9)),
                GET(sv_imag, term + SV16IDX(10)), GET(sv_imag, term + SV16IDX(11)),
                GET(sv_imag, term + SV16IDX(12)), GET(sv_imag, term + SV16IDX(13)),
                GET(sv_imag, term + SV16IDX(14)), GET(sv_imag, term + SV16IDX(15))};
#pragma unroll
            for (unsigned j = 0; j < 16; j++)
            {
                ValType res_real = 0;
                ValType res_imag = 0;
#pragma unroll
                for (unsigned k = 0; k < 16; k++)
                {
                    res_real += (el_real[k] * gm_real[j * 16 + k]) - (el_imag[k] * gm_imag[j * 16 + k]);
                    res_imag += (el_real[k] * gm_imag[j * 16 + k]) + (el_imag[k] * gm_real[j * 16 + k]);
                }
                PUT(sv_real, term + SV16IDX(j), res_real);
                PUT(sv_imag, term + SV16IDX(j), res_imag);
            }
        }
    }
#endif

    void M_GATE(const Simulation *sim, ValType *sv_real, ValType *sv_imag,
                ValType *gm_real, ValType *gm_imag, const IdxType qubit, const ValType rand)
    {
        IdxType mask = ((IdxType)1 << qubit);
        ValType prob_of_one = 0;
        for (IdxType i = 0; i < ((IdxType)1 << (sim->n_qubits)); i++)
        {
            if ((i & mask) != 0)
                prob_of_one += (sv_real[i] * sv_real[i]) + (sv_imag[i] * sv_imag[i]); // square
        }
        if (rand < prob_of_one)
        {
            ValType normalize_factor = sqrt(prob_of_one);
            for (IdxType i = 0; i < sim->dim; i++)
            {
                if ((i & mask) == 0)
                {
                    sv_real[i] = 0.;
                    sv_imag[i] = 0.;
                }
                else
                {
                    sv_real[i] /= normalize_factor;
                    sv_imag[i] /= normalize_factor;
                }
            }
        }
        else
        {
            for (IdxType i = 0; i < sim->dim; i++)
            {
                ValType normalize_factor = sqrt(1.0 - prob_of_one);
                if ((i & mask) == 0)
                {
                    sv_real[i] /= normalize_factor;
                    sv_imag[i] /= normalize_factor;
                }
                else
                {
                    sv_real[i] = 0;
                    sv_imag[i] = 0;
                }
            }
        }
        sim->results[0] = (rand < prob_of_one ? 1 : 0);
    }

    //============== MA Gate (Measure all qubits in Pauli-Z) ================
    inline void MA_GATE(const Simulation *sim, ValType *sv_real, ValType *sv_imag, const IdxType repetition)
    {
        ValType *m_real = sim->m_real;
        IdxType n_size = (IdxType)1 << (sim->n_qubits);
        m_real[0] = 0;
        for (IdxType i = 1; i < (((IdxType)1 << (sim->n_qubits)) + 1); i++)
            m_real[i] = m_real[i - 1] + ((sv_real[i - 1] * sv_real[i - 1]) + (sv_imag[i - 1] * sv_imag[i - 1]));
        ValType purity = fabs(m_real[((IdxType)1 << (sim->n_qubits))]);
        if (abs(purity - 1.0) > ERROR_BAR)
            printf("MA: Purity Check fails with %lf\n", purity);

        if (repetition < n_size)
        {
            for (IdxType j = 0; j < n_size; j++)
            {
                ValType lower = m_real[j];
                ValType upper = (j + 1 == n_size) ? 1 : m_real[j + 1];
                for (IdxType i = 0; i < repetition; i++)
                {
                    ValType r = sim->randoms[i];
                    if (lower <= r && r < upper)
                        sim->results[i] = j;
                }
            }
        }
        else
        {
            for (IdxType i = 0; i < repetition; i++)
            {
                IdxType lo = 0;
                IdxType hi = ((IdxType)1 << (sim->n_qubits));
                IdxType mid;
                ValType r = sim->randoms[i];
                while (hi - lo > 1)
                {
                    mid = lo + (hi - lo) / 2;
                    if (r >= m_real[mid])
                        lo = mid;
                    else
                        hi = mid;
                }
                sim->results[i] = lo;
            }
        }
    }

    //============== Reset ================
    inline void RESET_GATE(const Simulation *sim, ValType *sv_real, ValType *sv_imag, const IdxType qubit)
    {
        IdxType mask = ((IdxType)1 << qubit);
        ValType prob_of_one = 0;
        for (IdxType i = 0; i < ((IdxType)1 << (sim->n_qubits)); i++)
        {
            if ((i & mask) != 0)
                prob_of_one += (sv_real[i] * sv_real[i]) + (sv_imag[i] * sv_imag[i]); // square
        }
        assert(prob_of_one <= 1.0);

        if (prob_of_one < 1.0) // still possible to normalize
        {
            ValType normalize_factor = sqrt(1.0 - prob_of_one);
            for (IdxType i = 0; i < sim->dim; i++)
            {
                if ((i & mask) == 0)
                {
                    sv_real[i] /= normalize_factor;
                    sv_imag[i] /= normalize_factor;
                }
                else
                {
                    sv_real[i] = 0;
                    sv_imag[i] = 0;
                }
            }
        }
        else
        {
            for (IdxType i = 0; i < sim->dim; i++)
            {
                if ((i & mask) == 0)
                {
                    IdxType dual_i = i ^ mask;
                    sv_real[i] = sv_real[dual_i];
                    sv_imag[i] = sv_imag[dual_i];
                    sv_real[dual_i] = 0;
                    sv_imag[dual_i] = 0;
                }
            }
        }
    }

    //============== Purity Check  ================
    void Purity_Check(const Simulation *sim, const IdxType t, ValType *sv_real, ValType *sv_imag)
    {
        ValType purity = 0;
        for (IdxType i = 0; i < (((IdxType)1 << (sim->n_qubits))); i++)
            purity += ((sv_real[i] * sv_real[i]) + (sv_imag[i] * sv_imag[i]));
        if (fabs(purity - 1.0) > ERROR_BAR)
        {
            Gate *g = &sim->circuit_handle_cpu[t];
            printf("Purity Check fails after Gate-%lld=>%s(ctrl:%lld,qubit:%lld,theta:%lf) with %lf\n", t, OP_NAMES[g->op_name], g->ctrl, g->qubit, g->theta, purity);
        }
    }

    //=====================================================================================
    // Per-gate execution function
    void Gate::exe_op(Simulation *sim, ValType *sv_real, ValType *sv_imag)
    {
        if (op_name == RESET)
        {
            RESET_GATE(sim, sv_real, sv_imag, qubit);
        }
        else if (op_name == M)
        {
            M_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, qubit, theta);
        }
        else if (op_name == MA)
        {
            MA_GATE(sim, sv_real, sv_imag, ctrl);
        }
        else if (op_name == C1)
        {
            C1_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, qubit);
        }
        else if (op_name == C2)
        {
            C2_GATE(sim, sv_real, sv_imag, gm_real, gm_imag, ctrl, qubit);
        }
    }

}; // namespace NWQSim

#endif
