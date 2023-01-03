// ---------------------------------------------------------------------------
// NWQsim: Northwest Quantum Circuit Simulation Environment
// ---------------------------------------------------------------------------
// Ang Li, Senior Computer Scientist
// Pacific Northwest National Laboratory(PNNL), U.S.
// Homepage: http://www.angliphd.com
// GitHub repo: http://www.github.com/pnnl/SV-Sim
// PNNL-IPID: 32166, ECCN: EAR99, IR: PNNL-SA-161181
// BSD Lincese.
// ---------------------------------------------------------------------------
// File: gate.h
// DMSim Overall gate definition.
// ---------------------------------------------------------------------------
#ifndef GATE_H_
#define GATE_H_

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

namespace NWQSim
{
//==========================================================
// NWQSim supported gate set based on Q#,OpenQASM2 and OpenQASM3. 
// QIR: https://github.com/microsoft/qsharp-runtime/blob/
// main/src/Qir/Runtime/public/QSharpSimApi_I.hpp
// OpenQASM2: https://arxiv.org/pdf/1707.03429.pdf 
// OpenQASM3: https://arxiv.org/pdf/2104.14722.pdf
// Also find standard gates defined in Qiskit:
// https://qiskit.org/documentation/apidoc/circuit_library.html
// and summary of quantum operations:
// https://qiskit.org/documentation/tutorials/circuits/3_summary_of_quantum_operations.html
//==========================================================
enum OP 
{
    /******************************************
     * Pauli-X gate: bit-flip or NOT gate
     *  X = [0 1]
     *      [1 0]  
     ******************************************/ 
     X, 
    /******************************************
     * Pauli-Y gate: bit-flip and phase flip gate
     * Y = [0 -i]
     *     [i  0]  
     ******************************************/
     Y, 
    /******************************************
     * Pauli-Z gate: phase flip gate
     * Z = [1  0]
     *     [0 -1]  
     ******************************************/
     Z,
    /******************************************
     * Clifford gate: Hadamard gate 
     * H = 1/sqrt(2) * [1  1]
                       [1 -1]
     ******************************************/
     H,
    /******************************************
     * Clifford gate: sqrt(Z) phase gate
     * S = [1 0]
           [0 i]   
     ******************************************/
     S,
    /******************************************
     * Clifford gate: inverse of sqrt(Z)
     * SDG = [1  0]
     *       [0 -i]
     ******************************************/
     SDG,
    /******************************************
     * sqrt(S) phase gate or T gate
     * T = [1 0]
           [0 s2i+s2i*i]
     ******************************************/
     T,
    /******************************************
     * Inverse of sqrt(S) gate
     * TDG = [1 0]
     *       [0 s2i-s2i*i]
     ******************************************/
     TDG,
     /******************************************
     * Global phase gate, defined as ctrl@gphase(a)
     * in QASM3, which is U(0,0,a) or RI in Q#
     * RI = [1 0]     = [1 0]
     *      [0 e^(ia)]  [0 cos(a)+i*sin(a)]
     ******************************************/
     RI,
     /******************************************
     * Rotation around X axis
     * RX = [cos(a/2) -i*sin(a/2)]
     *      [-i*sin(a/2) cos(a/2)]
     ******************************************/
     RX,
     /******************************************
     * Rotation around Y axis
     * RY = [cos(a/2) -sin(a/2)]
     *      [sin(a/2)  cos(a/2)]
     ******************************************/
     RY,
     /******************************************
     * Rotation around Z axis
     * RZ = [cos(a/2)-i*sin(a/2)  0]
     *      [0  cos(a/2)+i*sin(a/2)]
     ******************************************/
     RZ,
    /******************************************
     * sqrt(X) gate, basis gate for IBM-Q
     * SX = 1/2 [1+i 1-i]
     *          [1-i 1+i]
     ******************************************/
     SX,
    /******************************************
     * Phase gate (not global phase gate, or phase shift gate
     * P = [1 0]
     *     [0 cos(a)+i*sin(a)]
     ******************************************/
     P,
    /******************************************
     * Unitary gate.
     * U(a,b,c) = [cos(a/2)        -e^(ic)*sin(a/2)]
     *          = [e^(ib)*sin(a/2) e^(i(b+c))*cos(a/2)]
     ******************************************/
     U,
     /******************************************
     * Controlled X gate (CNOT)
     * Apply X when the control qubit is 1
     ******************************************/
     CX,
     /******************************************
     * Controlled Y gate
     * Apply Y when the control qubit is 1
     ******************************************/
     CY, 
     /******************************************
     * Controlled Z gate
     * Apply Z when the control qubit is 1
     ******************************************/
     CZ,
     /******************************************
     * Controlled H gate
     * Apply H when the control qubit is 1
     ******************************************/
     CH,
     /******************************************
     * Controlled S gate
     * Apply S when the control qubit is 1
     ******************************************/
     CS,
     /******************************************
     * Controlled SDG gate
     * Apply SDG when the control qubit is 1
     ******************************************/
     CSDG,
     /******************************************
     * Controlled T gate
     * Apply T when the control qubit is 1
     ******************************************/
     CT,
     /******************************************
     * Controlled TDG gate
     * Apply TDG when the control qubit is 1
     ******************************************/
     CTDG,
     /******************************************
     * Controlled RI gate
     * Apply RI when the control qubit is 1 
     ******************************************/
     CRI,
     /******************************************
     * Controlled RX gate
     * Apply RX when the control qubit is 1
     ******************************************/
     CRX,
     /******************************************
     * Controlled RY gate
     * Apply RY when the control qubit is 1
     ******************************************/
     CRY,
     /******************************************
     * Controlled RZ gate
     * Apply RZ when the control qubit is 1
     ******************************************/
     CRZ,
    /******************************************
     * Controlled sqrt(X) gate
     * Apply SX when the control qubit is 1
     ******************************************/
     CSX,
    /******************************************
     * Controlled phase gate
     * Apply P when the control qubit is 1
     ******************************************/
     CP,
     /******************************************
     * Controlled U gate
     * Apply U(a,b,c) when the control qubit is 1
     ******************************************/
     CU,
    /******************************************
     * Identiy gate, this is meaningful
     * for noisy simulation
     ******************************************/
     ID,
     /******************************************
     * SWAP gate: swap the position of two qubits
     * SWAP = [1,0,0,0]
     *        [0,0,1,0]
     *        [0,1,0,0]
     *        [0,0,0,1]
     ******************************************/
     SWAP,
     /******************************************
     * Measure gate: it measure a single qubit
     * and depends on the value, normalize the 
     * other coefficients. It returns a single bit.
     ******************************************/
     M,
     /******************************************
     * Measure all gate: it measure all
     * the qubits at once and return a bit string.
     ******************************************/
     MA,
     /******************************************
     * Reset gate: it resets a qubit to |0> and
     * leave other qubits unchaged.
     ******************************************/
     RESET,
     /************
      * C1 Gate
      ***********/
     C1,
     /************
      * C2 Gate
      ***********/
     C2,
     /************
      * C4 Gate
      ***********/
     C4
};
//Name of the gate for tracing purpose
const char *OP_NAMES[] = {
    //Basic
    "X", 
    "Y", 
    "Z", 
    "H", 
    "S", 
    "SDG", 
    "T", 
    "TDG", 
    "RI", 
    "RX", 
    "RY", 
    "RZ", 
    "SX", 
    "P",
    "U", 
    //Controlled
    "CX", 
    "CY", 
    "CZ", 
    "CH", 
    "CS", 
    "CSDG", 
    "CT", 
    "CTDG", 
    "CRI", 
    "CRX", 
    "CRY", 
    "CRZ", 
    "CSX", 
    "CP",
    "CU", 
    //Other
    "ID",
    "SWAP", 
    "M", 
    "MA", 
    "RESET",
    "C1",
    "C2",
    "C4"
};
/*#ifdef USE_NVGPU*/
//Name of the gate for tracing on GPU side
const __device__ char *OP_NAMES_NVGPU[] = {
    //Basic
    "X", 
    "Y", 
    "Z", 
    "H", 
    "S", 
    "SDG", 
    "T", 
    "TDG", 
    "RI", 
    "RX", 
    "RY", 
    "RZ", 
    "SX", 
    "P",
    "U", 
    //Controlled
    "CX", 
    "CY", 
    "CZ", 
    "CH", 
    "CS", 
    "CSDG", 
    "CT", 
    "CTDG", 
    "CRI", 
    "CRX", 
    "CRY", 
    "CRZ", 
    "CSX",
    "CP",
    "CU", 
    //Other
    "ID",
    "SWAP", 
    "M", 
    "MA", 
    "RESET",
    "C1",
    "C2",
    "C4"
};

class Simulation;

/***********************************************
 * Gate Definition
 ***********************************************/
class Gate
{
public:
    Gate(enum OP _op_name, IdxType _qubit, IdxType _ctrl=-1, ValType _theta=0) : 
        op_name(_op_name), qubit(_qubit), ctrl(_ctrl), theta(_theta)
    {
        memset(gm_real, 0, sizeof(ValType)*16);
        memset(gm_imag, 0, sizeof(ValType)*16);
    }
    ~Gate() {}
    Gate(const Gate& g):op_name(g.op_name), qubit(g.qubit), ctrl(g.ctrl), theta(g.theta)
    {
        memcpy(gm_real, g.gm_real, 16*sizeof(ValType));
        memcpy(gm_imag, g.gm_imag, 16*sizeof(ValType));
    }
    //set gm
    void set_gm(ValType* real, ValType* imag, IdxType dim)
    {
        if (!(dim==2 || dim==4))
            throw std::logic_error("Dim should be 2 (1-qubit gate) or 4 (2-qubit gate)!");
        memcpy(gm_real, real, dim*dim*sizeof(ValType));
        memcpy(gm_imag, imag, dim*dim*sizeof(ValType));
    }
    //applying the embedded gate operation on GPU side
    __device__ void exe_op(Simulation* sim, ValType* sv_real, ValType* sv_imag);
    //for dumping the gate
    string gateToString()
    {
        std::stringstream ss;
        ss << OP_NAMES[op_name] << "(qubit:" << qubit << ", ctrl:" << ctrl << ", theta:" 
           << theta << ");" << std::endl;
        return ss.str();
    }
    //Gate name
    enum OP op_name;
    IdxType qubit;
    IdxType ctrl;
    ValType theta;
    //4-qubit gate parameters (after fusion)
    ValType gm_real[16];
    ValType gm_imag[16];
}; //end of Gate definition



}//end of namespace DMSim
#endif
