#ifndef CIRCUIT
#define CIRCUIT

#include "util.hpp"
#include "gate.hpp"

#include <string>
#include <vector>
#include <sstream>
#include <memory>

namespace NWQSim
{
    class Circuit
    {

    private:
        // number of qubits
        IdxType n_qubits = 0;
     
    public:
       // user input gate sequence
        std::shared_ptr<std::vector<std::shared_ptr<Gate>>> gates;

        
        Circuit();
        ~Circuit();

        IdxType num_qubits();
        IdxType num_gates();

        void clear();
        void reset();
        std::string circuitToString();

        // ===================== Standard Gates =========================

        void X(IdxType qubit);
        void Y(IdxType qubit);
        void Z(IdxType qubit);
        void H(IdxType qubit);
        void S(IdxType qubit);
        void SDG(IdxType qubit);
        void T(IdxType qubit);
        void TDG(IdxType qubit);
        void RI(ValType theta, IdxType qubit);
        void RX(ValType theta, IdxType qubit);
        void RY(ValType theta, IdxType qubit);
        void RZ(ValType theta, IdxType qubit);
        void P(ValType theta, IdxType qubit);
        void U(ValType theta, ValType phi, ValType lam, IdxType qubit);
        void CX(IdxType ctrl, IdxType qubit);
        void CY(IdxType ctrl, IdxType qubit);
        void CZ(IdxType ctrl, IdxType qubit);
        void CH(IdxType ctrl, IdxType qubit);
        void CS(IdxType ctrl, IdxType qubit);
        void CSDG(IdxType ctrl, IdxType qubit);
        void CT(IdxType ctrl, IdxType qubit);
        void CTDG(IdxType ctrl, IdxType qubit);
        void CRX(ValType theta, IdxType ctrl, IdxType qubit);
        void CRY(ValType theta, IdxType ctrl, IdxType qubit);
        void CRZ(ValType theta, IdxType ctrl, IdxType qubit);
        void CSX(IdxType ctrl, IdxType qubit);
        void CP(ValType theta, IdxType ctrl, IdxType qubit);
        void CU(ValType theta, ValType phi, ValType lam, ValType gamma,
                IdxType ctrl, IdxType qubit);
        void RXX(ValType theta, IdxType qubit0, IdxType qubit1);
        void RYY(ValType theta, IdxType qubit0, IdxType qubit1);
        void RZZ(ValType theta, IdxType qubit0, IdxType qubit1);
        void SX(IdxType qubit);
        void ID(IdxType qubit);
        void SWAP(IdxType ctrl, IdxType qubit);
        void M(IdxType qubit);
        void MA(IdxType repetition);
        void RESET(IdxType qubit);

        // ============================== Other Gate Definition ================================
        void U3(ValType theta, ValType phi, ValType lam, IdxType qubit);
        void U2(ValType phi, ValType lam, IdxType qubit);
        void U1(ValType lam, IdxType qubit);
        void CCX(IdxType qubit0, IdxType qubit1, IdxType qubit2);
        void CSWAP(IdxType qubit0, IdxType qubit1, IdxType qubit2);
    };

} // namespace NWQSim

#endif // CIRCUIT