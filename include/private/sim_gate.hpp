#pragma once

#include "../public/util.hpp"
#include "../public/gate.hpp"

namespace NWQSim
{

    enum SIM_TYPE
    {
        SV,
        DM
    };

    class SimGate
    {
    public:
        OP op_name;
        IdxType qubit;
        IdxType ctrl;

        SimGate(OP _op, IdxType _qubit, IdxType _ctrl)
        {
            op_name = _op;
            qubit = _qubit;
            ctrl = _ctrl;
        }
    };

    class SVGate : public SimGate
    {
    public:
        SVGate(OP _op, IdxType _qubit0, IdxType _qubit1) : SimGate(_op, _qubit0, _qubit1)
        {
            memset(gm_real, 0, sizeof(ValType) * 16);
            memset(gm_imag, 0, sizeof(ValType) * 16);
        }
        // gate parameters
        ValType gm_real[16];
        ValType gm_imag[16];
    };

    class DMSimGate : public SimGate
    {
    public:
        DMSimGate(OP _op, IdxType _qubit0, IdxType _qubit1) : SimGate(_op, _qubit0, _qubit1) {}
        // gate parameters
        ValType gm_real[256];
        ValType gm_imag[256];
    };

}