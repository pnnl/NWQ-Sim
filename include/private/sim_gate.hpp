#pragma once

#include "../public/util.hpp"
#include "../public/gate.hpp"
#include <complex>
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

    class DMGate : public SimGate
    {
    public:
        DMGate(OP _op, IdxType _qubit0, IdxType _qubit1) : SimGate(_op, _qubit0, _qubit1) {}

        void set_gm(std::complex<ValType> *gm, IdxType dim)
        {
            if (!(dim == 4 || dim == 16))
                throw std::logic_error("Dim should be 4 (1-qubit gate) or 16 (2-qubit gate)!");
            for (IdxType i = 0; i < dim; i++)
            {
                for (IdxType j = 0; j < dim; j++)
                {
                    gm_real[i * dim + j] = gm[i * dim + j].real();
                    gm_imag[i * dim + j] = gm[i * dim + j].imag();
                }
            }
        }

        // gate parameters
        ValType gm_real[256];
        ValType gm_imag[256];
    };
}