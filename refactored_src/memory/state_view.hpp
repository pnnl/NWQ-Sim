#pragma once

#include "../nwqsim_core.hpp"

namespace NWQSim
{
    namespace GateKernels
    {
        struct StateView
        {
            ValType *data_real = nullptr;
            ValType *data_imag = nullptr;
            ValType *buffer_real = nullptr;
            ValType *buffer_imag = nullptr;
            IdxType dim = 0;
            IdxType half_dim = 0;
            IdxType gpu_scale = 0;
            IdxType lg2_m_gpu = 0;
            IdxType m_gpu = 0;
            IdxType rank = 0;
        };
    } // namespace GateKernels
} // namespace NWQSim
