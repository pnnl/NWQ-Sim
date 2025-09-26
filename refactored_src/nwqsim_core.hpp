#pragma once

#if defined(__CUDACC__) || defined(__HIPCC__)
#define NWQSIM_HD __host__ __device__
#else
#define NWQSIM_HD
#endif

namespace NWQSim
{
    using IdxType = long long;
    using ValType = double;

    /***********************************************
     * Constant configuration:
     ***********************************************/
/* Constant value of PI */
#define PI 3.14159265358979323846
/* Constant value of sqrt(2) */
#define S2I 0.70710678118654752440
/* Constant value of 0.5 */
#define HALF 0.5

    enum class SimType
    {
        SV,
        DM
    };

    struct ZObservableTerm
    {
        IdxType mask;
        ValType coeff;
    };

    inline constexpr ValType MEASUREMENT_ERROR_BAR = 1e-3;

    struct StateSpan
    {
        ValType *real = nullptr;
        ValType *imag = nullptr;
        ValType *work_real = nullptr;
        ValType *work_imag = nullptr;
        IdxType size = 0;
        IdxType half_size = 0;

        NWQSIM_HD IdxType full_extent() const
        {
            return size ? size : (half_size << 1);
        }

        NWQSIM_HD IdxType half_extent() const
        {
            return half_size ? half_size : (size >> 1);
        }
    };

    NWQSIM_HD inline int popcount64(IdxType value)
    {
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
        return __popcll(static_cast<unsigned long long>(value));
#else
        return __builtin_popcountll(static_cast<unsigned long long>(value));
#endif
    }

    NWQSIM_HD inline ValType z_mask_sign(IdxType index, IdxType mask)
    {
        return (popcount64(index & mask) & 1) ? ValType(-1.0) : ValType(1.0);
    }
}

#undef NWQSIM_HD
