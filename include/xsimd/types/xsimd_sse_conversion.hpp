// ==========================================================================
// Copyright (C) 2016 by Wolf Vollprecht and Johan Mabille
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.
// ==========================================================================

#ifndef XSIMD_SSE_CONVERSION_HPP
#define XSIMD_SSE_CONVERSION_HPP

#include "xsimd_sse_double.hpp"
#include "xsimd_sse_float.hpp"
#include "xsimd_sse_int.hpp"

namespace xsimd
{
    /************************
     * conversion functions *
     ************************/

    batch<int32_t, 4> to_int(const batch<float, 4>& x);
    batch<int64_t, 2> to_int(const batch<double, 2>& x);

    batch<float, 4> to_float(const batch<int32_t, 4>& x);
    batch<double, 2> to_float(const batch<int64_t, 2>& x);

    /***********************
     * boolean conversions *
     ***********************/

    batch<float, 4> bitwise_cast(const batch<int32_t, 4>& x);
    batch<int32_t, 4> bitwise_cast(const batch<float, 4>& x);
    batch<double, 2> bitwise_cast(const batch<int64_t, 2>& x);
    batch<int64_t, 2> bitwise_cast(const batch<double, 2>& x);

    /**********************
     * conversion details *
     **********************/

    inline batch<int32_t, 4> to_int(const batch<float, 4>& x)
    {
        return _mm_cvttps_epi32(x);
    }

    inline batch<int64_t, 2> to_int(const batch<double, 2>& x)
    {
#if defined(__x86_64__)
        return _mm_cvttpd_epi64(x);
#else
        return batch<int64_t, 2>(static_cast<int64_t>(x[0]), static_cast<int64_t>(x[1]));
#endif
    }

    inline batch<float, 4> to_float(const batch<int32_t, 4>& x)
    {
        return _mm_cvtepi32_ps(x);
    }

    inline batch<double, 2> to_float(const batch<int64_t, 2>& x)
    {
#if defined(__x86_64__)
        return _mm_cvtepi64_pd(x);
#else
        return batch<double, 2>(static_cast<double>(x[0]), static_cast<double>(x[1]));
#endif
    }

    /*****************************
     * bitwise_cast details      *
     *****************************/

    inline batch<float, 4> bitwise_cast(const batch<int32_t, 4>& x)
    {
        return _mm_castsi128_ps(x);
    }

    inline batch<int32_t, 4> bitwise_cast(const batch<float, 4>& x)
    {
        return _mm_castps_si128(x);
    }

    inline batch<double, 2> bitwise_cast(const batch<int64_t, 2>& x)
    {
        return _mm_castsi128_pd(x);
    }

    inline batch<int64_t, 2> bitwise_cast(const batch<double, 2>& x)
    {
        return _mm_castpd_si128(x);
    }
}

#endif
