// ==========================================================================
// Copyright (C) 2016 by Wolf Vollprecht and Johan Mabille
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.
// ==========================================================================

#ifndef XSIMD_SSE_DOUBLE_HPP
#define XSIMD_SSE_DOUBLE_HPP

#include <cstdint>
#include <immintrin.h>

#include "xsimd_base.hpp"

namespace xsimd
{

    /**************************
     * batch_bool<double, 2> *
     **************************/

    template <>
    struct simd_batch_traits<batch_bool<double, 2>>
    {
        using value_type = double;
        static constexpr std::size_t size = 2;
        using batch_type = batch<double, 2>;
        static constexpr std::size_t align = 16;
    };

    template <>
    class batch_bool<double, 2> : public simd_batch_bool<batch_bool<double, 2>>
    {
    public:
        batch_bool() = default;
        using base_type::base_type;

        batch_bool(const __m128d& rhs)
            : m_value(rhs)
        {
        }

        batch_bool& operator=(const __m128d& rhs)
        {
            m_value = rhs;
            return *this;
        }

        operator __m128d() const
        {
            return m_value;
        }

    private:
        __m128d m_value;
    };

    namespace detail
    {
        template <>
        struct batch_bool_kernel<double, 2>
        {
            using batch_type = batch_bool<double, 2>;

            static batch_type bitwise_and(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_and_pd(lhs, rhs);
            }

            static batch_type bitwise_or(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_or_pd(lhs, rhs);
            }

            static batch_type bitwise_xor(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_xor_pd(lhs, rhs);
            }

            static batch_type bitwise_not(const batch_type& rhs)
            {
                return _mm_xor_pd(rhs, _mm_castsi128_pd(_mm_set1_epi32(-1)));
            }


            static batch_type bitwise_andnot(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_andnot_pd(lhs, rhs);
            }

            static batch_type equal(const batch_type& lhs, const batch_type& rhs)
            {
                return bitwise_not(bitwise_xor(lhs, rhs));
            }

            static batch_type not_equal(const batch_type& lhs, const batch_type& rhs)
            {
                return bitwise_xor(lhs, rhs);
            }

            static bool all(const batch_type& rhs)
            {
                return _mm_movemask_pd(rhs) == 3;
            }

            static bool any(const batch_type& rhs)
            {
                return _mm_movemask_pd(rhs) != 0;
            }
        };
    }

    /*********************
     * batch<double, 2> *
     *********************/

    template <>
    struct simd_batch_traits<batch<double, 2>>
    {
        using value_type = double;
        static constexpr std::size_t size = 2;
        using batch_bool_type = batch_bool<double, 2>;
        static constexpr std::size_t align = 16;
        using storage_type = __m128d;
    };

    template <>
    class batch<double, 2> : public simd_batch<batch<double, 2>>
    {
    public:
        using self_type = batch<double, 2>;
        using base_type = simd_batch<self_type>;

        batch() = default;
        using base_type::base_type;

        batch(const __m128d& rhs)
            : base_type(rhs)
        {
        }

        batch& operator=(const __m128d& rhs)
        {
            this->m_value = rhs;
            return *this;
        }

        operator __m128d() const
        {
            return this->m_value;
        }

        batch<double, 2>& load_aligned(const double* src)
        {
            this->m_value = _mm_load_pd(src);
            return *this;
        }

        batch<double, 2>& load_unaligned(const double* src)
        {
            this->m_value = _mm_loadu_pd(src);
            return *this;
        }

        void store_aligned(double* dst) const
        {
            _mm_store_pd(dst, this->m_value);
        }

        void store_unaligned(double* dst) const
        {
            _mm_storeu_pd(dst, this->m_value);
        }
    };

    namespace detail
    {
        template <>
        struct batch_kernel<double, 2>
        {
            using batch_type = batch<double, 2>;
            using batch_bool_type = batch_bool<double, 2>;
            static constexpr bool is_int = false;

            static batch_type neg(const batch_type& rhs)
            {
                return _mm_xor_pd(rhs, _mm_set1_pd(-0.));
            }

            static batch_type add(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_add_pd(lhs, rhs);
            }

            static batch_type sub(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_sub_pd(lhs, rhs);
            }

            static batch_type mul(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_mul_pd(lhs, rhs);
            }

            static batch_type div(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_div_pd(lhs, rhs);
            }

            static batch_bool_type eq(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_cmpeq_pd(lhs, rhs);
            }

            static batch_bool_type neq(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_cmpneq_pd(lhs, rhs);
            }

            static batch_bool_type lt(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_cmplt_pd(lhs, rhs);
            }

            static batch_bool_type lte(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_cmple_pd(lhs, rhs);
            }

            static batch_type bitwise_and(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_and_pd(lhs, rhs);
            }

            static batch_type bitwise_or(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_or_pd(lhs, rhs);
            }

            static batch_type bitwise_xor(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_xor_pd(lhs, rhs);
            }

            static batch_type bitwise_not(const batch_type& rhs)
            {
                return _mm_xor_pd(rhs, _mm_castsi128_pd(_mm_set1_epi32(-1)));
            }

            static batch_type bitwise_andnot(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_andnot_pd(lhs, rhs);
            }

            static batch_type min(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_min_pd(lhs, rhs);
            }

            static batch_type max(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_max_pd(lhs, rhs);
            }

            static batch_type fmin(const batch_type& lhs, const batch_type& rhs)
            {
                return min(lhs, rhs);
            }

            static batch_type fmax(const batch_type& lhs, const batch_type& rhs)
            {
                return max(lhs, rhs);
            }

            static batch_type abs(const batch_type& rhs)
            {
                __m128d sign_mask = _mm_set1_pd(-0.);
                return _mm_andnot_pd(sign_mask, rhs);
            }

            static batch_type fabs(const batch_type& rhs)
            {
                return abs(rhs);
            }

            static batch_type sqrt(const batch_type& rhs)
            {
                return _mm_sqrt_pd(rhs);
            }

            static batch_type fma(const batch_type& x, const batch_type& y, const batch_type& z)
            {
                return add(mul(x, y), z);
            }

            static batch_type fms(const batch_type& x, const batch_type& y, const batch_type& z)
            {
                return sub(mul(x, y), z);
            }

            static batch_type fnma(const batch_type& x, const batch_type& y, const batch_type& z)
            {
                return sub(z, mul(x, y));
            }

            static batch_type fnms(const batch_type& x, const batch_type& y, const batch_type& z)
            {
                return neg(add(mul(x, y), z));
            }

            static batch_type floor(const batch_type& rhs)
            {
#if XSIMD_X86_INSTR_SET >= 3
                return _mm_floor_pd(rhs);
#else
                return batch_type(std::floor(rhs[0]), std::floor(rhs[1]));
#endif
            }

            static batch_type ceil(const batch_type& rhs)
            {
#if XSIMD_X86_INSTR_SET >= 3
                return _mm_ceil_pd(rhs);
#else
                return batch_type(std::ceil(rhs[0]), std::ceil(rhs[1]));
#endif
            }

            static batch_type round(const batch_type& rhs)
            {
#if XSIMD_X86_INSTR_SET >= 3
                return _mm_round_pd(rhs, _MM_FROUND_TO_NEAREST_INT);
#else
                return batch_type(std::round(rhs[0]), std::round(rhs[1]));
#endif
            }

            static batch_type trunc(const batch_type& rhs)
            {
#if XSIMD_X86_INSTR_SET >= 3
                return _mm_round_pd(rhs, _MM_FROUND_TO_ZERO);
#else
                return batch_type(std::trunc(rhs[0]), std::trunc(rhs[1]));
#endif
            }

            static batch_type select(const batch_bool_type& cond, const batch_type& a, const batch_type& b)
            {
                return _mm_or_pd(_mm_and_pd(cond, a), _mm_andnot_pd(cond, b));
            }

            static double hadd(const batch_type& rhs)
            {
                __m128d tmp = _mm_add_pd(rhs, _mm_shuffle_pd(rhs, rhs, 1));
                return _mm_cvtsd_f64(tmp);
            }
        };
    }
}

#endif
