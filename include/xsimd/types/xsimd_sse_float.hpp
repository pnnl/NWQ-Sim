// ==========================================================================
// Copyright (C) 2016 by Wolf Vollprecht and Johan Mabille
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.
// ==========================================================================

#ifndef XSIMD_SSE_FLOAT_HPP
#define XSIMD_SSE_FLOAT_HPP

#include <cstdint>
#include <immintrin.h>

#include "xsimd_base.hpp"

namespace xsimd
{

    /*************************
     * batch_bool<float, 4> *
     *************************/

    template <>
    struct simd_batch_traits<batch_bool<float, 4>>
    {
        using value_type = float;
        static constexpr std::size_t size = 4;
        using batch_type = batch<float, 4>;
        static constexpr std::size_t align = 16;
    };

    template <>
    class batch_bool<float, 4> : public simd_batch_bool<batch_bool<float, 4>>
    {
    public:
        batch_bool() = default;
        using base_type::base_type;

        batch_bool(const __m128& rhs)
            : m_value(rhs)
        {
        }

        batch_bool& operator=(const __m128& rhs)
        {
            m_value = rhs;
            return *this;
        }

        operator __m128() const
        {
            return m_value;
        }

    private:
        __m128 m_value;
    };

    namespace detail
    {
        template <>
        struct batch_bool_kernel<float, 4>
        {
            using batch_type = batch_bool<float, 4>;

            static batch_type bitwise_and(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_and_ps(lhs, rhs);
            }

            static batch_type bitwise_or(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_or_ps(lhs, rhs);
            }

            static batch_type bitwise_xor(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_xor_ps(lhs, rhs);
            }

            static batch_type bitwise_not(const batch_type& rhs)
            {
                return _mm_xor_ps(rhs, _mm_castsi128_ps(_mm_set1_epi32(-1)));
            }

            static batch_type bitwise_andnot(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_andnot_ps(lhs, rhs);
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
                return _mm_movemask_ps(rhs) == 15;
            }

            static bool any(const batch_type& rhs)
            {
                return _mm_movemask_ps(rhs) != 0;
            }
        };
    }

    /********************
     * batch<float, 4> *
     ********************/

    template <>
    struct simd_batch_traits<batch<float, 4>>
    {
        using value_type = float;
        static constexpr std::size_t size = 4;
        using batch_bool_type = batch_bool<float, 4>;
        static constexpr std::size_t align = 16;
        using storage_type = __m128;
    };

    template <>
    class batch<float, 4> : public simd_batch<batch<float, 4>>
    {
    public:
        using self_type = batch<float, 4>;
        using base_type = simd_batch<self_type>;

        batch() = default;
        using base_type::base_type;

        batch(const __m128& rhs)
            : base_type(rhs)
        {
        }

        batch& operator=(const __m128& rhs)
        {
            this->m_value = rhs;
            return *this;
        }

        operator __m128() const
        {
            return this->m_value;
        }

        batch<float, 4>& load_aligned(const float* src)
        {
            this->m_value = _mm_load_ps(src);
            return *this;
        }

        batch<float, 4>& load_unaligned(const float* src)
        {
            this->m_value = _mm_loadu_ps(src);
            return *this;
        }

        void store_aligned(float* dst) const
        {
            _mm_store_ps(dst, this->m_value);
        }

        void store_unaligned(float* dst) const
        {
            _mm_storeu_ps(dst, this->m_value);
        }
    };

    namespace detail
    {
        template <>
        struct batch_kernel<float, 4>
        {
            using batch_type = batch<float, 4>;
            using batch_bool_type = batch_bool<float, 4>;
            static constexpr bool is_int = false;

            static batch_type neg(const batch_type& rhs)
            {
                return _mm_xor_ps(rhs, _mm_set1_ps(-0.f));
            }

            static batch_type add(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_add_ps(lhs, rhs);
            }

            static batch_type sub(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_sub_ps(lhs, rhs);
            }

            static batch_type mul(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_mul_ps(lhs, rhs);
            }

            static batch_type div(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_div_ps(lhs, rhs);
            }

            static batch_bool_type eq(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_cmpeq_ps(lhs, rhs);
            }

            static batch_bool_type neq(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_cmpneq_ps(lhs, rhs);
            }

            static batch_bool_type lt(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_cmplt_ps(lhs, rhs);
            }

            static batch_bool_type lte(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_cmple_ps(lhs, rhs);
            }

            static batch_type bitwise_and(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_and_ps(lhs, rhs);
            }

            static batch_type bitwise_or(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_or_ps(lhs, rhs);
            }

            static batch_type bitwise_xor(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_xor_ps(lhs, rhs);
            }

            static batch_type bitwise_not(const batch_type& rhs)
            {
                return _mm_xor_ps(rhs, _mm_castsi128_ps(_mm_set1_epi32(-1)));
            }

            static batch_type bitwise_andnot(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_andnot_ps(lhs, rhs);
            }

            static batch_type min(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_min_ps(lhs, rhs);
            }

            static batch_type max(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_max_ps(lhs, rhs);
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
                __m128 sign_mask = _mm_set1_ps(-0.f);
                return _mm_andnot_ps(sign_mask, rhs);
            }

            static batch_type fabs(const batch_type& rhs)
            {
                return abs(rhs);
            }

            static batch_type sqrt(const batch_type& rhs)
            {
                return _mm_sqrt_ps(rhs);
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
                return _mm_floor_ps(rhs);
#else
                return batch_type(std::floor(rhs[0]), std::floor(rhs[1]), std::floor(rhs[2]), std::floor(rhs[3]));
#endif
            }

            static batch_type ceil(const batch_type& rhs)
            {
#if XSIMD_X86_INSTR_SET >= 3
                return _mm_ceil_ps(rhs);
#else
                return batch_type(std::ceil(rhs[0]), std::ceil(rhs[1]), std::ceil(rhs[2]), std::ceil(rhs[3]));
#endif
            }

            static batch_type round(const batch_type& rhs)
            {
#if XSIMD_X86_INSTR_SET >= 3
                return _mm_round_ps(rhs, _MM_FROUND_TO_NEAREST_INT);
#else
                return batch_type(std::round(rhs[0]), std::round(rhs[1]), std::round(rhs[2]), std::round(rhs[3]));
#endif
            }

            static batch_type trunc(const batch_type& rhs)
            {
#if XSIMD_X86_INSTR_SET >= 3
                return _mm_round_ps(rhs, _MM_FROUND_TO_ZERO);
#else
                return batch_type(std::trunc(rhs[0]), std::trunc(rhs[1]), std::trunc(rhs[2]), std::trunc(rhs[3]));
#endif
            }

            static batch_type select(const batch_bool_type& cond, const batch_type& a, const batch_type& b)
            {
                return _mm_or_ps(_mm_and_ps(cond, a), _mm_andnot_ps(cond, b));
            }

            static float hadd(const batch_type& rhs)
            {
                __m128 tmp0 = _mm_add_ps(rhs, _mm_movehl_ps(rhs, rhs));
                __m128 tmp1 = _mm_add_ss(tmp0, _mm_shuffle_ps(tmp0, tmp0, 1));
                return _mm_cvtss_f32(tmp1);
            }
        };
    }
}

#endif
