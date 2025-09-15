// ==========================================================================
// Copyright (C) 2016 by Wolf Vollprecht and Johan Mabille
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.
// ==========================================================================

#ifndef XSIMD_SSE_INT_HPP
#define XSIMD_SSE_INT_HPP

#include <cstdint>
#include <immintrin.h>

#include "xsimd_base.hpp"

namespace xsimd
{

    /**************************
     * batch_bool<int32_t, 4> *
     **************************/

    template <>
    struct simd_batch_traits<batch_bool<int32_t, 4>>
    {
        using value_type = int32_t;
        static constexpr std::size_t size = 4;
        using batch_type = batch<int32_t, 4>;
        static constexpr std::size_t align = 16;
    };

    template <>
    class batch_bool<int32_t, 4> : public simd_batch_bool<batch_bool<int32_t, 4>>
    {
    public:
        batch_bool() = default;
        using base_type::base_type;

        batch_bool(const __m128i& rhs)
            : m_value(rhs)
        {
        }

        batch_bool& operator=(const __m128i& rhs)
        {
            m_value = rhs;
            return *this;
        }

        operator __m128i() const
        {
            return m_value;
        }

    private:
        __m128i m_value;
    };

    namespace detail
    {
        template <>
        struct batch_bool_kernel<int32_t, 4>
        {
            using batch_type = batch_bool<int32_t, 4>;

            static batch_type bitwise_and(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_and_si128(lhs, rhs);
            }

            static batch_type bitwise_or(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_or_si128(lhs, rhs);
            }

            static batch_type bitwise_xor(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_xor_si128(lhs, rhs);
            }

            static batch_type bitwise_not(const batch_type& rhs)
            {
                return _mm_xor_si128(rhs, _mm_set1_epi32(-1));
            }

            static batch_type bitwise_andnot(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_andnot_si128(lhs, rhs);
            }

            static batch_type equal(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_cmpeq_epi32(lhs, rhs);
            }

            static batch_type not_equal(const batch_type& lhs, const batch_type& rhs)
            {
                return bitwise_not(equal(lhs, rhs));
            }

            static bool all(const batch_type& rhs)
            {
                return _mm_movemask_epi8(rhs) == 0xFFFF;
            }

            static bool any(const batch_type& rhs)
            {
                return _mm_movemask_epi8(rhs) != 0;
            }
        };
    }

    /*********************
     * batch<int32_t, 4> *
     *********************/

    template <>
    struct simd_batch_traits<batch<int32_t, 4>>
    {
        using value_type = int32_t;
        static constexpr std::size_t size = 4;
        using batch_bool_type = batch_bool<int32_t, 4>;
        static constexpr std::size_t align = 16;
        using storage_type = __m128i;
    };

    template <>
    class batch<int32_t, 4> : public simd_batch<batch<int32_t, 4>>
    {
    public:
        using self_type = batch<int32_t, 4>;
        using base_type = simd_batch<self_type>;

        batch() = default;
        using base_type::base_type;

        batch(const __m128i& rhs)
            : base_type(rhs)
        {
        }

        batch& operator=(const __m128i& rhs)
        {
            this->m_value = rhs;
            return *this;
        }

        operator __m128i() const
        {
            return this->m_value;
        }

        batch<int32_t, 4>& load_aligned(const int32_t* src)
        {
            this->m_value = _mm_load_si128((__m128i const*)src);
            return *this;
        }

        batch<int32_t, 4>& load_unaligned(const int32_t* src)
        {
            this->m_value = _mm_loadu_si128((__m128i const*)src);
            return *this;
        }

        void store_aligned(int32_t* dst) const
        {
            _mm_store_si128((__m128i*)dst, this->m_value);
        }

        void store_unaligned(int32_t* dst) const
        {
            _mm_storeu_si128((__m128i*)dst, this->m_value);
        }
    };

    namespace detail
    {
        template <>
        struct batch_kernel<int32_t, 4>
        {
            using batch_type = batch<int32_t, 4>;
            using batch_bool_type = batch_bool<int32_t, 4>;
            static constexpr bool is_int = true;

            static batch_type neg(const batch_type& rhs)
            {
                return _mm_sub_epi32(_mm_setzero_si128(), rhs);
            }

            static batch_type add(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_add_epi32(lhs, rhs);
            }

            static batch_type sub(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_sub_epi32(lhs, rhs);
            }

            static batch_type mul(const batch_type& lhs, const batch_type& rhs)
            {
#if XSIMD_X86_INSTR_SET >= 3
                return _mm_mullo_epi32(lhs, rhs);
#else
                // SSE2 does not have a 32-bit integer multiplication instruction.
                // This is a workaround.
                __m128i tmp1 = _mm_mul_epu32(lhs, rhs);
                __m128i tmp2 = _mm_mul_epu32(_mm_srli_si128(lhs, 4), _mm_srli_si128(rhs, 4));
                return _mm_unpacklo_epi32(_mm_shuffle_epi32(tmp1, _MM_SHUFFLE(0, 0, 2, 0)),
                                          _mm_shuffle_epi32(tmp2, _MM_SHUFFLE(0, 0, 2, 0)));
#endif
            }

            static batch_type div(const batch_type& lhs, const batch_type& rhs)
            {
                return batch_type(lhs[0] / rhs[0], lhs[1] / rhs[1], lhs[2] / rhs[2], lhs[3] / rhs[3]);
            }

            static batch_bool_type eq(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_cmpeq_epi32(lhs, rhs);
            }

            static batch_bool_type neq(const batch_type& lhs, const batch_type& rhs)
            {
                return ~(eq(lhs, rhs));
            }

            static batch_bool_type lt(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_cmplt_epi32(lhs, rhs);
            }

            static batch_bool_type lte(const batch_type& lhs, const batch_type& rhs)
            {
                return ~(_mm_cmplt_epi32(rhs, lhs));
            }

            static batch_type bitwise_and(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_and_si128(lhs, rhs);
            }

            static batch_type bitwise_or(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_or_si128(lhs, rhs);
            }

            static batch_type bitwise_xor(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_xor_si128(lhs, rhs);
            }

            static batch_type bitwise_not(const batch_type& rhs)
            {
                return _mm_xor_si128(rhs, _mm_set1_epi32(-1));
            }

            static batch_type bitwise_andnot(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_andnot_si128(lhs, rhs);
            }

            static batch_type min(const batch_type& lhs, const batch_type& rhs)
            {
#if XSIMD_X86_INSTR_SET >= 3
                return _mm_min_epi32(lhs, rhs);
#else
                __m128i cond = _mm_cmplt_epi32(lhs, rhs);
                return select(cond, lhs, rhs);
#endif
            }

            static batch_type max(const batch_type& lhs, const batch_type& rhs)
            {
#if XSIMD_X86_INSTR_SET >= 3
                return _mm_max_epi32(lhs, rhs);
#else
                __m128i cond = _mm_cmpgt_epi32(lhs, rhs);
                return select(cond, lhs, rhs);
#endif
            }

            static batch_type abs(const batch_type& rhs)
            {
#if XSIMD_X86_INSTR_SET >= 2
                return _mm_abs_epi32(rhs);
#else
                __m128i sign = _mm_cmplt_epi32(rhs, _mm_setzero_si128());
                __m128i inv = _mm_xor_si128(rhs, sign);
                return _mm_sub_epi32(inv, sign);
#endif
            }

            static batch_type select(const batch_bool_type& cond, const batch_type& a, const batch_type& b)
            {
                return _mm_or_si128(_mm_and_si128(cond, a), _mm_andnot_si128(cond, b));
            }

            static int32_t hadd(const batch_type& rhs)
            {
                __m128i tmp1 = _mm_add_epi32(rhs, _mm_shuffle_epi32(rhs, 0x4E));
                __m128i tmp2 = _mm_add_epi32(tmp1, _mm_shuffle_epi32(tmp1, 0xB1));
                return _mm_cvtsi128_si32(tmp2);
            }
        };
    }

    /**************************
     * batch_bool<int64_t, 2> *
     **************************/

    template <>
    struct simd_batch_traits<batch_bool<int64_t, 2>>
    {
        using value_type = int64_t;
        static constexpr std::size_t size = 2;
        using batch_type = batch<int64_t, 2>;
        static constexpr std::size_t align = 16;
    };

    template <>
    class batch_bool<int64_t, 2> : public simd_batch_bool<batch_bool<int64_t, 2>>
    {
    public:
        batch_bool() = default;
        using base_type::base_type;

        batch_bool(const __m128i& rhs)
            : m_value(rhs)
        {
        }

        batch_bool& operator=(const __m128i& rhs)
        {
            m_value = rhs;
            return *this;
        }

        operator __m128i() const
        {
            return m_value;
        }

    private:
        __m128i m_value;
    };

    namespace detail
    {
        template <>
        struct batch_bool_kernel<int64_t, 2>
        {
            using batch_type = batch_bool<int64_t, 2>;

            static batch_type bitwise_and(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_and_si128(lhs, rhs);
            }

            static batch_type bitwise_or(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_or_si128(lhs, rhs);
            }

            static batch_type bitwise_xor(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_xor_si128(lhs, rhs);
            }

            static batch_type bitwise_not(const batch_type& rhs)
            {
                return _mm_xor_si128(rhs, _mm_set1_epi32(-1));
            }

            static batch_type bitwise_andnot(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_andnot_si128(lhs, rhs);
            }

            static batch_type equal(const batch_type& lhs, const batch_type& rhs)
            {
#if XSIMD_X86_INSTR_SET >= 3
                return _mm_cmpeq_epi64(lhs, rhs);
#else
                __m128i tmp = _mm_cmpeq_epi32(lhs, rhs);
                return _mm_and_si128(tmp, _mm_shuffle_epi32(tmp, 0xB1));
#endif
            }

            static batch_type not_equal(const batch_type& lhs, const batch_type& rhs)
            {
                return bitwise_not(equal(lhs, rhs));
            }

            static bool all(const batch_type& rhs)
            {
                return _mm_movemask_epi8(rhs) == 0xFFFF;
            }

            static bool any(const batch_type& rhs)
            {
                return _mm_movemask_epi8(rhs) != 0;
            }
        };
    }

    /*********************
     * batch<int64_t, 2> *
     *********************/

    template <>
    struct simd_batch_traits<batch<int64_t, 2>>
    {
        using value_type = int64_t;
        static constexpr std::size_t size = 2;
        using batch_bool_type = batch_bool<int64_t, 2>;
        static constexpr std::size_t align = 16;
        using storage_type = __m128i;
    };

    template <>
    class batch<int64_t, 2> : public simd_batch<batch<int64_t, 2>>
    {
    public:
        using self_type = batch<int64_t, 2>;
        using base_type = simd_batch<self_type>;

        batch() = default;
        using base_type::base_type;

        batch(const __m128i& rhs)
            : base_type(rhs)
        {
        }

        batch& operator=(const __m128i& rhs)
        {
            this->m_value = rhs;
            return *this;
        }

        operator __m128i() const
        {
            return this->m_value;
        }

        batch<int64_t, 2>& load_aligned(const int64_t* src)
        {
            this->m_value = _mm_load_si128((__m128i const*)src);
            return *this;
        }

        batch<int64_t, 2>& load_unaligned(const int64_t* src)
        {
            this->m_value = _mm_loadu_si128((__m128i const*)src);
            return *this;
        }

        void store_aligned(int64_t* dst) const
        {
            _mm_store_si128((__m128i*)dst, this->m_value);
        }

        void store_unaligned(int64_t* dst) const
        {
            _mm_storeu_si128((__m128i*)dst, this->m_value);
        }
    };

    namespace detail
    {
        template <>
        struct batch_kernel<int64_t, 2>
        {
            using batch_type = batch<int64_t, 2>;
            using batch_bool_type = batch_bool<int64_t, 2>;
            static constexpr bool is_int = true;

            static batch_type neg(const batch_type& rhs)
            {
                return _mm_sub_epi64(_mm_setzero_si128(), rhs);
            }

            static batch_type add(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_add_epi64(lhs, rhs);
            }

            static batch_type sub(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_sub_epi64(lhs, rhs);
            }

            static batch_type mul(const batch_type& lhs, const batch_type& rhs)
            {
                return batch_type(lhs[0] * rhs[0], lhs[1] * rhs[1]);
            }

            static batch_type div(const batch_type& lhs, const batch_type& rhs)
            {
                return batch_type(lhs[0] / rhs[0], lhs[1] / rhs[1]);
            }

            static batch_bool_type eq(const batch_type& lhs, const batch_type& rhs)
            {
#if XSIMD_X86_INSTR_SET >= 3
                return _mm_cmpeq_epi64(lhs, rhs);
#else
                __m128i tmp = _mm_cmpeq_epi32(lhs, rhs);
                return _mm_and_si128(tmp, _mm_shuffle_epi32(tmp, 0xB1));
#endif
            }

            static batch_bool_type neq(const batch_type& lhs, const batch_type& rhs)
            {
                return ~(eq(lhs, rhs));
            }

            static batch_bool_type lt(const batch_type& lhs, const batch_type& rhs)
            {
#if XSIMD_X86_INSTR_SET >= 3
                return _mm_cmpgt_epi64(rhs, lhs);
#else
                return batch_bool_type(lhs[0] < rhs[0], lhs[1] < rhs[1]);
#endif
            }

            static batch_bool_type lte(const batch_type& lhs, const batch_type& rhs)
            {
                return ~(lt(rhs, lhs));
            }

            static batch_type bitwise_and(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_and_si128(lhs, rhs);
            }

            static batch_type bitwise_or(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_or_si128(lhs, rhs);
            }

            static batch_type bitwise_xor(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_xor_si128(lhs, rhs);
            }

            static batch_type bitwise_not(const batch_type& rhs)
            {
                return _mm_xor_si128(rhs, _mm_set1_epi32(-1));
            }

            static batch_type bitwise_andnot(const batch_type& lhs, const batch_type& rhs)
            {
                return _mm_andnot_si128(lhs, rhs);
            }

            static batch_type min(const batch_type& lhs, const batch_type& rhs)
            {
                return select(lhs < rhs, lhs, rhs);
            }

            static batch_type max(const batch_type& lhs, const batch_type& rhs)
            {
                return select(lhs > rhs, lhs, rhs);
            }

            static batch_type abs(const batch_type& rhs)
            {
                return select(rhs > 0, rhs, -rhs);
            }

            static batch_type select(const batch_bool_type& cond, const batch_type& a, const batch_type& b)
            {
                return _mm_or_si128(_mm_and_si128(cond, a), _mm_andnot_si128(cond, b));
            }

            static int64_t hadd(const batch_type& rhs)
            {
                __m128i tmp = _mm_add_epi64(rhs, _mm_shuffle_epi32(rhs, 0x4E));
                return _mm_cvtsi128_si64(tmp);
            }
        };
    }
}

#endif
