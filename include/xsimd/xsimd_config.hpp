// ==========================================================================
// Copyright (C) 2016 by Wolf Vollprecht and Johan Mabille
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.
// ==========================================================================

#ifndef XSIMD_CONFIG_HPP
#define XSIMD_CONFIG_HPP

// __SSE__ and __SSE2__ are defined by default on amd64 architecture.
// See https://gcc.gnu.org/onlinedocs/gcc-4.3.4/gcc/X86-Built-in-Functions.html
//
// __FMA__ is defined when the processor supports FMA instructions,
// see https://www.felixcloutier.com/x86/FMA.html
//
// __AVX__ and __AVX2__ are defined when the processor supports AVX and AVX2
// instructions, see https://www.felixcloutier.com/x86/AVX.html
//
// __AVX512F__ is defined for AVX-512 Foundation, see
// https://www.felixcloutier.com/x86/AVX-512.html
//
// All these macros should be defined by the compiler, but it's possible
// to force the use of an instruction set with the following macros.
//
// XSIMD_FORCE_SSE2
// XSIMD_FORCE_SSE3
// XSIMD_FORCE_SSSE3
// XSIMD_FORCE_SSE4_1
// XSIMD_FORCE_SSE4_2
// XSIMD_FORCE_FMA3_SSE
// XSIMD_FORCE_FMA4
// XSIMD_FORCE_AVX
// XSIMD_FORCE_FMA3_AVX
// XSIMD_FORCE_AVX2
// XSIMD_FORCE_FMA3_AVX2
// XSIMD_FORCE_AVX512F

#if defined(XSIMD_FORCE_SSE2)
#define __SSE2__ 1
#endif

#if defined(XSIMD_FORCE_SSE3)
#define __SSE3__ 1
#endif

#if defined(XSIMD_FORCE_SSSE3)
#define __SSSE3__ 1
#endif

#if defined(XSIMD_FORCE_SSE4_1)
#define __SSE4_1__ 1
#endif

#if defined(XSIMD_FORCE_SSE4_2)
#define __SSE4_2__ 1
#endif

#if defined(XSIMD_FORCE_FMA3_SSE)
#define __FMA__ 1
#define __SSE3__ 1
#endif

#if defined(XSIMD_FORCE_FMA4)
#define __FMA4__ 1
#endif

#if defined(XSIMD_FORCE_AVX)
#define __AVX__ 1
#endif

#if defined(XSIMD_FORCE_FMA3_AVX)
#define __FMA__ 1
#define __AVX__ 1
#endif

#if defined(XSIMD_FORCE_AVX2)
#define __AVX2__ 1
#endif

#if defined(XSIMD_FORCE_FMA3_AVX2)
#define __FMA__ 1
#define __AVX2__ 1
#endif

#if defined(XSIMD_FORCE_AVX512F)
#define __AVX512F__ 1
#endif

#if defined(__AVX512F__)
#define XSIMD_X86_INSTR_SET 6
#elif defined(__AVX2__)
#define XSIMD_X86_INSTR_SET 5
#elif defined(__AVX__)
#define XSIMD_X86_INSTR_SET 4
#elif defined(__SSE4_2__)
#define XSIMD_X86_INSTR_SET 3
#elif defined(__SSE4_1__)
#define XSIMD_X86_INSTR_SET 3
#elif defined(__SSSE3__)
#define XSIMD_X86_INSTR_SET 2
#elif defined(__SSE3__)
#define XSIMD_X86_INSTR_SET 1
#elif defined(__SSE2__)
#define XSIMD_X86_INSTR_SET 0
#endif

#if defined(__FMA__)
#if XSIMD_X86_INSTR_SET >= 5
#define XSIMD_X86_FMA_INSTR_SET 2
#elif XSIMD_X86_INSTR_SET >= 4
#define XSIMD_X86_FMA_INSTR_SET 1
#elif XSIMD_X86_INSTR_SET >= 1
#define XSIMD_X86_FMA_INSTR_SET 0
#endif
#endif

#if defined(__FMA4__)
#define XSIMD_X86_FMA4_INSTR_SET 1
#endif

#if defined(__x86_64__)
#define XSIMD_X86_64_BIT_ABI 1
#else
#define XSIMD_X86_64_BIT_ABI 0
#endif

#if defined(_MSC_VER)
#define XSIMD_ALL_SITE_STATIC_INIT(SITE) \
    static int SITE##_init_variable = (SITE.init(), 0);
#else
#define XSIMD_ALL_SITE_STATIC_INIT(SITE) \
    [[gnu::unused]] static int SITE##_init_variable = (SITE.init(), 0);
#endif

#endif
