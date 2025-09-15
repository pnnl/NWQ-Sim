// ==========================================================================
// Copyright (C) 2016 by Wolf Vollprecht and Johan Mabille
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.
// ==========================================================================

#ifndef XSIMD_HPP
#define XSIMD_HPP

#include "types/xsimd_types_include.hpp"

namespace xsimd
{
    template <class T, class Arch>
    batch<T, Arch> abs(const batch<T, Arch>& rhs);

    template <class T, class Arch>
    batch<T, Arch> fabs(const batch<T, Arch>& rhs);

    template <class T, class Arch>
    batch<T, Arch> sqrt(const batch<T, Arch>& rhs);

    template <class T, class Arch>
    batch<T, Arch> fma(const batch<T, Arch>& x, const batch<T, Arch>& y, const batch<T, Arch>& z);

    template <class T, class Arch>
    batch<T, Arch> fms(const batch<T, Arch>& x, const batch<T, Arch>& y, const batch<T, Arch>& z);

    template <class T, class Arch>
    batch<T, Arch> fnma(const batch<T, Arch>& x, const batch<T, Arch>& y, const batch<T, Arch>& z);

    template <class T, class Arch>
    batch<T, Arch> fnms(const batch<T, Arch>& x, const batch<T, Arch>& y, const batch<T, Arch>& z);

    template <class T, class Arch>
    batch<T, Arch> min(const batch<T, Arch>& lhs, const batch<T, Arch>& rhs);

    template <class T, class Arch>
    batch<T, Arch> max(const batch<T, Arch>& lhs, const batch<T, Arch>& rhs);

    template <class T, class Arch>
    batch<T, Arch> fmin(const batch<T, Arch>& lhs, const batch<T, Arch>& rhs);

    template <class T, class Arch>
    batch<T, Arch> fmax(const batch<T, Arch>& lhs, const batch<T, Arch>& rhs);

    template <class T, class Arch>
    batch<T, Arch> select(const batch_bool<T, Arch>& cond, const batch<T, Arch>& a, const batch<T, Arch>& b);

    template <class T, class Arch>
    batch<T, Arch> load_aligned(const T* src);

    template <class T, class Arch>
    batch<T, Arch> load_unaligned(const T* src);

    template <class T, class Arch>
    void store_aligned(T* dst, const batch<T, Arch>& src);

    template <class T, class Arch>
    void store_unaligned(T* dst, const batch<T, Arch>& src);
}

#include "xsimd/xsimd_available_architectures.hpp"
#include "xsimd/xsimd_complex.hpp"
#include "xsimd/xsimd_config.hpp"
#include "xsimd/xsimd_math.hpp"
#include "xsimd/xsimd_relational.hpp"
#include "xsimd/xsimd_serial.hpp"
#include "xsimd/xsimd_utils.hpp"

#endif
