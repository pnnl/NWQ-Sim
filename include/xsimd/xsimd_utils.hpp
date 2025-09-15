// ==========================================================================
// Copyright (C) 2016 by Wolf Vollprecht and Johan Mabille
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.
// ==========================================================================

#ifndef XSIMD_UTILS_HPP
#define XSIMD_UTILS_HPP

#include <cstddef>
#include <ostream>

#include "xsimd_config.hpp"

namespace xsimd
{
    template <class T>
    inline T* aligned_malloc(std::size_t size, std::size_t alignment)
    {
        if (alignment & (alignment - 1))
        {
            // alignment must be a power of 2
            return nullptr;
        }
        void* ptr;
        if (posix_memalign(&ptr, alignment, size * sizeof(T)) != 0)
        {
            return nullptr;
        }
        return static_cast<T*>(ptr);
    }

    inline void aligned_free(void* ptr)
    {
        free(ptr);
    }

    template <class T, std::size_t N>
    class batch;

    template <class T, std::size_t N>
    std::ostream& operator<<(std::ostream& out, const batch<T, N>& rhs)
    {
        out << '(';
        for (std::size_t i = 0; i < N - 1; ++i)
        {
            out << rhs[i] << ", ";
        }
        out << rhs[N - 1] << ')';
        return out;
    }

    template <class T, std::size_t N>
    class batch_bool;

    template <class T, std::size_t N>
    std::ostream& operator<<(std::ostream& out, const batch_bool<T, N>& rhs)
    {
        out << '(';
        for (std::size_t i = 0; i < N - 1; ++i)
        {
            out << rhs[i] << ", ";
        }
        out << rhs[N - 1] << ')';
        return out;
    }
}

#endif
