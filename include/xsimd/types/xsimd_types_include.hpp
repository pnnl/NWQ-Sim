// ==========================================================================
// Copyright (C) 2016 by Wolf Vollprecht and Johan Mabille
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.
// ==========================================================================

#ifndef XSIMD_TYPES_HPP
#define XSIMD_TYPES_HPP

#include <complex>
#include <cstddef>
#include <functional>

#include "../xsimd_config.hpp"
#include "../xsimd_utils.hpp"

#if XSIMD_WITH_XTL
#include "xtl/xcomplex.hpp"
#endif

namespace xsimd
{

    /********************
     * batch<T, N>      *
     ********************/

    template <class T, std::size_t N>
    class batch
    {
    public:
        using value_type = T;
        static constexpr std::size_t size = N;
        using batch_bool_type = batch<bool, N>;

        // Constructors
        batch() = default;
        batch(const batch<T, N>&) = default;
        batch<T, N>& operator=(const batch<T, N>&) = default;
        batch(batch<T, N>&&) = default;
        batch<T, N>& operator=(batch<T, N>&&) = default;

        batch(T val);
        template <class... Args>
        batch(Args... args);
        batch(const T* buffer);

        // Store methods
        void store_unaligned(T* buffer) const;
        void store_aligned(T* buffer) const;

        // Array access
        T operator[](std::size_t index) const;

        // Unary operators
        batch_bool_type operator==(const batch<T, N>& other) const;
        batch_bool_type operator!=(const batch<T, N>& other) const;
        batch_bool_type operator<(const batch<T, N>& other) const;
        batch_bool_type operator<=(const batch<T, N>& other) const;

        batch<T, N> operator-() const;
        batch<T, N> operator+() const;

        // Arithmetic operators
        batch<T, N> operator+(const batch<T, N>& other) const;
        batch<T, N> operator-(const batch<T, N>& other) const;
        batch<T, N> operator*(const batch<T, N>& other) const;
        batch<T, N> operator/(const batch<T, N>& other) const;

        // Bitwise operators
        batch<T, N> operator&(const batch<T, N>& other) const;
        batch<T, N> operator|(const batch<T, N>& other) const;
        batch<T, N> operator^(const batch<T, N>& other) const;
        batch<T, N> operator~() const;

        // Assignment operators
        batch<T, N>& operator+=(const batch<T, N>& other);
        batch<T, N>& operator-=(const batch<T, N>& other);
        batch<T, N>& operator*=(const batch<T, N>& other);
        batch<T, N>& operator/=(const batch<T, N>& other);
        batch<T, N>& operator&=(const batch<T, N>& other);
        batch<T, N>& operator|=(const batch<T, N>& other);
        batch<T, N>& operator^=(const batch<T, N>& other);

    private:
        T m_value[N];
    };

    /********************
     * batch<bool, N>   *
     ********************/

    template <class T, std::size_t N>
    class batch_bool
    {
    public:
        using value_type = bool;
        static constexpr std::size_t size = N;

        // Constructors
        batch_bool() = default;
        batch_bool(const batch_bool<T, N>&) = default;
        batch_bool<T, N>& operator=(const batch_bool<T, N>&) = default;
        batch_bool(batch_bool<T, N>&&) = default;
        batch_bool<T, N>& operator=(batch_bool<T, N>&&) = default;

        batch_bool(bool val);
        template <class... Args>
        batch_bool(Args... args);
        batch_bool(const bool* buffer);

        // Store methods
        void store_unaligned(bool* buffer) const;

        // Array access
        bool operator[](std::size_t index) const;

        // Bitwise operators
        batch_bool<T, N> operator&(const batch_bool<T, N>& other) const;
        batch_bool<T, N> operator|(const batch_bool<T, N>& other) const;
        batch_bool<T, N> operator^(const batch_bool<T, N>& other) const;
        batch_bool<T, N> operator~() const;
        batch_bool<T, N> operator==(const batch_bool<T, N>& other) const;
        batch_bool<T, N> operator!=(const batch_bool<T, N>& other) const;

        // Assignment operators
        batch_bool<T, N>& operator&=(const batch_bool<T, N>& other);
        batch_bool<T, N>& operator|=(const batch_bool<T, N>& other);
        batch_bool<T, N>& operator^=(const batch_bool<T, N>& other);

    private:
        T m_value;
    };

    /*************************
     * batch<complex<T>, N>  *
     *************************/

    template <class T, std::size_t N>
    class batch<std::complex<T>, N>
    {
    public:
        using value_type = std::complex<T>;
        static constexpr std::size_t size = N;
        using real_batch = batch<T, N>;
        using batch_bool_type = batch_bool<T, N>;

        // Constructors
        batch() = default;
        batch(const batch<std::complex<T>, N>&) = default;
        batch<std::complex<T>, N>& operator=(const batch<std::complex<T>, N>&) = default;
        batch(batch<std::complex<T>, N>&&) = default;
        batch<std::complex<T>, N>& operator=(batch<std::complex<T>, N>&&) = default;

        batch(const std::complex<T>& val);
        batch(const real_batch& real, const real_batch& imag);
        batch(const std::complex<T>* buffer);

        // Store methods
        void store_unaligned(std::complex<T>* buffer) const;
        void store_aligned(std::complex<T>* buffer) const;

        // Array access
        std::complex<T> operator[](std::size_t index) const;

        // Unary operators
        batch_bool_type operator==(const batch<std::complex<T>, N>& other) const;
        batch_bool_type operator!=(const batch<std::complex<T>, N>& other) const;

        batch<std::complex<T>, N> operator-() const;
        batch<T, N> operator+() const;

        // Arithmetic operators
        batch<std::complex<T>, N> operator+(const batch<std::complex<T>, N>& other) const;
        batch<std::complex<T>, N> operator-(const batch<std::complex<T>, N>& other) const;
        batch<std::complex<T>, N> operator*(const batch<std::complex<T>, N>& other) const;
        batch<std::complex<T>, N> operator/(const batch<std::complex<T>, N>& other) const;

        // Assignment operators
        batch<std::complex<T>, N>& operator+=(const batch<std::complex<T>, N>& other);
        batch<std::complex<T>, N>& operator-=(const batch<std::complex<T>, N>& other);
        batch<std::complex<T>, N>& operator*=(const batch<std::complex<T>, N>& other);
        batch<std::complex<T>, N>& operator/=(const batch<std::complex<T>, N>& other);

        // Accessors
        real_batch real() const;
        real_batch imag() const;

    private:
        real_batch m_real;
        real_batch m_imag;
    };

#if XSIMD_WITH_XTL
    /*************************
     * batch<xtl::xcomplex<T, T, i>, N>  *
     *************************/

    template <class T, std::size_t N>
    class batch<xtl::xcomplex<T, T>, N>
    {
    public:
        using value_type = xtl::xcomplex<T, T>;
        static constexpr std::size_t size = N;
        using real_batch = batch<T, N>;
        using batch_bool_type = batch_bool<T, N>;

        // Constructors
        batch() = default;
        batch(const batch<xtl::xcomplex<T, T>, N>&) = default;
        batch<xtl::xcomplex<T, T>, N>& operator=(const batch<xtl::xcomplex<T, T>, N>&) = default;
        batch(batch<xtl::xcomplex<T, T>, N>&&) = default;
        batch<xtl::xcomplex<T, T>, N>& operator=(batch<xtl::xcomplex<T, T>, N>&&) = default;

        batch(const xtl::xcomplex<T, T>& val);
        batch(const real_batch& real, const real_batch& imag);
        batch(const xtl::xcomplex<T, T>* buffer);

        // Store methods
        void store_unaligned(xtl::xcomplex<T, T>* buffer) const;
        void store_aligned(xtl::xcomplex<T, T>* buffer) const;

        // Array access
        xtl::xcomplex<T, T> operator[](std::size_t index) const;

        // Unary operators
        batch_bool_type operator==(const batch<xtl::xcomplex<T, T>, N>& other) const;
        batch_bool_type operator!=(const batch<xtl::xcomplex<T, T>, N>& other) const;

        batch<xtl::xcomplex<T, T>, N> operator-() const;
        batch<T, N> operator+() const;

        // Arithmetic operators
        batch<xtl::xcomplex<T, T>, N> operator+(const batch<xtl::xcomplex<T, T>, N>& other) const;
        batch<xtl::xcomplex<T, T>, N> operator-(const batch<xtl::xcomplex<T, T>, N>& other) const;
        batch<xtl::xcomplex<T, T>, N> operator*(const batch<xtl::xcomplex<T, T>, N>& other) const;
        batch<xtl::xcomplex<T, T>, N> operator/(const batch<xtl::xcomplex<T, T>, N>& other) const;

        // Assignment operators
        batch<xtl::xcomplex<T, T>, N>& operator+=(const batch<xtl::xcomplex<T, T>, N>& other);
        batch<xtl::xcomplex<T, T>, N>& operator-=(const batch<xtl::xcomplex<T, T>, N>& other);
        batch<xtl::xcomplex<T, T>, N>& operator*=(const batch<xtl::xcomplex<T, T>, N>& other);
        batch<xtl::xcomplex<T, T>, N>& operator/=(const batch<xtl::xcomplex<T, T>, N>& other);

        // Accessors
        real_batch real() const;
        real_batch imag() const;

    private:
        real_batch m_real;
        real_batch m_imag;
    };
#endif

    /**************************
     * batch implementation   *
     **************************/

    template <class T, std::size_t N>
    inline batch<T, N>::batch(T val)
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            m_value[i] = val;
        }
    }

    template <class T, std::size_t N>
    template <class... Args>
    inline batch<T, N>::batch(Args... args)
        : m_value{static_cast<T>(args)...}
    {
    }

    template <class T, std::size_t N>
    inline batch<T, N>::batch(const T* buffer)
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            m_value[i] = buffer[i];
        }
    }

    template <class T, std::size_t N>
    inline void batch<T, N>::store_unaligned(T* buffer) const
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            buffer[i] = m_value[i];
        }
    }

    template <class T, std::size_t N>
    inline void batch<T, N>::store_aligned(T* buffer) const
    {
        store_unaligned(buffer);
    }

    template <class T, std::size_t N>
    inline T batch<T, N>::operator[](std::size_t index) const
    {
        return m_value[index];
    }

    template <class T, std::size_t N>
    inline batch<bool, N> batch<T, N>::operator==(const batch<T, N>& other) const
    {
        bool buffer[N];
        for (std::size_t i = 0; i < N; ++i)
        {
            buffer[i] = m_value[i] == other.m_value[i];
        }
        return batch<bool, N>(buffer);
    }

    template <class T, std::size_t N>
    inline batch<bool, N> batch<T, N>::operator!=(const batch<T, N>& other) const
    {
        bool buffer[N];
        for (std::size_t i = 0; i < N; ++i)
        {
            buffer[i] = m_value[i] != other.m_value[i];
        }
        return batch<bool, N>(buffer);
    }

    template <class T, std::size_t N>
    inline batch<bool, N> batch<T, N>::operator<(const batch<T, N>& other) const
    {
        bool buffer[N];
        for (std::size_t i = 0; i < N; ++i)
        {
            buffer[i] = m_value[i] < other.m_value[i];
        }
        return batch<bool, N>(buffer);
    }

    template <class T, std::size_t N>
    inline batch<bool, N> batch<T, N>::operator<=(const batch<T, N>& other) const
    {
        bool buffer[N];
        for (std::size_t i = 0; i < N; ++i)
        {
            buffer[i] = m_value[i] <= other.m_value[i];
        }
        return batch<bool, N>(buffer);
    }

    template <class T, std::size_t N>
    inline batch<T, N> batch<T, N>::operator-() const
    {
        T buffer[N];
        for (std::size_t i = 0; i < N; ++i)
        {
            buffer[i] = -m_value[i];
        }
        return batch<T, N>(buffer);
    }

    template <class T, std::size_t N>
    inline batch<T, N> batch<T, N>::operator+() const
    {
        return *this;
    }

    template <class T, std::size_t N>
    inline batch<T, N> batch<T, N>::operator+(const batch<T, N>& other) const
    {
        T buffer[N];
        for (std::size_t i = 0; i < N; ++i)
        {
            buffer[i] = m_value[i] + other.m_value[i];
        }
        return batch<T, N>(buffer);
    }

    template <class T, std::size_t N>
    inline batch<T, N> batch<T, N>::operator-(const batch<T, N>& other) const
    {
        T buffer[N];
        for (std::size_t i = 0; i < N; ++i)
        {
            buffer[i] = m_value[i] - other.m_value[i];
        }
        return batch<T, N>(buffer);
    }

    template <class T, std::size_t N>
    inline batch<T, N> batch<T, N>::operator*(const batch<T, N>& other) const
    {
        T buffer[N];
        for (std::size_t i = 0; i < N; ++i)
        {
            buffer[i] = m_value[i] * other.m_value[i];
        }
        return batch<T, N>(buffer);
    }

    template <class T, std::size_t N>
    inline batch<T, N> batch<T, N>::operator/(const batch<T, N>& other) const
    {
        T buffer[N];
        for (std::size_t i = 0; i < N; ++i)
        {
            buffer[i] = m_value[i] / other.m_value[i];
        }
        return batch<T, N>(buffer);
    }

    template <class T, std::size_t N>
    inline batch<T, N> batch<T, N>::operator&(const batch<T, N>& other) const
    {
        T buffer[N];
        for (std::size_t i = 0; i < N; ++i)
        {
            buffer[i] = m_value[i] & other.m_value[i];
        }
        return batch<T, N>(buffer);
    }

    template <class T, std::size_t N>
    inline batch<T, N> batch<T, N>::operator|(const batch<T, N>& other) const
    {
        T buffer[N];
        for (std::size_t i = 0; i < N; ++i)
        {
            buffer[i] = m_value[i] | other.m_value[i];
        }
        return batch<T, N>(buffer);
    }

    template <class T, std::size_t N>
    inline batch<T, N> batch<T, N>::operator^(const batch<T, N>& other) const
    {
        T buffer[N];
        for (std::size_t i = 0; i < N; ++i)
        {
            buffer[i] = m_value[i] ^ other.m_value[i];
        }
        return batch<T, N>(buffer);
    }

    template <class T, std::size_t N>
    inline batch<T, N> batch<T, N>::operator~() const
    {
        T buffer[N];
        for (std::size_t i = 0; i < N; ++i)
        {
            buffer[i] = ~m_value[i];
        }
        return batch<T, N>(buffer);
    }

    template <class T, std::size_t N>
    inline batch<T, N>& batch<T, N>::operator+=(const batch<T, N>& other)
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            m_value[i] += other.m_value[i];
        }
        return *this;
    }

    template <class T, std::size_t N>
    inline batch<T, N>& batch<T, N>::operator-=(const batch<T, N>& other)
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            m_value[i] -= other.m_value[i];
        }
        return *this;
    }

    template <class T, std::size_t N>
    inline batch<T, N>& batch<T, N>::operator*=(const batch<T, N>& other)
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            m_value[i] *= other.m_value[i];
        }
        return *this;
    }

    template <class T, std::size_t N>
    inline batch<T, N>& batch<T, N>::operator/=(const batch<T, N>& other)
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            m_value[i] /= other.m_value[i];
        }
        return *this;
    }

    template <class T, std::size_t N>
    inline batch<T, N>& batch<T, N>::operator&=(const batch<T, N>& other)
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            m_value[i] &= other.m_value[i];
        }
        return *this;
    }

    template <class T, std::size_t N>
    inline batch<T, N>& batch<T, N>::operator|=(const batch<T, N>& other)
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            m_value[i] |= other.m_value[i];
        }
        return *this;
    }

    template <class T, std::size_t N>
    inline batch<T, N>& batch<T, N>::operator^=(const batch<T, N>& other)
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            m_value[i] ^= other.m_value[i];
        }
        return *this;
    }

    /*****************************
     * batch_bool implementation *
     *****************************/

    template <class T, std::size_t N>
    inline batch_bool<T, N>::batch_bool(bool val)
        : m_value(val ? -1 : 0)
    {
    }

    template <class T, std::size_t N>
    template <class... Args>
    inline batch_bool<T, N>::batch_bool(Args... args)
    {
        static_assert(sizeof...(Args) == N, "consistent initialization");
        bool buffer[N] = {args...};
        for (std::size_t i = 0; i < N; ++i)
        {
            m_value[i] = buffer[i] ? -1 : 0;
        }
    }

    template <class T, std::size_t N>
    inline batch_bool<T, N>::batch_bool(const bool* buffer)
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            m_value[i] = buffer[i] ? -1 : 0;
        }
    }

    template <class T, std::size_t N>
    inline void batch_bool<T, N>::store_unaligned(bool* buffer) const
    {
        for (std::size_t i = 0; i < N; ++i)
        {
            buffer[i] = m_value[i] != 0;
        }
    }

    template <class T, std::size_t N>
    inline bool batch_bool<T, N>::operator[](std::size_t index) const
    {
        return m_value[index] != 0;
    }

    template <class T, std::size_t N>
    inline batch_bool<T, N> batch_bool<T, N>::operator&(const batch_bool<T, N>& other) const
    {
        batch_bool<T, N> res;
        res.m_value = m_value & other.m_value;
        return res;
    }

    template <class T, std::size_t N>
    inline batch_bool<T, N> batch_bool<T, N>::operator|(const batch_bool<T, N>& other) const
    {
        batch_bool<T, N> res;
        res.m_value = m_value | other.m_value;
        return res;
    }

    template <class T, std::size_t N>
    inline batch_bool<T, N> batch_bool<T, N>::operator^(const batch_bool<T, N>& other) const
    {
        batch_bool<T, N> res;
        res.m_value = m_value ^ other.m_value;
        return res;
    }

    template <class T, std::size_t N>
    inline batch_bool<T, N> batch_bool<T, N>::operator~() const
    {
        batch_bool<T, N> res;
        res.m_value = ~m_value;
        return res;
    }

    template <class T, std::size_t N>
    inline batch_bool<T, N> batch_bool<T, N>::operator==(const batch_bool<T, N>& other) const
    {
        batch_bool<T, N> res;
        res.m_value = ~(m_value ^ other.m_value);
        return res;
    }

    template <class T, std::size_t N>
    inline batch_bool<T, N> batch_bool<T, N>::operator!=(const batch_bool<T, N>& other) const
    {
        batch_bool<T, N> res;
        res.m_value = m_value ^ other.m_value;
        return res;
    }

    template <class T, std::size_t N>
    inline batch_bool<T, N>& batch_bool<T, N>::operator&=(const batch_bool<T, N>& other)
    {
        m_value &= other.m_value;
        return *this;
    }

    template <class T, std::size_t N>
    inline batch_bool<T, N>& batch_bool<T, N>::operator|=(const batch_bool<T, N>& other)
    {
        m_value |= other.m_value;
        return *this;
    }

    template <class T, std::size_t N>
    inline batch_bool<T, N>& batch_bool<T, N>::operator^=(const batch_bool<T, N>& other)
    {
        m_value ^= other.m_value;
        return *this;
    }

    /*********************************
     * batch<complex> implementation *
     *********************************/

    template <class T, std::size_t N>
    inline batch<std::complex<T>, N>::batch(const std::complex<T>& val)
        : m_real(val.real())
        , m_imag(val.imag())
    {
    }

    template <class T, std::size_t N>
    inline batch<std::complex<T>, N>::batch(const real_batch& real, const real_batch& imag)
        : m_real(real)
        , m_imag(imag)
    {
    }

    template <class T, std::size_t N>
    inline batch<std::complex<T>, N>::batch(const std::complex<T>* buffer)
    {
        T real_buffer[N];
        T imag_buffer[N];
        for (std::size_t i = 0; i < N; ++i)
        {
            real_buffer[i] = buffer[i].real();
            imag_buffer[i] = buffer[i].imag();
        }
        m_real = real_batch(real_buffer);
        m_imag = imag_batch(imag_buffer);
    }

    template <class T, std::size_t N>
    inline void batch<std::complex<T>, N>::store_unaligned(std::complex<T>* buffer) const
    {
        T real_buffer[N];
        T imag_buffer[N];
        m_real.store_unaligned(real_buffer);
        m_imag.store_unaligned(imag_buffer);
        for (std::size_t i = 0; i < N; ++i)
        {
            buffer[i] = std::complex<T>(real_buffer[i], imag_buffer[i]);
        }
    }

    template <class T, std::size_t N>
    inline void batch<std::complex<T>, N>::store_aligned(std::complex<T>* buffer) const
    {
        store_unaligned(buffer);
    }

    template <class T, std::size_t N>
    inline std::complex<T> batch<std::complex<T>, N>::operator[](std::size_t index) const
    {
        return std::complex<T>(m_real[index], m_imag[index]);
    }

    template <class T, std::size_t N>
    inline batch_bool<T, N> batch<std::complex<T>, N>::operator==(const batch<std::complex<T>, N>& other) const
    {
        return m_real == other.m_real && m_imag == other.m_imag;
    }

    template <class T, std::size_t N>
    inline batch_bool<T, N> batch<std::complex<T>, N>::operator!=(const batch<std::complex<T>, N>& other) const
    {
        return m_real != other.m_real || m_imag != other.m_imag;
    }

    template <class T, std::size_t N>
    inline batch<std::complex<T>, N> batch<std::complex<T>, N>::operator-() const
    {
        return batch<std::complex<T>, N>(-m_real, -m_imag);
    }

    template <class T, std::size_t N>
    inline batch<T, N> batch<std::complex<T>, N>::operator+() const
    {
        return *this;
    }

    template <class T, std::size_t N>
    inline batch<std::complex<T>, N> batch<std::complex<T>, N>::operator+(const batch<std::complex<T>, N>& other) const
    {
        return batch<std::complex<T>, N>(m_real + other.m_real, m_imag + other.m_imag);
    }

    template <class T, std::size_t N>
    inline batch<std::complex<T>, N> batch<std::complex<T>, N>::operator-(const batch<std::complex<T>, N>& other) const
    {
        return batch<std::complex<T>, N>(m_real - other.m_real, m_imag - other.m_imag);
    }

    template <class T, std::size_t N>
    inline batch<std::complex<T>, N> batch<std::complex<T>, N>::operator*(const batch<std::complex<T>, N>& other) const
    {
        return batch<std::complex<T>, N>(m_real * other.m_real - m_imag * other.m_imag,
                                         m_real * other.m_imag + m_imag * other.m_real);
    }

    template <class T, std::size_t N>
    inline batch<std::complex<T>, N> batch<std::complex<T>, N>::operator/(const batch<std::complex<T>, N>& other) const
    {
        real_batch tmp = other.m_real * other.m_real + other.m_imag * other.m_imag;
        return batch<std::complex<T>, N>((m_real * other.m_real + m_imag * other.m_imag) / tmp,
                                         (m_imag * other.m_real - m_real * other.m_imag) / tmp);
    }

    template <class T, std::size_t N>
    inline batch<std::complex<T>, N>& batch<std::complex<T>, N>::operator+=(const batch<std::complex<T>, N>& other)
    {
        m_real += other.m_real;
        m_imag += other.m_imag;
        return *this;
    }

    template <class T, std::size_t N>
    inline batch<std::complex<T>, N>& batch<std::complex<T>, N>::operator-=(const batch<std::complex<T>, N>& other)
    {
        m_real -= other.m_real;
        m_imag -= other.m_imag;
        return *this;
    }

    template <class T, std::size_t N>
    inline batch<std::complex<T>, N>& batch<std::complex<T>, N>::operator*=(const batch<std::complex<T>, N>& other)
    {
        real_batch tmp_real = m_real * other.m_real - m_imag * other.m_imag;
        m_imag = m_real * other.m_imag + m_imag * other.m_real;
        m_real = tmp_real;
        return *this;
    }

    template <class T, std::size_t N>
    inline batch<std::complex<T>, N>& batch<std::complex<T>, N>::operator/=(const batch<std::complex<T>, N>& other)
    {
        real_batch tmp = other.m_real * other.m_real + other.m_imag * other.m_imag;
        real_batch tmp_real = (m_real * other.m_real + m_imag * other.m_imag) / tmp;
        m_imag = (m_imag * other.m_real - m_real * other.m_imag) / tmp;
        m_real = tmp_real;
        return *this;
    }

    template <class T, std::size_t N>
    inline batch<T, N> batch<std::complex<T>, N>::real() const
    {
        return m_real;
    }

    template <class T, std::size_t N>
    inline batch<T, N> batch<std::complex<T>, N>::imag() const
    {
        return m_imag;
    }

#if XSIMD_WITH_XTL
    /*********************************
     * batch<xcomplex> implementation *
     *********************************/

    template <class T, std::size_t N>
    inline batch<xtl::xcomplex<T, T>, N>::batch(const xtl::xcomplex<T, T>& val)
        : m_real(val.real())
        , m_imag(val.imag())
    {
    }

    template <class T, std::size_t N>
    inline batch<xtl::xcomplex<T, T>, N>::batch(const real_batch& real, const real_batch& imag)
        : m_real(real)
        , m_imag(imag)
    {
    }

    template <class T, std::size_t N>
    inline batch<xtl::xcomplex<T, T>, N>::batch(const xtl::xcomplex<T, T>* buffer)
    {
        T real_buffer[N];
        T imag_buffer[N];
        for (std::size_t i = 0; i < N; ++i)
        {
            real_buffer[i] = buffer[i].real();
            imag_buffer[i] = buffer[i].imag();
        }
        m_real = real_batch(real_buffer);
        m_imag = imag_batch(imag_buffer);
    }

    template <class T, std::size_t N>
    inline void batch<xtl::xcomplex<T, T>, N>::store_unaligned(xtl::xcomplex<T, T>* buffer) const
    {
        T real_buffer[N];
        T imag_buffer[N];
        m_real.store_unaligned(real_buffer);
        m_imag.store_unaligned(imag_buffer);
        for (std::size_t i = 0; i < N; ++i)
        {
            buffer[i] = xtl::xcomplex<T, T>(real_buffer[i], imag_buffer[i]);
        }
    }

    template <class T, std::size_t N>
    inline void batch<xtl::xcomplex<T, T>, N>::store_aligned(xtl::xcomplex<T, T>* buffer) const
    {
        store_unaligned(buffer);
    }

    template <class T, std::size_t N>
    inline xtl::xcomplex<T, T> batch<xtl::xcomplex<T, T>, N>::operator[](std::size_t index) const
    {
        return xtl::xcomplex<T, T>(m_real[index], m_imag[index]);
    }

    template <class T, std::size_t N>
    inline batch_bool<T, N> batch<xtl::xcomplex<T, T>, N>::operator==(const batch<xtl::xcomplex<T, T>, N>& other) const
    {
        return m_real == other.m_real && m_imag == other.m_imag;
    }

    template <class T, std::size_t N>
    inline batch_bool<T, N> batch<xtl::xcomplex<T, T>, N>::operator!=(const batch<xtl::xcomplex<T, T>, N>& other) const
    {
        return m_real != other.m_real || m_imag != other.m_imag;
    }

    template <class T, std::size_t N>
    inline batch<xtl::xcomplex<T, T>, N> batch<xtl::xcomplex<T, T>, N>::operator-() const
    {
        return batch<xtl::xcomplex<T, T>, N>(-m_real, -m_imag);
    }

    template <class T, std::size_t N>
    inline batch<T, N> batch<xtl::xcomplex<T, T>, N>::operator+() const
    {
        return *this;
    }

    template <class T, std::size_t N>
    inline batch<xtl::xcomplex<T, T>, N> batch<xtl::xcomplex<T, T>, N>::operator+(const batch<xtl::xcomplex<T, T>, N>& other) const
    {
        return batch<xtl::xcomplex<T, T>, N>(m_real + other.m_real, m_imag + other.m_imag);
    }

    template <class T, std::size_t N>
    inline batch<xtl::xcomplex<T, T>, N> batch<xtl::xcomplex<T, T>, N>::operator-(const batch<xtl::xcomplex<T, T>, N>& other) const
    {
        return batch<xtl::xcomplex<T, T>, N>(m_real - other.m_real, m_imag - other.m_imag);
    }

    template <class T, std::size_t N>
    inline batch<xtl::xcomplex<T, T>, N> batch<xtl::xcomplex<T, T>, N>::operator*(const batch<xtl::xcomplex<T, T>, N>& other) const
    {
        return batch<xtl::xcomplex<T, T>, N>(m_real * other.m_real - m_imag * other.m_imag,
                                             m_real * other.m_imag + m_imag * other.m_real);
    }

    template <class T, std::size_t N>
    inline batch<xtl::xcomplex<T, T>, N> batch<xtl::xcomplex<T, T>, N>::operator/(const batch<xtl::xcomplex<T, T>, N>& other) const
    {
        real_batch tmp = other.m_real * other.m_real + other.m_imag * other.m_imag;
        return batch<xtl::xcomplex<T, T>, N>((m_real * other.m_real + m_imag * other.m_imag) / tmp,
                                             (m_imag * other.m_real - m_real * other.m_imag) / tmp);
    }

    template <class T, std::size_t N>
    inline batch<xtl::xcomplex<T, T>, N>& batch<xtl::xcomplex<T, T>, N>::operator+=(const batch<xtl::xcomplex<T, T>, N>& other)
    {
        m_real += other.m_real;
        m_imag += other.m_imag;
        return *this;
    }

    template <class T, std::size_t N>
    inline batch<xtl::xcomplex<T, T>, N>& batch<xtl::xcomplex<T, T>, N>::operator-=(const batch<xtl::xcomplex<T, T>, N>& other)
    {
        m_real -= other.m_real;
        m_imag -= other.m_imag;
        return *this;
    }

    template <class T, std::size_t N>
    inline batch<xtl::xcomplex<T, T>, N>& batch<xtl::xcomplex<T, T>, N>::operator*=(const batch<xtl::xcomplex<T, T>, N>& other)
    {
        real_batch tmp_real = m_real * other.m_real - m_imag * other.m_imag;
        m_imag = m_real * other.m_imag + m_imag * other.m_real;
        m_real = tmp_real;
        return *this;
    }

    template <class T, std::size_t N>
    inline batch<xtl::xcomplex<T, T>, N>& batch<xtl::xcomplex<T, T>, N>::operator/=(const batch<xtl::xcomplex<T, T>, N>& other)
    {
        real_batch tmp = other.m_real * other.m_real + other.m_imag * other.m_imag;
        real_batch tmp_real = (m_real * other.m_real + m_imag * other.m_imag) / tmp;
        m_imag = (m_imag * other.m_real - m_real * other.m_imag) / tmp;
        m_real = tmp_real;
        return *this;
    }

    template <class T, std::size_t N>
    inline batch<T, N> batch<xtl::xcomplex<T, T>, N>::real() const
    {
        return m_real;
    }

    template <class T, std::size_t N>
    inline batch<T, N> batch<xtl::xcomplex<T, T>, N>::imag() const
    {
        return m_imag;
    }
#endif

    /************************************
     * backends                         *
     ************************************/

    namespace arch
    {
        struct sse;
        struct avx;
        struct avx512;
        struct fma3;
        struct fma4;
    }

    template <class T, class Arch, std::size_t N>
    struct batch_constant;

    template <class T, class Arch>
    class batch;

    template <class T, class Arch>
    class batch_bool;

}

#include "xsimd_sse_conversion.hpp"
#include "xsimd_sse_double.hpp"
#include "xsimd_sse_float.hpp"
#include "xsimd_sse_int.hpp"

// #include "types/xsimd_avx_conversion.hpp"
// #include "types/xsimd_avx_double.hpp"
// #include "types/xsimd_avx_float.hpp"
// #include "types/xsimd_avx_int.hpp"

// #include "types/xsimd_avx512_conversion.hpp"
// #include "types/xsimd_avx512_double.hpp"
// #include "types/xsimd_avx512_float.hpp"
// #include "types/xsimd_avx512_int.hpp"

// #include "types/xsimd_fma_sse.hpp"
// #include "types/xsimd_fma_avx.hpp"

// #include "types/xsimd_complex_batch.hpp"
// #include "types/xsimd_complex_batch_bool.hpp"

// #include "types/xsimd_dispatch.hpp"

#endif
