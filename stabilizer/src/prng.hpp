/*
 * 32-bit PRNG utilities to equivocate random measurements in stab_cpu and stab_cuda.
 * Uses only 32-bit operations for deterministic CPU/GPU parity.
 */

#pragma once
#include <stdint.h>

#if !defined(__CUDACC__)
#ifndef __host__
#define __host__
#endif
#ifndef __device__
#define __device__
#endif
#ifndef __forceinline__
#define __forceinline__ inline
#endif
#endif

// Murmur3 fmix32: good 32-bit avalanche/mixing, 32-bit ops only.
static __host__ __device__ __forceinline__ uint32_t fmix32(uint32_t x) {
    x ^= x >> 16;
    x *= 0x85ebca6bu;
    x ^= x >> 13;
    x *= 0xc2b2ae35u;
    x ^= x >> 16;
    return x;
}

// Stateless 32-bit generator: mix(seed + phi * (ctr+1))
static __host__ __device__ __forceinline__ int32_t prng_32(int32_t seed, int32_t ctr) {
    // Work in unsigned for defined wraparound; keep API int32_t-based.
    uint32_t s = static_cast<uint32_t>(seed);
    uint32_t c = static_cast<uint32_t>(ctr);
    // 0x9E3779B9 is 2^32 * golden ratio fraction; good Weyl increment.
    uint32_t z = s + 0x9e3779b9u * (c + 1u);
    uint32_t m = fmix32(z);
    return static_cast<int32_t>(m);
}

static __host__ __device__ __forceinline__ int prng_bit(int32_t seed, int32_t ctr) {
    // Use LSB; mask after casting to unsigned for portability.
    uint32_t v = static_cast<uint32_t>(prng_32(seed, ctr));
    return static_cast<int>(v & 1u);
}

static __host__ __device__ __forceinline__ double prng_uniform01(int32_t seed, int32_t ctr) {
    // Map to [0,1) with 32-bit resolution. Convert via uint32_t to avoid negatives.
    uint32_t v = static_cast<uint32_t>(prng_32(seed, ctr));
    return (static_cast<double>(v) + 0.5) * (1.0 / 4294967296.0); // 1/2^32
}