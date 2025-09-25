#pragma once

#include <cassert>
#include <cooperative_groups.h>
#include <nvshmem.h>
#include <nvshmemx.h>
#include <cmath>
#include <cstdio>

#include "nwqsim_core.hpp"
#include "../memory/state_view.hpp"

namespace NWQSim
{
    namespace GateKernels
    {
        namespace MultiGPU
        {

            template <typename Group>
            __device__ __forceinline__ void nvshmem_grid_barrier(Group grid)
            {
                if (threadIdx.x == 0 && blockIdx.x == 0)
                    nvshmem_barrier_all();
                grid.sync();
            }

            __device__ __forceinline__ ValType symmetric_get(ValType *arr, IdxType idx, const StateView &state)
            {
                return nvshmem_double_g(arr + ((idx) & (state.m_gpu - 1)), idx >> state.lg2_m_gpu);
            }

            __device__ __forceinline__ void symmetric_put(ValType *arr, IdxType idx, const StateView &state, ValType value)
            {
                nvshmem_double_p(arr + ((idx) & (state.m_gpu - 1)), value, idx >> state.lg2_m_gpu);
            }

            __device__ __forceinline__ IdxType pow2i(IdxType power)
            {
                return static_cast<IdxType>(IdxType(1) << power);
            }

            __device__ __forceinline__ IdxType wrap_local(IdxType idx, IdxType mask)
            {
                return idx & mask;
            }

            __device__ __forceinline__ void apply_c1_gate_local(const ValType *gm_real,
                                                                const ValType *gm_imag,
                                                                IdxType qubit,
                                                                const StateView &state)
            {
                namespace cg = cooperative_groups;
                cg::grid_group grid = cg::this_grid();
                const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
                const IdxType mask = state.m_gpu - 1;
                const IdxType per_work = state.half_dim >> state.gpu_scale;

                for (IdxType i = state.rank * per_work + tid; i < (state.rank + 1) * per_work; i += blockDim.x * gridDim.x)
                {
                    IdxType outer = (i >> qubit);
                    IdxType inner = (i & (pow2i(qubit) - 1));
                    IdxType offset = (outer << (qubit + 1));
                    IdxType pos0 = wrap_local(offset + inner, mask);
                    IdxType pos1 = wrap_local(offset + inner + pow2i(qubit), mask);

                    const ValType el0_real = state.data_real[pos0];
                    const ValType el0_imag = state.data_imag[pos0];
                    const ValType el1_real = state.data_real[pos1];
                    const ValType el1_imag = state.data_imag[pos1];

                    ValType real0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag);
                    ValType imag0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real);
                    ValType real1 = (gm_real[2] * el0_real) - (gm_imag[2] * el0_imag) + (gm_real[3] * el1_real) - (gm_imag[3] * el1_imag);
                    ValType imag1 = (gm_real[2] * el0_imag) + (gm_imag[2] * el0_real) + (gm_real[3] * el1_imag) + (gm_imag[3] * el1_real);

                    state.data_real[pos0] = real0;
                    state.data_imag[pos0] = imag0;
                    state.data_real[pos1] = real1;
                    state.data_imag[pos1] = imag1;
                }
                grid.sync();
            }

            __device__ __forceinline__ void apply_c1_gate_v1(const ValType *gm_real,
                                                             const ValType *gm_imag,
                                                             IdxType qubit,
                                                             const StateView &state)
            {
                namespace cg = cooperative_groups;
                cg::grid_group grid = cg::this_grid();
                const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
                const IdxType mask = state.m_gpu - 1;
                const IdxType per_work = state.half_dim >> state.gpu_scale;

                if (qubit < state.lg2_m_gpu)
                {
                    apply_c1_gate_local(gm_real, gm_imag, qubit, state);
                    return;
                }

                IdxType partner = state.rank ^ pow2i(qubit - state.lg2_m_gpu);
                IdxType role = (state.rank < partner) ? 0 : 1;

                ValType *primary_real;
                ValType *primary_imag;
                ValType *secondary_real;
                ValType *secondary_imag;

                if (role == 0)
                {
                    primary_real = state.data_real;
                    primary_imag = state.data_imag;
                    secondary_real = state.buffer_real;
                    secondary_imag = state.buffer_imag;
                    if (tid == 0)
                    {
                        nvshmem_double_get(secondary_real, state.data_real, per_work, partner);
                        nvshmem_double_get(secondary_imag, state.data_imag, per_work, partner);
                    }
                }
                else
                {
                    primary_real = state.buffer_real;
                    primary_imag = state.buffer_imag;
                    secondary_real = state.data_real + per_work;
                    secondary_imag = state.data_imag + per_work;
                    if (tid == 0)
                    {
                        nvshmem_double_get(primary_real, state.data_real + per_work, per_work, partner);
                        nvshmem_double_get(primary_imag, state.data_imag + per_work, per_work, partner);
                    }
                }
                grid.sync();

                for (IdxType i = tid; i < per_work; i += blockDim.x * gridDim.x)
                {
                    const ValType el0_real = primary_real[i & mask];
                    const ValType el0_imag = primary_imag[i & mask];
                    const ValType el1_real = secondary_real[i & mask];
                    const ValType el1_imag = secondary_imag[i & mask];

                    ValType real0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag);
                    ValType imag0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real);
                    ValType real1 = (gm_real[2] * el0_real) - (gm_imag[2] * el0_imag) + (gm_real[3] * el1_real) - (gm_imag[3] * el1_imag);
                    ValType imag1 = (gm_real[2] * el0_imag) + (gm_imag[2] * el0_real) + (gm_real[3] * el1_imag) + (gm_imag[3] * el1_real);

                    primary_real[i & mask] = real0;
                    primary_imag[i & mask] = imag0;
                    secondary_real[i & mask] = real1;
                    secondary_imag[i & mask] = imag1;
                }
                grid.sync();

                if (tid == 0)
                {
                    if (role == 0)
                    {
                        nvshmem_double_put(state.data_real, secondary_real, per_work, partner);
                        nvshmem_double_put(state.data_imag, secondary_imag, per_work, partner);
                    }
                    else
                    {
                        nvshmem_double_put(state.data_real + per_work, primary_real, per_work, partner);
                        nvshmem_double_put(state.data_imag + per_work, primary_imag, per_work, partner);
                    }
                }
            }

            __device__ __forceinline__ void apply_c1_gate_v2(const ValType *gm_real,
                                                             const ValType *gm_imag,
                                                             IdxType qubit,
                                                             const StateView &state)
            {
                namespace cg = cooperative_groups;
                cg::grid_group grid = cg::this_grid();
                const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
                const IdxType mask = state.m_gpu - 1;
                const IdxType per_work = state.dim >> state.gpu_scale;

                if (qubit < state.lg2_m_gpu)
                {
                    apply_c1_gate_local(gm_real, gm_imag, qubit, state);
                    return;
                }

                IdxType partner = state.rank ^ pow2i(qubit - state.lg2_m_gpu);
                if (state.rank > partner)
                    return; // idle half for unidirectional option

                if (tid == 0)
                {
                    nvshmem_double_get(state.buffer_real, state.data_real, per_work, partner);
                    nvshmem_double_get(state.buffer_imag, state.data_imag, per_work, partner);
                }
                grid.sync();

                for (IdxType i = tid; i < per_work; i += blockDim.x * gridDim.x)
                {
                    const ValType el0_real = state.data_real[i & mask];
                    const ValType el0_imag = state.data_imag[i & mask];
                    const ValType el1_real = state.buffer_real[i & mask];
                    const ValType el1_imag = state.buffer_imag[i & mask];

                    ValType real0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag);
                    ValType imag0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real);
                    ValType real1 = (gm_real[2] * el0_real) - (gm_imag[2] * el0_imag) + (gm_real[3] * el1_real) - (gm_imag[3] * el1_imag);
                    ValType imag1 = (gm_real[2] * el0_imag) + (gm_imag[2] * el0_real) + (gm_real[3] * el1_imag) + (gm_imag[3] * el1_real);

                    state.data_real[i & mask] = real0;
                    state.data_imag[i & mask] = imag0;
                    state.buffer_real[i & mask] = real1;
                    state.buffer_imag[i & mask] = imag1;
                }
                grid.sync();

                if (tid == 0)
                {
                    nvshmem_double_put(state.data_real, state.buffer_real, per_work, partner);
                    nvshmem_double_put(state.data_imag, state.buffer_imag, per_work, partner);
                }
            }

            __device__ __forceinline__ void apply_c2_gate_local(const ValType *gm_real,
                                                                const ValType *gm_imag,
                                                                IdxType qubit0,
                                                                IdxType qubit1,
                                                                const StateView &state)
            {
                namespace cg = cooperative_groups;
                cg::grid_group grid = cg::this_grid();
                const int tid = blockDim.x * blockIdx.x + threadIdx.x;
                const IdxType mask = state.m_gpu - 1;
                const IdxType per_work = state.dim >> (state.gpu_scale + 2);

                const IdxType q0dim = pow2i(max(qubit0, qubit1));
                const IdxType q1dim = pow2i(min(qubit0, qubit1));
                const IdxType mid_factor = (q0dim + q1dim + q1dim - 1) >> (min(qubit0, qubit1) + 1);
                const IdxType inner_factor = q1dim;
                const IdxType qubit0_dim = pow2i(qubit0);
                const IdxType qubit1_dim = pow2i(qubit1);

                for (IdxType i = state.rank * per_work + tid; i < (state.rank + 1) * per_work; i += blockDim.x * gridDim.x)
                {
                    IdxType outer = ((i / inner_factor) / (mid_factor)) * (q0dim + q0dim);
                    IdxType mid = ((i / inner_factor) % (mid_factor)) * (q1dim + q1dim);
                    IdxType inner = i % inner_factor;
                    IdxType pos0 = wrap_local(outer + mid + inner, mask);
                    IdxType pos1 = wrap_local(outer + mid + inner + qubit1_dim, mask);
                    IdxType pos2 = wrap_local(outer + mid + inner + qubit0_dim, mask);
                    IdxType pos3 = wrap_local(outer + mid + inner + q0dim + q1dim, mask);

                    const ValType el0_real = state.data_real[pos0];
                    const ValType el0_imag = state.data_imag[pos0];
                    const ValType el1_real = state.data_real[pos1];
                    const ValType el1_imag = state.data_imag[pos1];
                    const ValType el2_real = state.data_real[pos2];
                    const ValType el2_imag = state.data_imag[pos2];
                    const ValType el3_real = state.data_real[pos3];
                    const ValType el3_imag = state.data_imag[pos3];

                    ValType data_real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag) + (gm_real[2] * el2_real) - (gm_imag[2] * el2_imag) + (gm_real[3] * el3_real) - (gm_imag[3] * el3_imag);
                    ValType data_real_pos1 = (gm_real[4] * el0_real) - (gm_imag[4] * el0_imag) + (gm_real[5] * el1_real) - (gm_imag[5] * el1_imag) + (gm_real[6] * el2_real) - (gm_imag[6] * el2_imag) + (gm_real[7] * el3_real) - (gm_imag[7] * el3_imag);
                    ValType data_real_pos2 = (gm_real[8] * el0_real) - (gm_imag[8] * el0_imag) + (gm_real[9] * el1_real) - (gm_imag[9] * el1_imag) + (gm_real[10] * el2_real) - (gm_imag[10] * el2_imag) + (gm_real[11] * el3_real) - (gm_imag[11] * el3_imag);
                    ValType data_real_pos3 = (gm_real[12] * el0_real) - (gm_imag[12] * el0_imag) + (gm_real[13] * el1_real) - (gm_imag[13] * el1_imag) + (gm_real[14] * el2_real) - (gm_imag[14] * el2_imag) + (gm_real[15] * el3_real) - (gm_imag[15] * el3_imag);

                    ValType data_imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real) + (gm_real[2] * el2_imag) + (gm_imag[2] * el2_real) + (gm_real[3] * el3_imag) + (gm_imag[3] * el3_real);
                    ValType data_imag_pos1 = (gm_real[4] * el0_imag) + (gm_imag[4] * el0_real) + (gm_real[5] * el1_imag) + (gm_imag[5] * el1_real) + (gm_real[6] * el2_imag) + (gm_imag[6] * el2_real) + (gm_real[7] * el3_imag) + (gm_imag[7] * el3_real);
                    ValType data_imag_pos2 = (gm_real[8] * el0_imag) + (gm_imag[8] * el0_real) + (gm_real[9] * el1_imag) + (gm_imag[9] * el1_real) + (gm_real[10] * el2_imag) + (gm_imag[10] * el2_real) + (gm_real[11] * el3_imag) + (gm_imag[11] * el3_real);
                    ValType data_imag_pos3 = (gm_real[12] * el0_imag) + (gm_imag[12] * el0_real) + (gm_real[13] * el1_imag) + (gm_imag[13] * el1_real) + (gm_real[14] * el2_imag) + (gm_imag[14] * el2_real) + (gm_real[15] * el3_imag) + (gm_imag[15] * el3_real);

                    state.data_real[pos0] = data_real_pos0;
                    state.data_real[pos1] = data_real_pos1;
                    state.data_real[pos2] = data_real_pos2;
                    state.data_real[pos3] = data_real_pos3;

                    state.data_imag[pos0] = data_imag_pos0;
                    state.data_imag[pos1] = data_imag_pos1;
                    state.data_imag[pos2] = data_imag_pos2;
                    state.data_imag[pos3] = data_imag_pos3;
                }
                grid.sync();
            }

            __device__ __forceinline__ void apply_c2_gate_v1(const ValType *gm_real,
                                                             const ValType *gm_imag,
                                                             IdxType qubit0,
                                                             IdxType qubit1,
                                                             const StateView &state)
            {
                namespace cg = cooperative_groups;
                cg::grid_group grid = cg::this_grid();
                const int tid = blockDim.x * blockIdx.x + threadIdx.x;
                const IdxType mask = state.m_gpu - 1;

                if (qubit0 < state.lg2_m_gpu && qubit1 < state.lg2_m_gpu)
                {
                    apply_c2_gate_local(gm_real, gm_imag, qubit0, qubit1, state);
                    return;
                }

                const IdxType per_work = state.dim >> (state.gpu_scale + 1);
                const IdxType per_num = state.dim >> state.gpu_scale;
                const IdxType p = min(qubit0, qubit1);
                const IdxType q = max(qubit0, qubit1);

                IdxType partner = state.rank ^ pow2i(q - state.lg2_m_gpu);
                if (state.rank > partner)
                    return;

                if (tid == 0)
                {
                    nvshmem_double_get(state.buffer_real, state.data_real, per_num, partner);
                    nvshmem_double_get(state.buffer_imag, state.data_imag, per_num, partner);
                }
                grid.sync();

                IdxType index = (state.rank >> (q - state.lg2_m_gpu + 1)) << (q - state.lg2_m_gpu);
                index |= state.rank & (pow2i(q - state.lg2_m_gpu) - 1);

                for (IdxType i = index * per_work + tid; i < (index + 1) * per_work; i += blockDim.x * gridDim.x)
                {
                    ValType el_real[4];
                    ValType el_imag[4];
                    ValType res_real[4] = {0};
                    ValType res_imag[4] = {0};

                    const IdxType term0 = i & (pow2i(p) - 1);
                    const IdxType term1 = ((i >> p) & (pow2i(q - p - 1) - 1)) * pow2i(p + 1);
                    const IdxType term2 = (i >> p >> (q - p - 1)) * pow2i(q + 1);
                    const IdxType term = term2 + term1 + term0;

                    el_real[0] = state.data_real[wrap_local(term, mask)];
                    el_imag[0] = state.data_imag[wrap_local(term, mask)];
                    el_real[3] = state.buffer_real[wrap_local(term + pow2i(qubit0) + pow2i(qubit1), mask)];
                    el_imag[3] = state.buffer_imag[wrap_local(term + pow2i(qubit0) + pow2i(qubit1), mask)];

                    if (qubit0 == q)
                    {
                        el_real[1] = state.data_real[wrap_local(term + pow2i(qubit1), mask)];
                        el_imag[1] = state.data_imag[wrap_local(term + pow2i(qubit1), mask)];
                        el_real[2] = state.buffer_real[wrap_local(term + pow2i(qubit0), mask)];
                        el_imag[2] = state.buffer_imag[wrap_local(term + pow2i(qubit0), mask)];
                    }
                    else
                    {
                        el_real[1] = state.buffer_real[wrap_local(term + pow2i(qubit1), mask)];
                        el_imag[1] = state.buffer_imag[wrap_local(term + pow2i(qubit1), mask)];
                        el_real[2] = state.data_real[wrap_local(term + pow2i(qubit0), mask)];
                        el_imag[2] = state.data_imag[wrap_local(term + pow2i(qubit0), mask)];
                    }

                    for (unsigned row = 0; row < 4; row++)
                    {
                        for (unsigned col = 0; col < 4; col++)
                        {
                            res_real[row] += (el_real[col] * gm_real[row * 4 + col]) - (el_imag[col] * gm_imag[row * 4 + col]);
                            res_imag[row] += (el_real[col] * gm_imag[row * 4 + col]) + (el_imag[col] * gm_real[row * 4 + col]);
                        }
                    }

                    state.data_real[wrap_local(term, mask)] = res_real[0];
                    state.data_imag[wrap_local(term, mask)] = res_imag[0];
                    state.buffer_real[wrap_local(term + pow2i(qubit0) + pow2i(qubit1), mask)] = res_real[3];
                    state.buffer_imag[wrap_local(term + pow2i(qubit0) + pow2i(qubit1), mask)] = res_imag[3];

                    if (qubit0 == q)
                    {
                        state.data_real[wrap_local(term + pow2i(qubit1), mask)] = res_real[1];
                        state.data_imag[wrap_local(term + pow2i(qubit1), mask)] = res_imag[1];
                        state.buffer_real[wrap_local(term + pow2i(qubit0), mask)] = res_real[2];
                        state.buffer_imag[wrap_local(term + pow2i(qubit0), mask)] = res_imag[2];
                    }
                    else
                    {
                        state.buffer_real[wrap_local(term + pow2i(qubit1), mask)] = res_real[1];
                        state.buffer_imag[wrap_local(term + pow2i(qubit1), mask)] = res_imag[1];
                        state.data_real[wrap_local(term + pow2i(qubit0), mask)] = res_real[2];
                        state.data_imag[wrap_local(term + pow2i(qubit0), mask)] = res_imag[2];
                    }
                }
                grid.sync();

                if (tid == 0)
                {
                    nvshmem_double_put(state.data_real, state.buffer_real, per_num, partner);
                    nvshmem_double_put(state.data_imag, state.buffer_imag, per_num, partner);
                }
            }

            __device__ __forceinline__ void measure_qubit(const StateView &state,
                                                          IdxType n_qubits,
                                                          IdxType qubit,
                                                          ValType *probabilities,
                                                          const ValType *random_values,
                                                          IdxType random_index,
                                                          IdxType *results,
                                                          IdxType result_index)
            {
                namespace cg = cooperative_groups;
                cg::grid_group grid = cg::this_grid();
                const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
                const IdxType per_work = state.dim >> state.gpu_scale;
                const IdxType local_mask = state.m_gpu - 1;
                const IdxType mask = pow2i(qubit);
                ValType rand = random_values[random_index];

                for (IdxType i = tid; i < state.m_gpu; i += blockDim.x * gridDim.x)
                {
                    IdxType global_idx = state.rank * per_work + i;
                    if ((global_idx & mask) == 0)
                        probabilities[i] = 0;
                    else
                        probabilities[i] = state.data_real[i] * state.data_real[i] + state.data_imag[i] * state.data_imag[i];
                }
                nvshmem_grid_barrier(grid);

                if (n_qubits > 0)
                {
                    for (IdxType k = pow2i(n_qubits - 1); k > 0; k >>= 1)
                    {
                        for (IdxType i = tid; i < k; i += blockDim.x * gridDim.x)
                        {
                            IdxType owner = i >> state.lg2_m_gpu;
                            if (owner == state.rank)
                            {
                                IdxType local_idx = i & local_mask;
                                probabilities[local_idx] += symmetric_get(probabilities, i + k, state);
                            }
                        }
                        nvshmem_grid_barrier(grid);
                    }
                }

                if (tid == 0 && state.rank != 0)
                    probabilities[0] = symmetric_get(probabilities, 0, state);
                grid.sync();
                ValType prob_of_one = probabilities[0];
                grid.sync();

                if (rand < prob_of_one)
                {
                    ValType factor = 1. / sqrt(prob_of_one);
                    for (IdxType i = tid; i < state.m_gpu; i += blockDim.x * gridDim.x)
                    {
                        IdxType global_idx = state.rank * per_work + i;
                        if ((global_idx & mask) == 0)
                        {
                            state.data_real[i] = 0;
                            state.data_imag[i] = 0;
                        }
                        else
                        {
                            state.data_real[i] *= factor;
                            state.data_imag[i] *= factor;
                        }
                    }
                }
                else
                {
                    ValType factor = 1. / sqrt(1. - prob_of_one);
                    for (IdxType i = tid; i < state.m_gpu; i += blockDim.x * gridDim.x)
                    {
                        IdxType global_idx = state.rank * per_work + i;
                        if ((global_idx & mask) == 0)
                        {
                            state.data_real[i] *= factor;
                            state.data_imag[i] *= factor;
                        }
                        else
                        {
                            state.data_real[i] = 0;
                            state.data_imag[i] = 0;
                        }
                    }
                }
                if (tid == 0)
                    results[result_index] = (rand <= prob_of_one ? 1 : 0);
                nvshmem_grid_barrier(grid);
            }

            __device__ __forceinline__ void measure_all(const StateView &state,
                                                        IdxType n_qubits,
                                                        IdxType repetition,
                                                        ValType *probabilities,
                                                        IdxType *results,
                                                        IdxType result_index,
                                                        const ValType *random_values,
                                                        IdxType random_index)
            {
                namespace cg = cooperative_groups;
                cg::grid_group grid = cg::this_grid();
                const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
                const IdxType local_mask = state.m_gpu - 1;
                const IdxType n_size = pow2i(n_qubits);
                const IdxType total_gpus = pow2i(state.gpu_scale);

                for (IdxType i = tid; i < state.m_gpu; i += blockDim.x * gridDim.x)
                {
                    probabilities[i] = state.data_real[i] * state.data_real[i] + state.data_imag[i] * state.data_imag[i];
                }
                nvshmem_grid_barrier(grid);

                for (IdxType d = 0; d < n_qubits; d++)
                {
                    IdxType step = pow2i(d + 1);
                    for (IdxType k = tid * step; k < n_size; k += step * blockDim.x * gridDim.x)
                    {
                        IdxType idx = k + step - 1;
                        IdxType owner = idx >> state.lg2_m_gpu;
                        if (owner == state.rank)
                        {
                            IdxType local_idx = idx & local_mask;
                            probabilities[local_idx] += symmetric_get(probabilities, k + pow2i(d) - 1, state);
                        }
                    }
                    nvshmem_grid_barrier(grid);
                }

                if (state.rank == total_gpus - 1 && tid == 0)
                {
                    IdxType local_idx = (n_size - 1) & local_mask;
                    ValType val = probabilities[local_idx];
                    probabilities[local_idx] = 0;
                    ValType purity = fabs(val);
                    if (fabs(purity - 1.0) > MEASUREMENT_ERROR_BAR)
                        std::printf("MA: Purity Check fails with %lf\n", purity);
                }
                nvshmem_grid_barrier(grid);

                if (n_qubits > 0)
                {
                    for (IdxType d = n_qubits - 1;; d--)
                    {
                        IdxType step = pow2i(d + 1);
                        for (IdxType k = tid * step; k < n_size - 1; k += step * blockDim.x * gridDim.x)
                        {
                            IdxType idx = k + step - 1;
                            IdxType owner = idx >> state.lg2_m_gpu;
                            if (owner == state.rank)
                            {
                                IdxType local_idx = idx & local_mask;
                                IdxType remote_idx = k + pow2i(d) - 1;
                                ValType tmp = symmetric_get(probabilities, remote_idx, state);
                                ValType tmp2 = probabilities[local_idx];
                                symmetric_put(probabilities, remote_idx, state, tmp2);
                                probabilities[local_idx] = tmp + tmp2;
                            }
                        }
                        nvshmem_grid_barrier(grid);
                        if (d == 0)
                            break;
                    }
                }

                for (IdxType j = tid; j < n_size; j += blockDim.x * gridDim.x)
                {
                    IdxType owner = j >> state.lg2_m_gpu;
                    if (owner == state.rank)
                    {
                        IdxType local_idx = j & local_mask;
                        ValType lower = probabilities[local_idx];
                        ValType upper = (j + 1 == n_size) ? 1 : symmetric_get(probabilities, j + 1, state);
                        for (IdxType shot = 0; shot < repetition; shot++)
                        {
                            ValType r = random_values[random_index + shot];
                            if (lower <= r && r < upper)
                                nvshmem_longlong_p(&results[result_index + shot], j, 0);
                        }
                    }
                }
                nvshmem_grid_barrier(grid);

                if (state.rank != 0 && tid == 0)
                    nvshmem_longlong_get(results + result_index, results + result_index, repetition, 0);
                nvshmem_grid_barrier(grid);
            }

            __device__ __forceinline__ void reset_qubit(const StateView &state,
                                                        IdxType n_qubits,
                                                        IdxType qubit,
                                                        ValType *probabilities)
            {
                namespace cg = cooperative_groups;
                cg::grid_group grid = cg::this_grid();
                const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
                const IdxType total_threads = blockDim.x * gridDim.x;
                const IdxType local_mask = state.m_gpu - 1;
                const IdxType per_work = state.dim >> state.gpu_scale;
                IdxType mask = pow2i(qubit);

                for (IdxType local = tid; local < state.m_gpu; local += total_threads)
                {
                    IdxType global_idx = state.rank * per_work + local;
                    if ((global_idx & mask) == 0)
                        probabilities[local] = 0;
                    else
                        probabilities[local] = state.data_real[local] * state.data_real[local] +
                                               state.data_imag[local] * state.data_imag[local];
                }
                nvshmem_grid_barrier(grid);

                if (n_qubits > 0)
                {
                    for (IdxType k = pow2i(n_qubits - 1); k > 0; k >>= 1)
                    {
                        for (IdxType i = tid; i < k; i += total_threads)
                        {
                            IdxType owner = i >> state.lg2_m_gpu;
                            if (owner == state.rank)
                            {
                                IdxType local_idx = i & local_mask;
                                probabilities[local_idx] += symmetric_get(probabilities, i + k, state);
                            }
                        }
                        nvshmem_grid_barrier(grid);
                    }
                }

                if (tid == 0 && state.rank != 0)
                    probabilities[0] = symmetric_get(probabilities, 0, state);
                grid.sync();
                ValType prob_of_one = probabilities[0];
                grid.sync();

                if (prob_of_one < 1.0)
                {
                    ValType factor = 1.0 / sqrt(1.0 - prob_of_one);
                    for (IdxType local = tid; local < state.m_gpu; local += total_threads)
                    {
                        IdxType global_idx = state.rank * per_work + local;
                        if ((global_idx & mask) == 0)
                        {
                            state.data_real[local] *= factor;
                            state.data_imag[local] *= factor;
                        }
                        else
                        {
                            state.data_real[local] = 0;
                            state.data_imag[local] = 0;
                        }
                    }
                }
                else
                {
                    if (qubit >= state.lg2_m_gpu)
                    {
                        IdxType partner = state.rank ^ pow2i(qubit - state.lg2_m_gpu);
                        if (tid == 0)
                        {
                            nvshmem_double_get(state.buffer_real, state.data_real, state.m_gpu, partner);
                            nvshmem_double_get(state.buffer_imag, state.data_imag, state.m_gpu, partner);
                        }
                        nvshmem_grid_barrier(grid);
                        for (IdxType local = tid; local < state.m_gpu; local += total_threads)
                        {
                            state.data_real[local] = state.buffer_real[local];
                            state.data_imag[local] = state.buffer_imag[local];
                        }
                    }
                    else
                    {
                        for (IdxType local = tid; local < state.m_gpu; local += total_threads)
                        {
                            IdxType global_idx = state.rank * per_work + local;
                            if ((global_idx & mask) == 0)
                            {
                                IdxType partner_idx = local ^ mask;
                                state.data_real[local] = state.data_real[partner_idx];
                                state.data_imag[local] = state.data_imag[partner_idx];
                                state.data_real[partner_idx] = 0;
                                state.data_imag[partner_idx] = 0;
                            }
                        }
                    }
                }
                nvshmem_grid_barrier(grid);
            }

            __device__ __forceinline__ void expectation_mask(const StateView &state,
                                                             IdxType n_qubits,
                                                             IdxType mask,
                                                             ValType coeff,
                                                             ValType *workspace)
            {
                namespace cg = cooperative_groups;
                cg::grid_group grid = cg::this_grid();
                const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
                const IdxType total_threads = blockDim.x * gridDim.x;
                const IdxType per_work = state.dim >> state.gpu_scale;
                const IdxType local_span = state.m_gpu;

                ValType local_sum = 0;
                for (IdxType local = tid; local < local_span; local += total_threads)
                {
                    IdxType global_idx = state.rank * per_work + local;
                    ValType prob = state.data_real[local] * state.data_real[local] +
                                   state.data_imag[local] * state.data_imag[local];
                    local_sum += z_mask_sign(global_idx, mask) * prob;
                }

                if (tid < local_span)
                    workspace[tid] += coeff * local_sum;

                grid.sync();
            }

            __device__ __forceinline__ void expectation_observable(const StateView &state,
                                                                   IdxType n_qubits,
                                                                   const ZObservableTerm *terms,
                                                                   IdxType term_count,
                                                                   ValType *workspace,
                                                                   ValType *result_ptr)
            {
                namespace cg = cooperative_groups;
                cg::grid_group grid = cg::this_grid();
                const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
                const IdxType local_span = state.m_gpu;

                if (tid < local_span)
                    workspace[tid] = 0;
                nvshmem_grid_barrier(grid);

                for (IdxType t = 0; t < term_count; t++)
                {
                    expectation_mask(state, n_qubits, terms[t].mask, terms[t].coeff, workspace);
                    grid.sync();
                }

                nvshmem_grid_barrier(grid);

                IdxType reduce_limit = 1;
                while (reduce_limit < local_span)
                    reduce_limit <<= 1;

                for (IdxType offset = reduce_limit >> 1; offset > 0; offset >>= 1)
                {
                    if (tid < offset && tid + offset < local_span)
                        workspace[tid] += workspace[tid + offset];
                    grid.sync();
                }

                if (tid == 0)
                    workspace[0] = workspace[0];
                nvshmem_grid_barrier(grid);

                if (state.rank == 0 && tid == 0)
                    *result_ptr = 0;
                nvshmem_grid_barrier(grid);
                if (tid == 0)
                    nvshmem_double_atomic_add(result_ptr, workspace[0], 0);
                nvshmem_grid_barrier(grid);
                ValType global = nvshmem_double_g(result_ptr, 0);
                grid.sync();
                if (tid == 0)
                    *result_ptr = global;
                nvshmem_grid_barrier(grid);
            }

        } // namespace MultiGPU
    } // namespace GateKernels
} // namespace NWQSim
