#pragma once

#include <cassert>
#include <cooperative_groups.h>
#include <cmath>
#include <cstdio>

#include "nwqsim_core.hpp"
#include "../memory/state_view.hpp"

namespace NWQSim
{
    namespace GateKernels
    {

        namespace GPU
        {

            __device__ __forceinline__ IdxType pow2i(IdxType power)
            {
                return static_cast<IdxType>(IdxType(1) << power);
            }

            __device__ __forceinline__ IdxType mod_pow2(IdxType value, IdxType power)
            {
                return value & (pow2i(power) - IdxType(1));
            }

            __device__ __forceinline__ IdxType div_pow2(IdxType value, IdxType power)
            {
                return value >> power;
            }

            __device__ __forceinline__ IdxType state16_offset(unsigned pattern, IdxType q0, IdxType q1, IdxType q2, IdxType q3)
            {
                return ((pattern >> 3) & 1 ? pow2i(q0) : IdxType(0)) +
                       ((pattern >> 2) & 1 ? pow2i(q1) : IdxType(0)) +
                       ((pattern >> 1) & 1 ? pow2i(q2) : IdxType(0)) +
                       ((pattern) & 1 ? pow2i(q3) : IdxType(0));
            }

            __device__ __forceinline__ void apply_c1_gate(const ValType *gm_real,
                                                          const ValType *gm_imag,
                                                          IdxType qubit,
                                                          const StateView &state)
            {
                namespace cg = cooperative_groups;
                cg::grid_group grid = cg::this_grid();
                const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
                ValType *state_real = state.data_real;
                ValType *state_imag = state.data_imag;
                const IdxType half_dim = state.half_dim ? state.half_dim : (state.dim >> 1);
                for (IdxType i = tid; i < half_dim; i += blockDim.x * gridDim.x)
                {
                    IdxType outer = (i >> qubit);
                    IdxType inner = (i & (pow2i(qubit) - 1));
                    IdxType offset = (outer << (qubit + 1));
                    IdxType pos0 = (offset + inner);
                    IdxType pos1 = (offset + inner + pow2i(qubit));

                    const ValType el0_real = state_real[pos0];
                    const ValType el0_imag = state_imag[pos0];
                    const ValType el1_real = state_real[pos1];
                    const ValType el1_imag = state_imag[pos1];

                    ValType state_real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag);
                    ValType state_imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real);
                    ValType state_real_pos1 = (gm_real[2] * el0_real) - (gm_imag[2] * el0_imag) + (gm_real[3] * el1_real) - (gm_imag[3] * el1_imag);
                    ValType state_imag_pos1 = (gm_real[2] * el0_imag) + (gm_imag[2] * el0_real) + (gm_real[3] * el1_imag) + (gm_imag[3] * el1_real);

                    state_real[pos0] = state_real_pos0;
                    state_imag[pos0] = state_imag_pos0;
                    state_real[pos1] = state_real_pos1;
                    state_imag[pos1] = state_imag_pos1;
                }
                grid.sync();
            }

            __device__ __forceinline__ void apply_c2_gate(const ValType *gm_real,
                                                          const ValType *gm_imag,
                                                          IdxType qubit0,
                                                          IdxType qubit1,
                                                          const StateView &state)
            {
                namespace cg = cooperative_groups;
                cg::grid_group grid = cg::this_grid();
                const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
                ValType *state_real = state.data_real;
                ValType *state_imag = state.data_imag;
                const IdxType dim = state.dim ? state.dim : (state.half_dim << 1);
                const IdxType per_work = (dim >> 2);

                const IdxType q0dim = pow2i(max(qubit0, qubit1));
                const IdxType q1dim = pow2i(min(qubit0, qubit1));
                const IdxType mid_factor = (q0dim + q1dim + q1dim - 1) >> (min(qubit0, qubit1) + 1);
                const IdxType inner_factor = q1dim;
                const IdxType qubit0_dim = pow2i(qubit0);
                const IdxType qubit1_dim = pow2i(qubit1);

                for (IdxType i = tid; i < per_work; i += blockDim.x * gridDim.x)
                {
                    IdxType outer = ((i / inner_factor) / (mid_factor)) * (q0dim + q0dim);
                    IdxType mid = ((i / inner_factor) % (mid_factor)) * (q1dim + q1dim);
                    IdxType inner = i % inner_factor;
                    IdxType pos0 = outer + mid + inner;
                    IdxType pos1 = outer + mid + inner + qubit1_dim;
                    IdxType pos2 = outer + mid + inner + qubit0_dim;
                    IdxType pos3 = outer + mid + inner + q0dim + q1dim;

                    const ValType el0_real = state_real[pos0];
                    const ValType el0_imag = state_imag[pos0];
                    const ValType el1_real = state_real[pos1];
                    const ValType el1_imag = state_imag[pos1];
                    const ValType el2_real = state_real[pos2];
                    const ValType el2_imag = state_imag[pos2];
                    const ValType el3_real = state_real[pos3];
                    const ValType el3_imag = state_imag[pos3];

                    ValType state_real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag) + (gm_real[2] * el2_real) - (gm_imag[2] * el2_imag) + (gm_real[3] * el3_real) - (gm_imag[3] * el3_imag);
                    ValType state_real_pos1 = (gm_real[4] * el0_real) - (gm_imag[4] * el0_imag) + (gm_real[5] * el1_real) - (gm_imag[5] * el1_imag) + (gm_real[6] * el2_real) - (gm_imag[6] * el2_imag) + (gm_real[7] * el3_real) - (gm_imag[7] * el3_imag);
                    ValType state_real_pos2 = (gm_real[8] * el0_real) - (gm_imag[8] * el0_imag) + (gm_real[9] * el1_real) - (gm_imag[9] * el1_imag) + (gm_real[10] * el2_real) - (gm_imag[10] * el2_imag) + (gm_real[11] * el3_real) - (gm_imag[11] * el3_imag);
                    ValType state_real_pos3 = (gm_real[12] * el0_real) - (gm_imag[12] * el0_imag) + (gm_real[13] * el1_real) - (gm_imag[13] * el1_imag) + (gm_real[14] * el2_real) - (gm_imag[14] * el2_imag) + (gm_real[15] * el3_real) - (gm_imag[15] * el3_imag);

                    ValType state_imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real) + (gm_real[2] * el2_imag) + (gm_imag[2] * el2_real) + (gm_real[3] * el3_imag) + (gm_imag[3] * el3_real);
                    ValType state_imag_pos1 = (gm_real[4] * el0_imag) + (gm_imag[4] * el0_real) + (gm_real[5] * el1_imag) + (gm_imag[5] * el1_real) + (gm_real[6] * el2_imag) + (gm_imag[6] * el2_real) + (gm_real[7] * el3_imag) + (gm_imag[7] * el3_real);
                    ValType state_imag_pos2 = (gm_real[8] * el0_imag) + (gm_imag[8] * el0_real) + (gm_real[9] * el1_imag) + (gm_imag[9] * el1_real) + (gm_real[10] * el2_imag) + (gm_imag[10] * el2_real) + (gm_real[11] * el3_imag) + (gm_imag[11] * el3_real);
                    ValType state_imag_pos3 = (gm_real[12] * el0_imag) + (gm_imag[12] * el0_real) + (gm_real[13] * el1_imag) + (gm_imag[13] * el1_real) + (gm_real[14] * el2_imag) + (gm_imag[14] * el2_real) + (gm_real[15] * el3_imag) + (gm_imag[15] * el3_real);

                    state_real[pos0] = state_real_pos0;
                    state_real[pos1] = state_real_pos1;
                    state_real[pos2] = state_real_pos2;
                    state_real[pos3] = state_real_pos3;

                    state_imag[pos0] = state_imag_pos0;
                    state_imag[pos1] = state_imag_pos1;
                    state_imag[pos2] = state_imag_pos2;
                    state_imag[pos3] = state_imag_pos3;
                }
                grid.sync();
            }

            __device__ __forceinline__ void apply_c4_gate(const ValType *gm_real,
                                                          const ValType *gm_imag,
                                                          IdxType qubit0,
                                                          IdxType qubit1,
                                                          IdxType qubit2,
                                                          IdxType qubit3,
                                                          const StateView &state)
            {
                namespace cg = cooperative_groups;
                cg::grid_group grid = cg::this_grid();
                const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;

                assert(qubit0 != qubit1);
                assert(qubit0 != qubit2);
                assert(qubit0 != qubit3);
                assert(qubit1 != qubit2);
                assert(qubit1 != qubit3);
                assert(qubit2 != qubit3);

                ValType *state_real = state.data_real;
                ValType *state_imag = state.data_imag;
                const IdxType dim = state.dim ? state.dim : (state.half_dim << 1);
                const IdxType v0 = min(qubit0, qubit1);
                const IdxType v1 = min(qubit2, qubit3);
                const IdxType v2 = max(qubit0, qubit1);
                const IdxType v3 = max(qubit2, qubit3);
                const IdxType p = min(v0, v1);
                const IdxType q = min(min(v2, v3), max(v0, v1));
                const IdxType r = max(min(v2, v3), max(v0, v1));
                const IdxType s = max(v2, v3);

                for (IdxType i = tid; i < (dim >> 4); i += blockDim.x * gridDim.x)
                {
                    const IdxType term0 = mod_pow2(i, p);
                    const IdxType term1 = mod_pow2(div_pow2(i, p), q - p - 1) * pow2i(p + 1);
                    const IdxType term2 = mod_pow2(div_pow2(div_pow2(i, p), q - p - 1), r - q - 1) * pow2i(q + 1);
                    const IdxType term3 = mod_pow2(div_pow2(div_pow2(div_pow2(i, p), q - p - 1), r - q - 1), s - r - 1) * pow2i(r + 1);
                    const IdxType term4 = div_pow2(div_pow2(div_pow2(div_pow2(i, p), q - p - 1), r - q - 1), s - r - 1) * pow2i(s + 1);
                    const IdxType term = term4 + term3 + term2 + term1 + term0;

                    ValType el_real[16];
                    ValType el_imag[16];
                    for (unsigned pattern = 0; pattern < 16; ++pattern)
                    {
                        const IdxType offset = state16_offset(pattern, qubit0, qubit1, qubit2, qubit3);
                        el_real[pattern] = state_real[term + offset];
                        el_imag[pattern] = state_imag[term + offset];
                    }

                    for (unsigned row = 0; row < 16; row++)
                    {
                        ValType res_real = 0;
                        ValType res_imag = 0;
                        const unsigned row_offset = row * 16;
                        for (unsigned col = 0; col < 16; col++)
                        {
                            res_real += (el_real[col] * gm_real[row_offset + col]) - (el_imag[col] * gm_imag[row_offset + col]);
                            res_imag += (el_real[col] * gm_imag[row_offset + col]) + (el_imag[col] * gm_real[row_offset + col]);
                        }
                        const IdxType offset = state16_offset(row, qubit0, qubit1, qubit2, qubit3);
                        state_real[term + offset] = res_real;
                        state_imag[term + offset] = res_imag;
                    }
                }
                grid.sync();
            }

            __device__ __forceinline__ void measure_qubit(const StateView &state,
                                                          IdxType n_qubits,
                                                          IdxType qubit,
                                                          ValType *prob_buffer,
                                                          const ValType *random_values,
                                                          IdxType random_index,
                                                          IdxType *results,
                                                          IdxType result_index)
            {
                namespace cg = cooperative_groups;
                cg::grid_group grid = cg::this_grid();
                const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
                ValType *state_real = state.data_real;
                ValType *state_imag = state.data_imag;
                IdxType dim = state.dim ? state.dim : pow2i(n_qubits);
                IdxType mask = pow2i(qubit);
                ValType rand = random_values[random_index];

                for (IdxType i = tid; i < dim; i += blockDim.x * gridDim.x)
                {
                    if ((i & mask) == 0)
                        prob_buffer[i] = 0;
                    else
                        prob_buffer[i] = state_real[i] * state_real[i] + state_imag[i] * state_imag[i];
                }
                grid.sync();

                if (n_qubits > 0)
                {
                    for (IdxType k = pow2i(n_qubits - 1); k > 0; k >>= 1)
                    {
                        for (IdxType i = tid; i < k; i += blockDim.x * gridDim.x)
                        {
                            prob_buffer[i] += prob_buffer[i + k];
                        }
                        grid.sync();
                    }
                }

                ValType prob_of_one = prob_buffer[0];
                grid.sync();

                if (rand < prob_of_one)
                {
                    ValType factor = 1. / sqrt(prob_of_one);
                    for (IdxType i = tid; i < dim; i += blockDim.x * gridDim.x)
                    {
                        if ((i & mask) == 0)
                        {
                            state_real[i] = 0;
                            state_imag[i] = 0;
                        }
                        else
                        {
                            state_real[i] *= factor;
                            state_imag[i] *= factor;
                        }
                    }
                }
                else
                {
                    ValType factor = 1. / sqrt(1. - prob_of_one);
                    for (IdxType i = tid; i < dim; i += blockDim.x * gridDim.x)
                    {
                        if ((i & mask) == 0)
                        {
                            state_real[i] *= factor;
                            state_imag[i] *= factor;
                        }
                        else
                        {
                            state_real[i] = 0;
                            state_imag[i] = 0;
                        }
                    }
                }
                if (tid == 0)
                    results[result_index] = (rand <= prob_of_one ? 1 : 0);
                grid.sync();
            }

            __device__ __forceinline__ void measure_all(const StateView &state,
                                                        IdxType n_qubits,
                                                        IdxType repetition,
                                                        ValType *prefix_buffer,
                                                        const ValType *random_values,
                                                        IdxType random_index,
                                                        IdxType *results)
            {
                namespace cg = cooperative_groups;
                cg::grid_group grid = cg::this_grid();
                const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
                ValType *state_real = state.data_real;
                ValType *state_imag = state.data_imag;
                IdxType dim = state.dim ? state.dim : pow2i(n_qubits);
                IdxType n_size = dim;

                for (IdxType i = tid; i < dim; i += blockDim.x * gridDim.x)
                {
                    prefix_buffer[i] = state_real[i] * state_real[i] + state_imag[i] * state_imag[i];
                }
                grid.sync();

                for (IdxType d = 0; d < n_qubits; d++)
                {
                    IdxType step = pow2i(d + 1);
                    for (IdxType k = tid * step; k < n_size; k += step * blockDim.x * gridDim.x)
                    {
                        IdxType target = k + step - 1;
                        prefix_buffer[target] += prefix_buffer[k + pow2i(d) - 1];
                    }
                    grid.sync();
                }

                if (tid == 0)
                {
                    ValType total = prefix_buffer[n_size - 1];
                    prefix_buffer[n_size] = total;
                    prefix_buffer[n_size - 1] = 0;
                    ValType purity = fabs(total);
                    if (fabs(purity - 1.0) > MEASUREMENT_ERROR_BAR)
                    printf("MA: Purity Check fails with %lf
    ", purity);
                }
                grid.sync();

                if (n_qubits > 0)
                {
                    for (IdxType d = n_qubits - 1;; d--)
                    {
                        IdxType step = pow2i(d + 1);
                        for (IdxType k = tid * step; k < n_size - 1; k += step * blockDim.x * gridDim.x)
                        {
                            ValType tmp = prefix_buffer[k + pow2i(d) - 1];
                            ValType tmp2 = prefix_buffer[k + step - 1];
                            prefix_buffer[k + pow2i(d) - 1] = tmp2;
                            prefix_buffer[k + step - 1] = tmp + tmp2;
                        }
                        grid.sync();
                        if (d == 0)
                            break;
                    }
                }

                for (IdxType j = tid; j < n_size; j += blockDim.x * gridDim.x)
                {
                    ValType lower = prefix_buffer[j];
                    ValType upper = (j + 1 == n_size) ? 1 : prefix_buffer[j + 1];
                    for (IdxType shot = 0; shot < repetition; shot++)
                    {
                        ValType r = random_values[random_index + shot];
                        if (lower <= r && r < upper)
                            results[result_index + shot] = j;
                    }
                }
                grid.sync();
            }

            __device__ __forceinline__ void reset_qubit(const StateView &state,
                                                        IdxType n_qubits,
                                                        IdxType qubit,
                                                        ValType *workspace)
            {
                namespace cg = cooperative_groups;
                cg::grid_group grid = cg::this_grid();
                const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
                ValType *state_real = state.data_real;
                ValType *state_imag = state.data_imag;
                IdxType dim = state.dim ? state.dim : pow2i(n_qubits);
                IdxType mask = pow2i(qubit);

                for (IdxType i = tid; i < dim; i += blockDim.x * gridDim.x)
                {
                    if ((i & mask) == 0)
                        workspace[i] = 0;
                    else
                        workspace[i] = state_real[i] * state_real[i] + state_imag[i] * state_imag[i];
                }
                grid.sync();

                for (IdxType k = pow2i(n_qubits - 1); k > 0; k >>= 1)
                {
                    for (IdxType i = tid; i < k; i += blockDim.x * gridDim.x)
                    {
                        workspace[i] += workspace[i + k];
                    }
                    grid.sync();
                }

                ValType prob_of_one = workspace[0];
                grid.sync();

                if (prob_of_one < 1.0)
                {
                    ValType factor = 1.0 / sqrt(1.0 - prob_of_one);
                    for (IdxType i = tid; i < dim; i += blockDim.x * gridDim.x)
                    {
                        if ((i & mask) == 0)
                        {
                            state_real[i] *= factor;
                            state_imag[i] *= factor;
                        }
                        else
                        {
                            state_real[i] = 0;
                            state_imag[i] = 0;
                        }
                    }
                }
                else
                {
                    for (IdxType i = tid; i < dim; i += blockDim.x * gridDim.x)
                    {
                        if ((i & mask) == 0)
                        {
                            IdxType partner = i ^ mask;
                            ValType real_val = state_real[partner];
                            ValType imag_val = state_imag[partner];
                            state_real[i] = real_val;
                            state_imag[i] = imag_val;
                            state_real[partner] = 0;
                            state_imag[partner] = 0;
                        }
                    }
                }
                grid.sync();
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
                const ValType *state_real = state.data_real;
                const ValType *state_imag = state.data_imag;
                IdxType dim = state.dim ? state.dim : pow2i(n_qubits);

                ValType local = 0;
                for (IdxType i = tid; i < dim; i += total_threads)
                {
                    ValType prob = state_real[i] * state_real[i] + state_imag[i] * state_imag[i];
                    local += z_mask_sign(i, mask) * prob;
                }

                workspace[tid] += coeff * local;
                grid.sync();
            }

            __device__ __forceinline__ ValType compute_expectation_observable(const StateView &state,
                                                                              IdxType n_qubits,
                                                                              const ZObservableTerm *terms,
                                                                              IdxType term_count,
                                                                              ValType *workspace)
            {
                namespace cg = cooperative_groups;
                cg::grid_group grid = cg::this_grid();
                const IdxType tid = blockDim.x * blockIdx.x + threadIdx.x;
                const IdxType total_threads = blockDim.x * gridDim.x;

                if (tid < total_threads)
                    workspace[tid] = 0;
                grid.sync();

                for (IdxType t = 0; t < term_count; t++)
                {
                    expectation_mask(state, n_qubits, terms[t].mask, terms[t].coeff, workspace);
                    grid.sync();
                }

                IdxType reduce_limit = 1;
                while (reduce_limit < total_threads)
                    reduce_limit <<= 1;
                if (reduce_limit > total_threads)
                    reduce_limit = total_threads;

                for (IdxType offset = reduce_limit >> 1; offset > 0; offset >>= 1)
                {
                    if (tid < offset && tid + offset < total_threads)
                        workspace[tid] += workspace[tid + offset];
                    grid.sync();
                }

                ValType expectation = workspace[0];
                grid.sync();
                return expectation;
            }

            __device__ __forceinline__ void expectation_observable(const StateView &state,
                                                                   IdxType n_qubits,
                                                                   const ZObservableTerm *terms,
                                                                   IdxType term_count,
                                                                   ValType *workspace,
                                                                   ValType *result_ptr)
            {
                ValType value = compute_expectation_observable(state, n_qubits, terms, term_count, workspace);
                if (blockIdx.x == 0 && threadIdx.x == 0)
                    *result_ptr = value;
            }

        } // namespace GPU
    } // namespace GateKernels
} // namespace NWQSim
