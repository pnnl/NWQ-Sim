#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <omp.h>

#include "nwqsim_core.hpp"

namespace NWQSim
{
    namespace GateKernels
    {
        namespace OMP
        {
            inline IdxType pow2i(IdxType power)
            {
                return static_cast<IdxType>(IdxType(1) << power);
            }

            inline IdxType mod_pow2(IdxType value, IdxType power)
            {
                return value & (pow2i(power) - IdxType(1));
            }

            inline IdxType div_pow2(IdxType value, IdxType power)
            {
                return value >> power;
            }

            inline IdxType state16_offset(unsigned pattern, IdxType q0, IdxType q1, IdxType q2, IdxType q3)
            {
                return ((pattern >> 3) & 1 ? pow2i(q0) : IdxType(0)) +
                       ((pattern >> 2) & 1 ? pow2i(q1) : IdxType(0)) +
                       ((pattern >> 1) & 1 ? pow2i(q2) : IdxType(0)) +
                       ((pattern) & 1 ? pow2i(q3) : IdxType(0));
            }

            inline void apply_c1_gate(const ValType *gm_real,
                                      const ValType *gm_imag,
                                      IdxType qubit,
                                      const StateSpan &state)
            {
                ValType *state_real = state.real;
                ValType *state_imag = state.imag;
                const IdxType half_dim = state.half_extent();
#pragma omp for schedule(auto)
                for (IdxType i = 0; i < half_dim; i++)
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
            }

            inline void apply_c2_gate(const ValType *gm_real,
                                      const ValType *gm_imag,
                                      IdxType qubit0,
                                      IdxType qubit1,
                                      const StateSpan &state)
            {
                ValType *state_real = state.real;
                ValType *state_imag = state.imag;
                const IdxType dim = state.full_extent();
                const IdxType per_work = (dim >> 2);
                const IdxType q0dim = pow2i(std::max(qubit0, qubit1));
                const IdxType q1dim = pow2i(std::min(qubit0, qubit1));
                const IdxType mid_factor = (q0dim + q1dim + q1dim - 1) >> (std::min(qubit0, qubit1) + 1);
                const IdxType inner_factor = q1dim;
                const IdxType qubit0_dim = pow2i(qubit0);
                const IdxType qubit1_dim = pow2i(qubit1);

#pragma omp for schedule(auto)
                for (IdxType i = 0; i < per_work; i++)
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
            }

            inline void apply_c4_gate(const ValType *gm_real,
                                      const ValType *gm_imag,
                                      IdxType qubit0,
                                      IdxType qubit1,
                                      IdxType qubit2,
                                      IdxType qubit3,
                                      const StateSpan &state)
            {
                assert(qubit0 != qubit1);
                assert(qubit0 != qubit2);
                assert(qubit0 != qubit3);
                assert(qubit1 != qubit2);
                assert(qubit1 != qubit3);
                assert(qubit2 != qubit3);

                ValType *state_real = state.real;
                ValType *state_imag = state.imag;
                const IdxType dim = state.full_extent();
                const IdxType v0 = std::min(qubit0, qubit1);
                const IdxType v1 = std::min(qubit2, qubit3);
                const IdxType v2 = std::max(qubit0, qubit1);
                const IdxType v3 = std::max(qubit2, qubit3);
                const IdxType p = std::min(v0, v1);
                const IdxType q = std::min(std::min(v2, v3), std::max(v0, v1));
                const IdxType r = std::max(std::min(v2, v3), std::max(v0, v1));
                const IdxType s = std::max(v2, v3);

#pragma omp for schedule(auto)
                for (IdxType i = 0; i < (dim >> 4); i++)
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
            }

            inline IdxType measure_qubit(const StateSpan &state,
                                         IdxType n_qubits,
                                         IdxType qubit,
                                         ValType *scratch,
                                         ValType random_value)
            {
                ValType *state_real = state.real;
                ValType *state_imag = state.imag;
                IdxType dim = state.size ? state.size : pow2i(n_qubits);
                IdxType half_dim = dim >> 1;
                IdxType mask = pow2i(qubit);

#pragma omp for schedule(auto)
                for (IdxType i = 0; i < dim; i++)
                {
                    if ((i & mask) == 0)
                        scratch[i] = 0;
                    else
                        scratch[i] = state_real[i] * state_real[i] + state_imag[i] * state_imag[i];
                }
#pragma omp barrier

                for (IdxType k = half_dim; k > 0; k >>= 1)
                {
#pragma omp for schedule(auto)
                    for (IdxType i = 0; i < k; i++)
                    {
                        scratch[i] += scratch[i + k];
                    }
#pragma omp barrier
                }

                ValType prob_of_one = scratch[0];
#pragma omp barrier
                if (random_value < prob_of_one)
                {
                    ValType normalize_factor = std::sqrt(prob_of_one);
#pragma omp for schedule(auto)
                    for (IdxType i = 0; i < dim; i++)
                    {
                        if ((i & mask) == 0)
                        {
                            state_real[i] = 0;
                            state_imag[i] = 0;
                        }
                        else
                        {
                            state_real[i] /= normalize_factor;
                            state_imag[i] /= normalize_factor;
                        }
                    }
                }
                else
                {
#pragma omp for schedule(auto)
                    for (IdxType i = 0; i < dim; i++)
                    {
                        ValType normalize_factor = std::sqrt(1.0 - prob_of_one);
                        if ((i & mask) == 0)
                        {
                            state_real[i] /= normalize_factor;
                            state_imag[i] /= normalize_factor;
                        }
                        else
                        {
                            state_real[i] = 0;
                            state_imag[i] = 0;
                        }
                    }
                }
#pragma omp barrier
                return (random_value < prob_of_one) ? 1 : 0;
            }

            inline void measure_all(const StateSpan &state,
                                    IdxType n_qubits,
                                    IdxType repetition,
                                    ValType *scratch,
                                    const ValType *random_values,
                                    IdxType *results_out)
            {
                ValType *state_real = state.real;
                ValType *state_imag = state.imag;
                IdxType n_size = state.size ? state.size : pow2i(n_qubits);
#pragma omp for schedule(auto)
                for (IdxType i = 0; i < n_size; i++)
                    scratch[i] = state_real[i] * state_real[i] + state_imag[i] * state_imag[i];
#pragma omp barrier

                for (IdxType d = 0; d < n_qubits; d++)
                {
                    IdxType step = pow2i(d + 1);
#pragma omp for schedule(auto)
                    for (IdxType k = 0; k < n_size; k += step)
                        scratch[k + step - 1] = scratch[k + pow2i(d) - 1] + scratch[k + step - 1];
#pragma omp barrier
                }

                if (omp_get_thread_num() == 0)
                {
                    ValType purity = std::fabs(scratch[n_size - 1]);
                    scratch[n_size - 1] = 0;
                    if (std::fabs(purity - 1.0) > MEASUREMENT_ERROR_BAR)
                        std::printf("MA: Purity Check fails with %lf\n", purity);
                }
#pragma omp barrier

                if (n_qubits > 0)
                {
                    for (IdxType d = n_qubits - 1;; d--)
                    {
                        IdxType step = pow2i(d + 1);
#pragma omp for schedule(auto)
                        for (IdxType k = 0; k < n_size - 1; k += step)
                        {
                            ValType tmp = scratch[k + pow2i(d) - 1];
                            ValType tmp2 = scratch[k + step - 1];
                            scratch[k + pow2i(d) - 1] = tmp2;
                            scratch[k + step - 1] = tmp + tmp2;
                        }
#pragma omp barrier
                        if (d == 0)
                            break;
                    }
                }

#pragma omp for schedule(auto)
                for (IdxType shot = 0; shot < repetition; shot++)
                {
                    IdxType lo = 0;
                    IdxType hi = n_size;
                    ValType r = random_values[shot];
                    while (hi - lo > 1)
                    {
                        IdxType mid = lo + (hi - lo) / 2;
                        if (r >= scratch[mid])
                            lo = mid;
                        else
                            hi = mid;
                    }
                    results_out[shot] = lo;
                }
#pragma omp barrier
            }

            inline void reset_qubit(const StateSpan &state,
                                    IdxType n_qubits,
                                    IdxType qubit,
                                    ValType *scratch)
            {
                ValType *state_real = state.real;
                ValType *state_imag = state.imag;
                IdxType dim = state.size ? state.size : pow2i(n_qubits);
                IdxType half_dim = dim >> 1;
                IdxType mask = pow2i(qubit);

#pragma omp for schedule(auto)
                for (IdxType i = 0; i < dim; i++)
                {
                    if ((i & mask) == 0)
                        scratch[i] = 0;
                    else
                        scratch[i] = state_real[i] * state_real[i] + state_imag[i] * state_imag[i];
                }
#pragma omp barrier

                for (IdxType k = half_dim; k > 0; k >>= 1)
                {
#pragma omp for schedule(auto)
                    for (IdxType i = 0; i < k; i++)
                        scratch[i] += scratch[i + k];
#pragma omp barrier
                }

                ValType prob_of_one = scratch[0];
#pragma omp barrier
                if (prob_of_one < 1.0)
                {
                    ValType factor = 1.0 / std::sqrt(1.0 - prob_of_one);
#pragma omp for schedule(auto)
                    for (IdxType i = 0; i < dim; i++)
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
#pragma omp for schedule(auto)
                    for (IdxType i = 0; i < dim; i++)
                    {
                        if ((i & mask) == 0)
                        {
                            IdxType partner = i ^ mask;
                            state_real[i] = state_real[partner];
                            state_imag[i] = state_imag[partner];
                            state_real[partner] = 0;
                            state_imag[partner] = 0;
                        }
                    }
                }
#pragma omp barrier
            }

            inline ValType expectation_mask(const StateSpan &state,
                                            IdxType n_qubits,
                                            IdxType mask)
            {
                const ValType *state_real = state.real;
                const ValType *state_imag = state.imag;
                IdxType dim = state.size ? state.size : pow2i(n_qubits);
                ValType total = 0;
#pragma omp for reduction(+ : total) schedule(auto)
                for (IdxType i = 0; i < dim; i++)
                {
                    ValType prob = state_real[i] * state_real[i] + state_imag[i] * state_imag[i];
                    total += z_mask_sign(i, mask) * prob;
                }
#pragma omp barrier
                return total;
            }

            inline ValType expectation_observable(const StateSpan &state,
                                                  IdxType n_qubits,
                                                  const ZObservableTerm *terms,
                                                  IdxType term_count)
            {
                ValType result = 0;
                for (IdxType t = 0; t < term_count; t++)
                {
                    result += terms[t].coeff * expectation_mask(state, n_qubits, terms[t].mask);
                }
                return result;
            }

        } // namespace OMP
    } // namespace GateKernels
} // namespace NWQSim
