#ifndef SVSIM_CPU_HPP
#define SVSIM_CPU_HPP

// #include "../../include/NWQSim.hpp"
#include "../util.hpp"
#include "../gate.hpp"
#include "../circuit.hpp"

#include "../fusion.hpp"

#include <random>
#include <cstring>
#include <algorithm>

#ifdef USE_OMP
#include "svsim_omp.hpp"
#endif

#define PRINT_SIM_TRACE

namespace NWQSim
{
    class SVSIM_CPU
    {
    public:
        // n_qubits is the number of qubits
        IdxType n_qubits;
        IdxType sv_size;
        IdxType dim;
        IdxType half_dim;

        // CPU arrays
        ValType *sv_real = NULL;
        ValType *sv_imag = NULL;
        ValType *m_real = NULL;
        // For measurement randoms
        ValType *randoms = NULL;
        // For measurement result
        IdxType *results = NULL;

        // Random
        std::mt19937 rng;
        std::uniform_real_distribution<ValType> uni_dist;

        // CPU memory usage
        ValType cpu_mem;

        SVSIM_CPU(IdxType _n_qubits)
        {
            n_qubits = _n_qubits;
            dim = (IdxType)1 << (n_qubits);
            half_dim = (IdxType)1 << (n_qubits - 1);

            sv_size = dim * (IdxType)sizeof(ValType);
            // CPU side initialization
            SAFE_ALOC_HOST(sv_real, sv_size);
            SAFE_ALOC_HOST(sv_imag, sv_size);
            memset(sv_real, 0, sv_size);
            memset(sv_imag, 0, sv_size);
            // State-vector initial state [0..0] = 1
            sv_real[0] = 1.;
            cpu_mem += sv_size * 4;

            SAFE_ALOC_HOST(m_real, sv_size + sizeof(ValType));
            memset(m_real, 0, sv_size + sizeof(ValType));
            rng.seed(time(0));
#ifdef PRINT_SIM_TRACE
            printf("SVSim_cpu is initialized!\n");
#endif
        }

        ~SVSIM_CPU()
        {
            // Release for CPU side
            SAFE_FREE_HOST(sv_real);
            SAFE_FREE_HOST(sv_imag);
            SAFE_FREE_HOST(randoms);
            SAFE_FREE_HOST(results);
            SAFE_FREE_HOST(m_real);
#ifdef PRINT_SIM_TRACE
            printf("SVSim_cpu is finalized!\n\n");
#endif
        }

        void reset_sim()
        {
            // Reset CPU input & output
            memset(sv_real, 0, sv_size);
            memset(sv_imag, 0, sv_size);
            memset(m_real, 0, sv_size + sizeof(ValType));
            // State Vector initial state [0..0] = 1
            sv_real[0] = 1.;
        }

        void set_seed(IdxType seed)
        {
            rng.seed(seed);
        }

        void simulation_kernel(std::shared_ptr<std::vector<Gate>> gates)
        {

            //=========================================
            for (auto g : *gates)
            {
                if (g.op_name == OP::RESET)
                {
#ifdef USE_OMP
                    RESET_GATE_OMP(sv_real, sv_imag, m_real,
                                   g.qubit, n_qubits, dim, half_dim);
#else
                    RESET_GATE(g.qubit);
#endif
                }

                else if (g.op_name == OP::M)
                {
                    SAFE_FREE_HOST(results);
                    SAFE_ALOC_HOST(results, sizeof(IdxType));
                    memset(results, 0, sizeof(IdxType));
                    ValType rand = uni_dist(rng);

#ifdef USE_OMP
                    M_GATE_OMP(sv_real, sv_imag, m_real,
                               results, g.qubit, rand, n_qubits, dim, half_dim);
#else
                    M_GATE(g.qubit, rand);
#endif
                }
                else if (g.op_name == OP::MA)
                {
                    auto repetition = g.qubit;
                    SAFE_FREE_HOST(results);
                    SAFE_ALOC_HOST(results, sizeof(IdxType) * repetition);
                    memset(results, 0, sizeof(IdxType) * repetition);
                    SAFE_FREE_HOST(randoms);
                    SAFE_ALOC_HOST(randoms, sizeof(ValType) * repetition);
                    for (IdxType i = 0; i < repetition; i++)
                        randoms[i] = uni_dist(rng);
#ifdef USE_OMP
                    MA_GATE_OMP(sv_real, sv_imag, m_real,
                                results, randoms,
                                n_qubits, g.qubit);
#else
                    MA_GATE(repetition);
#endif
                }
                else if (g.n_qubits == 1)
                {

#ifdef USE_OMP
                    C1_GATE_OMP(sv_real, sv_imag,
                                g.gm_real, g.gm_imag,
                                g.qubit, half_dim);
#else
                    C1_GATE(g.gm_real, g.gm_imag, g.qubit);
#endif
                }
                else if (g.n_qubits == 2)
                {

#ifdef USE_OMP
                    C2_GATE_OMP(sv_real, sv_imag,
                                g.gm_real, g.gm_imag,
                                g.ctrl, g.qubit,
                                dim);
#else
                    C2_GATE(g.gm_real, g.gm_imag, g.ctrl, g.qubit);
#endif
                }
                else
                {
                    std::cout << "unrecognized gates" << std::endl;
                }
            }
        }

        void sim(Circuit &circuit)
        {
            fuse_circuit(circuit);
            auto gates = circuit.gates;
            IdxType n_gates = gates->size();
            assert(circuit.num_qubits() == n_qubits);

            auto ncpus = 1;

#ifdef PRINT_SIM_TRACE
            double sim_time;
            cpu_timer sim_timer;
            sim_timer.start_timer();
#endif
#ifdef USE_OMP
            ncpus = omp_get_max_threads();
#pragma omp parallel num_threads(ncpus)
            {
#endif
                simulation_kernel(gates);
#ifdef USE_OMP
            }
#endif

#ifdef PRINT_SIM_TRACE
            sim_timer.stop_timer();
            sim_time = sim_timer.measure();
            printf("\n============== SV-Sim ===============\n");
            printf("nqubits:%lld, ngates:%lld, ncpus:%d, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms, mem:%.3lf MB, mem_per_cpu:%.3lf MB\n",
                   n_qubits, n_gates, ncpus, sim_time, 0.,
                   sim_time, cpu_mem / 1024 / 1024, cpu_mem / 1024 / 1024);
            printf("=====================================\n");
            fflush(stdout);
#endif
            //=========================================
        }

        IdxType *measure_all(IdxType repetition = DEFAULT_REPETITIONS)
        {
            SAFE_FREE_HOST(results);
            SAFE_ALOC_HOST(results, sizeof(IdxType) * repetition);
            memset(results, 0, sizeof(IdxType) * repetition);
            SAFE_FREE_HOST(randoms);
            SAFE_ALOC_HOST(randoms, sizeof(ValType) * repetition);
            for (IdxType i = 0; i < repetition; i++)
                randoms[i] = uni_dist(rng);

            MA_GATE(repetition);
            return results;
        }

        IdxType measure(IdxType qubit)
        {
            SAFE_FREE_HOST(results);
            SAFE_ALOC_HOST(results, sizeof(IdxType));
            memset(results, 0, sizeof(IdxType));

            ValType rand = uni_dist(rng);
            M_GATE(qubit, rand);
            return results[0];
        }

        /***********************************************
         * Key Macros
         ***********************************************/
#define PUT(arr, i, val) (arr[(i)] = (val))
#define GET(arr, i) (arr[(i)])
#define BARR  \
    while (0) \
    {         \
    };

// For C2 and C4 gates
#define DIV2E(x, y) ((x) >> (y))
#define MOD2E(x, y) ((x) & (((IdxType)1 << (y)) - (IdxType)1))
#define EXP2E(x) ((IdxType)1 << (x))
#define SV16IDX(x) (((x >> 3) & 1) * EXP2E(qubit0) + ((x >> 2) & 1) * EXP2E(qubit1) + ((x >> 1) & 1) * EXP2E(qubit2) + ((x & 1) * EXP2E(qubit3)))

        //============== C1 Gate ================
        // Arbitrary 1-qubit gate
        void C1_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit)
        {
            for (IdxType i = 0; i < half_dim; i++)
            {
                IdxType outer = (i >> qubit);
                IdxType inner = (i & (((IdxType)1 << qubit) - 1));
                IdxType offset = (outer << (qubit + 1));
                IdxType pos0 = (offset + inner);
                IdxType pos1 = (offset + inner + ((IdxType)1 << qubit));

                const ValType el0_real = sv_real[pos0];
                const ValType el0_imag = sv_imag[pos0];
                const ValType el1_real = sv_real[pos1];
                const ValType el1_imag = sv_imag[pos1];

                ValType sv_real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag);
                ValType sv_imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real);
                ValType sv_real_pos1 = (gm_real[2] * el0_real) - (gm_imag[2] * el0_imag) + (gm_real[3] * el1_real) - (gm_imag[3] * el1_imag);
                ValType sv_imag_pos1 = (gm_real[2] * el0_imag) + (gm_imag[2] * el0_real) + (gm_real[3] * el1_imag) + (gm_imag[3] * el1_real);

                sv_real[pos0] = sv_real_pos0;
                sv_imag[pos0] = sv_imag_pos0;
                sv_real[pos1] = sv_real_pos1;
                sv_imag[pos1] = sv_imag_pos1;
            }
        }

        //============== C2 Gate ================
        // Arbitrary 2-qubit gate
        inline void C2_GATE(const ValType *gm_real, const ValType *gm_imag,
                            const IdxType qubit0, const IdxType qubit1)
        {

            const IdxType per_pe_work = (dim >> 2);
            assert(qubit0 != qubit1); // Non-cloning

            const IdxType q0dim = ((IdxType)1 << std::max(qubit0, qubit1));
            const IdxType q1dim = ((IdxType)1 << std::min(qubit0, qubit1));
            // const IdxType outer_factor = (dim + q0dim + q0dim - 1) >> (std::max(qubit0, qubit1) + 1);
            const IdxType mider_factor = (q0dim + q1dim + q1dim - 1) >> (std::min(qubit0, qubit1) + 1);
            const IdxType inner_factor = q1dim;
            const IdxType qubit0_dim = ((IdxType)1 << qubit0);
            const IdxType qubit1_dim = ((IdxType)1 << qubit1);

            for (IdxType i = 0; i < per_pe_work; i++)
            {
                IdxType outer = ((i / inner_factor) / (mider_factor)) * (q0dim + q0dim);
                IdxType mider = ((i / inner_factor) % (mider_factor)) * (q1dim + q1dim);
                IdxType inner = i % inner_factor;
                IdxType pos0 = outer + mider + inner;
                IdxType pos1 = outer + mider + inner + qubit1_dim;
                IdxType pos2 = outer + mider + inner + qubit0_dim;
                IdxType pos3 = outer + mider + inner + q0dim + q1dim;

                const ValType el0_real = GET(sv_real, pos0);
                const ValType el0_imag = GET(sv_imag, pos0);
                const ValType el1_real = GET(sv_real, pos1);
                const ValType el1_imag = GET(sv_imag, pos1);
                const ValType el2_real = GET(sv_real, pos2);
                const ValType el2_imag = GET(sv_imag, pos2);
                const ValType el3_real = GET(sv_real, pos3);
                const ValType el3_imag = GET(sv_imag, pos3);
                // Real part
                ValType sv_real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag) + (gm_real[2] * el2_real) - (gm_imag[2] * el2_imag) + (gm_real[3] * el3_real) - (gm_imag[3] * el3_imag);
                ValType sv_real_pos1 = (gm_real[4] * el0_real) - (gm_imag[4] * el0_imag) + (gm_real[5] * el1_real) - (gm_imag[5] * el1_imag) + (gm_real[6] * el2_real) - (gm_imag[6] * el2_imag) + (gm_real[7] * el3_real) - (gm_imag[7] * el3_imag);
                ValType sv_real_pos2 = (gm_real[8] * el0_real) - (gm_imag[8] * el0_imag) + (gm_real[9] * el1_real) - (gm_imag[9] * el1_imag) + (gm_real[10] * el2_real) - (gm_imag[10] * el2_imag) + (gm_real[11] * el3_real) - (gm_imag[11] * el3_imag);
                ValType sv_real_pos3 = (gm_real[12] * el0_real) - (gm_imag[12] * el0_imag) + (gm_real[13] * el1_real) - (gm_imag[13] * el1_imag) + (gm_real[14] * el2_real) - (gm_imag[14] * el2_imag) + (gm_real[15] * el3_real) - (gm_imag[15] * el3_imag);
                // Imag part
                ValType sv_imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real) + (gm_real[2] * el2_imag) + (gm_imag[2] * el2_real) + (gm_real[3] * el3_imag) + (gm_imag[3] * el3_real);
                ValType sv_imag_pos1 = (gm_real[4] * el0_imag) + (gm_imag[4] * el0_real) + (gm_real[5] * el1_imag) + (gm_imag[5] * el1_real) + (gm_real[6] * el2_imag) + (gm_imag[6] * el2_real) + (gm_real[7] * el3_imag) + (gm_imag[7] * el3_real);
                ValType sv_imag_pos2 = (gm_real[8] * el0_imag) + (gm_imag[8] * el0_real) + (gm_real[9] * el1_imag) + (gm_imag[9] * el1_real) + (gm_real[10] * el2_imag) + (gm_imag[10] * el2_real) + (gm_real[11] * el3_imag) + (gm_imag[11] * el3_real);
                ValType sv_imag_pos3 = (gm_real[12] * el0_imag) + (gm_imag[12] * el0_real) + (gm_real[13] * el1_imag) + (gm_imag[13] * el1_real) + (gm_real[14] * el2_imag) + (gm_imag[14] * el2_real) + (gm_real[15] * el3_imag) + (gm_imag[15] * el3_real);

                PUT(sv_real, pos0, sv_real_pos0);
                PUT(sv_real, pos1, sv_real_pos1);
                PUT(sv_real, pos2, sv_real_pos2);
                PUT(sv_real, pos3, sv_real_pos3);

                PUT(sv_imag, pos0, sv_imag_pos0);
                PUT(sv_imag, pos1, sv_imag_pos1);
                PUT(sv_imag, pos2, sv_imag_pos2);
                PUT(sv_imag, pos3, sv_imag_pos3);
            }
        }

        //============== C4 Gate ================
        // Arbitrary 4-qubit gate
        inline void C4_GATE(ValType *sv_real, ValType *sv_imag,
                            const ValType *gm_real, const ValType *gm_imag,
                            const IdxType qubit0, const IdxType qubit1,
                            const IdxType qubit2, const IdxType qubit3)
        {
            assert(qubit0 != qubit1); // Non-cloning
            assert(qubit0 != qubit2); // Non-cloning
            assert(qubit0 != qubit3); // Non-cloning
            assert(qubit1 != qubit2); // Non-cloning
            assert(qubit1 != qubit3); // Non-cloning
            assert(qubit2 != qubit3); // Non-cloning
            // need to sort qubits: min->max: p, q, r, s
            const IdxType v0 = std::min(qubit0, qubit1);
            const IdxType v1 = std::min(qubit2, qubit3);
            const IdxType v2 = std::max(qubit0, qubit1);
            const IdxType v3 = std::max(qubit2, qubit3);
            const IdxType p = std::min(v0, v1);
            const IdxType q = std::min(std::min(v2, v3), std::max(v0, v1));
            const IdxType r = std::max(std::min(v2, v3), std::max(v0, v1));
            const IdxType s = std::max(v2, v3);
            for (IdxType i = 0; i < (dim >> 4); i++)
            {
                const IdxType term0 = MOD2E(i, p);
                const IdxType term1 = MOD2E(DIV2E(i, p), q - p - 1) * EXP2E(p + 1);
                const IdxType term2 = MOD2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1) * EXP2E(q + 1);
                const IdxType term3 = MOD2E(DIV2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(r + 1);
                const IdxType term4 = DIV2E(DIV2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(s + 1);
                const IdxType term = term4 + term3 + term2 + term1 + term0;
                const ValType el_real[16] = {
                    GET(sv_real, term + SV16IDX(0)), GET(sv_real, term + SV16IDX(1)),
                    GET(sv_real, term + SV16IDX(2)), GET(sv_real, term + SV16IDX(3)),
                    GET(sv_real, term + SV16IDX(4)), GET(sv_real, term + SV16IDX(5)),
                    GET(sv_real, term + SV16IDX(6)), GET(sv_real, term + SV16IDX(7)),
                    GET(sv_real, term + SV16IDX(8)), GET(sv_real, term + SV16IDX(9)),
                    GET(sv_real, term + SV16IDX(10)), GET(sv_real, term + SV16IDX(11)),
                    GET(sv_real, term + SV16IDX(12)), GET(sv_real, term + SV16IDX(13)),
                    GET(sv_real, term + SV16IDX(14)), GET(sv_real, term + SV16IDX(15))};
                const ValType el_imag[16] = {
                    GET(sv_imag, term + SV16IDX(0)), GET(sv_imag, term + SV16IDX(1)),
                    GET(sv_imag, term + SV16IDX(2)), GET(sv_imag, term + SV16IDX(3)),
                    GET(sv_imag, term + SV16IDX(4)), GET(sv_imag, term + SV16IDX(5)),
                    GET(sv_imag, term + SV16IDX(6)), GET(sv_imag, term + SV16IDX(7)),
                    GET(sv_imag, term + SV16IDX(8)), GET(sv_imag, term + SV16IDX(9)),
                    GET(sv_imag, term + SV16IDX(10)), GET(sv_imag, term + SV16IDX(11)),
                    GET(sv_imag, term + SV16IDX(12)), GET(sv_imag, term + SV16IDX(13)),
                    GET(sv_imag, term + SV16IDX(14)), GET(sv_imag, term + SV16IDX(15))};
#pragma unroll
                for (unsigned j = 0; j < 16; j++)
                {
                    ValType res_real = 0;
                    ValType res_imag = 0;
#pragma unroll
                    for (unsigned k = 0; k < 16; k++)
                    {
                        res_real += (el_real[k] * gm_real[j * 16 + k]) - (el_imag[k] * gm_imag[j * 16 + k]);
                        res_imag += (el_real[k] * gm_imag[j * 16 + k]) + (el_imag[k] * gm_real[j * 16 + k]);
                    }
                    PUT(sv_real, term + SV16IDX(j), res_real);
                    PUT(sv_imag, term + SV16IDX(j), res_imag);
                }
            }
        }

        void M_GATE(const IdxType qubit, const ValType rand)
        {
            IdxType mask = ((IdxType)1 << qubit);
            ValType prob_of_one = 0;
            for (IdxType i = 0; i < ((IdxType)1 << n_qubits); i++)
            {
                if ((i & mask) != 0)
                    prob_of_one += (sv_real[i] * sv_real[i]) + (sv_imag[i] * sv_imag[i]); // square
            }
            if (rand < prob_of_one)
            {
                ValType normalize_factor = sqrt(prob_of_one);
                for (IdxType i = 0; i < dim; i++)
                {
                    if ((i & mask) == 0)
                    {
                        sv_real[i] = 0.;
                        sv_imag[i] = 0.;
                    }
                    else
                    {
                        sv_real[i] /= normalize_factor;
                        sv_imag[i] /= normalize_factor;
                    }
                }
            }
            else
            {
                for (IdxType i = 0; i < dim; i++)
                {
                    ValType normalize_factor = sqrt(1.0 - prob_of_one);
                    if ((i & mask) == 0)
                    {
                        sv_real[i] /= normalize_factor;
                        sv_imag[i] /= normalize_factor;
                    }
                    else
                    {
                        sv_real[i] = 0;
                        sv_imag[i] = 0;
                    }
                }
            }
            results[0] = (rand < prob_of_one ? 1 : 0);
        }

        //============== MA Gate (Measure all qubits in Pauli-Z) ================
        inline void MA_GATE(const IdxType repetition)
        {
            IdxType n_size = (IdxType)1 << n_qubits;
            m_real[0] = 0;
            for (IdxType i = 1; i < (((IdxType)1 << n_qubits) + 1); i++)
                m_real[i] = m_real[i - 1] + ((sv_real[i - 1] * sv_real[i - 1]) + (sv_imag[i - 1] * sv_imag[i - 1]));
            ValType purity = fabs(m_real[((IdxType)1 << n_qubits)]);
            if (abs(purity - 1.0) > ERROR_BAR)
                printf("MA: Purity Check fails with %lf\n", purity);

            if (repetition < n_size)
            {
                for (IdxType j = 0; j < n_size; j++)
                {
                    ValType lower = m_real[j];
                    ValType upper = (j + 1 == n_size) ? 1 : m_real[j + 1];
                    for (IdxType i = 0; i < repetition; i++)
                    {
                        ValType r = randoms[i];
                        if (lower <= r && r < upper)
                            results[i] = j;
                    }
                }
            }
            else
            {
                for (IdxType i = 0; i < repetition; i++)
                {
                    IdxType lo = 0;
                    IdxType hi = ((IdxType)1 << n_qubits);
                    IdxType mid;
                    ValType r = randoms[i];
                    while (hi - lo > 1)
                    {
                        mid = lo + (hi - lo) / 2;
                        if (r >= m_real[mid])
                            lo = mid;
                        else
                            hi = mid;
                    }
                    results[i] = lo;
                }
            }
        }

        //============== Reset ================
        inline void RESET_GATE(const IdxType qubit)
        {
            IdxType mask = ((IdxType)1 << qubit);
            ValType prob_of_one = 0;
            for (IdxType i = 0; i < ((IdxType)1 << n_qubits); i++)
            {
                if ((i & mask) != 0)
                    prob_of_one += (sv_real[i] * sv_real[i]) + (sv_imag[i] * sv_imag[i]); // square
            }
            assert(prob_of_one <= 1.0);

            if (prob_of_one < 1.0) // still possible to normalize
            {
                ValType normalize_factor = sqrt(1.0 - prob_of_one);
                for (IdxType i = 0; i < dim; i++)
                {
                    if ((i & mask) == 0)
                    {
                        sv_real[i] /= normalize_factor;
                        sv_imag[i] /= normalize_factor;
                    }
                    else
                    {
                        sv_real[i] = 0;
                        sv_imag[i] = 0;
                    }
                }
            }
            else
            {
                for (IdxType i = 0; i < dim; i++)
                {
                    if ((i & mask) == 0)
                    {
                        IdxType dual_i = i ^ mask;
                        sv_real[i] = sv_real[dual_i];
                        sv_imag[i] = sv_imag[dual_i];
                        sv_real[dual_i] = 0;
                        sv_imag[dual_i] = 0;
                    }
                }
            }
        }

        //============== Purity Check  ================
        void Purity_Check(Gate *g, const IdxType t, ValType *sv_real, ValType *sv_imag)
        {
            ValType purity = 0;
            for (IdxType i = 0; i < (((IdxType)1 << n_qubits)); i++)
                purity += ((sv_real[i] * sv_real[i]) + (sv_imag[i] * sv_imag[i]));
            if (fabs(purity - 1.0) > ERROR_BAR)
            {
                printf("Purity Check fails after Gate-%lld=>%s(ctrl:%lld,qubit:%lld,theta:%lf) with %lf\n", t, OP_NAMES[g->op_name], g->ctrl, g->qubit, g->theta, purity);
            }
        }
    };

} // namespace NWQSim

#endif // SVSIM_CPU