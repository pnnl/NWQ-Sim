#pragma once

#include "../state.hpp"

#include "../nwq_util.hpp"
#include "../gate.hpp"
#include "../circuit.hpp"
#include "../private/config.hpp"

#include "../circuit_pass/fusion.hpp"
#include "../private/macros.hpp"
#include "../private/sim_gate.hpp"

#include "../private/gate_factory/dm_gates.hpp"

#include <random>
#include <cstring>
#include <algorithm>
#include <vector>
#include <string>

namespace NWQSim
{
    class DM_CPU : public QuantumState
    {

    public:
        DM_CPU(IdxType _n_qubits) : QuantumState(_n_qubits)
        {
            // Initialize CPU side
            n_qubits = _n_qubits;

            if (Config::ENABLE_NOISE)
            {
                IdxType device_qubits = Config::backend_config["num_qubits"];
                if (n_qubits > device_qubits)
                {
                    std::string msg = "Error: Circuit uses " + std::to_string(n_qubits) + " qubits, more than " + std::to_string(device_qubits) + "qubits in the device!!\n";
                    throw std::logic_error(msg.c_str());
                }
            }

            dim = (IdxType)1 << (2 * n_qubits);
            half_dim = (IdxType)1 << (2 * n_qubits - 1);
            dm_size = dim * (IdxType)sizeof(ValType);
            n_cpu = 1;

            diag_num = ((IdxType)1 << n_qubits) + 1;
            diag_size = diag_num * sizeof(ValType);

            // CPU side initialization
            SAFE_ALOC_HOST(dm_real, dm_size);
            SAFE_ALOC_HOST(dm_imag, dm_size);

            SAFE_ALOC_HOST(m_real, diag_size);
            memset(dm_real, 0, dm_size);
            memset(dm_imag, 0, dm_size);
            memset(m_real, 0, diag_size);

            // State-vector initial state [0..0] = 1
            dm_real[0] = 1.;
            cpu_mem += dm_size * 2 + diag_size;

            rng.seed(time(0));
        }

        ~DM_CPU()
        {
            // Release for CPU side
            SAFE_FREE_HOST(dm_real);
            SAFE_FREE_HOST(dm_imag);
            SAFE_FREE_HOST(m_real);
            SAFE_FREE_HOST(results);
        }

        void reset_state() override
        {
            // Reset CPU input & output
            memset(dm_real, 0, dm_size);
            memset(dm_imag, 0, dm_size);
            memset(m_real, 0, diag_size);
            // State Vector initial state [0..0] = 1
            dm_real[0] = 1.;
        }

        void set_seed(IdxType seed) override
        {
            rng.seed(seed);
        }

        void sim(std::shared_ptr<NWQSim::Circuit> circuit) override
        {
            IdxType origional_gates = circuit->num_gates();

            std::vector<DMGate> gates = fuse_circuit_dm(circuit);

            // getDMGates(circuit->get_gates(), circuit->num_qubits());

            IdxType n_gates = gates.size();

            assert(circuit->num_qubits() == n_qubits);

            double sim_time;
            cpu_timer sim_timer;
            sim_timer.start_timer();

            simulation_kernel(gates);

            sim_timer.stop_timer();
            sim_time = sim_timer.measure();

            if (Config::PRINT_SIM_TRACE)
            {
                printf("\n============== DM-Sim ===============\n");
                printf("n_qubits:%lld, n_gates:%lld, sim_gates:%lld, ncpus:%lld, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms, mem:%.3lf MB, mem_per_cpu:%.3lf MB\n",
                       n_qubits, origional_gates, n_gates, n_cpu, sim_time, 0.,
                       sim_time, cpu_mem / 1024 / 1024, cpu_mem / 1024 / 1024);
                printf("=====================================\n");
            }

            //=========================================
        }

        IdxType *get_results() override
        {
            return results;
        }

        IdxType measure(IdxType qubit) override
        {
            M_GATE(qubit);
            return results[0];
        }

        IdxType *measure_all(IdxType repetition) override
        {
            MA_GATE(repetition);
            return results;
        }

        virtual ValType get_exp_z(const std::vector<size_t> &in_bits) override
        {
            throw std::logic_error("get_exp_z Not implemented (DM_CPU)");
        }

        virtual ValType get_exp_z() override
        {
            throw std::logic_error("get_exp_z Not implemented (DM_CPU)");
        }

        void print_res_state() override
        {
            IdxType num = ((IdxType)1 << n_qubits);
            printf("----- Real DM diag ------\n");
            for (IdxType i = 0; i < num; i++)
            {
                printf("%lf ", dm_real[i * num + i]);
                if ((i + 1) % 8 == 0)
                    printf("\n");
            }
            printf("\n");
        }

    protected:
        // n_qubits is the number of qubits
        IdxType n_qubits;
        IdxType dm_size;
        IdxType dim;
        IdxType half_dim;
        IdxType n_cpu;

        // for measurement along the dm diagnol
        IdxType diag_num;
        IdxType diag_size;
        // CPU arrays
        ValType *dm_real;
        ValType *dm_imag;
        ValType *m_real;

        IdxType *results = NULL;

        // Random
        std::mt19937 rng;
        std::uniform_real_distribution<ValType> uni_dist;

        // CPU memory usage
        ValType cpu_mem;

        virtual void simulation_kernel(const std::vector<DMGate> &gates)
        {
            auto start = std::chrono::steady_clock::now();
            int n_gates = gates.size();
            for (int i = 0; i < n_gates; i++)
            {
                if (Config::PRINT_SIM_TRACE)
                    printProgressBar(i, n_gates, start);

                auto g = gates[i];

                if (g.op_name == OP::C2)
                {
                    C2_GATE(g.gm_real, g.gm_imag, g.qubit, g.qubit + n_qubits);
                }
                else if (g.op_name == OP::C4)
                {
                    C4_GATE(g.gm_real, g.gm_imag, g.ctrl, g.qubit, g.ctrl + n_qubits, g.qubit + n_qubits);
                }
                else if (g.op_name == OP::RESET)
                {
                    RESET_GATE(g.qubit);
                }
                else if (g.op_name == OP::M)
                {
                    M_GATE(g.qubit);
                }
                else if (g.op_name == OP::MA)
                {
                    MA_GATE(g.qubit);
                }
                else
                {
                    std::cout << "Unrecognized gates" << std::endl
                              << OP_NAMES[g.op_name] << std::endl;
                    std::logic_error("Invalid gate type");
                }
            }
            std::cout << std::endl;
        }

        //============== Local 2-qubit Gate  ================
        void C2_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit0, const IdxType qubit1)
        {
            const IdxType per_pe_work = ((dim) >> 2);
            assert(qubit0 != qubit1); // Non-cloning
            const IdxType q0dim = ((IdxType)1 << std::max(qubit0, qubit1));
            const IdxType q1dim = ((IdxType)1 << std::min(qubit0, qubit1));
            const IdxType outer_factor = ((dim) + q0dim + q0dim - 1) >> (std::max(qubit0, qubit1) + 1);
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

                const ValType el0_real = GET(dm_real, pos0);
                const ValType el0_imag = GET(dm_imag, pos0);
                const ValType el1_real = GET(dm_real, pos1);
                const ValType el1_imag = GET(dm_imag, pos1);
                const ValType el2_real = GET(dm_real, pos2);
                const ValType el2_imag = GET(dm_imag, pos2);
                const ValType el3_real = GET(dm_real, pos3);
                const ValType el3_imag = GET(dm_imag, pos3);

                // Real part
                ValType dm_real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag) + (gm_real[2] * el2_real) - (gm_imag[2] * el2_imag) + (gm_real[3] * el3_real) - (gm_imag[3] * el3_imag);
                ValType dm_real_pos1 = (gm_real[4] * el0_real) - (gm_imag[4] * el0_imag) + (gm_real[5] * el1_real) - (gm_imag[5] * el1_imag) + (gm_real[6] * el2_real) - (gm_imag[6] * el2_imag) + (gm_real[7] * el3_real) - (gm_imag[7] * el3_imag);
                ValType dm_real_pos2 = (gm_real[8] * el0_real) - (gm_imag[8] * el0_imag) + (gm_real[9] * el1_real) - (gm_imag[9] * el1_imag) + (gm_real[10] * el2_real) - (gm_imag[10] * el2_imag) + (gm_real[11] * el3_real) - (gm_imag[11] * el3_imag);
                ValType dm_real_pos3 = (gm_real[12] * el0_real) - (gm_imag[12] * el0_imag) + (gm_real[13] * el1_real) - (gm_imag[13] * el1_imag) + (gm_real[14] * el2_real) - (gm_imag[14] * el2_imag) + (gm_real[15] * el3_real) - (gm_imag[15] * el3_imag);

                // Imag part
                ValType dm_imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real) + (gm_real[2] * el2_imag) + (gm_imag[2] * el2_real) + (gm_real[3] * el3_imag) + (gm_imag[3] * el3_real);
                ValType dm_imag_pos1 = (gm_real[4] * el0_imag) + (gm_imag[4] * el0_real) + (gm_real[5] * el1_imag) + (gm_imag[5] * el1_real) + (gm_real[6] * el2_imag) + (gm_imag[6] * el2_real) + (gm_real[7] * el3_imag) + (gm_imag[7] * el3_real);
                ValType dm_imag_pos2 = (gm_real[8] * el0_imag) + (gm_imag[8] * el0_real) + (gm_real[9] * el1_imag) + (gm_imag[9] * el1_real) + (gm_real[10] * el2_imag) + (gm_imag[10] * el2_real) + (gm_real[11] * el3_imag) + (gm_imag[11] * el3_real);
                ValType dm_imag_pos3 = (gm_real[12] * el0_imag) + (gm_imag[12] * el0_real) + (gm_real[13] * el1_imag) + (gm_imag[13] * el1_real) + (gm_real[14] * el2_imag) + (gm_imag[14] * el2_real) + (gm_real[15] * el3_imag) + (gm_imag[15] * el3_real);

                PUT(dm_real, pos0, dm_real_pos0);
                PUT(dm_real, pos1, dm_real_pos1);
                PUT(dm_real, pos2, dm_real_pos2);
                PUT(dm_real, pos3, dm_real_pos3);

                PUT(dm_imag, pos0, dm_imag_pos0);
                PUT(dm_imag, pos1, dm_imag_pos1);
                PUT(dm_imag, pos2, dm_imag_pos2);
                PUT(dm_imag, pos3, dm_imag_pos3);
            }
            // BARR;
        }

        //============== Local 4-qubit Gate ================
        void C4_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit0, const IdxType qubit1,
                     const IdxType qubit2, const IdxType qubit3)
        {
            const IdxType per_pe_work = ((dim) >> 4);
            assert(qubit0 != qubit1); // Non-cloning
            assert(qubit0 != qubit2); // Non-cloning
            assert(qubit0 != qubit3); // Non-cloning
            assert(qubit1 != qubit2); // Non-cloning
            assert(qubit1 != qubit3); // Non-cloning
            assert(qubit2 != qubit3); // Non-cloning

            // need to sort qubits: std::min->std::max: p, q, r, s
            const IdxType v0 = std::min(qubit0, qubit1);
            const IdxType v1 = std::min(qubit2, qubit3);
            const IdxType v2 = std::max(qubit0, qubit1);
            const IdxType v3 = std::max(qubit2, qubit3);
            const IdxType p = std::min(v0, v1);
            const IdxType q = std::min(std::min(v2, v3), std::max(v0, v1));
            const IdxType r = std::max(std::min(v2, v3), std::max(v0, v1));
            const IdxType s = std::max(v2, v3);

            for (IdxType i = 0; i < per_pe_work; i++)
            {
                const IdxType term0 = MOD2E(i, p);
                const IdxType term1 = MOD2E(DIV2E(i, p), q - p - 1) * EXP2E(p + 1);
                const IdxType term2 = MOD2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1) * EXP2E(q + 1);
                const IdxType term3 = MOD2E(DIV2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(r + 1);
                const IdxType term4 = DIV2E(DIV2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(s + 1);
                const IdxType term = term4 + term3 + term2 + term1 + term0;
                const ValType el_real[16] = {
                    GET(dm_real, term + SV16IDX(0)), GET(dm_real, term + SV16IDX(1)),
                    GET(dm_real, term + SV16IDX(2)), GET(dm_real, term + SV16IDX(3)),
                    GET(dm_real, term + SV16IDX(4)), GET(dm_real, term + SV16IDX(5)),
                    GET(dm_real, term + SV16IDX(6)), GET(dm_real, term + SV16IDX(7)),
                    GET(dm_real, term + SV16IDX(8)), GET(dm_real, term + SV16IDX(9)),
                    GET(dm_real, term + SV16IDX(10)), GET(dm_real, term + SV16IDX(11)),
                    GET(dm_real, term + SV16IDX(12)), GET(dm_real, term + SV16IDX(13)),
                    GET(dm_real, term + SV16IDX(14)), GET(dm_real, term + SV16IDX(15))};
                const ValType el_imag[16] = {
                    GET(dm_imag, term + SV16IDX(0)), GET(dm_imag, term + SV16IDX(1)),
                    GET(dm_imag, term + SV16IDX(2)), GET(dm_imag, term + SV16IDX(3)),
                    GET(dm_imag, term + SV16IDX(4)), GET(dm_imag, term + SV16IDX(5)),
                    GET(dm_imag, term + SV16IDX(6)), GET(dm_imag, term + SV16IDX(7)),
                    GET(dm_imag, term + SV16IDX(8)), GET(dm_imag, term + SV16IDX(9)),
                    GET(dm_imag, term + SV16IDX(10)), GET(dm_imag, term + SV16IDX(11)),
                    GET(dm_imag, term + SV16IDX(12)), GET(dm_imag, term + SV16IDX(13)),
                    GET(dm_imag, term + SV16IDX(14)), GET(dm_imag, term + SV16IDX(15))};
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
                    PUT(dm_real, term + SV16IDX(j), res_real);
                    PUT(dm_imag, term + SV16IDX(j), res_imag);
                }
            }
            // BARR;
        }

        void M_GATE(const IdxType qubit)
        {
            SAFE_FREE_HOST(results);
            SAFE_ALOC_HOST(results, sizeof(IdxType));
            memset(results, 0, sizeof(IdxType));

            ValType random = uni_dist(rng);
            IdxType mask = ((IdxType)1 << qubit);
            ValType prob_of_one = 0;
            for (IdxType i = 0; i < ((IdxType)1 << (n_qubits)); i++)
            {
                if ((i & mask) != 0)
                    prob_of_one += fabs(dm_real[(i << (n_qubits)) + i]);
            }

            ValType gm_real[16];
            ValType gm_imag[16];
            if (random <= prob_of_one)
                gm_real[15] = 1.0 / prob_of_one;
            else
                gm_real[0] = 1.0 / (1.0 - prob_of_one);
            BARR;
            C2_GATE(gm_real, gm_imag, qubit, qubit + n_qubits);
            BARR;
            results[0] = (random <= prob_of_one ? 1 : 0);
        }

        void MA_GATE(const IdxType repetition)
        {
            SAFE_FREE_HOST(results);
            SAFE_ALOC_HOST(results, sizeof(IdxType) * repetition);
            memset(results, 0, sizeof(IdxType) * repetition);

            const IdxType n_size = (IdxType)1 << (n_qubits);
            m_real[0] = 0;
            for (IdxType i = 1; i <= (n_size); i++)
                m_real[i] = m_real[i - 1] + fabs(dm_real[((i - 1) << (n_qubits)) + i - 1]);
            ValType purity = fabs(m_real[n_size]);

            if (fabs(purity - 1.0) > ERROR_BAR)
                printf("MA: Purity Check fails with %lf\n", purity);

            for (IdxType i = 0; i < repetition; i++)
            {
                IdxType lo = 0;
                IdxType hi = ((IdxType)1 << (n_qubits));
                IdxType mid;
                ValType r = uni_dist(rng);
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

        void RESET_GATE(const IdxType qubit)
        {
            IdxType mask = ((IdxType)1 << qubit);
            ValType prob_of_one = 0;
            for (IdxType i = 0; i < ((IdxType)1 << (n_qubits)); i++)
            {
                if ((i & mask) != 0)
                    prob_of_one += fabs(dm_real[(i << (n_qubits)) + i]);
            }
            assert(prob_of_one <= 1.0);

            if (prob_of_one < 1.0) // still possible to normalize
            {
                ValType factor = 1.0 / (1.0 - prob_of_one);
                for (IdxType i = 0; i < dim; i++)
                {
                    if ((i & mask) == 0)
                    {
                        dm_real[i] *= factor;
                        dm_imag[i] *= factor;
                    }
                    else
                    {
                        dm_real[i] = 0;
                        dm_imag[i] = 0;
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
                        dm_real[i] = dm_real[dual_i];
                        dm_imag[i] = dm_imag[dual_i];
                        dm_real[dual_i] = 0;
                        dm_imag[dual_i] = 0;
                    }
                }
            }
        }
    };

} // namespace NWQSim
