#pragma once

#include "../state.hpp"

#include "../util.hpp"
#include "../gate.hpp"
#include "../circuit.hpp"

#include "../circuit_pass/fusion.hpp"
#include "../macros.hpp"

#include <random>
#include <cstring>
#include <algorithm>

#include <mpi.h>

#define PRINT_SIM_TRACE

namespace NWQSim
{
    class SV_MPI : public QuantumState
    {

    public:
        SV_MPI(IdxType _n_qubits) : QuantumState(_n_qubits)
        {
            // SV initialization
            n_qubits = _n_qubits;
            dim = (IdxType)1 << (n_qubits);
            half_dim = (IdxType)1 << (n_qubits - 1);
            sv_size = dim * (IdxType)sizeof(ValType);

            // MPI initialization
            int rank, size;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &size);
            comm_global = MPI_COMM_WORLD;

            i_cpu = (IdxType)rank;
            n_cpus = (IdxType)size;

            cpu_scale = floor(log((double)n_cpus + 0.5) / log(2.0));
            lg2_m_cpu = n_qubits - cpu_scale;
            m_cpu = ((IdxType)1 << (lg2_m_cpu));

            sv_size_per_cpu = sv_size / n_cpus;

            // CPU side initialization
            if (!is_power_of_2(n_cpus))
                throw std::logic_error("Error: Number of CPU nodes should be power of 2.");
            if (dim % n_cpus != 0)
                throw std::logic_error("Error: Number of CPU nodes is too large or too small.");
            if (lg2_m_cpu < 2)
                throw std::logic_error("Error: Each CPU should have at least 2 qubits for multi-node version. Please increase qubits or reduce CPUs");

            // CPU side initialization
            SAFE_ALOC_HOST(sv_real, sv_size_per_cpu);
            SAFE_ALOC_HOST(sv_imag, sv_size_per_cpu);
            SAFE_ALOC_HOST(m_real, sv_size_per_cpu);
            SAFE_ALOC_HOST(m_imag, sv_size_per_cpu);
            memset(sv_real, 0, sv_size_per_cpu);
            memset(sv_imag, 0, sv_size_per_cpu);
            memset(m_real, 0, sv_size_per_cpu);
            memset(m_imag, 0, sv_size_per_cpu);

            // State-vector initial state [0..0] = 1
            if (i_cpu == 0)
                sv_real[0] = 1.;
            cpu_mem = sv_size_per_cpu * 4;

            rng.seed(time(0));

#ifdef PRINT_SIM_TRACE
            if (i_cpu == 0)
                printf("SVSIM MPI is initialized!\n");
#endif
        }

        ~SV_MPI()
        {
            // Release for CPU side
            SAFE_FREE_HOST(sv_real);
            SAFE_FREE_HOST(sv_imag);
            SAFE_FREE_HOST(m_real);
            SAFE_FREE_HOST(m_imag);

            SAFE_FREE_HOST(results);
            SAFE_FREE_HOST(results_local);
        }

        void reset_state() override
        {
            // Reset CPU input & output
            memset(sv_real, 0, sv_size_per_cpu);
            memset(sv_imag, 0, sv_size_per_cpu);
            memset(m_real, 0, sv_size_per_cpu);
            memset(m_imag, 0, sv_size_per_cpu);
            // State Vector initial state [0..0] = 1
            if (i_cpu == 0)
                sv_real[0] = 1.;
        }

        void set_seed(IdxType seed) override
        {
            rng.seed(seed);
        }

        void sim(Circuit *circuit) override
        {
            fuse_circuit(circuit);
            auto gates = circuit->gates;
            IdxType n_gates = gates->size();
            assert(circuit->num_qubits() == n_qubits);

#ifdef PRINT_SIM_TRACE
            double *sim_times;
            double sim_time;
            cpu_timer sim_timer;
            if (i_cpu == 0)
            {
                SAFE_ALOC_HOST(sim_times, sizeof(double) * n_cpus);
                memset(sim_times, 0, sizeof(double) * n_cpus);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            sim_timer.start_timer();
#endif
            // Run simulation
            simulation_kernel(gates);

#ifdef PRINT_SIM_TRACE
            sim_timer.stop_timer();
            sim_time = sim_timer.measure();

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Gather(&sim_time, 1, MPI_DOUBLE,
                       &sim_times[i_cpu], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            if (i_cpu == 0)
            {
                double avg_sim_time = 0;
                for (unsigned d = 0; d < n_cpus; d++)
                {
                    avg_sim_time += sim_times[d];
                }
                avg_sim_time /= (double)n_cpus;
                printf("\n============== SV-Sim ===============\n");
                printf("nqubits:%lld, sim_gates:%lld, n_nodes:%lld, sim:%.3lf ms, mem_per_node:%.3lf MB, total_mem:%.3lf MB, \n",
                       n_qubits, n_gates, n_cpus,
                       avg_sim_time, cpu_mem / 1024 / 1024,
                       n_cpus * cpu_mem / 1024 / 1024);
                printf("=====================================\n");
                SAFE_FREE_HOST(sim_times);
            }
#endif
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

        IdxType *measure_all(IdxType repetition = DEFAULT_REPETITIONS) override
        {
            MA_GATE(repetition);
            return results;
        }

        void print_res_sv() override
        {
            ValType *sv_diag_real = NULL;
            ValType *sv_diag_imag = NULL;
            if (i_cpu == 0)
                SAFE_ALOC_HOST(sv_diag_real, dim * sizeof(ValType));
            if (i_cpu == 0)
                SAFE_ALOC_HOST(sv_diag_imag, dim * sizeof(ValType));

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Gather(sv_real, m_cpu, MPI_DOUBLE,
                       &sv_diag_real[i_cpu * m_cpu], m_cpu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gather(sv_imag, m_cpu, MPI_DOUBLE,
                       &sv_diag_imag[i_cpu * m_cpu], m_cpu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            if (i_cpu == 0)
            {
                IdxType num = ((IdxType)1 << n_qubits);
                printf("----- SVSim ------\n");
                for (IdxType i = 0; i < num; i++)
                {
                    printf("(%.3lf,%.3lfj) ", sv_diag_real[i], sv_diag_imag[i]);
                    if ((i + 1) % 8 == 0)
                        printf("\n");
                }
                printf("\n");
                SAFE_FREE_HOST(sv_diag_real);
                SAFE_FREE_HOST(sv_diag_imag);
            }
        }

    protected:
        // SV variables
        IdxType n_qubits;
        IdxType sv_size_per_cpu;
        IdxType sv_size;
        IdxType dim;
        IdxType half_dim;

        // MPI variables
        IdxType i_cpu;
        MPI_Comm comm_global;
        ValType cpu_mem;
        IdxType cpu_scale;
        IdxType n_cpus;
        IdxType lg2_m_cpu;
        IdxType m_cpu;

        // CPU arrays
        ValType *sv_real;
        ValType *sv_imag;

        // For joint measurement
        ValType *m_real;
        ValType *m_imag;
        // For measurement randoms
        ValType *randoms = NULL;
        // For measurement result
        IdxType *results = NULL;
        // Local measurement result for MA
        IdxType *results_local = NULL;

        // Random
        std::mt19937 rng;
        std::uniform_real_distribution<ValType> uni_dist;

        void simulation_kernel(std::shared_ptr<std::vector<Gate>> gates)
        {

            //=========================================
            for (auto g : *gates)
            {
                // only need sync when operating on remote qubits
                if ((g.ctrl >= lg2_m_cpu) || (g.qubit >= lg2_m_cpu))
                {
                    BARR_MPI;
                }

                if (g.op_name == OP::RESET)
                {
                    RESET_GATE(g.qubit);
                }
                else if (g.op_name == OP::M)
                {
                    M_GATE(g.qubit);
                }
                else if (g.op_name == OP::MA)
                {
                    MA_GATE(g.repetation);
                }
                else if (g.n_qubits == 1)
                {
                    C1V2_GATE(g.gm_real, g.gm_imag, g.qubit);
                }
                else if (g.n_qubits == 2)
                {
                    if ((g.ctrl >= lg2_m_cpu) && (g.qubit >= lg2_m_cpu))
                    {
                        SWAP_GATE(0, g.ctrl);
                        BARR_MPI;
                        C2V1_GATE(g.gm_real, g.gm_imag, 0, g.qubit);
                        BARR_MPI;
                        SWAP_GATE(0, g.ctrl);
                    }
                    else
                    {
                        C2V1_GATE(g.gm_real, g.gm_imag, g.ctrl, g.qubit);
                    }
                }
                else
                {
                    if (i_cpu == 0)
                    {
                        {
                            std::cout << "Unrecognized gates" << std::endl
                                      << g.gateToString() << std::endl;
                            std::logic_error("Invalid gate type");
                        }
                    }
                }

                // only need sync when operating on remote qubits
                if ((g.ctrl >= lg2_m_cpu) || (g.qubit >= lg2_m_cpu))
                {
                    BARR_MPI;
                }
            }
        }

        //============== Local Unified 1-qubit Gate ================
        void C1_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit)
        {
            const IdxType per_pe_work = ((half_dim) >> (cpu_scale));
            for (IdxType i = (i_cpu)*per_pe_work; i < (i_cpu + 1) * per_pe_work; i++)
            {
                IdxType outer = (i >> qubit);
                IdxType inner = (i & (((IdxType)1 << qubit) - 1));
                IdxType offset = (outer << (qubit + 1));
                IdxType pos0 = ((offset + inner) & (m_cpu - 1));
                IdxType pos1 = ((offset + inner + ((IdxType)1 << qubit)) & (m_cpu - 1));
                const ValType el0_real = LOCAL_G(sv_real, pos0);
                const ValType el0_imag = LOCAL_G(sv_imag, pos0);
                const ValType el1_real = LOCAL_G(sv_real, pos1);
                const ValType el1_imag = LOCAL_G(sv_imag, pos1);
                ValType sv_real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag);
                ValType sv_imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real);
                ValType sv_real_pos1 = (gm_real[2] * el0_real) - (gm_imag[2] * el0_imag) + (gm_real[3] * el1_real) - (gm_imag[3] * el1_imag);
                ValType sv_imag_pos1 = (gm_real[2] * el0_imag) + (gm_imag[2] * el0_real) + (gm_real[3] * el1_imag) + (gm_imag[3] * el1_real);
                LOCAL_P(sv_real, pos0, sv_real_pos0);
                LOCAL_P(sv_imag, pos0, sv_imag_pos0);
                LOCAL_P(sv_real, pos1, sv_real_pos1);
                LOCAL_P(sv_imag, pos1, sv_imag_pos1);
            }
        }

        void C1V2_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit)
        {
            const IdxType per_pe_work = ((dim) >> (cpu_scale));
            if (qubit < lg2_m_cpu)
            {
                C1_GATE(gm_real, gm_imag, qubit);
            }
            else
            {
                IdxType pair_cpu = (i_cpu) ^ ((IdxType)1 << (qubit - (lg2_m_cpu)));
                assert(pair_cpu != i_cpu);
                // these nodes send their sv to the pairs
                if (i_cpu > pair_cpu)
                {
                    // Send own partial statevector to remote nodes
                    MPI_Send(sv_real, per_pe_work, MPI_DOUBLE, pair_cpu, 0, MPI_COMM_WORLD);
                    MPI_Send(sv_imag, per_pe_work, MPI_DOUBLE, pair_cpu, 1, MPI_COMM_WORLD);
                    // Recevive partial statevector back
                    MPI_Recv(sv_real, per_pe_work, MPI_DOUBLE, pair_cpu, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(sv_imag, per_pe_work, MPI_DOUBLE, pair_cpu, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                else
                {
                    ValType *sv_real_remote = m_real;
                    ValType *sv_imag_remote = m_imag;

                    MPI_Recv(sv_real_remote, per_pe_work, MPI_DOUBLE, pair_cpu, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(sv_imag_remote, per_pe_work, MPI_DOUBLE, pair_cpu, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                    for (IdxType i = 0; i < per_pe_work; i++)
                    {
                        const ValType el0_real = LOCAL_G(sv_real, i);
                        const ValType el0_imag = LOCAL_G(sv_imag, i);
                        const ValType el1_real = LOCAL_G(sv_real_remote, i);
                        const ValType el1_imag = LOCAL_G(sv_imag_remote, i);

                        ValType real_pos0 = (gm_real[0] * el0_real) - (gm_imag[0] * el0_imag) + (gm_real[1] * el1_real) - (gm_imag[1] * el1_imag);
                        ValType imag_pos0 = (gm_real[0] * el0_imag) + (gm_imag[0] * el0_real) + (gm_real[1] * el1_imag) + (gm_imag[1] * el1_real);
                        ValType real_pos1 = (gm_real[2] * el0_real) - (gm_imag[2] * el0_imag) + (gm_real[3] * el1_real) - (gm_imag[3] * el1_imag);
                        ValType imag_pos1 = (gm_real[2] * el0_imag) + (gm_imag[2] * el0_real) + (gm_real[3] * el1_imag) + (gm_imag[3] * el1_real);
                        LOCAL_P(sv_real, i, real_pos0);
                        LOCAL_P(sv_imag, i, imag_pos0);
                        LOCAL_P(sv_real_remote, i, real_pos1);
                        LOCAL_P(sv_imag_remote, i, imag_pos1);
                    }
                    MPI_Send(sv_real_remote, per_pe_work, MPI_DOUBLE, pair_cpu, 2, MPI_COMM_WORLD);
                    MPI_Send(sv_imag_remote, per_pe_work, MPI_DOUBLE, pair_cpu, 3, MPI_COMM_WORLD);
                }
            }
        }

        //============== Local 2-qubit Gate  ================
        void C2_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit0, const IdxType qubit1)
        {
            const IdxType per_pe_work = ((dim) >> (cpu_scale + 2));
            assert(qubit0 != qubit1); // Non-cloning
            const IdxType q0dim = ((IdxType)1 << std::max(qubit0, qubit1));
            const IdxType q1dim = ((IdxType)1 << std::min(qubit0, qubit1));
            // const IdxType outer_factor = ((dim) + q0dim + q0dim - 1) >> (std::max(qubit0, qubit1) + 1);
            const IdxType mider_factor = (q0dim + q1dim + q1dim - 1) >> (std::min(qubit0, qubit1) + 1);
            const IdxType inner_factor = q1dim;
            const IdxType qubit0_dim = ((IdxType)1 << qubit0);
            const IdxType qubit1_dim = ((IdxType)1 << qubit1);

            for (IdxType i = (i_cpu)*per_pe_work; i < (i_cpu + 1) * per_pe_work; i++)
            {
                IdxType outer = ((i / inner_factor) / (mider_factor)) * (q0dim + q0dim);
                IdxType mider = ((i / inner_factor) % (mider_factor)) * (q1dim + q1dim);
                IdxType inner = i % inner_factor;
                IdxType pos0 = outer + mider + inner;
                IdxType pos1 = outer + mider + inner + qubit1_dim;
                IdxType pos2 = outer + mider + inner + qubit0_dim;
                IdxType pos3 = outer + mider + inner + q0dim + q1dim;

                const ValType el0_real = LOCAL_G(sv_real, pos0);
                const ValType el0_imag = LOCAL_G(sv_imag, pos0);
                const ValType el1_real = LOCAL_G(sv_real, pos1);
                const ValType el1_imag = LOCAL_G(sv_imag, pos1);
                const ValType el2_real = LOCAL_G(sv_real, pos2);
                const ValType el2_imag = LOCAL_G(sv_imag, pos2);
                const ValType el3_real = LOCAL_G(sv_real, pos3);
                const ValType el3_imag = LOCAL_G(sv_imag, pos3);

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

                LOCAL_P(sv_real, pos0, sv_real_pos0);
                LOCAL_P(sv_real, pos1, sv_real_pos1);
                LOCAL_P(sv_real, pos2, sv_real_pos2);
                LOCAL_P(sv_real, pos3, sv_real_pos3);

                LOCAL_P(sv_imag, pos0, sv_imag_pos0);
                LOCAL_P(sv_imag, pos1, sv_imag_pos1);
                LOCAL_P(sv_imag, pos2, sv_imag_pos2);
                LOCAL_P(sv_imag, pos3, sv_imag_pos3);
            }
            // BARR_MPI;
        }

        //============== Unified 2-qubit Gate ================
        // Perform communication optimization here
        void C2V1_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit0, const IdxType qubit1)
        {
            assert(qubit0 != qubit1); // Non-cloning

            if (qubit0 < lg2_m_cpu && qubit1 < lg2_m_cpu)
            {
                C2_GATE(gm_real, gm_imag, qubit0, qubit1);
            }
            else
            {
                const IdxType per_pe_work = ((dim) >> (cpu_scale + 1));
                const IdxType per_pe_num = ((dim) >> (cpu_scale));
                const IdxType p = std::min(qubit0, qubit1);
                const IdxType q = std::max(qubit0, qubit1);

                // load data from pair node
                IdxType pair_cpu = (i_cpu) ^ ((IdxType)1 << (q - (lg2_m_cpu)));
                assert(pair_cpu != i_cpu);

                if (i_cpu > pair_cpu)
                {
                    // Send own partial statevector to remote nodes
                    MPI_Send(sv_real, per_pe_num, MPI_DOUBLE, pair_cpu, 0, MPI_COMM_WORLD);
                    MPI_Send(sv_imag, per_pe_num, MPI_DOUBLE, pair_cpu, 1, MPI_COMM_WORLD);
                    // Recevive partial statevector back
                    MPI_Recv(sv_real, per_pe_num, MPI_DOUBLE, pair_cpu, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(sv_imag, per_pe_num, MPI_DOUBLE, pair_cpu, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
                else
                {
                    ValType *sv_real_remote = m_real;
                    ValType *sv_imag_remote = m_imag;
                    MPI_Recv(sv_real_remote, per_pe_num, MPI_DOUBLE, pair_cpu, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(sv_imag_remote, per_pe_num, MPI_DOUBLE, pair_cpu, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    for (IdxType i = (i_cpu)*per_pe_work; i < (i_cpu + 1) * per_pe_work; i++)
                    {
                        ValType el_real[4];
                        ValType el_imag[4];
                        ValType res_real[4] = {0};
                        ValType res_imag[4] = {0};
                        const IdxType term0 = MOD2E(i, p);
                        const IdxType term1 = MOD2E(DIV2E(i, p), q - p - 1) * EXP2E(p + 1);
                        const IdxType term2 = DIV2E(DIV2E(i, p), q - p - 1) * EXP2E(q + 1);
                        const IdxType term = term2 + term1 + term0;
                        el_real[0] = LOCAL_G(sv_real, term + SV4IDX(0));
                        el_imag[0] = LOCAL_G(sv_imag, term + SV4IDX(0));
                        el_real[3] = LOCAL_G(sv_real_remote, term + SV4IDX(3));
                        el_imag[3] = LOCAL_G(sv_imag_remote, term + SV4IDX(3));

                        if (qubit0 == q) // qubit0 is the remote qubit
                        {
                            el_real[1] = LOCAL_G(sv_real, term + SV4IDX(1));
                            el_imag[1] = LOCAL_G(sv_imag, term + SV4IDX(1));
                            el_real[2] = LOCAL_G(sv_real_remote, term + SV4IDX(2));
                            el_imag[2] = LOCAL_G(sv_imag_remote, term + SV4IDX(2));
                        }
                        else // qubit1 is the remote qubit
                        {
                            el_real[1] = LOCAL_G(sv_real_remote, term + SV4IDX(1));
                            el_imag[1] = LOCAL_G(sv_imag_remote, term + SV4IDX(1));
                            el_real[2] = LOCAL_G(sv_real, term + SV4IDX(2));
                            el_imag[2] = LOCAL_G(sv_imag, term + SV4IDX(2));
                        }
                        // #pragma unroll
                        for (unsigned j = 0; j < 4; j++)
                        {
                            // #pragma unroll
                            for (unsigned k = 0; k < 4; k++)
                            {
                                res_real[j] += (el_real[k] * gm_real[j * 4 + k]) - (el_imag[k] * gm_imag[j * 4 + k]);
                                res_imag[j] += (el_real[k] * gm_imag[j * 4 + k]) + (el_imag[k] * gm_real[j * 4 + k]);
                            }
                        }
                        LOCAL_P(sv_real, term + SV4IDX(0), res_real[0]);
                        LOCAL_P(sv_imag, term + SV4IDX(0), res_imag[0]);
                        LOCAL_P(sv_real_remote, term + SV4IDX(3), res_real[3]);
                        LOCAL_P(sv_imag_remote, term + SV4IDX(3), res_imag[3]);

                        if (qubit0 == q) // qubit0 is the remote qubit
                        {
                            LOCAL_P(sv_real, term + SV4IDX(1), res_real[1]);
                            LOCAL_P(sv_imag, term + SV4IDX(1), res_imag[1]);
                            LOCAL_P(sv_real_remote, term + SV4IDX(2), res_real[2]);
                            LOCAL_P(sv_imag_remote, term + SV4IDX(2), res_imag[2]);
                        }
                        else // qubit1 is the remote qubit
                        {
                            LOCAL_P(sv_real_remote, term + SV4IDX(1), res_real[1]);
                            LOCAL_P(sv_imag_remote, term + SV4IDX(1), res_imag[1]);
                            LOCAL_P(sv_real, term + SV4IDX(2), res_real[2]);
                            LOCAL_P(sv_imag, term + SV4IDX(2), res_imag[2]);
                        }
                    }
                    MPI_Send(sv_real_remote, per_pe_num, MPI_DOUBLE, pair_cpu, 2, MPI_COMM_WORLD);
                    MPI_Send(sv_imag_remote, per_pe_num, MPI_DOUBLE, pair_cpu, 3, MPI_COMM_WORLD);
                }
            }
        }

        void M_GATE(const IdxType qubit)
        {
            SAFE_FREE_HOST(results);
            SAFE_ALOC_HOST(results, sizeof(IdxType));
            memset(results, 0, sizeof(IdxType));
            auto rand = uni_dist(rng);

            const IdxType per_pe_work = ((dim) >> (cpu_scale));
            // ValType *m_real = this->m_real;
            IdxType mask = ((IdxType)1 << qubit);

            ValType sum = 0;
            ValType prob_of_one = 0;
            for (IdxType i = 0; i < m_cpu; i++)
            {
                IdxType idx = (i_cpu)*per_pe_work + i;
                if ((idx & mask) != 0)
                    sum += sv_real[i] * sv_real[i] + sv_imag[i] * sv_imag[i];
            }
            BARR_MPI;
            MPI_Allreduce(&sum, &prob_of_one, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            BARR_MPI;
            if (rand < prob_of_one)
            {
                ValType factor = 1. / sqrt(prob_of_one);
                for (IdxType i = 0; i < m_cpu; i++)
                {
                    IdxType idx = (i_cpu)*per_pe_work + i;
                    if ((idx & mask) == 0)
                    {
                        sv_real[i] = 0;
                        sv_imag[i] = 0;
                    }
                    else
                    {
                        sv_real[i] *= factor;
                        sv_imag[i] *= factor;
                    }
                }
            }
            else
            {
                ValType factor = 1. / sqrt(1. - prob_of_one);
                for (IdxType i = 0; i < m_cpu; i++)
                {
                    IdxType idx = (i_cpu)*per_pe_work + i;
                    if ((idx & mask) == 0)
                    {
                        sv_real[i] *= factor;
                        sv_imag[i] *= factor;
                    }
                    else
                    {
                        sv_real[i] = 0;
                        sv_imag[i] = 0;
                    }
                }
            }
            results[0] = (rand <= prob_of_one ? 1 : 0);
            BARR_MPI;
        }

        // We first do a local reduction, then we do a MPI scan, then update local vector
        void MA_GATE(const IdxType repetition)
        {
            SAFE_FREE_HOST(results);
            SAFE_ALOC_HOST(results, sizeof(IdxType) * repetition);
            memset(results, 0, sizeof(IdxType) * repetition);
            SAFE_FREE_HOST(results_local);
            SAFE_ALOC_HOST(results_local, sizeof(IdxType) * repetition);
            memset(results_local, 0, sizeof(IdxType) * repetition);

            SAFE_FREE_HOST(randoms);
            SAFE_ALOC_HOST(randoms, sizeof(ValType) * repetition);
            for (IdxType i = 0; i < repetition; i++)
                randoms[i] = uni_dist(rng);

            const IdxType n_size = (IdxType)1 << (n_qubits);
            ValType *m_real = this->m_real;

            ValType reduce = 0;
            ValType partial = 0;
            for (IdxType i = 0; i < m_cpu; i++)
                reduce += sv_real[i] * sv_real[i] + sv_imag[i] * sv_imag[i];
            BARR_MPI;
            MPI_Scan(&reduce, &partial, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            if (i_cpu == (n_cpus - 1)) // last node
            {
                ValType purity = fabs(partial);
                if (fabs(purity - 1.0) > ERROR_BAR)
                    printf("MA: Purity Check fails with %lf\n", purity);
            }

            BARR_MPI;
            m_real[0] = (partial - reduce); // per-node incremental val
            for (IdxType i = 1; i < m_cpu; i++)
                m_real[i] = m_real[i - 1] + ((sv_real[i - 1] * sv_real[i - 1]) + (sv_imag[i - 1] * sv_imag[i - 1]));
            BARR_MPI;

            for (IdxType j = 0; j < n_size; j++)
            {
                IdxType local_cpu = j >> (lg2_m_cpu);
                IdxType local_j = j & (m_cpu - 1);
                if (local_cpu == i_cpu)
                {
                    ValType lower = LOCAL_G(m_real, j);
                    ValType upper = 0;
                    if (j + 1 == n_size)
                        upper = 1.0; // last element
                    else
                        upper = (local_j + 1 == m_cpu) ? partial : m_real[local_j + 1]; // last element per node

                    for (IdxType i = 0; i < repetition; i++)
                    {
                        ValType r = randoms[i];
                        // all nodes store partial results locally. since other entires are all
                        // zeros, we can do all-reduce to merge
                        if (lower <= r && r < upper)
                            results_local[i] = j;
                    }
                }
            }
            BARR_MPI;
            MPI_Allreduce(results_local, results, repetition, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        }

        void RESET_GATE(const IdxType qubit)
        {
            const IdxType per_pe_work = ((dim) >> (cpu_scale));
            ValType *m_real = this->m_real;
            IdxType mask = ((IdxType)1 << qubit);

            ValType sum = 0;
            ValType prob_of_one = 0;
            for (IdxType i = 0; i < m_cpu; i++)
            {
                IdxType idx = (i_cpu)*per_pe_work + i;
                if ((idx & mask) != 0)
                    m_real[i] = sv_real[i] * sv_real[i] + sv_imag[i] * sv_imag[i];
            }
            BARR_MPI;
            MPI_Allreduce(&sum, &prob_of_one, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            BARR_MPI;
            if (prob_of_one < 1.0) // still possible to normalize
            {
                ValType factor = 1.0 / sqrt(1.0 - prob_of_one);
                for (IdxType i = 0; i < m_cpu; i++)
                {
                    IdxType idx = (i_cpu)*per_pe_work + i;
                    if ((idx & mask) == 0)
                    {
                        sv_real[i] *= factor;
                        sv_imag[i] *= factor;
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
                if ((qubit + n_qubits) >= lg2_m_cpu) // remote qubit, need switch
                {
                    IdxType pair_cpu = (i_cpu) ^ ((IdxType)1 << (qubit - (lg2_m_cpu)));
                    assert(pair_cpu != i_cpu);
                    ValType *sv_real_remote = m_real;
                    ValType *sv_imag_remote = m_imag;

                    if (i_cpu > pair_cpu)
                    {
                        MPI_Send(sv_real, per_pe_work, MPI_DOUBLE, pair_cpu, 0, MPI_COMM_WORLD);
                        MPI_Send(sv_imag, per_pe_work, MPI_DOUBLE, pair_cpu, 1, MPI_COMM_WORLD);
                        MPI_Recv(sv_real_remote, per_pe_work, MPI_DOUBLE, pair_cpu, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(sv_imag_remote, per_pe_work, MPI_DOUBLE, pair_cpu, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                    else
                    {
                        MPI_Recv(sv_real_remote, per_pe_work, MPI_DOUBLE, pair_cpu, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        MPI_Recv(sv_imag_remote, per_pe_work, MPI_DOUBLE, pair_cpu, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                        MPI_Send(sv_real, per_pe_work, MPI_DOUBLE, pair_cpu, 2, MPI_COMM_WORLD);
                        MPI_Send(sv_imag, per_pe_work, MPI_DOUBLE, pair_cpu, 3, MPI_COMM_WORLD);
                    }
                    BARR_MPI;
                    memcpy(sv_real, sv_real_remote, per_pe_work * sizeof(ValType));
                    memcpy(sv_imag, sv_imag_remote, per_pe_work * sizeof(ValType));
                }
                else // local (ensuring dual_i is local)
                {
                    for (IdxType i = 0; i < per_pe_work; i++)
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
            BARR_MPI;
        }

        //============== SWAP Gate ================
        // This gate is for internal usage. It is used
        // when ctrl and target qubits are remote qubis, we then
        // swap one of them to a local qubit,
        // perform the C2 gate, and then swap back
        // It is assumed qubit0 is local, qubit1 is remote
        void SWAP_GATE(const IdxType qubit0, const IdxType qubit1)
        {
            assert(qubit0 != qubit1); // Non-cloning
            // For processing a remote qubit, half of the nodes will be idle for
            // better communication efficiency. Depending on qubit position,
            // for a node pair, we only use nodes with smaller ids.
            // Consequently, each node should take double of the workload than before
            // Therefore, here it is cpu_scale+1 not cpu_scale+2
            const IdxType per_pe_work = ((dim) >> (cpu_scale + 1));
            const IdxType per_pe_num = ((dim) >> (cpu_scale));
            const IdxType p = std::min(qubit0, qubit1);
            const IdxType q = std::max(qubit0, qubit1);

            // load data from pair node
            IdxType pair_cpu = (i_cpu) ^ ((IdxType)1 << (q - (lg2_m_cpu)));
            assert(pair_cpu != i_cpu);

            if (i_cpu > pair_cpu)
            {
                // Send own partial statevector to remote nodes
                MPI_Send(sv_real, per_pe_num, MPI_DOUBLE, pair_cpu, 0, MPI_COMM_WORLD);
                MPI_Send(sv_imag, per_pe_num, MPI_DOUBLE, pair_cpu, 1, MPI_COMM_WORLD);
                // Recevive partial statevector back
                MPI_Recv(sv_real, per_pe_num, MPI_DOUBLE, pair_cpu, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(sv_imag, per_pe_num, MPI_DOUBLE, pair_cpu, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            else
            {
                ValType *sv_real_remote = m_real;
                ValType *sv_imag_remote = m_imag;
                MPI_Recv(sv_real_remote, per_pe_num, MPI_DOUBLE, pair_cpu, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(sv_imag_remote, per_pe_num, MPI_DOUBLE, pair_cpu, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for (IdxType i = (i_cpu)*per_pe_work; i < (i_cpu + 1) * per_pe_work; i++)
                {
                    ValType el_real[4];
                    ValType el_imag[4];
                    ValType res_real[4] = {0};
                    ValType res_imag[4] = {0};
                    const IdxType term0 = MOD2E(i, p);
                    const IdxType term1 = MOD2E(DIV2E(i, p), q - p - 1) * EXP2E(p + 1);
                    const IdxType term2 = DIV2E(DIV2E(i, p), q - p - 1) * EXP2E(q + 1);
                    const IdxType term = term2 + term1 + term0;

                    el_real[1] = LOCAL_G(sv_real_remote, term + SV4IDX(1));
                    el_imag[1] = LOCAL_G(sv_imag_remote, term + SV4IDX(1));
                    el_real[2] = LOCAL_G(sv_real, term + SV4IDX(2));
                    el_imag[2] = LOCAL_G(sv_imag, term + SV4IDX(2));

                    res_real[1] = el_real[2];
                    res_imag[1] = el_imag[2];
                    res_real[2] = el_real[1];
                    res_imag[2] = el_imag[1];

                    LOCAL_P(sv_real_remote, term + SV4IDX(1), res_real[1]);
                    LOCAL_P(sv_imag_remote, term + SV4IDX(1), res_imag[1]);
                    LOCAL_P(sv_real, term + SV4IDX(2), res_real[2]);
                    LOCAL_P(sv_imag, term + SV4IDX(2), res_imag[2]);
                }
                MPI_Send(sv_real_remote, per_pe_num, MPI_DOUBLE, pair_cpu, 2, MPI_COMM_WORLD);
                MPI_Send(sv_imag_remote, per_pe_num, MPI_DOUBLE, pair_cpu, 3, MPI_COMM_WORLD);
            }
            // BARR_MPI;
        }
    };

} // namespace NWQSim
