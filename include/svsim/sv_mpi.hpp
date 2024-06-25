#pragma once

#include "../state.hpp"

#include "../nwq_util.hpp"
#include "../gate.hpp"
#include "../circuit.hpp"

#include "../circuit_pass/fusion.hpp"
#include "../private/macros.hpp"
#include "../private/sim_gate.hpp"
#include "../private/config.hpp"
#include "../private/exp_gate_declarations_host.hpp"
#include <random>
#include <cstring>
#include <algorithm>

#include <mpi.h>

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

            i_proc = (IdxType)rank;
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
            if (i_proc == 0)
                sv_real[0] = 1.;
            cpu_mem = sv_size_per_cpu * 4;

            rng.seed(time(0));

            if (i_proc == 0 && Config::PRINT_SIM_TRACE)
                printf("SVSIM MPI is initialized!\n");
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
            if (i_proc == 0)
                sv_real[0] = 1.;
        }

        void set_seed(IdxType seed) override
        {
            rng.seed(seed);
        }

        void sim(std::shared_ptr<NWQSim::Circuit> circuit) override
        {
            IdxType origional_gates = circuit->num_gates();

            std::vector<SVGate> gates = fuse_circuit_sv(circuit);

            IdxType n_gates = gates.size();
            assert(circuit->num_qubits() == n_qubits);

            double *sim_times;
            double sim_time;
            cpu_timer sim_timer;

            if (Config::PRINT_SIM_TRACE)
            {
                if (i_proc == 0)
                {
                    SAFE_ALOC_HOST(sim_times, sizeof(double) * n_cpus);
                    memset(sim_times, 0, sizeof(double) * n_cpus);
                }
                MPI_Barrier(MPI_COMM_WORLD);
                sim_timer.start_timer();
            }

            // Run simulation
            simulation_kernel(gates);

            if (Config::PRINT_SIM_TRACE)
            {
                sim_timer.stop_timer();
                sim_time = sim_timer.measure();

                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Gather(&sim_time, 1, MPI_DOUBLE,
                           &sim_times[i_proc], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                if (i_proc == 0)
                {
                    double avg_sim_time = 0;
                    for (unsigned d = 0; d < n_cpus; d++)
                    {
                        avg_sim_time += sim_times[d];
                    }
                    avg_sim_time /= (double)n_cpus;
                    printf("\n============== SV-Sim ===============\n");
                    printf("nqubits:%lld, n_gates:%lld, sim_gates:%lld, n_nodes:%lld, sim:%.3lf ms, mem_per_node:%.3lf MB, total_mem:%.3lf MB, \n",
                           n_qubits, origional_gates, n_gates, n_cpus,
                           avg_sim_time, cpu_mem / 1024 / 1024,
                           n_cpus * cpu_mem / 1024 / 1024);
                    printf("=====================================\n");
                    SAFE_FREE_HOST(sim_times);
                }
            }
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

        // This function calculates the expectation value for specified qubits in the MPI environment
        ValType get_exp_z(const std::vector<size_t> &in_bits) override
        {
            // Determine the start and end index for this node
            IdxType offset = (dim / n_cpus) * i_proc;

            // Calculate the expectation value for this node
            double local_result = 0.0;
            for (IdxType i = 0; i < m_cpu; ++i)
            {
                local_result += (hasEvenParity(i + offset, in_bits) ? 1.0 : -1.0) *
                                (sv_real[i] * sv_real[i] + sv_imag[i] * sv_imag[i]);
            }

            // Gather the results from all nodes at the root
            double result = 0.0;

            MPI_Reduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            // Return the final result
            return result;
        }

        // This function calculates the expectation value for all qubits in the MPI environment
        ValType get_exp_z() override
        {

            // Determine the start and end index for this node
            IdxType offset = (dim / n_cpus) * i_proc;

            // Calculate the expectation value for this node
            double local_result = 0.0;
            for (IdxType i = 0; i < m_cpu; ++i)
            {
                bool parity = __builtin_parity(i + offset);
                local_result += (parity ? -1.0 : 1.0) * (sv_real[i] * sv_real[i] + sv_imag[i] * sv_imag[i]);
            }

            // Gather the results from all nodes at the root
            double result = 0.0;

            MPI_Reduce(&local_result, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            // Return the final result
            return result;
        }

        void print_res_state() override
        {
            ValType *sv_diag_real = NULL;
            ValType *sv_diag_imag = NULL;
            if (i_proc == 0)
                SAFE_ALOC_HOST(sv_diag_real, dim * sizeof(ValType));
            if (i_proc == 0)
                SAFE_ALOC_HOST(sv_diag_imag, dim * sizeof(ValType));

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Gather(sv_real, m_cpu, MPI_DOUBLE,
                       &sv_diag_real[i_proc * m_cpu], m_cpu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Gather(sv_imag, m_cpu, MPI_DOUBLE,
                       &sv_diag_imag[i_proc * m_cpu], m_cpu, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            if (i_proc == 0)
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

        void simulation_kernel(std::vector<SVGate> gates)
        {

            auto start = std::chrono::steady_clock::now();
            int n_gates = gates.size();
            for (int i = 0; i < n_gates; i++)
            {

                if (Config::PRINT_SIM_TRACE && i_proc == 0)
                    printProgressBar(i, n_gates, start);
            
                auto g = gates[i];

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
                    MA_GATE(g.qubit);
                }
                else if (g.op_name == OP::C1)
                {
                    C1V2_GATE(g.gm_real, g.gm_imag, g.qubit);
                }
                else if (g.op_name == OP::C2)
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
                else if (g.op_name == OP::EXPECT)
                {
                    BARR_MPI;
                    ObservableList o = *(ObservableList*)(g.data);
                    IdxType* xinds = o.x_indices;
                    for (IdxType obs_ind = 0; obs_ind < o.numterms; obs_ind++) {
                        EXPECT_GATE(xinds, 
                                    o.x_index_sizes[obs_ind],
                                    o.xmasks[obs_ind],
                                    o.zmasks[obs_ind],
                                    o.exp_output,
                                    obs_ind);
                        xinds += o.x_index_sizes[obs_ind];
                    }
                }
                else
                {
                    if (i_proc == 0)
                    {
                        {
                            // std::cout << "Unrecognized gates" << std::endl
                            //           << g.gateToString() << std::endl;
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
            for (IdxType i = (i_proc)*per_pe_work; i < (i_proc + 1) * per_pe_work; i++)
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
                IdxType pair_cpu = (i_proc) ^ ((IdxType)1 << (qubit - (lg2_m_cpu)));
                assert(pair_cpu != i_proc);
                // these nodes send their sv to the pairs
                if (i_proc > pair_cpu)
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

            for (IdxType i = (i_proc)*per_pe_work; i < (i_proc + 1) * per_pe_work; i++)
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
                IdxType pair_cpu = (i_proc) ^ ((IdxType)1 << (q - (lg2_m_cpu)));
                assert(pair_cpu != i_proc);

                if (i_proc > pair_cpu)
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
                    
                    IdxType index = (i_proc >> (q - (lg2_m_cpu) + 1)) << (q - (lg2_m_cpu));
                    index |= i_proc & ((1 << (q - (lg2_m_cpu))) - 1);
                    for (IdxType i = (index)*per_pe_work; i < (index + 1) * per_pe_work; i++)
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


        virtual double EXPECT_C4_GATE(const ValType* gm_real, const ValType* gm_imag, IdxType qubit0, IdxType qubit1, IdxType qubit2, IdxType qubit3, IdxType mask) {
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
            ValType exp_val = 0.0;
            const IdxType per_pe_work = ((dim) >> (cpu_scale + 4));
            
            for (IdxType i = (i_proc)*per_pe_work; i < (i_proc + 1) * per_pe_work; i++)
            {
                const IdxType term0 = MOD2E(i, p);
                const IdxType term1 = MOD2E(DIV2E(i, p), q - p - 1) * EXP2E(p + 1);
                const IdxType term2 = MOD2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1) * EXP2E(q + 1);
                const IdxType term3 = MOD2E(DIV2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(r + 1);
                const IdxType term4 = DIV2E(DIV2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(s + 1);
                const IdxType term = term4 + term3 + term2 + term1 + term0;
                const ValType el_real[16] = {
                    LOCAL_G(sv_real, term + SV16IDX(0)), LOCAL_G(sv_real, term + SV16IDX(1)),
                    LOCAL_G(sv_real, term + SV16IDX(2)), LOCAL_G(sv_real, term + SV16IDX(3)),
                    LOCAL_G(sv_real, term + SV16IDX(4)), LOCAL_G(sv_real, term + SV16IDX(5)),
                    LOCAL_G(sv_real, term + SV16IDX(6)), LOCAL_G(sv_real, term + SV16IDX(7)),
                    LOCAL_G(sv_real, term + SV16IDX(8)), LOCAL_G(sv_real, term + SV16IDX(9)),
                    LOCAL_G(sv_real, term + SV16IDX(10)), LOCAL_G(sv_real, term + SV16IDX(11)),
                    LOCAL_G(sv_real, term + SV16IDX(12)), LOCAL_G(sv_real, term + SV16IDX(13)),
                    LOCAL_G(sv_real, term + SV16IDX(14)), LOCAL_G(sv_real, term + SV16IDX(15))};
                const ValType el_imag[16] = {
                    LOCAL_G(sv_imag, term + SV16IDX(0)), LOCAL_G(sv_imag, term + SV16IDX(1)),
                    LOCAL_G(sv_imag, term + SV16IDX(2)), LOCAL_G(sv_imag, term + SV16IDX(3)),
                    LOCAL_G(sv_imag, term + SV16IDX(4)), LOCAL_G(sv_imag, term + SV16IDX(5)),
                    LOCAL_G(sv_imag, term + SV16IDX(6)), LOCAL_G(sv_imag, term + SV16IDX(7)),
                    LOCAL_G(sv_imag, term + SV16IDX(8)), LOCAL_G(sv_imag, term + SV16IDX(9)),
                    LOCAL_G(sv_imag, term + SV16IDX(10)), LOCAL_G(sv_imag, term + SV16IDX(11)),
                    LOCAL_G(sv_imag, term + SV16IDX(12)), LOCAL_G(sv_imag, term + SV16IDX(13)),
                    LOCAL_G(sv_imag, term + SV16IDX(14)), LOCAL_G(sv_imag, term + SV16IDX(15))};
                // #pragma unroll
                for (unsigned j = 0; j < 16; j++)
                {
                    ValType res_real = 0;
                    ValType res_imag = 0;
                    // #pragma unroll
                    for (unsigned k = 0; k < 16; k++)
                    {
                        res_real += (el_real[k] * gm_real[j * 16 + k]) - (el_imag[k] * gm_imag[j * 16 + k]);
                        res_imag += (el_real[k] * gm_imag[j * 16 + k]) + (el_imag[k] * gm_real[j * 16 + k]);
                    }

                    ValType val = res_real * res_real + res_imag * res_imag;
                
                    exp_val += hasEvenParity((term + SV16IDX(j)) & mask, n_qubits) ? val : -val;
                    // printf("val %f %lld\n", val, term + SV16IDX(j));
            
                }
            }
            return exp_val;
        }
        ValType EXPECT_C4V1_GATE(const ValType* gm_real, const ValType* gm_imag, IdxType qubit0, IdxType qubit1, IdxType qubit2, IdxType qubit3, IdxType mask) {
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
            const IdxType per_pe_work = ((dim) >> (cpu_scale + 3));
            const IdxType per_pe_num = ((dim) >> (cpu_scale));
            if (s < lg2_m_cpu)
            {
                return EXPECT_C4_GATE(gm_real, gm_imag, qubit0, qubit1, qubit2, qubit3, mask);
                
            } else {
                ValType exp_val = 0.0;

                // load data from pair node
                IdxType pair_cpu = (i_proc) ^ ((IdxType)1 << (s - (lg2_m_cpu)));
                assert(pair_cpu != i_proc);
                if (i_proc > pair_cpu)
                {
                    // Send own partial statevector to remote nodes
                    MPI_Send(sv_real, per_pe_num, MPI_DOUBLE, pair_cpu, 0, MPI_COMM_WORLD);
                    MPI_Send(sv_imag, per_pe_num, MPI_DOUBLE, pair_cpu, 1, MPI_COMM_WORLD);
                }
                else
                {

                    IdxType index = (i_proc >> (s - (lg2_m_cpu) + 1)) << (s - (lg2_m_cpu));
                    index |= i_proc & ((1 << (s - (lg2_m_cpu))) - 1);
                    // printf("%d %d %d %d\n", i_proc, (index)*per_pe_work, (index + 1) * per_pe_work, dim);
                    ValType *sv_real_remote = m_real;
                    ValType *sv_imag_remote = m_imag;
                    MPI_Recv(sv_real_remote, per_pe_num, MPI_DOUBLE, pair_cpu, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(sv_imag_remote, per_pe_num, MPI_DOUBLE, pair_cpu, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    std::vector<bool> markers;
                    if (i_proc == 0) {
                        markers.resize(per_pe_num);
                    }
                    
                    for (IdxType i = (index)*per_pe_work; i < (index + 1) * per_pe_work; i++)
                    {
                        ValType el_real[16];
                        ValType el_imag[16];
                        const IdxType term0 = MOD2E(i, p);
                        const IdxType term1 = MOD2E(DIV2E(i, p), q - p - 1) * EXP2E(p + 1);
                        const IdxType term2 = MOD2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1) * EXP2E(q + 1);
                        const IdxType term3 = MOD2E(DIV2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(r + 1);
                        const IdxType term4 = DIV2E(DIV2E(DIV2E(DIV2E(i, p), q - p - 1), r - q - 1), s - r - 1) * EXP2E(s + 1);
                        const IdxType term = term4 + term3 + term2 + term1 + term0;
                        if (qubit3 == s) // qubit3 is remote qubit
                        {
                            
                            // printf("Qubit 3");
                            el_real[0] = LOCAL_G(sv_real, term + SV16IDX(0));
                            el_real[1] = LOCAL_G(sv_real_remote, term + SV16IDX(1));
                            el_real[2] = LOCAL_G(sv_real, term + SV16IDX(2));
                            el_real[3] = LOCAL_G(sv_real_remote, term + SV16IDX(3));
                            el_real[4] = LOCAL_G(sv_real, term + SV16IDX(4));
                            el_real[5] = LOCAL_G(sv_real_remote, term + SV16IDX(5));
                            el_real[6] = LOCAL_G(sv_real, term + SV16IDX(6));
                            el_real[7] = LOCAL_G(sv_real_remote, term + SV16IDX(7));
                            el_real[8] = LOCAL_G(sv_real, term + SV16IDX(8));
                            el_real[9] = LOCAL_G(sv_real_remote, term + SV16IDX(9));
                            el_real[10] = LOCAL_G(sv_real, term + SV16IDX(10));
                            el_real[11] = LOCAL_G(sv_real_remote, term + SV16IDX(11));
                            el_real[12] = LOCAL_G(sv_real, term + SV16IDX(12));
                            el_real[13] = LOCAL_G(sv_real_remote, term + SV16IDX(13));
                            el_real[14] = LOCAL_G(sv_real, term + SV16IDX(14));
                            el_real[15] = LOCAL_G(sv_real_remote, term + SV16IDX(15));

                            el_imag[0] = LOCAL_G(sv_imag, term + SV16IDX(0));
                            el_imag[1] = LOCAL_G(sv_imag_remote, term + SV16IDX(1));
                            el_imag[2] = LOCAL_G(sv_imag, term + SV16IDX(2));
                            el_imag[3] = LOCAL_G(sv_imag_remote, term + SV16IDX(3));
                            el_imag[4] = LOCAL_G(sv_imag, term + SV16IDX(4));
                            el_imag[5] = LOCAL_G(sv_imag_remote, term + SV16IDX(5));
                            el_imag[6] = LOCAL_G(sv_imag, term + SV16IDX(6));
                            el_imag[7] = LOCAL_G(sv_imag_remote, term + SV16IDX(7));
                            el_imag[8] = LOCAL_G(sv_imag, term + SV16IDX(8));
                            el_imag[9] = LOCAL_G(sv_imag_remote, term + SV16IDX(9));
                            el_imag[10] = LOCAL_G(sv_imag, term + SV16IDX(10));
                            el_imag[11] = LOCAL_G(sv_imag_remote, term + SV16IDX(11));
                            el_imag[12] = LOCAL_G(sv_imag, term + SV16IDX(12));
                            el_imag[13] = LOCAL_G(sv_imag_remote, term + SV16IDX(13));
                            el_imag[14] = LOCAL_G(sv_imag, term + SV16IDX(14));
                            el_imag[15] = LOCAL_G(sv_imag_remote, term + SV16IDX(15));
                        }
                        else // qubit2 is remote (not possible qubit0 or 1 is remote)
                        {
                            // if (( (laneid>>1)&1)!=0) el_real_s[j*16+laneid] = LOCAL_G(sv_real_remote,addr);
                            // else el_real_s[j*16+laneid] = LOCAL_G(sv_real, addr);
                            // laneid = 2,3,6,7,10,11,14,15;

                            el_real[0] = LOCAL_G(sv_real, term + SV16IDX(0));
                            el_real[1] = LOCAL_G(sv_real, term + SV16IDX(1));
                            el_real[2] = LOCAL_G(sv_real_remote, term + SV16IDX(2));
                            el_real[3] = LOCAL_G(sv_real_remote, term + SV16IDX(3));
                            el_real[4] = LOCAL_G(sv_real, term + SV16IDX(4));
                            el_real[5] = LOCAL_G(sv_real, term + SV16IDX(5));
                            el_real[6] = LOCAL_G(sv_real_remote, term + SV16IDX(6));
                            el_real[7] = LOCAL_G(sv_real_remote, term + SV16IDX(7));
                            el_real[8] = LOCAL_G(sv_real, term + SV16IDX(8));
                            el_real[9] = LOCAL_G(sv_real, term + SV16IDX(9));
                            el_real[10] = LOCAL_G(sv_real_remote, term + SV16IDX(10));
                            el_real[11] = LOCAL_G(sv_real_remote, term + SV16IDX(11));
                            el_real[12] = LOCAL_G(sv_real, term + SV16IDX(12));
                            el_real[13] = LOCAL_G(sv_real, term + SV16IDX(13));
                            el_real[14] = LOCAL_G(sv_real_remote, term + SV16IDX(14));
                            el_real[15] = LOCAL_G(sv_real_remote, term + SV16IDX(15));

                            el_imag[0] = LOCAL_G(sv_imag, term + SV16IDX(0));
                            el_imag[1] = LOCAL_G(sv_imag, term + SV16IDX(1));
                            el_imag[2] = LOCAL_G(sv_imag_remote, term + SV16IDX(2));
                            el_imag[3] = LOCAL_G(sv_imag_remote, term + SV16IDX(3));
                            el_imag[4] = LOCAL_G(sv_imag, term + SV16IDX(4));
                            el_imag[5] = LOCAL_G(sv_imag, term + SV16IDX(5));
                            el_imag[6] = LOCAL_G(sv_imag_remote, term + SV16IDX(6));
                            el_imag[7] = LOCAL_G(sv_imag_remote, term + SV16IDX(7));
                            el_imag[8] = LOCAL_G(sv_imag, term + SV16IDX(8));
                            el_imag[9] = LOCAL_G(sv_imag, term + SV16IDX(9));
                            el_imag[10] = LOCAL_G(sv_imag_remote, term + SV16IDX(10));
                            el_imag[11] = LOCAL_G(sv_imag_remote, term + SV16IDX(11));
                            el_imag[12] = LOCAL_G(sv_imag, term + SV16IDX(12));
                            el_imag[13] = LOCAL_G(sv_imag, term + SV16IDX(13));
                            el_imag[14] = LOCAL_G(sv_imag_remote, term + SV16IDX(14));
                            el_imag[15] = LOCAL_G(sv_imag_remote, term + SV16IDX(15));
                        }
                        
                        // #pragma unroll
                        for (unsigned j = 0; j < 16; j++)
                        {
                            ValType res_real = 0.0;
                            ValType res_imag = 0.0;
                            // #pragma unroll
                            for (unsigned k = 0; k < 16; k++)
                            {
                                res_real += (el_real[k] * gm_real[j * 16 + k]) - (el_imag[k] * gm_imag[j * 16 + k]);
                                res_imag += (el_real[k] * gm_imag[j * 16 + k]) + (el_imag[k] * gm_real[j * 16 + k]);
                            }
                            ValType val = res_real * res_real + res_imag * res_imag;

                            exp_val += hasEvenParity((term + SV16IDX(j)) & mask, n_qubits) ? val : -val;
                            
                        }
                       
                    }
                }
                BARR_MPI;
                return exp_val;
            }
        }


        //============== Local 2-qubit Expectation Gate  ================
        ValType EXPECT_C2_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit0, const IdxType qubit1, const IdxType mask)
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
            const IdxType p = qubit0;
            const IdxType q = qubit1;
            ValType exp_val = 0.0;
            for (IdxType i = (i_proc)*per_pe_work; i < (i_proc + 1) * per_pe_work; i++)
            {
                // IdxType outer = ((i / inner_factor) / (mider_factor)) * (q0dim + q0dim);
                // IdxType mider = ((i / inner_factor) % (mider_factor)) * (q1dim + q1dim);
                // IdxType inner = i % inner_factor;
                // IdxType pos0 = outer + mider + inner;
                // IdxType pos1 = outer + mider + inner + qubit1_dim;
                // IdxType pos2 = outer + mider + inner + qubit0_dim;
                // IdxType pos3 = outer + mider + inner + q0dim + q1dim;
                ValType el_real[4];
                ValType el_imag[4];
                ValType res_real[4] = {0};
                ValType res_imag[4] = {0};
                const IdxType term0 = MOD2E(i, p);
                const IdxType term1 = MOD2E(DIV2E(i, p), q - p - 1) * EXP2E(p + 1);
                const IdxType term2 = DIV2E(DIV2E(i, p), q - p - 1) * EXP2E(q + 1);
                const IdxType term = term2 + term1 + term0;

                const ValType el0_real = LOCAL_G(sv_real, term + SV4IDX(0));
                const ValType el0_imag = LOCAL_G(sv_imag, term + SV4IDX(0));
                const ValType el1_real = LOCAL_G(sv_real, term + SV4IDX(1));
                const ValType el1_imag = LOCAL_G(sv_imag, term + SV4IDX(1));
                const ValType el2_real = LOCAL_G(sv_real, term + SV4IDX(2));
                const ValType el2_imag = LOCAL_G(sv_imag, term + SV4IDX(2));
                const ValType el3_real = LOCAL_G(sv_real, term + SV4IDX(3));
                const ValType el3_imag = LOCAL_G(sv_imag, term + SV4IDX(4));

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

                ValType v0 = sv_real_pos0 * sv_real_pos0 + sv_imag_pos0 * sv_imag_pos0;
                ValType v1 = sv_real_pos1 * sv_real_pos1 + sv_imag_pos1 * sv_imag_pos1;
                ValType v2 = sv_real_pos2 * sv_real_pos2 + sv_imag_pos2 * sv_imag_pos2;
                ValType v3 = sv_real_pos3 * sv_real_pos3 + sv_imag_pos3 * sv_imag_pos3;

                exp_val += hasEvenParity((term + SV4IDX(0)) & mask, n_qubits) ? v0 : -v0;
                exp_val += hasEvenParity((term + SV4IDX(1)) & mask, n_qubits) ? v1 : -v1;
                exp_val += hasEvenParity((term + SV4IDX(2)) & mask, n_qubits) ? v2 : -v2;
                exp_val += hasEvenParity((term + SV4IDX(3)) & mask, n_qubits) ? v3 : -v3;
            }
            // BARR_MPI;
            return exp_val;
        }

        //============== Unified 2-qubit Gate ================
        // Perform communication optimization here
        ValType EXPECT_C2V1_GATE(const ValType *gm_real, const ValType *gm_imag, const IdxType qubit0, const IdxType qubit1, IdxType mask)
        {
            assert(qubit0 != qubit1); // Non-cloning

            if (qubit0 < lg2_m_cpu && qubit1 < lg2_m_cpu)
            {
                return EXPECT_C2_GATE(gm_real, gm_imag, qubit0, qubit1, mask);
            }
            else
            {
                ValType exp_val = 0.0;
                const IdxType per_pe_work = ((dim) >> (cpu_scale + 1));
                const IdxType per_pe_num = ((dim) >> (cpu_scale));
                const IdxType p = std::min(qubit0, qubit1);
                const IdxType q = std::max(qubit0, qubit1);

                // load data from pair node
                IdxType pair_cpu = (i_proc) ^ ((IdxType)1 << (q - (lg2_m_cpu)));
                assert(pair_cpu != i_proc);
                if (i_proc > pair_cpu)
                {
                    // Send own partial statevector to remote nodes
                    MPI_Send(sv_real, per_pe_num, MPI_DOUBLE, pair_cpu, 0, MPI_COMM_WORLD);
                    MPI_Send(sv_imag, per_pe_num, MPI_DOUBLE, pair_cpu, 1, MPI_COMM_WORLD);
                    // Recieve partial statevector back
                    return 0.0;
                }
                else
                {
                    IdxType index = (i_proc >> (q - (lg2_m_cpu) + 1)) << (q - (lg2_m_cpu));
                    index |= i_proc & ((1 << (q - (lg2_m_cpu))) - 1);
                    // printf("%d %d %d %d\n", i_proc, (index)*per_pe_work, (index + 1) * per_pe_work, dim);
                    ValType *sv_real_remote = m_real;
                    ValType *sv_imag_remote = m_imag;
                    MPI_Recv(sv_real_remote, per_pe_num, MPI_DOUBLE, pair_cpu, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(sv_imag_remote, per_pe_num, MPI_DOUBLE, pair_cpu, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    
                    for (IdxType i = index * per_pe_work; i < (index + 1) * per_pe_work; i++)
                    {
                        
                        ValType el_real[4];
                        ValType el_imag[4];
                        ValType res_real[4] = {0};
                        ValType res_imag[4] = {0};
                        const IdxType term0 = MOD2E(i, p);
                        const IdxType term1 = MOD2E(DIV2E(i, p), q - p - 1) * EXP2E(p + 1);
                        const IdxType term2 = DIV2E(DIV2E(i, p), q - p - 1) * EXP2E(q + 1);
                        const IdxType term = term2 + term1 + term0;
                        IdxType indices[4] = {
                            term + SV4IDX(0),
                            term + SV4IDX(1),
                            term + SV4IDX(2),
                            term + SV4IDX(3)
                        };
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
                            ValType r_real = 0.0;
                            ValType r_imag = 0.0;
                            // #pragma unroll
                            for (unsigned k = 0; k < 4; k++)
                            {
                                r_real += (el_real[k] * gm_real[j * 4 + k]) - (el_imag[k] * gm_imag[j * 4 + k]);
                                r_imag += (el_real[k] * gm_imag[j * 4 + k]) + (el_imag[k] * gm_real[j * 4 + k]);
                            }
                            ValType val = r_real * r_real + r_imag * r_imag;
                            exp_val += hasEvenParity((term + SV4IDX(j)) & mask, n_qubits) ? val: -val;
                        }
                    }
                    return exp_val;
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
                IdxType idx = (i_proc)*per_pe_work + i;
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
                    IdxType idx = (i_proc)*per_pe_work + i;
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
                    IdxType idx = (i_proc)*per_pe_work + i;
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
        double EXPECT_C0_GATE(IdxType zmask) {
            const IdxType per_pe_work = ((dim) >> (cpu_scale));
            ValType exp_val = 0;
            for (IdxType i = (i_proc)*per_pe_work; i < (i_proc + 1) * per_pe_work; i++)
            {
                ValType val = LOCAL_G(sv_real, i) * LOCAL_G(sv_real, i) + LOCAL_G(sv_imag, i) * LOCAL_G(sv_imag, i);
                exp_val += hasEvenParity(i & zmask, n_qubits) ? val : -val;
            }
            return exp_val;
            
        }
        void EXPECT_GATE(IdxType* x_indices, 
                                 IdxType num_x_indices, 
                                 IdxType xmask, 
                                 IdxType zmask, 
                                 ValType* output,
                                 IdxType output_index)  {
            ValType result = 0.0;
            BARR_MPI;
            if (num_x_indices == 2) {
                IdxType q0 = x_indices[0];
                IdxType q1 = x_indices[1];
                IdxType q0t = q0;

                if (q0 >= lg2_m_cpu) {
                    SWAP_GATE(0, q0);
                    q0t = 0;
                    xmask = swapBits(xmask, q0t, q0);
                    zmask = swapBits(zmask, q0t, q0);
                    BARR_MPI;
                    assert(q0t != q1);
                }
                IdxType zind0 = ((zmask & (1 << q0t)) >> q0t);
                IdxType zind1 = ((zmask & (1 << q1)) >> q1) << 1;
                const ValType* gm_real = exp_gate_perms_2q[zind0 + zind1];
                const ValType* gm_imag = exp_gate_perms_2q[zind0 + zind1] + 16;
                result = EXPECT_C2V1_GATE(gm_real, gm_imag, q0t, q1, xmask | zmask);
                if (q0 >= lg2_m_cpu) {
                    BARR_MPI;
                    SWAP_GATE(0, q0);
                    BARR_MPI;
                }
            } else if (num_x_indices == 4) {
                IdxType q0 = x_indices[0];
                IdxType q1 = x_indices[1];
                IdxType q2 = x_indices[2];
                IdxType q3 = x_indices[3];
                IdxType q0t = q0;
                IdxType q1t = q1;
                IdxType q2t = q2;
                
                IdxType local_index = 0;
                // assume the indices are sorted
                if (q0 >= lg2_m_cpu) {
                    q0t = local_index++;
                    while((q0t == q1t || q0t == q2t) && q0t < lg2_m_cpu) {
                        q0t++;
                    }
                    SWAP_GATE(q0t, q0);
                    zmask = (IdxType)swapBits(zmask, (uint64_t)q0t, (uint64_t)q0);
                    xmask = (IdxType)swapBits(xmask, (uint64_t)q0t, (uint64_t)q0);
                    BARR_MPI;
                }
                if (q1 >= lg2_m_cpu) {

                    q1t = local_index++;
                    while((q1t == q0t || q1t == q2t) && q1t < lg2_m_cpu) {
                        q1t++;
                    }
                    SWAP_GATE(q1t, q1);
                    zmask = (IdxType)swapBits(zmask, (uint64_t)q1t, (uint64_t)q1);
                    xmask = (IdxType)swapBits(xmask, (uint64_t)q1t, (uint64_t)q1);
                    BARR_MPI;
                }
                if (q2 >= lg2_m_cpu) {
                    q2t = local_index++;

                    while((q2t == q0t || q2t == q1t) && q2t < lg2_m_cpu) {
                        q2t++;
                    }
                    SWAP_GATE(q2t, q2);
                    zmask = (IdxType)swapBits(zmask, (uint64_t)q2t, (uint64_t)q2);
                    xmask = (IdxType)swapBits(xmask, (uint64_t)q2t, (uint64_t)q2);
                    BARR_MPI;
                }
                assert (q0t < lg2_m_cpu && q1t < lg2_m_cpu && q2t < lg2_m_cpu);

                const IdxType v0 = std::min(q0t, q1t);
                const IdxType v1 = std::min(q2t, q3);
                const IdxType v2 = std::max(q0t, q1t);
                const IdxType v3 = std::max(q2t, q3);
                const IdxType p = std::min(v0, v1);
                const IdxType q = std::min(std::min(v2, v3), std::max(v0, v1));
                const IdxType r = std::max(std::min(v2, v3), std::max(v0, v1));
                const IdxType s = std::max(v2, v3);                
                IdxType zind0 = ((zmask & (1 << p)) >> p);
                IdxType zind1 = ((zmask & (1 << q)) >> q) << 1;
                IdxType zind2 = ((zmask & (1 << r)) >> r) << 2;
                IdxType zind3 = ((zmask & (1 << s)) >> s) << 3;
                const ValType* gm_real = exp_gate_perms_4q[zind0 + zind1 + zind2 + zind3];
                const ValType* gm_imag = exp_gate_perms_4q[zind0 + zind1 + zind2 + zind3] + 256;
                result = EXPECT_C4V1_GATE(gm_real, gm_imag, p, q, r, s, xmask | zmask);
                BARR_MPI;
                if (q2 >= lg2_m_cpu) {
                    SWAP_GATE(q2t, q2);
                    BARR_MPI;
                }
                if (q1 >= lg2_m_cpu) {
                    SWAP_GATE(q1t, q1);
                    BARR_MPI;
                }
                if (q0 >= lg2_m_cpu) {
                    SWAP_GATE(q0t, q0);
                    BARR_MPI;
                }
            } else if (num_x_indices == 0) {
                result = EXPECT_C0_GATE(zmask);
            }
            ValType expect = 0;
            // printf("%lld %f\n", i_proc, result);
            BARR_MPI;
            MPI_Reduce(&result, &expect, 1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            BARR_MPI;
            if (i_proc == 0) {
                output[output_index] = expect;
            }
           
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

            if (i_proc == (n_cpus - 1)) // last node
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
                if (local_cpu == i_proc)
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
                IdxType idx = (i_proc)*per_pe_work + i;
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
                    IdxType idx = (i_proc)*per_pe_work + i;
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
                if (qubit >= lg2_m_cpu) // remote qubit, need switch
                {
                    IdxType pair_cpu = (i_proc) ^ ((IdxType)1 << (qubit - (lg2_m_cpu)));
                    assert(pair_cpu != i_proc);
                    ValType *sv_real_remote = m_real;
                    ValType *sv_imag_remote = m_imag;

                    if (i_proc > pair_cpu)
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
            IdxType pair_cpu = (i_proc) ^ ((IdxType)1 << (q - (lg2_m_cpu)));
            assert(pair_cpu != i_proc);

            if (i_proc > pair_cpu)
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

                IdxType index = (i_proc >> (q - (lg2_m_cpu) + 1)) << (q - (lg2_m_cpu));
                index |= i_proc & ((1 << (q - (lg2_m_cpu))) - 1);
                for (IdxType i = (index)*per_pe_work; i < (index + 1) * per_pe_work; i++)
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
