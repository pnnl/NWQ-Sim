#pragma once

#include "state.hpp"

#include "nwq_util.hpp"
#include "gate.hpp"
#include "circuit.hpp"
#include "config.hpp"
#include "private/exp_gate_declarations_host.hpp"

#include "circuit_pass/fusion.hpp"
#include "private/macros.hpp"
#include "private/sim_gate.hpp"

#include "prng.hpp"

#include <random>
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <set>
// #include "../xsimd/xsimd.hpp"

namespace NWQSim
{
    class STAB_CPU : public QuantumState
    {

    public:
        // Default identity constructor
        bool has_destabilizers;

        STAB_CPU(IdxType _n_qubits) : QuantumState(SimType::STAB)
        {
            n = _n_qubits;
            rows = 2 * n + 1;
            cols = n;
            stabCounts = n;
            x.resize(rows, std::vector<int>(cols, 0)); // first 2n+1 x n block. first n represents destabilizers
                                                       // second n represents stabilizers + 1 extra row
            z.resize(rows, std::vector<int>(cols, 0)); // second 2n+1 x n block to form the 2n+1 x 2n sized tableau
            r.resize(rows, 0);                         // column on the right with 2n+1 rows
            // The 2n+1 th row is scratch space

            // Intialize the identity tableau
            for (int i = 0; i < n; i++)
            {
                x[i][i] = 1;
                z[i + n][i] = 1;
            }

            std::random_device rd;
            rng.seed(rd());
            dist = std::uniform_int_distribution<int>(0, 1);
            random_float = std::uniform_real_distribution<double>(0.0, 1.0);
            measurement_count = 0; // initialize
            seed = 0;

            has_destabilizers = true;
        }

        ~STAB_CPU()
        {
            // Release for CPU side
            // SAFE_FREE_HOST(totalResults);
        }

        // resets the tableau to a full identity
        void reset_state() override
        {
            rows = 2 * n + 1;
            cols = n;
            stabCounts = n;
            for (auto &row : x)
            {
                std::fill(row.begin(), row.end(), 0);
            }
            for (auto &row : z)
            {
                std::fill(row.begin(), row.end(), 0);
            }
            std::fill(r.begin(), r.end(), 0);
            // The 2n+1 th row is scratch space

            // Intialize the identity tableau
            for (int i = 0; i < n; i++)
            {
                x[i][i] = 1;
                z[i + n][i] = 1;
            }

            seed = 0;
            measurement_count = 0;

            m_results.clear();
            m_results.reserve(est_measurements);

            has_destabilizers = true;
        }

        void allocate_measurement_buffers(IdxType measurements)
        {
            est_measurements = measurements;
            m_results.reserve(est_measurements);
        }

        void set_seed(IdxType s) override
        {
            rng.seed(s);
            seed = s;
            measurement_count = 0;
        }

        // Prints the tableau including phase
        void print_res_state() override
        {
            for (int i = 0; i < rows; i++)
            {
                if (has_destabilizers && ((i == (rows / 2)) || (i == rows - 1)))
                {
                    for (int j = -5; j < (n * 2); j++)
                    {
                        std::cout << "-";
                    }
                    std::cout << std::endl;
                }
                for (int j = 0; j < cols; j++)
                {
                    std::cout << x[i][j];
                }
                std::cout << " | ";
                for (int j = 0; j < cols; j++)
                {
                    std::cout << z[i][j];
                }
                std::cout << "|" << r[i] << std::endl;
            }
            std::cout << std::endl;
        }

        // Prints a 2n x m tableau matrix w/o phase bits
        void print_table(std::vector<std::vector<int>> &M)
        {
            for (int i = 0; i < M[0].size() + 1; i++)
            {
                std::cout << "- ";
            }
            std::cout << std::endl;
            for (int i = 0; i < M.size(); i++)
            {
                for (int j = 0; j < M[i].size() / 2; j++)
                {
                    std::cout << M[i][j] << " ";
                }
                std::cout << "| ";
                for (int j = M[i].size() / 2; j < M[i].size(); j++)
                {
                    std::cout << M[i][j] << " ";
                }
                std::cout << std::endl;
            }
        }

        // Sub-process in measurement gates
        void rowsum(int h, int i)
        {
            int sum = 0;
            for (int j = 0; j < n; j++)
            {
                // Sum every column in the row
                if (x[i][j])
                {
                    if (z[i][j])
                        sum += z[h][j] - x[h][j];
                    else
                        sum += z[h][j] * (2 * x[h][j] - 1);
                }
                else if (z[i][j])
                    sum += x[h][j] * (1 - 2 * z[h][j]);

                // XOR x's and z's
                x[h][j] = x[i][j] ^ x[h][j];
                z[h][j] = z[i][j] ^ z[h][j];
            }

            // if((sum % 4 != 0) && (abs(sum % 4) != 2))
            // {
            //     printf("m%d tot: %d\n", measurement_count, sum);
            // }
            // printf("m%d sum: %d\n", measurement_count, sum);
            // printf("r[i]%d: %d\n", measurement_count, r[i]);

            sum += 2 * r[i];

            sum += 2 * r[h];
            //  printf("m%d sum: %d\n", *d_measurement_idx_counter, total);

            // printf("m%d tot: %d\n", measurement_count, sum);

            r[h] = (sum % 4 == 0) ? 0 : 1;

            // printf("r[h]%d: %d\n", measurement_count, r[h]);

            // std::cout << "r[" << h << "]" <<  " after if = " <<  r[h] << std::endl;

        } // End rowsum

        void i_rowsum(int h, int i)
        {
            int sum = 0;
            for (int j = 0; j < n; j++)
            {
                // Sum every column in the row
                if (x[i][j])
                {
                    if (z[i][j])
                        sum += z[h][j] - x[h][j];
                    else
                        sum += z[h][j] * (2 * x[h][j] - 1);
                }
                else if (z[i][j])
                    sum += x[h][j] * (1 - 2 * z[h][j]);

                // XOR x's and z's
                x[h][j] = x[i][j] ^ x[h][j];
                z[h][j] = z[i][j] ^ z[h][j];
            }
            sum += 1 + 2 * r[h] + 2 * r[i];

            if (sum % 4 == 0)
                r[h] = 0;
            else
                r[h] = 1;
        } // End rowsum

        // Simulate the gates from a circuit in the tableau
        void sim(std::shared_ptr<NWQSim::Circuit> circuit) override
        {
            std::vector<Gate> gates = circuit->get_gates();
            IdxType n_gates = gates.size();

            // Check that the circuit object matches the tableau
            assert((circuit->num_qubits() == n) && "Circuit does not match the qubit number of the state!");

            // Start a timer
            cpu_timer sim_timer;
            sim_timer.start_timer();

            // Perform the tableau simulation
            simulation_kernel(gates);

            // End timer
            sim_timer.stop_timer();
            double sim_time = sim_timer.measure();

            if (Config::PRINT_SIM_TRACE)
            {
                printf("\n============== STAB-Sim ===============\n");
                printf("n_qubits:%lld, n_gates:%lld, ncpus:%d, comp:%.3lf ms, comm:%.3lf ms, sim:%.3lf ms\n",
                       n, n_gates, 1, sim_time, 0.,
                       sim_time);
                printf("=====================================\n");
            }
            //=========================================
        }

        IdxType *get_results() override
        {
            static std::vector<IdxType> converted;
            converted.assign(m_results.begin(), m_results.end());
            return converted.data();
        }

        IdxType measure(IdxType qubit) override
        {
            throw std::logic_error("measure Not supported (STAB_CPU)");
        }

        IdxType *measure_all(IdxType repetition) override
        {
            throw std::logic_error("measure_all Not supported (STAB_CPU)");
        }

        void set_initial(std::string fpath, std::string format) override
        {
            throw std::logic_error("set_initial Not supported (STAB_CPU)");
        }

        void dump_res_state(std::string outfile) override
        {
            throw std::logic_error("dump_res_state Not supported (STAB_CPU)");
        }

        ValType *get_real() const override
        {
            throw std::logic_error("get_real Not supported in STAB_CPU");
        }

        ValType *get_imag() const override
        {
            throw std::logic_error("get_imag Not supported in STAB_CPU");
        }

        ValType get_exp_z(const std::vector<size_t> &in_bits) override
        {
            throw std::logic_error("get_exp_z Not supported in STAB_CPU");
        }

        ValType get_exp_z() override
        {
            throw std::logic_error("get_exp_z Not supported in STAB_CPU");
        }

    protected:
        IdxType n;
        int stabCounts;
        IdxType rows;
        IdxType cols;
        std::vector<std::vector<int>> x;
        std::vector<std::vector<int>> z;
        std::vector<int> r;

        std::vector<int32_t> m_results;

        IdxType seed;              // used to align random measurements with GPU
        IdxType measurement_count; // used to align random measurements with GPU
        IdxType est_measurements = 0;
        std::mt19937 rng;
        std::uniform_int_distribution<int> dist;
        std::uniform_real_distribution<double> random_float;

        void apply_X(int q)
        {
            for (int i = 0; i < rows - 1; i++)
            {
                r[i] ^= z[i][q];
            }
        }
        void apply_Y(int q)
        {
            for (int i = 0; i < rows - 1; i++)
            {
                r[i] ^= (x[i][q] ^ z[i][q]);
            }
        }
        void apply_Z(int q)
        {
            for (int i = 0; i < rows - 1; i++)
            {
                r[i] ^= x[i][q];
            }
        }

        // For shot based measurement
        void tempRowsum(int h, int i, std::vector<std::vector<int>> &temp_x, std::vector<std::vector<int>> &temp_z, std::vector<int> &temp_r)
        {
            int sum = 0;
            for (int j = 0; j < n; j++)
            {
                // Sum every column in the row
                if (temp_x[i][j])
                {
                    if (z[i][j])
                        sum += temp_z[h][j] - temp_x[h][j];
                    else
                        sum += temp_z[h][j] * (2 * temp_x[h][j] - 1);
                }
                else if (temp_z[i][j])
                    sum += temp_x[h][j] * (1 - 2 * temp_z[h][j]);

                // XOR x's and z's
                temp_x[h][j] = temp_x[i][j] ^ temp_x[h][j];
                temp_z[h][j] = temp_z[i][j] ^ temp_z[h][j];
            }
            sum = sum + 2 * temp_r[h] + 2 * temp_r[i];

            if (sum % 4 == 0)
                temp_r[h] = 0;
            else
                temp_r[h] = 1;
        } // End rowsum

        // Provides a bit/shot measurement output without affecting the original tableau
        void apply_ma(IdxType shots)
        {
            m_results.clear();
            m_results.reserve(shots);

            int half_rows = rows >> 1;

            for (int shot = 0; shot < shots; shot++)
            {
                // Make a copy of the class being measured so many shots can be performed
                std::vector<std::vector<int>> temp_x = x;
                std::vector<std::vector<int>> temp_z = z;
                std::vector<int> temp_r = r;
                for (int a = 0; a < n; a++)
                {
                    int p = -1;
                    for (int p_index = half_rows; p_index < rows - 1; p_index++)
                    {
                        // std::cout << "x at [" << p_index << "][" << a << "] = " << x[p_index][a] << std::endl;
                        if (temp_x[p_index][a])
                        {
                            p = p_index;
                            break;
                        }
                    }
                    // A p such that x[p][a] = 1 exists
                    if (p > -1)
                    {
                        for (int i = 0; i < rows - 1; i++)
                        {
                            if ((x[i][a]) && (i != p))
                            {
                                tempRowsum(i, p, temp_x, temp_z, temp_r);
                            }
                        }
                        temp_x[p - half_rows] = temp_x[p];
                        temp_z[p - half_rows] = temp_z[p];
                        // Change all the columns in row p to be 0
                        for (int i = 0; i < n; i++)
                        {
                            temp_x[p][i] = 0;
                            temp_z[p][i] = 0;
                        }

                        int randomBit = dist(rng);

                        if (randomBit)
                        {
                            // std::cout << "Random result of 1" << std::endl;
                            temp_r[p] = 1;
                        }
                        else
                        {
                            // std::cout << "Random result of 0" << std::endl;
                            temp_r[p] = 0;
                        }
                        temp_z[p][a] = 1;

                        m_results.push_back(temp_r[p] << a);
                    }
                    else
                    {
                        // Set the scratch space row to be 0
                        // i is the column indexer in this case
                        for (int i = 0; i < n; i++)
                        {
                            temp_x[rows - 1][i] = 0;
                            temp_z[rows - 1][i] = 0;
                        }
                        temp_r[rows - 1] = 0;

                        // Run rowsum subroutine
                        for (int i = 0; i < half_rows; i++)
                        {
                            if (temp_x[i][a] == 1)
                            {
                                // std::cout << "Perform rowsum at " << i << " + n" << std::endl;
                                tempRowsum(rows - 1, i + half_rows, temp_x, temp_z, temp_r);
                            }
                        }

                        m_results.push_back(temp_r[rows - 1] << a);
                    } // End if else
                } // End single shot for all qubits

            } // End shots
        } // End measure_all

        void simulation_kernel(std::vector<Gate> &gates)
        {
            int g = gates.size();

            // For swapping rows
            int tempVal;
            int half_rows = rows >> 1;
            // Loop over every gate in the circuit and apply them
            for (int k = 0; k < g; k++)
            {
                auto gate = gates[k];
                int a = gate.qubit;
                // std::cout << "Damping inactivated rrrrrr " << gate.gamma << std::endl;

                switch (gate.op_name)
                {
                case OP::CX:
                {
                    int a = gate.ctrl;
                    int b = gate.qubit;
                    for (int i = 0; i < rows - 1; i++)
                    {
                        // Phase
                        r[i] = r[i] ^ (x[i][a] & z[i][b] & (x[i][b] ^ z[i][a] ^ 1));

                        // Entry
                        x[i][b] ^= x[i][a];
                        z[i][a] ^= z[i][b];
                    }
                    break;
                }

                case OP::H:
                    for (int i = 0; i < rows - 1; i++)
                    {
                        // Phase
                        r[i] ^= (x[i][a] & z[i][a]);
                        // Entry -- swap x and z bits
                        tempVal = x[i][a];
                        x[i][a] = z[i][a];
                        z[i][a] = tempVal;
                    }
                    break;

                case OP::S:
                    for (int i = 0; i < rows - 1; i++)
                    {
                        // Phase
                        r[i] ^= (x[i][a] & z[i][a]);

                        // Entry
                        z[i][a] ^= x[i][a];
                    }
                    break;

                case OP::SDG:
                    for (int i = 0; i < rows - 1; i++)
                    {
                        // Phase -- Equal to Z S or x & !z
                        r[i] ^= x[i][a] ^ (x[i][a] & z[i][a]);

                        // Entry
                        z[i][a] ^= x[i][a];
                    }
                    break;

                case OP::RX:
                    // H SDG H
                    if (gate.theta == PI / 2)
                    {
                        for (int i = 0; i < rows - 1; i++)
                        {
                            // Phase
                            r[i] ^= (x[i][a] & z[i][a]);
                            // Entry -- swap x and z bits
                            tempVal = x[i][a];
                            x[i][a] = z[i][a];
                            z[i][a] = tempVal;

                            // Phase -- Equal to Z S or x & !z
                            r[i] ^= x[i][a] ^ (x[i][a] & z[i][a]);

                            // Entry
                            z[i][a] ^= x[i][a];

                            r[i] ^= (x[i][a] & z[i][a]);
                            // Entry -- swap x and z bits
                            tempVal = x[i][a];
                            x[i][a] = z[i][a];
                            z[i][a] = tempVal;
                        }
                    }
                    // H S
                    else if (gate.theta == -PI / 2)
                    {
                        for (int i = 0; i < rows - 1; i++)
                        {
                            // Phase
                            r[i] ^= (x[i][a] & z[i][a]);
                            // Entry -- swap x and z bits
                            tempVal = x[i][a];
                            x[i][a] = z[i][a];
                            z[i][a] = tempVal;

                            // S
                            // Phase
                            r[i] ^= (x[i][a] & z[i][a]);

                            // Entry
                            z[i][a] ^= x[i][a];

                            // Phase
                            r[i] ^= (x[i][a] & z[i][a]);
                            // Entry -- swap x and z bits
                            tempVal = x[i][a];
                            x[i][a] = z[i][a];
                            z[i][a] = tempVal;
                        }
                    }
                    else
                    {
                        std::cout << "Non-Clifford angle in RX gate! "
                                  << OP_NAMES[gate.op_name] << "(" << gate.theta << ")" << std::endl;
                        std::logic_error("Invalid gate type");
                    }
                    break;

                case OP::RY:
                    // H X -- X : r[i] ^= z[i][a]
                    if (gate.theta == PI / 2)
                    {
                        for (int i = 0; i < rows - 1; i++)
                        {
                            // Phase -- Equal to Z S or x & !z
                            r[i] ^= x[i][a] ^ (x[i][a] & z[i][a]);

                            // Entry
                            z[i][a] ^= x[i][a];

                            // Phase
                            r[i] ^= (x[i][a] & z[i][a]);
                            // Entry -- swap x and z bits
                            tempVal = x[i][a];
                            x[i][a] = z[i][a];
                            z[i][a] = tempVal;

                            // Phase
                            r[i] ^= (x[i][a] & z[i][a]);

                            // Entry
                            z[i][a] ^= x[i][a];
                        }
                    }
                    // X H
                    else if (gate.theta == -PI / 2)
                    {
                        for (int i = 0; i < rows - 1; i++)
                        {
                            // Phase -- Equal to Z S or x & !z
                            r[i] ^= (x[i][a] & z[i][a]);

                            // Entry
                            z[i][a] ^= x[i][a];

                            // Phase
                            r[i] ^= (x[i][a] & z[i][a]);
                            // Entry -- swap x and z bits
                            tempVal = x[i][a];
                            x[i][a] = z[i][a];
                            z[i][a] = tempVal;

                            // Phase
                            r[i] ^= x[i][a] ^ (x[i][a] & z[i][a]);

                            // Entry
                            z[i][a] ^= x[i][a];
                        }
                    }
                    else
                    {
                        std::cout << "Non-Clifford angle in RY gate! "
                                  << OP_NAMES[gate.op_name] << "(" << gate.phi << ")" << std::endl;
                        std::logic_error("Invalid gate type");
                    }
                    break;
                case OP::CY:
                {
                    int a = gate.ctrl;
                    int b = gate.qubit;
                    for (int i = 0; i < rows - 1; i++)
                    {
                        // Phase -- Equal to Z S or x & !z
                        r[i] ^= x[i][b] ^ (x[i][b] & z[i][b]);

                        // Entry
                        z[i][b] ^= x[i][b];
                    }
                    for (int i = 0; i < rows - 1; i++)
                    {
                        // Phase
                        r[i] ^= ((x[i][a] & z[i][b]) & (x[i][b] ^ z[i][a] ^ 1));

                        // Entry
                        x[i][b] ^= x[i][a];
                        z[i][a] ^= z[i][b];
                    }
                    for (int i = 0; i < rows - 1; i++)
                    {
                        // Phase
                        r[i] ^= (x[i][b] & z[i][b]);

                        // Entry
                        z[i][b] ^= x[i][b];
                    }
                    break;
                }
                case OP::CZ:
                {
                    int a = gate.ctrl;
                    int b = gate.qubit;

                    for (int i = 0; i < rows - 1; i++)
                    {
                        // Phase
                        r[i] ^= (x[i][b] & z[i][b]);
                        // Entry -- swap x and z bits
                        tempVal = x[i][b];
                        x[i][b] = z[i][b];
                        z[i][b] = tempVal;
                    }
                    for (int i = 0; i < rows - 1; i++)
                    {
                        // Phase
                        r[i] ^= ((x[i][a] & z[i][b]) & (x[i][b] ^ z[i][a] ^ 1));

                        // Entry
                        x[i][b] ^= x[i][a];
                        z[i][a] ^= z[i][b];
                    }
                    for (int i = 0; i < rows - 1; i++)
                    {
                        // Phase
                        r[i] ^= (x[i][b] & z[i][b]);
                        // Entry -- swap x and z bits
                        tempVal = x[i][b];
                        x[i][b] = z[i][b];
                        z[i][b] = tempVal;
                    }
                    break;
                }

                case OP::M:
                {
                    int p = -1;
                    for (int p_index = half_rows; p_index < rows - 1; p_index++)
                    {
                        // printf("x at p_index = %d, qubit=%d ", x[p_index][a], a);
                        // std::cout << "x at [" << p_index << "][" << a << "] = " << x[p_index][a] << std::endl;
                        if (x[p_index][a])
                        {
                            p = p_index;
                            break;
                        }
                    }
                    // std::cout << "p = " << p << std::endl;
                    // A p such that x[p][a] = 1 exists
                    // Random
                    if (p > -1)
                    {
                        // std::cout << "CPU p = " << p << std::endl;
                        for (int i = 0; i < rows - 1; i++)
                        {
                            // std::cout << "x = " << x[i][a] << std::endl;
                            if ((x[i][a]) && (i != p))
                            {
                                rowsum(i, p);
                                // printf("rand r%d row%d\n", r[i], i);
                            }
                        }

                        x[p - half_rows] = x[p];
                        z[p - half_rows] = z[p];
                        // Change all the columns in row p to be 0
                        for (int i = 0; i < n; i++)
                        {
                            x[p][i] = 0;
                            z[p][i] = 0;
                        }

                        int randomBit = prng_bit(seed, measurement_count);
                        // printf("rand m%d, res%d\n",  measurement_count, randomBit);

                        // std::cout << "Seed for measurement " << measurement_count << ": " << seed << std::endl;
                        r[p] = randomBit;
                        z[p][a] = 1;

                        // std::cout << "Random measurement " << measurement_count << ": " << randomBit << std::endl;
                        m_results.push_back(randomBit);
                        measurement_count++;
                    }
                    // Deterministic
                    else
                    {
                        // Set the scratch space row to be 0
                        // i is the column indexer in this case
                        for (int i = 0; i < n; i++)
                        {
                            x[rows - 1][i] = 0;
                            z[rows - 1][i] = 0;
                        }
                        r[rows - 1] = 0;

                        // Run rowsum subroutine
                        for (int i = 0; i < half_rows; i++)
                        {
                            if (x[i][a] == 1)
                            {
                                // printf("r[row]%d ", r[i+half_rows]);
                                rowsum(rows - 1, i + half_rows);
                                // printf("determ r%d row%lld\n", r[rows-1], i+half_rows);
                            }
                        }
                        // std::cout << "Determ measurement " << measurement_count << ": " << r[rows-1] << std::endl;
                        m_results.push_back(r[rows - 1]);
                        // printf("determ m%d, res%d\n",  measurement_count, r[rows-1]);
                        measurement_count++;
                    }
                    break;
                }

                case OP::RESET:
                {
                    int p = -1;
                    int temp_result = 0;
                    int half_rows = rows / 2;

                    for (int p_index = half_rows; p_index < rows - 1; p_index++)
                    {
                        // std::cout << "x at [" << p_index << "][" << a << "] = " << x[p_index][a] << std::endl;
                        if (x[p_index][a])
                        {
                            p = p_index;
                            break;
                        }
                    }
                    // std::cout << "p = " << p << std::endl;
                    // A p such that x[p][a] = 1 exists
                    // Random
                    if (p > -1)
                    {
                        for (int i = 0; i < rows - 1; i++)
                        {
                            // std::cout << "x = " << x[i][a] << std::endl;
                            if ((x[i][a]) && (i != p))
                            {
                                rowsum(i, p);
                            }
                        }

                        x[p - half_rows] = x[p];
                        z[p - half_rows] = z[p];
                        // Change all the columns in row p to be 0
                        for (int i = 0; i < n; i++)
                        {
                            x[p][i] = 0;
                            z[p][i] = 0;
                        }

                        int randomBit = prng_bit(seed, measurement_count);
                        // std::cout << "Seed for measurement " << measurement_count << ": " << seed << std::endl;
                        r[p] = randomBit;
                        z[p][a] = 1;

                        temp_result = r[p];
                        // std::cout << "Random measurement at qubit " << a << " value: " << (r[p] << a) << std::endl;
                    }
                    // Deterministic
                    else
                    {
                        // Set the scratch space row to be 0
                        // i is the column indexer in this case
                        for (int i = 0; i < n; i++)
                        {
                            x[rows - 1][i] = 0;
                            z[rows - 1][i] = 0;
                        }
                        r[rows - 1] = 0;

                        // Run rowsum subroutine
                        for (int i = 0; i < half_rows; i++)
                        {
                            if (x[i][a] == 1)
                            {
                                rowsum(rows - 1, i + half_rows);
                            }
                        }
                        // std::cout << "Deterministc measurement at qubit " << a << " value: " << (r[rows-1] << a) << std::endl;
                        temp_result = r[rows - 1];
                    }
                    // measurement_results.push_back(temp_result);
                    if (temp_result == 1) // Apply X to flip back to 0
                    {
                        apply_X(a);
                    }
                    break;
                }

                case OP::X:
                    // equiv to H S S H or H Z H
                    apply_X(a);
                    break;

                case OP::Y:
                    // equiv to Z X
                    apply_Y(a);
                    break;

                case OP::Z:
                    // equiv to S S gates
                    apply_Z(a);
                    break;

                case OP::MA:
                    apply_ma(gate.repetition);
                    break;

                default:
                    std::cerr << "Non-Clifford or unrecognized gate: "
                              << OP_NAMES[gate.op_name] << std::endl;
                    std::logic_error("Invalid gate type");

                } // End switch
            } // End gates for loop
        } // End simulate
    }; // End tableau class
} // namespace NWQSim

// #endif