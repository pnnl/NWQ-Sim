#include "stabsim/stab_cpu.hpp"
#include "xsimd/xsimd.hpp"
#include <vector>

namespace NWQSim
{
    // SIMD-optimized rowsum function
    void STAB_CPU::rowsum_simd(int h, int i)
    {
        using batch_type = xsimd::batch<int>;
        const std::size_t batch_size = batch_type::size;
        
        batch_type sum_batch(0);

        std::size_t j = 0;
        for (; j + batch_size <= n; j += batch_size)
        {
            batch_type x_i = batch_type::load_unaligned(&x[i][j]);
            batch_type z_i = batch_type::load_unaligned(&z[i][j]);
            batch_type x_h = batch_type::load_unaligned(&x[h][j]);
            batch_type z_h = batch_type::load_unaligned(&z[h][j]);

            // Vectorized sum calculation based on conditions
            batch_type term1 = (z_h - x_h) & (x_i & z_i);
            batch_type term2 = (z_h * (2 * x_h - 1)) & (x_i & ~z_i);
            batch_type term3 = (x_h * (1 - 2 * z_h)) & (~x_i & z_i);
            sum_batch += term1 + term2 + term3;

            // Vectorized XOR operations
            x_h ^= x_i;
            z_h ^= z_i;

            x_h.store_unaligned(&x[h][j]);
            z_h.store_unaligned(&z[h][j]);
        }

        // Scalar part for remaining elements
        int sum = xsimd::hadd(sum_batch);
        for (; j < n; ++j)
        {
            if (x[i][j] == 1 && z[i][j] == 1)
                sum += z[h][j] - x[h][j];
            else if (x[i][j] == 1 && z[i][j] == 0)
                sum += z[h][j] * (2 * x[h][j] - 1);
            else if (x[i][j] == 0 && z[i][j] == 1)
                sum += x[h][j] * (1 - 2 * z[h][j]);

            x[h][j] = x[i][j] ^ x[h][j];
            z[h][j] = z[i][j] ^ z[h][j];
        }

        sum += 2 * r[i];
        sum += 2 * r[h];
        r[h] = (sum % 4 == 0) ? 0 : 1;
    }

    // SIMD-optimized simulation kernel
    void STAB_CPU::simulation_kernel_simd(std::vector<Gate>& gates)
    {
        int g = gates.size();
        int half_rows = rows >> 1;

        for (int k = 0; k < g; k++)
        {
            auto gate = gates[k];
            int a = gate.qubit;

            switch (gate.op_name)
            {
                case OP::S:
                {
                    using batch_type = xsimd::batch<int>;
                    const std::size_t batch_size = batch_type::size;
                    for (std::size_t i = 0; i < rows - 1; i += batch_size)
                    {
                        auto r_batch = batch_type::load_unaligned(&r[i]);
                        auto x_batch = batch_type::load_unaligned(&x[i][a]);
                        auto z_batch = batch_type::load_unaligned(&z[i][a]);

                        r_batch ^= (x_batch & z_batch);
                        z_batch ^= x_batch;

                        r_batch.store_unaligned(&r[i]);
                        z_batch.store_unaligned(&z[i][a]);
                    }
                    break;
                }

                case OP::H:
                {
                    using batch_type = xsimd::batch<int>;
                    const std::size_t batch_size = batch_type::size;
                    for (std::size_t i = 0; i < rows - 1; i += batch_size)
                    {
                        auto r_batch = batch_type::load_unaligned(&r[i]);
                        auto x_batch = batch_type::load_unaligned(&x[i][a]);
                        auto z_batch = batch_type::load_unaligned(&z[i][a]);

                        r_batch ^= (x_batch & z_batch);
                        
                        r_batch.store_unaligned(&r[i]);
                        z_batch.store_unaligned(&x[i][a]); // swap
                        x_batch.store_unaligned(&z[i][a]); // swap
                    }
                    break;
                }

                case OP::CX:
                {
                    int b = gate.ctrl;
                    using batch_type = xsimd::batch<int>;
                    const std::size_t batch_size = batch_type::size;
                    for (std::size_t i = 0; i < rows - 1; i += batch_size)
                    {
                        auto r_batch = batch_type::load_unaligned(&r[i]);
                        auto x_a_batch = batch_type::load_unaligned(&x[i][a]);
                        auto z_a_batch = batch_type::load_unaligned(&z[i][a]);
                        auto x_b_batch = batch_type::load_unaligned(&x[i][b]);
                        auto z_b_batch = batch_type::load_unaligned(&z[i][b]);

                        r_batch ^= (x_a_batch & z_b_batch & (x_b_batch ^ z_a_batch ^ 1));
                        x_b_batch ^= x_a_batch;
                        z_a_batch ^= z_b_batch;

                        r_batch.store_unaligned(&r[i]);
                        x_b_batch.store_unaligned(&x[i][b]);
                        z_a_batch.store_unaligned(&z[i][a]);
                    }
                    break;
                }

                case OP::M:
                {
                    int p = -1;
                    for (int p_index = half_rows; p_index < rows - 1; p_index++)
                    {
                        if (x[p_index][a])
                        {
                            p = p_index;
                            break;
                        }
                    }

                    if (p > -1) // Random case
                    {
                        for (int i = 0; i < rows - 1; i++)
                        {
                            if ((x[i][a]) && (i != p))
                            {
                                rowsum_simd(i, p);
                            }
                        }
                        x[p - half_rows] = x[p];
                        z[p - half_rows] = z[p];
                        std::fill(x[p].begin(), x[p].end(), 0);
                        std::fill(z[p].begin(), z[p].end(), 0);

                        int randomBit = prng_bit(seed, measurement_count);
                        r[p] = randomBit;
                        z[p][a] = 1;
                        m_results.push_back(randomBit);
                        measurement_count++;
                    }
                    else // Deterministic case
                    {
                        std::fill(x[rows - 1].begin(), x[rows - 1].end(), 0);
                        std::fill(z[rows - 1].begin(), z[rows - 1].end(), 0);
                        r[rows - 1] = 0;

                        for (int i = 0; i < half_rows; i++)
                        {
                            if (x[i][a] == 1)
                            {
                                rowsum_simd(rows - 1, i + half_rows);
                            }
                        }
                        m_results.push_back(r[rows - 1]);
                        measurement_count++;
                    }
                    break;
                }
                default:
                    // Other gates are ignored as per the request
                    break;
            }
        }
    }

    // SIMD-optimized sim function
    void STAB_CPU::sim_simd(std::shared_ptr<NWQSim::Circuit> circuit, double &sim_time)
    {
        std::vector<Gate> gates = circuit->get_gates();
        assert((circuit->num_qubits() == n) && "Circuit does not match the qubit number of the state!");

        cpu_timer sim_timer;
        sim_timer.start_timer();

        simulation_kernel_simd(gates);

        sim_timer.stop_timer();
        sim_time = sim_timer.measure();
    }
}
