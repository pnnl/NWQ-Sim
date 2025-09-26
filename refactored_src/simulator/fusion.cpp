#include "simulator/fusion.hpp"

#include <algorithm>
#include <optional>
#include <vector>

namespace NWQSim
{
    namespace
    {
        inline void multiply_in_place(DeviceGate &accum, const DeviceGate &incoming, int dim)
        {
            const size_t elements = static_cast<size_t>(dim) * static_cast<size_t>(dim);
            std::vector<ValType> result_real(elements, 0.0);
            std::vector<ValType> result_imag(elements, 0.0);

            for (int row = 0; row < dim; ++row)
            {
                for (int col = 0; col < dim; ++col)
                {
                    ValType real = 0.0;
                    ValType imag = 0.0;
                    for (int k = 0; k < dim; ++k)
                    {
                        const ValType ar = incoming.gm_real[row * dim + k];
                        const ValType ai = incoming.gm_imag[row * dim + k];
                        const ValType br = accum.gm_real[k * dim + col];
                        const ValType bi = accum.gm_imag[k * dim + col];

                        real += ar * br - ai * bi;
                        imag += ar * bi + ai * br;
                    }
                    const size_t idx = static_cast<size_t>(row) * static_cast<size_t>(dim) + static_cast<size_t>(col);
                    result_real[idx] = real;
                    result_imag[idx] = imag;
                }
            }

            accum.gm_real.swap(result_real);
            accum.gm_imag.swap(result_imag);
        }

        inline void clear_pairs_for_qubit(std::vector<std::vector<std::optional<size_t>>> &last_two,
                                          IdxType qubit,
                                          IdxType keep = -1)
        {
            if (qubit < 0)
                return;
            const size_t q = static_cast<size_t>(qubit);
            if (q >= last_two.size())
                return;
            for (size_t other = 0; other < last_two.size(); ++other)
            {
                if (static_cast<IdxType>(other) == keep)
                    continue;
                last_two[q][other].reset();
                last_two[other][q].reset();
            }
        }

        inline void clear_all_pairs(std::vector<std::vector<std::optional<size_t>>> &last_two)
        {
            for (auto &row : last_two)
                for (auto &entry : row)
                    entry.reset();
        }

        inline void clear_all_singles(std::vector<std::optional<size_t>> &last_single)
        {
            for (auto &entry : last_single)
                entry.reset();
        }

        inline void clear_single(std::vector<std::optional<size_t>> &last_single, IdxType qubit)
        {
            if (qubit < 0)
                return;
            const size_t q = static_cast<size_t>(qubit);
            if (q < last_single.size())
                last_single[q].reset();
        }
    } // namespace

    std::vector<DeviceGate> fuse_device_gates(const std::vector<DeviceGate> &input, IdxType n_qubits)
    {
        if (input.empty())
            return {};

        std::vector<DeviceGate> output;
        output.reserve(input.size());

        const size_t qubits = static_cast<size_t>(std::max<IdxType>(0, n_qubits));
        std::vector<std::optional<size_t>> last_single(qubits);
        std::vector<std::vector<std::optional<size_t>>> last_two(qubits, std::vector<std::optional<size_t>>(qubits));

        auto flush_qubit = [&](IdxType qubit)
        {
            clear_single(last_single, qubit);
            clear_pairs_for_qubit(last_two, qubit);
        };

        auto record_two = [&](IdxType q0, IdxType q1, size_t index)
        {
            if (q0 < 0 || q1 < 0)
                return;
            const size_t s0 = static_cast<size_t>(q0);
            const size_t s1 = static_cast<size_t>(q1);
            if (s0 < qubits && s1 < qubits)
            {
                last_two[s0][s1] = index;
                last_two[s1][s0] = index;
            }
        };

        auto fetch_two = [&](IdxType q0, IdxType q1) -> std::optional<size_t>
        {
            if (q0 < 0 || q1 < 0)
                return std::nullopt;
            const size_t s0 = static_cast<size_t>(q0);
            const size_t s1 = static_cast<size_t>(q1);
            if (s0 < qubits && s1 < qubits)
                return last_two[s0][s1];
            return std::nullopt;
        };

        for (const auto &gate : input)
        {
            if (gate.type != DeviceGate::Type::Matrix)
            {
                if (gate.type == DeviceGate::Type::Measure || gate.type == DeviceGate::Type::Reset)
                {
                    flush_qubit(gate.qubit0);
                }
                else if (gate.type == DeviceGate::Type::MeasureAll)
                {
                    clear_all_singles(last_single);
                    clear_all_pairs(last_two);
                }

                output.push_back(gate);
                continue;
            }

            if (gate.matrix_dim <= 2)
            {
                const IdxType qubit = gate.qubit0;
                clear_pairs_for_qubit(last_two, qubit);

                if (qubit < 0 || static_cast<size_t>(qubit) >= qubits)
                {
                    output.push_back(gate);
                    continue;
                }

                auto &slot = last_single[static_cast<size_t>(qubit)];
                if (slot)
                {
                    multiply_in_place(output[*slot], gate, 2);
                }
                else
                {
                    output.push_back(gate);
                    slot = output.size() - 1;
                }
            }
            else if (gate.matrix_dim == 4)
            {
                const IdxType q0 = gate.qubit0;
                const IdxType q1 = gate.qubit1;

                clear_single(last_single, q0);
                clear_single(last_single, q1);
                clear_pairs_for_qubit(last_two, q0, q1);
                clear_pairs_for_qubit(last_two, q1, q0);

                auto slot = fetch_two(q0, q1);
                if (slot)
                {
                    multiply_in_place(output[*slot], gate, 4);
                }
                else
                {
                    output.push_back(gate);
                    const size_t index = output.size() - 1;
                    record_two(q0, q1, index);
                }
            }
            else
            {
                // Unsupported matrix dimension, flush state and append as-is.
                clear_all_singles(last_single);
                clear_all_pairs(last_two);
                output.push_back(gate);
            }
        }

        return output;
    }
}
