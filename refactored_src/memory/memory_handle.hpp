#pragma once

#include "../nwqsim_core.hpp"

#include <cstddef>
#include <cstring>
#include <stdexcept>
#include <utility>
#include <vector>

namespace NWQSim
{
    struct DeviceGate
    {
        enum class Type
        {
            Matrix,
            Measure,
            MeasureAll,
            Reset
        };

        Type type = Type::Matrix;
        IdxType qubit0 = 0;
        IdxType qubit1 = -1;
        IdxType repetition = 1;
        std::vector<ValType> gm_real;
        std::vector<ValType> gm_imag;
        IdxType matrix_dim = 0;

        void resize_matrix(IdxType dim)
        {
            matrix_dim = dim;
            const std::size_t entries = static_cast<std::size_t>(dim) * static_cast<std::size_t>(dim);
            gm_real.assign(entries, 0.0);
            gm_imag.assign(entries, 0.0);
        }
    };

    namespace Memory
    {
        enum class Backend
        {
            CPU
        };

        struct MemoryConfig
        {
            Backend backend = Backend::CPU;
            SimType method = SimType::SV;
            IdxType n_qubits = 0;
        };

        class MemoryHandle
        {
        public:
            StateSpan view{};
            Backend backend = Backend::CPU;
            SimType method = SimType::SV;
            IdxType n_qubits = 0;

            MemoryHandle() = default;
            MemoryHandle(MemoryHandle &&) noexcept = default;
            MemoryHandle &operator=(MemoryHandle &&) noexcept = default;

        private:
            std::vector<ValType> state_real_;
            std::vector<ValType> state_imag_;
            std::vector<ValType> buffer_;

            void initialise_state()
            {
                IdxType dim_log = (method == SimType::SV) ? n_qubits : 2 * n_qubits;
                const IdxType dim = dim_log > 0 ? (IdxType(1) << dim_log) : 1;
                const IdxType half_dim = (method == SimType::SV && dim_log > 0) ? (dim >> 1) : 0;

                state_real_.assign(static_cast<std::size_t>(dim), 0.0);
                state_imag_.assign(static_cast<std::size_t>(dim), 0.0);
                if (method == SimType::SV)
                {
                    buffer_.assign(static_cast<std::size_t>(dim), 0.0);
                }
                else
                {
                    buffer_.assign(static_cast<std::size_t>((IdxType(1) << n_qubits) + 1), 0.0);
                }

                state_real_[0] = 1.0;

                view.real = state_real_.data();
                view.imag = state_imag_.data();
                view.work_real = buffer_.data();
                view.work_imag = nullptr;
                view.size = dim;
                view.half_size = half_dim;
            }

            friend MemoryHandle create_memory(const MemoryConfig &config);
            friend void reset_state(MemoryHandle &handle);
        };

        inline MemoryHandle create_memory(const MemoryConfig &config)
        {
            if (config.backend != Backend::CPU)
                throw std::runtime_error("Only CPU backend is supported in the refactored build");

            MemoryHandle handle;
            handle.backend = Backend::CPU;
            handle.method = config.method;
            handle.n_qubits = config.n_qubits;
            handle.initialise_state();
            return handle;
        }

        inline void reset_state(MemoryHandle &handle)
        {
            handle.initialise_state();
        }
    } // namespace Memory
} // namespace NWQSim
