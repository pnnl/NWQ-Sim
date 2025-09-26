#pragma once

#include "simulator/device_gate_builder.hpp"
#include "simulator/sim_types.hpp"

#include "gate_kernels/cpu_gates.hpp"
#include "memory/memory_handle.hpp"

#include <random>

namespace NWQSim
{
    class SVSim
    {
    public:
        explicit SVSim(IdxType n_qubits,
                       Memory::Backend backend = Memory::Backend::CPU,
                       uint64_t seed = 0);

        void set_seed(uint64_t seed);
        void reset_state();
        void set_fusion_enabled(bool enabled);

        SimulationResult simulate(const Circuit &circuit);

        const StateSpan &state() const { return memory_.view; }

    private:
        SimulationResult run_cpu(const SimulationTask &task);

        Memory::MemoryHandle memory_;
        DeviceGateBuilder builder_{};
        std::mt19937_64 rng_;
        std::uniform_real_distribution<ValType> dist_{0.0, 1.0};
        bool fusion_enabled_ = true;
    };
}
