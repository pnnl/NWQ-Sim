#include "simulator/sv_sim.hpp"

#include <stdexcept>
#include <vector>

namespace NWQSim
{
    SVSim::SVSim(IdxType n_qubits,
                 Memory::Backend backend,
                 uint64_t seed)
        : rng_(seed)
    {
        Memory::MemoryConfig config{};
        config.method = SimType::SV;
        config.n_qubits = n_qubits;
        memory_ = Memory::create_memory(config);
    }

    void SVSim::set_seed(uint64_t seed)
    {
        rng_.seed(seed);
    }

    void SVSim::reset_state()
    {
        Memory::reset_state(memory_);
    }

    void SVSim::set_fusion_enabled(bool enabled)
    {
        fusion_enabled_ = enabled;
    }

    SimulationResult SVSim::simulate(const Circuit &circuit)
    {
        SimulationTask task = builder_.build(circuit, fusion_enabled_);

        if (task.metadata.n_qubits <= 0)
            throw std::runtime_error("Circuit reports zero qubits");

        if (task.metadata.n_qubits != memory_.n_qubits)
        {
            Memory::MemoryConfig config{};
            config.backend = memory_.backend;
            config.method = SimType::SV;
            config.n_qubits = task.metadata.n_qubits;
            memory_ = Memory::create_memory(config);
        }

        reset_state();

        return run_cpu(task);
    }

    SimulationResult SVSim::run_cpu(const SimulationTask &task)
    {
        SimulationResult result;
        result.metadata = task.metadata;
        result.single_qubit_measurements.reserve(task.metadata.measurements);
        result.measure_all_results.reserve(task.metadata.measure_all_ops);

        for (const auto &gate : task.device_gates)
        {
            switch (gate.type)
            {
            case DeviceGate::Type::Matrix:
                if (gate.matrix_dim <= 2)
                {
                    GateKernels::CPU::apply_c1_gate(gate.gm_real.data(),
                                                    gate.gm_imag.data(),
                                                    gate.qubit0,
                                                    memory_.view);
                }
                else if (gate.matrix_dim == 4)
                {
                    GateKernels::CPU::apply_c2_gate(gate.gm_real.data(),
                                                    gate.gm_imag.data(),
                                                    gate.qubit0,
                                                    gate.qubit1,
                                                    memory_.view);
                }
                else
                {
                    throw std::runtime_error("Unsupported matrix dimension in DeviceGate");
                }
                break;

            case DeviceGate::Type::Measure:
            {
                ValType random_value = dist_(rng_);
                IdxType outcome = GateKernels::CPU::measure_qubit(memory_.view,
                                                                  memory_.n_qubits,
                                                                  gate.qubit0,
                                                                  random_value);
                result.single_qubit_measurements.push_back(outcome);
                break;
            }

            case DeviceGate::Type::MeasureAll:
            {
                MeasureAllResult record;
                record.repetitions = gate.repetition;
                record.outcomes.resize(static_cast<size_t>(gate.repetition));

                std::vector<ValType> randoms(static_cast<size_t>(gate.repetition));
                for (auto &value : randoms)
                    value = dist_(rng_);

                GateKernels::CPU::measure_all(memory_.view,
                                              memory_.n_qubits,
                                              gate.repetition,
                                              memory_.view.work_real,
                                              randoms.data(),
                                              record.outcomes.data());

                result.measure_all_results.push_back(std::move(record));
                break;
            }

            case DeviceGate::Type::Reset:
                GateKernels::CPU::reset_qubit(memory_.view,
                                              memory_.n_qubits,
                                              gate.qubit0);
                break;
            }
        }

        return result;
    }
}
