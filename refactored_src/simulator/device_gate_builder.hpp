#pragma once

#include "circuit.hpp"
#include "gate.hpp"
#include "memory/memory_handle.hpp"
#include "simulator/fusion.hpp"
#include "simulator/sim_types.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace NWQSim
{
    struct SimulationTask
    {
        std::vector<DeviceGate> device_gates;
        SimulationMetadata metadata;
    };

    class DeviceGateBuilder
    {
    public:
        SimulationTask build(const Circuit &circuit, bool enable_fusion) const
        {
            SimulationTask task;
            task.metadata.n_qubits = circuit.num_qubits();

            std::vector<DeviceGate> gates;
            gates.reserve(circuit.gates().size());

            for (const auto &gate : circuit.gates())
            {
                if (auto classical = convert_classical_gate(gate))
                {
                    gates.emplace_back(std::move(*classical));
                    continue;
                }

                auto unitary = convert_unitary_gate(gate, circuit);
                if (!unitary)
                    throw std::runtime_error(std::string("Unsupported gate in circuit: ") + gate_kind_name(gate.kind));
                gates.emplace_back(std::move(*unitary));
            }

            if (enable_fusion)
            {
                gates = fuse_device_gates(gates, task.metadata.n_qubits);
            }

            task.metadata = summarize(gates, task.metadata.n_qubits);
            task.device_gates = std::move(gates);
            return task;
        }

    private:
        static SimulationMetadata summarize(const std::vector<DeviceGate> &gates, IdxType n_qubits)
        {
            SimulationMetadata metadata{};
            metadata.n_qubits = n_qubits;
            for (const auto &gate : gates)
            {
                metadata.total_device_gates++;
                switch (gate.type)
                {
                case DeviceGate::Type::Matrix:
                    if (gate.matrix_dim <= 2)
                        metadata.single_qubit_gates++;
                    else
                        metadata.two_qubit_gates++;
                    break;
                case DeviceGate::Type::Measure:
                    metadata.measurements++;
                    break;
                case DeviceGate::Type::MeasureAll:
                    metadata.measure_all_ops++;
                    metadata.measure_all_shots += gate.repetition;
                    break;
                case DeviceGate::Type::Reset:
                    metadata.resets++;
                    break;
                }
            }
            return metadata;
        }

        static DeviceGate make_single_qubit_gate(const Gate &gate,
                                                 const std::array<ValType, 4> &real,
                                                 const std::array<ValType, 4> &imag)
        {
            DeviceGate device_gate;
            device_gate.type = DeviceGate::Type::Matrix;
            device_gate.qubit0 = gate.target;
            device_gate.qubit1 = -1;
            device_gate.resize_matrix(2);
            std::copy(real.begin(), real.end(), device_gate.gm_real.begin());
            std::copy(imag.begin(), imag.end(), device_gate.gm_imag.begin());
            return device_gate;
        }

        static DeviceGate make_two_qubit_gate(const Gate &gate,
                                              const std::array<ValType, 16> &real,
                                              const std::array<ValType, 16> &imag)
        {
            DeviceGate device_gate;
            device_gate.type = DeviceGate::Type::Matrix;
            device_gate.qubit0 = gate.control;
            device_gate.qubit1 = gate.target;
            device_gate.resize_matrix(4);
            std::copy(real.begin(), real.end(), device_gate.gm_real.begin());
            std::copy(imag.begin(), imag.end(), device_gate.gm_imag.begin());
            return device_gate;
        }

        static std::optional<DeviceGate> convert_classical_gate(const Gate &gate)
        {
            switch (gate.kind)
            {
            case GateKind::Measure:
            {
                DeviceGate device_gate;
                device_gate.type = DeviceGate::Type::Measure;
                device_gate.qubit0 = gate.target;
                return device_gate;
            }
            case GateKind::MeasureAll:
            {
                DeviceGate device_gate;
                device_gate.type = DeviceGate::Type::MeasureAll;
                device_gate.repetition = gate.repetition;
                return device_gate;
            }
            case GateKind::Reset:
            {
                DeviceGate device_gate;
                device_gate.type = DeviceGate::Type::Reset;
                device_gate.qubit0 = gate.target;
                return device_gate;
            }
            default:
                return std::nullopt;
            }
        }

        static std::optional<DeviceGate> convert_unitary_gate(const Gate &gate, const Circuit &circuit)
        {
            ensure_supported_gate(gate.kind);

            switch (gate.kind)
            {
            case GateKind::H:
                return make_single_qubit_gate(gate,
                                              {S2I, S2I, S2I, -S2I},
                                              {0.0, 0.0, 0.0, 0.0});
            case GateKind::RX:
            {
                ValType angle = gate.uses_param ? circuit.parameter_value(gate.param_index) : gate.angle;
                ValType half = HALF * angle;
                return make_single_qubit_gate(gate,
                                              {std::cos(half), 0.0, 0.0, std::cos(half)},
                                              {0.0, -std::sin(half), -std::sin(half), 0.0});
            }
            case GateKind::RZ:
            {
                ValType angle = gate.uses_param ? circuit.parameter_value(gate.param_index) : gate.angle;
                ValType half = HALF * angle;
                return make_single_qubit_gate(gate,
                                              {std::cos(half), 0.0, 0.0, std::cos(half)},
                                              {-std::sin(half), 0.0, 0.0, std::sin(half)});
            }
            case GateKind::CX:
                return make_two_qubit_gate(gate,
                                           {1, 0, 0, 0,
                                            0, 1, 0, 0,
                                            0, 0, 0, 1,
                                            0, 0, 1, 0},
                                           {0});
            default:
                return std::nullopt;
            }
        }

        static void ensure_supported_gate(GateKind kind)
        {
            switch (kind)
            {
            case GateKind::RX:
            case GateKind::RZ:
            case GateKind::H:
            case GateKind::CX:
                return;
            default:
                throw std::runtime_error(std::string("Unsupported unitary gate: ") + gate_kind_name(kind));
            }
        }
    };
}

