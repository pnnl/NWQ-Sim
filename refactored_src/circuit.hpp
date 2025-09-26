#pragma once

#include "gate.hpp"

#include <stdexcept>
#include <utility>
#include <vector>

namespace NWQSim
{
    class Circuit
    {
    public:
        explicit Circuit(IdxType n_qubits) : n_qubits_(n_qubits) {}

        IdxType num_qubits() const { return n_qubits_; }
        const std::vector<Gate> &gates() const { return gates_; }
        std::vector<Gate> &gates() { return gates_; }

        size_t declare_parameter(ValType initial = 0.0)
        {
            parameters_.push_back(initial);
            return parameters_.size() - 1;
        }

        void set_parameter(size_t index, ValType value)
        {
            ValType &param = parameter_at(index);
            param = value;
        }

        void set_parameters(const std::vector<ValType> &values)
        {
            if (values.size() != parameters_.size())
                throw std::invalid_argument("Parameter vector length mismatch");
            for (size_t i = 0; i < values.size(); ++i)
                set_parameter(i, values[i]);
        }

        ValType parameter_value(size_t index) const
        {
            const ValType &param = parameter_at(index);
            return param;
        }

        void clear()
        {
            gates_.clear();
        }

        Gate &h(IdxType target)
        {
            Gate gate{};
            gate.kind = GateKind::H;
            gate.target = target;
            return append(std::move(gate));
        }

        Gate &rx(IdxType target, ValType angle)
        {
            Gate gate{};
            gate.kind = GateKind::RX;
            gate.target = target;
            gate.set_angle(angle);
            return append(std::move(gate));
        }

        Gate &rx_param(IdxType target, size_t param_index)
        {
            parameter_at(param_index);
            Gate gate{};
            gate.kind = GateKind::RX;
            gate.target = target;
            gate.set_parameter_index(param_index);
            return append(std::move(gate));
        }

        Gate &rz(IdxType target, ValType angle)
        {
            Gate gate{};
            gate.kind = GateKind::RZ;
            gate.target = target;
            gate.set_angle(angle);
            return append(std::move(gate));
        }

        Gate &rz_param(IdxType target, size_t param_index)
        {
            parameter_at(param_index);
            Gate gate{};
            gate.kind = GateKind::RZ;
            gate.target = target;
            gate.set_parameter_index(param_index);
            return append(std::move(gate));
        }

        Gate &cx(IdxType control, IdxType target)
        {
            Gate gate{};
            gate.kind = GateKind::CX;
            gate.control = control;
            gate.target = target;
            return append(std::move(gate));
        }

        Gate &measure(IdxType target)
        {
            Gate gate{};
            gate.kind = GateKind::Measure;
            gate.target = target;
            return append(std::move(gate));
        }

        Gate &measure_all(IdxType shots)
        {
            Gate gate{};
            gate.kind = GateKind::MeasureAll;
            gate.repetition = shots;
            return append(std::move(gate));
        }

        Gate &reset(IdxType target)
        {
            Gate gate{};
            gate.kind = GateKind::Reset;
            gate.target = target;
            return append(std::move(gate));
        }

    private:
        ValType &parameter_at(size_t index)
        {
            if (index >= parameters_.size())
                throw std::out_of_range("Parameter index out of range");
            return parameters_[index];
        }

        const ValType &parameter_at(size_t index) const
        {
            if (index >= parameters_.size())
                throw std::out_of_range("Parameter index out of range");
            return parameters_[index];
        }

        Gate &append(Gate &&gate)
        {
            gates_.push_back(std::move(gate));
            return gates_.back();
        }

        IdxType n_qubits_ = 0;
        std::vector<Gate> gates_;
        std::vector<ValType> parameters_;
    };
}
