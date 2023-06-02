#pragma once

#include "util.hpp"
#include "circuit.hpp"
#include <vector>

namespace NWQSim
{
    // Base class
    class QuantumState
    {
    public:
        QuantumState(IdxType _n_qubits) {} // constructor
        virtual ~QuantumState() {}         // virtual destructor

        virtual void reset_state() = 0;
        virtual void set_seed(IdxType seed) = 0;

        virtual void sim(Circuit *circuit) = 0;
        virtual IdxType *get_results() = 0;
        virtual IdxType measure(IdxType qubit) = 0;
        virtual IdxType *measure_all(IdxType repetition) = 0;

        virtual ValType get_exp_z(const std::vector<size_t> &in_bits) = 0;

        virtual void print_res_sv() = 0;
    };

} // namespace NWQSim