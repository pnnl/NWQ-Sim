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

        virtual void reset_sim() = 0;
        virtual void set_seed(IdxType seed) = 0;

        virtual std::vector<IdxType> sim(Circuit *circuit) = 0;
        virtual std::vector<IdxType> measure(IdxType qubit) = 0;
        virtual std::vector<IdxType> measure_all(IdxType repetition) = 0;
    };

} // namespace NWQSim