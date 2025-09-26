#pragma once

#include "nwqsim_core.hpp"

#include <vector>

namespace NWQSim
{
    struct SimulationMetadata
    {
        IdxType n_qubits = 0;
        IdxType total_device_gates = 0;
        IdxType single_qubit_gates = 0;
        IdxType two_qubit_gates = 0;
        IdxType measurements = 0;
        IdxType measure_all_ops = 0;
        IdxType measure_all_shots = 0;
        IdxType resets = 0;
    };

    struct MeasureAllResult
    {
        IdxType repetitions = 0;
        std::vector<IdxType> outcomes;
    };

struct SimulationResult
{
    SimulationMetadata metadata;
    std::vector<IdxType> single_qubit_measurements;
    std::vector<MeasureAllResult> measure_all_results;
};

}
