#pragma once

#include <vector>

#include "backend/sim_gate.hpp"
#include "core/circuit.hpp"

namespace vqe::backend {

void build_simulation_gates(const circuit& circ, std::vector<sim_gate>& out);

std::vector<sim_gate> build_simulation_gates(const circuit& circ);

}  // namespace vqe::backend
