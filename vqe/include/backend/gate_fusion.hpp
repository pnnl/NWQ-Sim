#pragma once

#include <cstddef>
#include <vector>

#include "backend/sim_gate.hpp"

namespace vqe::backend {

void fuse_simulation_gates(const std::vector<sim_gate>& gates,
                           std::size_t n_qubits,
                           std::vector<sim_gate>& buffer,
                           std::vector<sim_gate>& tmp1,
                           std::vector<sim_gate>& tmp2,
                           std::vector<sim_gate>& tmp3,
                           std::vector<sim_gate>& chunk,
                           std::vector<sim_gate>& fused);

std::vector<sim_gate> fuse_simulation_gates(const std::vector<sim_gate>& gates,
                                            std::size_t n_qubits);

}  // namespace vqe::backend
