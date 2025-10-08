#pragma once

#include <vector>

#include "fermion_term.hpp"
#include "pauli_term.hpp"

namespace vqe {

std::vector<pauli_term> jordan_wigner_transform(const hamiltonian_data& data);

}  // namespace vqe
