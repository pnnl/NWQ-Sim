#pragma once

#include <string>

#include "fermion_term.hpp"

namespace vqe {

hamiltonian_data read_hamiltonian_file(const std::string& path);

}  // namespace vqe
