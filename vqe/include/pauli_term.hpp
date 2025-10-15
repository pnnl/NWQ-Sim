#pragma once

#include <complex>
#include <cstdint>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace vqe {

struct pauli_term {
  std::uint64_t x_mask = 0;
  std::uint64_t z_mask = 0;
  std::complex<double> coefficient{0.0};
};

using pauli_accumulator = std::map<std::pair<std::uint64_t, std::uint64_t>, std::complex<double>>;

std::string pauli_to_string(std::uint64_t x_mask, std::uint64_t z_mask, std::size_t n_qubits);

std::vector<pauli_term> normalize_terms(const pauli_accumulator& accum, std::size_t n_qubits, double cutoff = 1e-12);

}  // namespace vqe
