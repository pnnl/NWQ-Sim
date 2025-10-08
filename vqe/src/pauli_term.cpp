#include "pauli_term.hpp"

#include <algorithm>
#include <cmath>

namespace vqe {
namespace {

constexpr char symbol_from_masks(bool x, bool z) {
  if (!x && !z) {
    return 'I';
  }
  if (x && !z) {
    return 'X';
  }
  if (!x && z) {
    return 'Z';
  }
  return 'Y';
}

}  // namespace

std::string pauli_to_string(std::uint64_t x_mask, std::uint64_t z_mask, std::size_t n_qubits) {
  std::string label;
  label.reserve(n_qubits ? n_qubits : 1);
  for (std::size_t idx = n_qubits; idx > 0; --idx) {
    const std::size_t qubit = idx - 1;
    const std::uint64_t bit = static_cast<std::uint64_t>(1) << qubit;
    const bool x = (x_mask & bit) != 0;
    const bool z = (z_mask & bit) != 0;
    label.push_back(symbol_from_masks(x, z));
  }
  if (label.empty()) {
    label.push_back('I');
  }
  return label;
}

std::vector<pauli_term> normalize_terms(const pauli_accumulator& accum,
                                        std::size_t n_qubits,
                                        double cutoff) {
  std::vector<pauli_term> result;
  result.reserve(accum.size());
  for (const auto& entry : accum) {
    const auto& coeff = entry.second;
    if (std::abs(coeff.real()) < cutoff && std::abs(coeff.imag()) < cutoff) {
      continue;
    }
    result.push_back(pauli_term{entry.first.first, entry.first.second, coeff});
  }
  std::sort(result.begin(), result.end(), [](const pauli_term& lhs, const pauli_term& rhs) {
    if (lhs.x_mask == rhs.x_mask) {
      if (lhs.z_mask == rhs.z_mask) {
        if (lhs.coefficient.real() == rhs.coefficient.real()) {
          return lhs.coefficient.imag() < rhs.coefficient.imag();
        }
        return lhs.coefficient.real() < rhs.coefficient.real();
      }
      return lhs.z_mask < rhs.z_mask;
    }
    return lhs.x_mask < rhs.x_mask;
  });
  return result;
}

}  // namespace vqe
