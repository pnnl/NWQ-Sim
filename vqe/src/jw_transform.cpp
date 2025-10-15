#include "jw_transform.hpp"

#include <array>
#include <stdexcept>

namespace vqe {
namespace {

constexpr int sign_relations[4][4] = {
    {1, 1, 1, 1},
    {1, 1, -1, 1},
    {1, 1, 1, -1},
    {1, -1, 1, 1}};

pauli_term multiply_pauli(const pauli_term& lhs,
                          const pauli_term& rhs,
                          std::size_t n_qubits) {
  pauli_term result{};
  result.x_mask = lhs.x_mask ^ rhs.x_mask;
  result.z_mask = lhs.z_mask ^ rhs.z_mask;
  result.coefficient = lhs.coefficient * rhs.coefficient;

  const std::uint64_t anticommute_mask = (lhs.x_mask & rhs.z_mask) ^ (rhs.x_mask & lhs.z_mask);
  for (std::size_t qubit = 0; qubit < n_qubits; ++qubit) {
    const std::uint64_t bit = static_cast<std::uint64_t>(1) << qubit;
    if ((anticommute_mask & bit) == 0) {
      continue;
    }
    const int x1 = (lhs.x_mask & bit) ? 1 : 0;
    const int z1 = (lhs.z_mask & bit) ? 1 : 0;
    const int x2 = (rhs.x_mask & bit) ? 1 : 0;
    const int z2 = (rhs.z_mask & bit) ? 1 : 0;
    result.coefficient *= std::complex<double>(0.0, sign_relations[x1 + 2 * z1][x2 + 2 * z2]);
  }
  return result;
}

std::array<pauli_term, 2> single_op_terms(const fermion_op& op) {
  const std::uint64_t bit = static_cast<std::uint64_t>(1) << op.index;
  const std::uint64_t z_lower = (bit == 0) ? 0 : (bit - 1);
  const std::uint64_t z_with_self = z_lower | bit;
  const int sign = (op.kind == fermion_op_kind::annihilation) ? 1 : -1;

  return {pauli_term{bit, z_lower, 0.5},
          pauli_term{bit, z_with_self, std::complex<double>(0.0, 0.5 * sign)}};
}

void expand_terms(const std::vector<std::array<pauli_term, 2>>& per_op_terms,
                  std::size_t index,
                  std::size_t n_qubits,
                  const pauli_term& current,
                  std::vector<pauli_term>& out) {
  if (index == per_op_terms.size()) {
    out.push_back(current);
    return;
  }
  for (const auto& candidate : per_op_terms[index]) {
    auto next = multiply_pauli(current, candidate, n_qubits);
    expand_terms(per_op_terms, index + 1, n_qubits, next, out);
  }
}

}  // namespace

std::vector<pauli_term> jordan_wigner_transform(const hamiltonian_data& data) {
  if (data.num_qubits() > 63) {
    throw std::runtime_error("jordan-wigner mapping currently supports up to 63 qubits");
  }

  const std::size_t n_qubits = data.num_qubits();
  pauli_accumulator accumulator;

  if (std::abs(data.constant.real()) > 0.0 || std::abs(data.constant.imag()) > 0.0) {
    accumulator[{0, 0}] += data.constant;
  }

  for (const auto& term : data.terms) {
    if (term.operators.empty()) {
      accumulator[{0, 0}] += term.coefficient;
      continue;
    }

    std::vector<std::array<pauli_term, 2>> per_op_terms;
    per_op_terms.reserve(term.operators.size());
    for (const auto& op : term.operators) {
      per_op_terms.push_back(single_op_terms(op));
    }

    std::vector<pauli_term> expanded;
    pauli_term identity{};
    identity.coefficient = 1.0;
    expand_terms(per_op_terms, 0, n_qubits, identity, expanded);

    for (auto entry : expanded) {
      entry.coefficient *= term.coefficient;
      accumulator[{entry.x_mask, entry.z_mask}] += entry.coefficient;
    }
  }

  return normalize_terms(accumulator, n_qubits);
}

}  // namespace vqe
