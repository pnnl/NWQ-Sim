#pragma once

#include <cstddef>

namespace vqe {

enum class spin { up, down };

enum class orbital_kind { occupied, virtual_orbital };

enum class operator_kind { creation, annihilation };

struct molecular_environment {
  std::size_t n_spatial = 0;
  std::size_t n_electrons = 0;
  bool xacc_indexing = false;
  double constant_energy = 0.0;

  [[nodiscard]] std::size_t occupied_orbitals() const {
    return n_electrons / 2;
  }

  [[nodiscard]] std::size_t virtual_orbitals() const {
    return n_spatial - occupied_orbitals();
  }

  [[nodiscard]] std::size_t total_qubits() const {
    return 2 * n_spatial;
  }

  [[nodiscard]] std::size_t qubit_index(std::size_t orbital_index,
                                        orbital_kind kind,
                                        spin spin_state) const;
};

struct fermion_operator {
  std::size_t orbital_index = 0;
  orbital_kind orbital = orbital_kind::occupied;
  spin spin_state = spin::up;
  operator_kind op = operator_kind::annihilation;

  fermion_operator() = default;

  fermion_operator(std::size_t orbital_idx,
                   orbital_kind orbital_kind_value,
                   spin spin_value,
                   operator_kind op_kind)
      : orbital_index(orbital_idx),
        orbital(orbital_kind_value),
        spin_state(spin_value),
        op(op_kind) {}
};

}  // namespace vqe

