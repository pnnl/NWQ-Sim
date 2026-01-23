#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace vqe {

enum class fermion_op_kind {
  creation,
  annihilation
};

struct fermion_op {
  std::size_t index;
  fermion_op_kind kind;
};

struct fermion_term {
  std::complex<double> coefficient{};
  std::vector<fermion_op> operators;
};

struct hamiltonian_data {
  std::complex<double> constant{};
  std::vector<fermion_term> terms;
  std::size_t max_index = 0;

  [[nodiscard]] std::size_t num_qubits() const {
    return max_index + 1;
  }
};

}  // namespace vqe
