#pragma once

#include <array>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "pauli_term.hpp"

namespace vqe::backend {

struct sim_gate {
  enum class kind : std::uint8_t { single = 0, two = 1, pauli = 2 };

  kind op = kind::single;
  std::uint8_t reserved[7] = {};
  std::uint32_t target = 0;
  std::uint32_t control = invalid_index;
  std::array<double, 16> real{};
  std::array<double, 16> imag{};
  pauli_term pauli{};
  double angle = 0.0;
  std::array<std::uint32_t, 32> pauli_qubits{};
  std::array<std::uint8_t, 32> pauli_basis{};
  std::uint32_t pauli_qubit_count = 0;

  [[nodiscard]] bool is_single() const { return op == kind::single; }
  [[nodiscard]] bool is_two_qubit() const { return op == kind::two; }
  [[nodiscard]] std::size_t matrix_size() const {
    if (op == kind::single) {
      return 4;
    }
    if (op == kind::two) {
      return 16;
    }
    return 0;
  }

  static constexpr std::uint32_t invalid_index = std::numeric_limits<std::uint32_t>::max();
};

static_assert(std::is_trivially_copyable_v<sim_gate>, "sim_gate must be trivially copyable");

}  // namespace vqe::backend
