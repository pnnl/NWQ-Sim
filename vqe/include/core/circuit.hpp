#pragma once

#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "pauli_term.hpp"

namespace vqe {

enum class gate_kind {
  x,
  h,
  cnot,
  rz,
  pauli_rotation
};

struct parameter_term {
  std::size_t parameter_index = 0;
  double coefficient = 0.0;
};

struct parameter_expression {
  std::vector<parameter_term> terms;
  double scale = 1.0;

  [[nodiscard]] bool empty() const {
    return terms.empty();
  }
};

struct gate {
  gate_kind kind{};
  std::vector<std::size_t> qubits;
  double angle = 0.0;
  bool parametrized = false;
  parameter_expression expression{};
  pauli_term pauli{};
};

class circuit {
 public:
  explicit circuit(std::size_t qubit_count = 0);

  [[nodiscard]] std::size_t num_qubits() const;
  void set_num_qubits(std::size_t qubits);

  [[nodiscard]] std::size_t add_parameter(double initial_value = 0.0);
  void set_parameter(std::size_t index, double value);
  [[nodiscard]] double parameter(std::size_t index) const;
  [[nodiscard]] const std::vector<double>& parameters() const;

  void add_x(std::size_t qubit);
  void add_h(std::size_t qubit);
  void add_cnot(std::size_t control, std::size_t target);
  void add_fixed_rz(std::size_t qubit, double angle);
  void add_parametric_rz(std::size_t qubit, const parameter_expression& expr);
  void add_pauli_rotation(const pauli_term& term, const parameter_expression& expr, double global_scale);

  [[nodiscard]] const std::vector<gate>& gates() const;

  [[nodiscard]] double evaluate(const parameter_expression& expr) const;

 private:
  void validate_qubit(std::size_t qubit) const;
  void validate_parameter(std::size_t index) const;

  std::size_t qubit_count_ = 0;
  std::vector<gate> gates_;
  std::vector<double> parameters_;
};

std::vector<std::size_t> non_trivial_qubits(const pauli_term& term);

}  // namespace vqe

