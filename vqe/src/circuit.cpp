#include "core/circuit.hpp"

#include <algorithm>
#include <cmath>

namespace vqe {

namespace {

double evaluate_terms(const std::vector<double>& values, const parameter_expression& expr) {
  double result = 0.0;
  for (const auto& term : expr.terms) {
    if (term.parameter_index >= values.size()) {
      throw std::out_of_range("parameter index out of range while evaluating expression");
    }
    result += values[term.parameter_index] * term.coefficient;
  }
  return expr.scale * result;
}

}  // namespace

circuit::circuit(std::size_t qubit_count) : qubit_count_(qubit_count) {}

std::size_t circuit::num_qubits() const {
  return qubit_count_;
}

void circuit::set_num_qubits(std::size_t qubits) {
  qubit_count_ = qubits;
}

std::size_t circuit::add_parameter(double initial_value) {
  const std::size_t index = parameters_.size();
  parameters_.push_back(initial_value);
  return index;
}

void circuit::set_parameter(std::size_t index, double value) {
  validate_parameter(index);
  parameters_[index] = value;
}

double circuit::parameter(std::size_t index) const {
  validate_parameter(index);
  return parameters_[index];
}

const std::vector<double>& circuit::parameters() const {
  return parameters_;
}

void circuit::add_x(std::size_t qubit) {
  validate_qubit(qubit);
  gate g{};
  g.kind = gate_kind::x;
  g.qubits = {qubit};
  gates_.push_back(std::move(g));
}

void circuit::add_h(std::size_t qubit) {
  validate_qubit(qubit);
  gate g{};
  g.kind = gate_kind::h;
  g.qubits = {qubit};
  gates_.push_back(std::move(g));
}

void circuit::add_cnot(std::size_t control, std::size_t target) {
  validate_qubit(control);
  validate_qubit(target);
  if (control == target) {
    throw std::invalid_argument("control and target cannot be the same for CNOT");
  }
  gate g{};
  g.kind = gate_kind::cnot;
  g.qubits = {control, target};
  gates_.push_back(std::move(g));
}

void circuit::add_fixed_rz(std::size_t qubit, double angle) {
  validate_qubit(qubit);
  gate g{};
  g.kind = gate_kind::rz;
  g.qubits = {qubit};
  g.angle = angle;
  gates_.push_back(std::move(g));
}

void circuit::add_parametric_rz(std::size_t qubit, const parameter_expression& expr) {
  validate_qubit(qubit);
  if (expr.empty()) {
    throw std::invalid_argument("parameter expression cannot be empty for parametric RZ gate");
  }
  gate g{};
  g.kind = gate_kind::rz;
  g.qubits = {qubit};
  g.parametrized = true;
  g.expression = expr;
  gates_.push_back(std::move(g));
}

void circuit::add_pauli_rotation(const pauli_term& term,
                                 const parameter_expression& expr,
                                 double global_scale) {
  if (expr.empty()) {
    throw std::invalid_argument("parameter expression cannot be empty for Pauli rotation");
  }
  gate g{};
  g.kind = gate_kind::pauli_rotation;
  g.pauli = term;
  g.parametrized = true;
  g.expression = expr;
  g.expression.scale *= global_scale;
  g.qubits = non_trivial_qubits(term);
  gates_.push_back(std::move(g));
}

const std::vector<gate>& circuit::gates() const {
  return gates_;
}

double circuit::evaluate(const parameter_expression& expr) const {
  return evaluate_terms(parameters_, expr);
}

void circuit::validate_qubit(std::size_t qubit) const {
  if (qubit >= qubit_count_) {
    throw std::out_of_range("qubit index out of range");
  }
}

void circuit::validate_parameter(std::size_t index) const {
  if (index >= parameters_.size()) {
    throw std::out_of_range("parameter index out of range");
  }
}

std::vector<std::size_t> non_trivial_qubits(const pauli_term& term) {
  std::vector<std::size_t> result;
  const std::size_t max_bits = 8 * sizeof(std::uint64_t);
  for (std::size_t idx = 0; idx < max_bits; ++idx) {
    const std::uint64_t bit = static_cast<std::uint64_t>(1) << idx;
    if ((term.x_mask & bit) != 0 || (term.z_mask & bit) != 0) {
      result.push_back(idx);
    }
  }
  return result;
}

}  // namespace vqe
