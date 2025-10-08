#pragma once

#include <complex>
#include <cstddef>
#include <vector>

#include "core/circuit.hpp"
#include "backend/sim_gate.hpp"

namespace vqe::backend {

class statevector_cpu {
 public:
  explicit statevector_cpu(std::size_t num_qubits);

  void reset();
  void apply(const circuit& circ);

  [[nodiscard]] std::complex<double> expectation(const pauli_term& term) const;
  [[nodiscard]] std::complex<double> expectation(const std::vector<pauli_term>& terms) const;

  [[nodiscard]] const std::vector<std::complex<double>>& amplitudes() const;
  void load_state(const std::vector<std::complex<double>>& amplitudes);

  [[nodiscard]] std::size_t num_qubits() const { return num_qubits_; }

  void copy_state_from(const statevector_cpu& other);

 private:
  void apply_gate(const sim_gate& gate);
  void apply_single(const sim_gate& gate);
  void apply_two(const sim_gate& gate);
  void apply_pauli(const sim_gate& gate);

  void apply_h_internal(std::size_t qubit);
  void apply_rz_internal(std::size_t qubit, double angle);
  void apply_cnot_internal(std::size_t control, std::size_t target);
  void apply_y_basis_forward(std::size_t qubit);
  void apply_y_basis_backward(std::size_t qubit);

  std::size_t num_qubits_ = 0;
  std::size_t dimension_ = 0;
  std::size_t half_dimension_ = 0;
  std::vector<double> state_real_;
  std::vector<double> state_imag_;
  mutable std::vector<std::complex<double>> amplitude_cache_;
  mutable bool amplitudes_valid_ = false;

  std::vector<sim_gate> gate_buffer_;
  std::vector<sim_gate> fuse_tmp1_;
  std::vector<sim_gate> fuse_tmp2_;
  std::vector<sim_gate> fuse_tmp3_;
  std::vector<sim_gate> fuse_buffer_;
  std::vector<sim_gate> fuse_chunk_;
  std::vector<sim_gate> fused_buffer_;

};

}  // namespace vqe::backend
