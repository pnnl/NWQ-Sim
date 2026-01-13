#pragma once

#include <complex>
#include <cstddef>
#include <memory>
#include <vector>

#include "core/circuit.hpp"
#include "backend/sim_gate.hpp"

namespace vqe::backend {

class statevector_gpu {
 public:
  explicit statevector_gpu(std::size_t num_qubits);
  ~statevector_gpu();

  void reset();
  void apply(const circuit& circ);

  [[nodiscard]] std::complex<double> expectation(const pauli_term& term) const;
  [[nodiscard]] std::complex<double> expectation(const std::vector<pauli_term>& terms) const;

  [[nodiscard]] const std::vector<std::complex<double>>& amplitudes() const;
  void load_state(const std::vector<std::complex<double>>& amplitudes);
  void copy_state_from(const statevector_gpu& other);

  [[nodiscard]] std::size_t num_qubits() const;

 private:
#if defined(VQE_ENABLE_CUDA) || defined(VQE_ENABLE_HIP)
  struct impl;
  std::unique_ptr<impl> impl_;
#else
  std::size_t num_qubits_ = 0;
#endif
};

}  // namespace vqe::backend
