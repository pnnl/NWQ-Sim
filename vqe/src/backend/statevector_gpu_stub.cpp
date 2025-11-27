#include "backend/statevector_gpu.hpp"

#include <stdexcept>

namespace vqe::backend {

statevector_gpu::statevector_gpu(std::size_t num_qubits) : num_qubits_(num_qubits) {
  throw std::runtime_error("vqe built without CUDA/HIP support; GPU backend unavailable");
}

statevector_gpu::~statevector_gpu() = default;

void statevector_gpu::reset() {
  throw std::runtime_error("vqe built without CUDA/HIP support; GPU backend unavailable");
}

void statevector_gpu::apply(const circuit&) {
  throw std::runtime_error("vqe built without CUDA/HIP support; GPU backend unavailable");
}

std::complex<double> statevector_gpu::expectation(const pauli_term&) const {
  throw std::runtime_error("vqe built without CUDA/HIP support; GPU backend unavailable");
}

std::complex<double> statevector_gpu::expectation(const std::vector<pauli_term>&) const {
  throw std::runtime_error("vqe built without CUDA/HIP support; GPU backend unavailable");
}

const std::vector<std::complex<double>>& statevector_gpu::amplitudes() const {
  throw std::runtime_error("vqe built without CUDA/HIP support; GPU backend unavailable");
}

void statevector_gpu::load_state(const std::vector<std::complex<double>>&) {
  throw std::runtime_error("vqe built without CUDA/HIP support; GPU backend unavailable");
}

void statevector_gpu::copy_state_from(const statevector_gpu&) {
  throw std::runtime_error("vqe built without CUDA/HIP support; GPU backend unavailable");
}

std::size_t statevector_gpu::num_qubits() const {
  return num_qubits_;
}

}  // namespace vqe::backend
