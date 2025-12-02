#include "backend/statevector_cpu.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

#include "backend/gate_builder.hpp"
#include "backend/gate_fusion.hpp"

namespace vqe::backend {
namespace {
constexpr double kHalfPi = 1.57079632679489661923132169163975144;
constexpr double kSqrtHalf = 0.707106781186547524400844362104849039;
constexpr double kHalf = 0.5;
constexpr std::size_t kParallelThreshold = static_cast<std::size_t>(1) << 12;  // 4096 amplitudes
}  // namespace

statevector_cpu::statevector_cpu(std::size_t num_qubits) : num_qubits_(num_qubits) {
  dimension_ = num_qubits_ == 0 ? 1 : (static_cast<std::size_t>(1) << num_qubits_);
  half_dimension_ = num_qubits_ < 1 ? 0 : (dimension_ >> 1);
  state_real_.assign(dimension_, 0.0);
  state_imag_.assign(dimension_, 0.0);
  if (!state_real_.empty()) {
    state_real_[0] = 1.0;
  }
  amplitudes_valid_ = false;
}

void statevector_cpu::reset() {
  std::fill(state_real_.begin(), state_real_.end(), 0.0);
  std::fill(state_imag_.begin(), state_imag_.end(), 0.0);
  if (!state_real_.empty()) {
    state_real_[0] = 1.0;
  }
  amplitudes_valid_ = false;
}

void statevector_cpu::apply(const circuit& circ) {
  if (circ.num_qubits() != num_qubits_) {
    throw std::invalid_argument("circuit qubit count does not match statevector");
  }

  build_simulation_gates(circ, gate_buffer_);
  fuse_simulation_gates(gate_buffer_, num_qubits_, fuse_buffer_, fuse_tmp1_, fuse_tmp2_, fuse_tmp3_, fuse_chunk_, fused_buffer_);

  for (const auto& gate : fused_buffer_) {
    apply_gate(gate);
  }

  amplitudes_valid_ = false;
}

void statevector_cpu::apply_gate(const sim_gate& gate) {
  switch (gate.op) {
    case sim_gate::kind::single:
      apply_single(gate);
      break;
    case sim_gate::kind::two:
      apply_two(gate);
      break;
    case sim_gate::kind::pauli:
      apply_pauli(gate);
      break;
  }
}

void statevector_cpu::apply_single(const sim_gate& gate) {
  const std::size_t qubit = static_cast<std::size_t>(gate.target);
  const std::size_t block = static_cast<std::size_t>(1) << qubit;
  const std::size_t stride = block << 1;

  const double m00r = gate.real[0];
  const double m00i = gate.imag[0];
  const double m01r = gate.real[1];
  const double m01i = gate.imag[1];
  const double m10r = gate.real[2];
  const double m10i = gate.imag[2];
  const double m11r = gate.real[3];
  const double m11i = gate.imag[3];

  #ifdef NWQSIM_USE_OPENMP
  #pragma omp parallel for schedule(static) if (dimension_ >= kParallelThreshold)
  #endif
  for (std::size_t base = 0; base < dimension_; base += stride) {
    #ifdef NWQSIM_USE_OPENMP
    #pragma omp simd
    #endif
    for (std::size_t offset = 0; offset < block; ++offset) {
      const std::size_t i0 = base + offset;
      const std::size_t i1 = i0 + block;

      const double a_r = state_real_[i0];
      const double a_i = state_imag_[i0];
      const double b_r = state_real_[i1];
      const double b_i = state_imag_[i1];

      const double out0_r = m00r * a_r - m00i * a_i + m01r * b_r - m01i * b_i;
      const double out0_i = m00r * a_i + m00i * a_r + m01r * b_i + m01i * b_r;

      const double out1_r = m10r * a_r - m10i * a_i + m11r * b_r - m11i * b_i;
      const double out1_i = m10r * a_i + m10i * a_r + m11r * b_i + m11i * b_r;

      state_real_[i0] = out0_r;
      state_imag_[i0] = out0_i;
      state_real_[i1] = out1_r;
      state_imag_[i1] = out1_i;
    }
  }
}

void statevector_cpu::apply_two(const sim_gate& gate) {
  if (gate.control == sim_gate::invalid_index || gate.target == sim_gate::invalid_index) {
    throw std::invalid_argument("two-qubit gate requires valid indices");
  }

  const std::size_t qubit0 = static_cast<std::size_t>(gate.control);
  const std::size_t qubit1 = static_cast<std::size_t>(gate.target);
  if (qubit0 == qubit1) {
    throw std::invalid_argument("control and target qubits must differ for two-qubit gates");
  }

  const std::size_t per_work = dimension_ >> 2;
  const std::size_t q0dim = static_cast<std::size_t>(1) << std::max(qubit0, qubit1);
  const std::size_t q1dim = static_cast<std::size_t>(1) << std::min(qubit0, qubit1);
  const std::size_t mider_factor = (q0dim + q1dim + q1dim - 1) >> (std::min(qubit0, qubit1) + 1);
  const std::size_t inner_factor = q1dim;
  const std::size_t qubit0_dim = static_cast<std::size_t>(1) << qubit0;
  const std::size_t qubit1_dim = static_cast<std::size_t>(1) << qubit1;

  #ifdef NWQSIM_USE_OPENMP
  #pragma omp parallel for schedule(static) if (dimension_ >= kParallelThreshold)
  #endif
  for (std::size_t i = 0; i < per_work; ++i) {
    const std::size_t outer = ((i / inner_factor) / mider_factor) * (q0dim + q0dim);
    const std::size_t mider = ((i / inner_factor) % mider_factor) * (q1dim + q1dim);
    const std::size_t inner = i % inner_factor;

    const std::size_t pos0 = outer + mider + inner;
    const std::size_t pos1 = pos0 + qubit1_dim;
    const std::size_t pos2 = pos0 + qubit0_dim;
    const std::size_t pos3 = pos0 + q0dim + q1dim;

    const double vec_real[4] = {state_real_[pos0], state_real_[pos1], state_real_[pos2], state_real_[pos3]};
    const double vec_imag[4] = {state_imag_[pos0], state_imag_[pos1], state_imag_[pos2], state_imag_[pos3]};

    double out_real[4] = {0.0};
    double out_imag[4] = {0.0};

    for (std::size_t row = 0; row < 4; ++row) {
      double r = 0.0;
      double im = 0.0;
      for (std::size_t col = 0; col < 4; ++col) {
        const double mr = gate.real[row * 4 + col];
        const double mi = gate.imag[row * 4 + col];
        const double vr = vec_real[col];
        const double vi = vec_imag[col];
        r += mr * vr - mi * vi;
        im += mr * vi + mi * vr;
      }
      out_real[row] = r;
      out_imag[row] = im;
    }

    state_real_[pos0] = out_real[0];
    state_imag_[pos0] = out_imag[0];
    state_real_[pos1] = out_real[1];
    state_imag_[pos1] = out_imag[1];
    state_real_[pos2] = out_real[2];
    state_imag_[pos2] = out_imag[2];
    state_real_[pos3] = out_real[3];
    state_imag_[pos3] = out_imag[3];
  }
}

void statevector_cpu::apply_h_internal(std::size_t qubit) {
  const std::size_t block = static_cast<std::size_t>(1) << qubit;
  const std::size_t stride = block << 1;
  #ifdef NWQSIM_USE_OPENMP
  #pragma omp parallel for schedule(static) if (dimension_ >= kParallelThreshold)
  #endif
  for (std::size_t base = 0; base < dimension_; base += stride) {
    #ifdef NWQSIM_USE_OPENMP
    #pragma omp simd
    #endif
    for (std::size_t offset = 0; offset < block; ++offset) {
      const std::size_t i0 = base + offset;
      const std::size_t i1 = i0 + block;
      const double a_r = state_real_[i0];
      const double a_i = state_imag_[i0];
      const double b_r = state_real_[i1];
      const double b_i = state_imag_[i1];
      state_real_[i0] = kSqrtHalf * (a_r + b_r);
      state_imag_[i0] = kSqrtHalf * (a_i + b_i);
      state_real_[i1] = kSqrtHalf * (a_r - b_r);
      state_imag_[i1] = kSqrtHalf * (a_i - b_i);
    }
  }
}

void statevector_cpu::apply_rz_internal(std::size_t qubit, double angle) {
  const double half = angle * 0.5;
  const double c = std::cos(half);
  const double s = std::sin(half);
  const std::size_t block = static_cast<std::size_t>(1) << qubit;
  const std::size_t stride = block << 1;
  #ifdef NWQSIM_USE_OPENMP
  #pragma omp parallel for schedule(static) if (dimension_ >= kParallelThreshold)
  #endif
  for (std::size_t base = 0; base < dimension_; base += stride) {
    #ifdef NWQSIM_USE_OPENMP
    #pragma omp simd
    #endif
    for (std::size_t offset = 0; offset < block; ++offset) {
      const std::size_t i0 = base + offset;
      const std::size_t i1 = i0 + block;

      const double a_r = state_real_[i0];
      const double a_i = state_imag_[i0];
      const double b_r = state_real_[i1];
      const double b_i = state_imag_[i1];

      state_real_[i0] = c * a_r + s * a_i;
      state_imag_[i0] = c * a_i - s * a_r;
      state_real_[i1] = c * b_r - s * b_i;
      state_imag_[i1] = c * b_i + s * b_r;
    }
  }
}

void statevector_cpu::apply_cnot_internal(std::size_t control, std::size_t target) {
  if (control == target) {
    throw std::invalid_argument("control and target qubits must differ for CNOT");
  }
  const std::size_t control_mask = static_cast<std::size_t>(1) << control;
  const std::size_t target_mask = static_cast<std::size_t>(1) << target;
  #ifdef NWQSIM_USE_OPENMP
  #pragma omp parallel for schedule(static) if (dimension_ >= kParallelThreshold)
  #endif
  for (std::size_t index = 0; index < dimension_; ++index) {
    if ((index & control_mask) == 0 || (index & target_mask) != 0) {
      continue;
    }
    const std::size_t partner = index | target_mask;
    std::swap(state_real_[index], state_real_[partner]);
    std::swap(state_imag_[index], state_imag_[partner]);
  }
}

void statevector_cpu::apply_y_basis_forward(std::size_t qubit) {
  const std::size_t block = static_cast<std::size_t>(1) << qubit;
  const std::size_t stride = block << 1;
  const double c0r = kHalf;
  const double c0i = -kHalf;
  const double c1r = kHalf;
  const double c1i = kHalf;

  #ifdef NWQSIM_USE_OPENMP
  #pragma omp parallel for schedule(static) if (dimension_ >= kParallelThreshold)
  #endif
  for (std::size_t base = 0; base < dimension_; base += stride) {
    #ifdef NWQSIM_USE_OPENMP
    #pragma omp simd
    #endif
    for (std::size_t offset = 0; offset < block; ++offset) {
      const std::size_t i0 = base + offset;
      const std::size_t i1 = i0 + block;

      const double a_r = state_real_[i0];
      const double a_i = state_imag_[i0];
      const double b_r = state_real_[i1];
      const double b_i = state_imag_[i1];

      const double a0_r = c0r * a_r - c0i * a_i;
      const double a0_i = c0r * a_i + c0i * a_r;
      const double b1_r = c1r * b_r - c1i * b_i;
      const double b1_i = c1r * b_i + c1i * b_r;

      state_real_[i0] = a0_r + b1_r;
      state_imag_[i0] = a0_i + b1_i;
      state_real_[i1] = a0_r - b1_r;
      state_imag_[i1] = a0_i - b1_i;
    }
  }
}

void statevector_cpu::apply_y_basis_backward(std::size_t qubit) {
  const std::size_t block = static_cast<std::size_t>(1) << qubit;
  const std::size_t stride = block << 1;
  const double c0r = kHalf;
  const double c0i = kHalf;
  const double c1r = kHalf;
  const double c1i = -kHalf;

  #ifdef NWQSIM_USE_OPENMP
  #pragma omp parallel for schedule(static) if (dimension_ >= kParallelThreshold)
  #endif
  for (std::size_t base = 0; base < dimension_; base += stride) {
    #ifdef NWQSIM_USE_OPENMP
    #pragma omp simd
    #endif
    for (std::size_t offset = 0; offset < block; ++offset) {
      const std::size_t i0 = base + offset;
      const std::size_t i1 = i0 + block;

      const double a_r = state_real_[i0];
      const double a_i = state_imag_[i0];
      const double b_r = state_real_[i1];
      const double b_i = state_imag_[i1];

      const double sum_r = a_r + b_r;
      const double sum_i = a_i + b_i;
      const double diff_r = a_r - b_r;
      const double diff_i = a_i - b_i;

      state_real_[i0] = c0r * sum_r - c0i * sum_i;
      state_imag_[i0] = c0r * sum_i + c0i * sum_r;
      state_real_[i1] = c1r * diff_r - c1i * diff_i;
      state_imag_[i1] = c1r * diff_i + c1i * diff_r;
    }
  }
}

void statevector_cpu::apply_pauli(const sim_gate& gate) {
  const std::uint64_t support = gate.pauli.x_mask | gate.pauli.z_mask;
  if (support == 0) {
    const double half = gate.angle * 0.5;
    const double c = std::cos(half);
    const double s = std::sin(half);
    #ifdef NWQSIM_USE_OPENMP
    #pragma omp parallel for schedule(static) if (dimension_ >= kParallelThreshold)
    #endif
    for (std::size_t idx = 0; idx < dimension_; ++idx) {
      const double r = state_real_[idx];
      const double i = state_imag_[idx];
      state_real_[idx] = c * r + s * i;
      state_imag_[idx] = -s * r + c * i;
    }
    return;
  }

  const std::uint32_t count = gate.pauli_qubit_count;
  if (count == 0) {
    return;
  }

  for (std::uint32_t idx = 0; idx < count; ++idx) {
    const auto basis = gate.pauli_basis[idx];
    const auto qubit = static_cast<std::size_t>(gate.pauli_qubits[idx]);
    switch (basis) {
      case 1:
        apply_h_internal(qubit);
        break;
      case 2:
        apply_y_basis_forward(qubit);
        break;
      default:
        break;
    }
  }

  for (std::uint32_t idx = count; idx > 1; --idx) {
    apply_cnot_internal(gate.pauli_qubits[idx - 1], gate.pauli_qubits[idx - 2]);
  }

  apply_rz_internal(gate.pauli_qubits[0], gate.angle);

  for (std::uint32_t idx = 1; idx < count; ++idx) {
    apply_cnot_internal(gate.pauli_qubits[idx], gate.pauli_qubits[idx - 1]);
  }

  for (std::uint32_t idx = 0; idx < count; ++idx) {
    const auto basis = gate.pauli_basis[idx];
    const auto qubit = static_cast<std::size_t>(gate.pauli_qubits[idx]);
    switch (basis) {
      case 1:
        apply_h_internal(qubit);
        break;
      case 2:
        apply_y_basis_backward(qubit);
        break;
      default:
        break;
    }
  }
}

std::complex<double> statevector_cpu::expectation(const pauli_term& term) const {
  const std::uint64_t x_mask = term.x_mask;
  const std::uint64_t z_mask = term.z_mask;
  const std::uint64_t y_mask = x_mask & z_mask;
  const std::size_t y_count = static_cast<std::size_t>(__builtin_popcountll(y_mask));
  const double sign = (y_count / 2) % 2 ? -1.0 : 1.0;
  const bool odd_y = (y_count & 1) != 0;

  double expectation_value = 0.0;

  if (x_mask == 0) {
    #ifdef NWQSIM_USE_OPENMP
    #pragma omp parallel for reduction(+:expectation_value) schedule(static) if (dimension_ >= kParallelThreshold)
    #endif
    for (std::size_t idx = 0; idx < dimension_; ++idx) {
      double contrib = state_real_[idx] * state_real_[idx] + state_imag_[idx] * state_imag_[idx];
      if (__builtin_popcountll(idx & z_mask) & 1) {
        contrib = -contrib;
      }
      expectation_value += contrib;
    }
    return term.coefficient * expectation_value;
  }

  const unsigned max_x = 63U - static_cast<unsigned>(__builtin_clzll(x_mask));
  const std::size_t mask_lower = (static_cast<std::size_t>(1) << max_x) - 1;
  const std::size_t mask_upper = ~((static_cast<std::size_t>(1) << (max_x + 1)) - 1);

  #ifdef NWQSIM_USE_OPENMP
  #pragma omp parallel for reduction(+:expectation_value) schedule(static) if (dimension_ >= kParallelThreshold)
  #endif
  for (std::size_t i = 0; i < half_dimension_; ++i) {
    const std::size_t idx0 = ((i << 1) & mask_upper) | (i & mask_lower);
    const std::size_t idx1 = idx0 ^ x_mask;

    const double a0_r = state_real_[idx0];
    const double a0_i = state_imag_[idx0];
    const double a1_r = state_real_[idx1];
    const double a1_i = state_imag_[idx1];

    double v0;
    double v1;

    if (odd_y) {
      v0 = -a1_i * a0_r + a0_i * a1_r;
      v1 = -a0_i * a1_r + a1_i * a0_r;
    } else {
      v0 = a1_r * a0_r + a1_i * a0_i;
      v1 = a1_r * a0_r + a0_i * a1_i;
    }

    v0 *= sign;
    v1 *= sign;

    double contrib = v0;
    if (__builtin_popcountll(idx0 & z_mask) & 1) {
      contrib = -contrib;
    }

    if (__builtin_popcountll(idx1 & z_mask) & 1) {
      contrib -= v1;
    } else {
      contrib += v1;
    }

    expectation_value += contrib;
  }

  return term.coefficient * expectation_value;
}

std::complex<double> statevector_cpu::expectation(const std::vector<pauli_term>& terms) const {
  std::complex<double> total{0.0, 0.0};
  for (const auto& term : terms) {
    total += expectation(term);
  }
  return total;
}

const std::vector<std::complex<double>>& statevector_cpu::amplitudes() const {
  if (!amplitudes_valid_) {
    amplitude_cache_.resize(dimension_);
    #ifdef NWQSIM_USE_OPENMP
    #pragma omp parallel for schedule(static) if (dimension_ >= kParallelThreshold)
    #endif
    for (std::size_t idx = 0; idx < dimension_; ++idx) {
      amplitude_cache_[idx] = {state_real_[idx], state_imag_[idx]};
    }
    amplitudes_valid_ = true;
  }
  return amplitude_cache_;
}

void statevector_cpu::load_state(const std::vector<std::complex<double>>& amplitudes) {
  if (amplitudes.size() != dimension_) {
    throw std::invalid_argument("amplitude vector size mismatch");
  }
  #ifdef NWQSIM_USE_OPENMP
  #pragma omp parallel for schedule(static) if (dimension_ >= kParallelThreshold)
  #endif
  for (std::size_t idx = 0; idx < dimension_; ++idx) {
    state_real_[idx] = amplitudes[idx].real();
    state_imag_[idx] = amplitudes[idx].imag();
  }
  amplitudes_valid_ = false;
}

void statevector_cpu::copy_state_from(const statevector_cpu& other) {
  if (other.dimension_ != dimension_) {
    throw std::invalid_argument("state dimension mismatch in copy_state_from");
  }
  std::copy(other.state_real_.begin(), other.state_real_.end(), state_real_.begin());
  std::copy(other.state_imag_.begin(), other.state_imag_.end(), state_imag_.begin());
  amplitudes_valid_ = false;
}

}  // namespace vqe::backend
