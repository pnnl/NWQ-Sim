#include "backend/gate_builder.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace vqe::backend {
namespace {

inline std::uint32_t checked_index(std::size_t value) {
  if (value >= static_cast<std::size_t>(sim_gate::invalid_index)) {
    throw std::out_of_range("qubit index exceeds sim_gate capacity");
  }
  return static_cast<std::uint32_t>(value);
}

inline double evaluate_expression(const circuit& circ, const parameter_expression& expr) {
  if (expr.empty()) {
    return 0.0;
  }
  const auto& parameters = circ.parameters();
  double value = 0.0;
  for (const auto& term : expr.terms) {
    if (term.parameter_index >= parameters.size()) {
      throw std::out_of_range("parameter index out of range while evaluating expression");
    }
    value += parameters[term.parameter_index] * term.coefficient;
  }
  return expr.scale * value;
}

inline sim_gate make_single(std::uint32_t qubit, const double (&real)[4], const double (&imag)[4]) {
  sim_gate gate{};
  gate.op = sim_gate::kind::single;
  gate.target = qubit;
  std::copy(real, real + 4, gate.real.begin());
  std::copy(imag, imag + 4, gate.imag.begin());
  return gate;
}

inline sim_gate make_two(std::uint32_t control, std::uint32_t target,
                         const double (&real)[16], const double (&imag)[16]) {
  sim_gate gate{};
  gate.op = sim_gate::kind::two;
  gate.control = control;
  gate.target = target;
  std::copy(real, real + 16, gate.real.begin());
  std::copy(imag, imag + 16, gate.imag.begin());
  return gate;
}

inline sim_gate make_x(std::uint32_t qubit) {
  constexpr double real[4] = {0.0, 1.0, 1.0, 0.0};
  constexpr double imag[4] = {0.0, 0.0, 0.0, 0.0};
  return make_single(qubit, real, imag);
}

inline sim_gate make_h(std::uint32_t qubit) {
  const double real[4] = {0.707106781186547524400844362104849039, 0.707106781186547524400844362104849039,
                          0.707106781186547524400844362104849039, -0.707106781186547524400844362104849039};
  constexpr double imag[4] = {0.0, 0.0, 0.0, 0.0};
  return make_single(qubit, real, imag);
}

inline sim_gate make_rz(std::uint32_t qubit, double angle) {
  const double half = 0.5 * angle;
  const double c = std::cos(half);
  const double s = std::sin(half);
  const double real[4] = {c, 0.0, 0.0, c};
  const double imag[4] = {-s, 0.0, 0.0, s};
  return make_single(qubit, real, imag);
}

inline sim_gate make_cnot(std::uint32_t control, std::uint32_t target) {
  constexpr double real[16] = {
      1.0, 0.0, 0.0, 0.0,
      0.0, 1.0, 0.0, 0.0,
      0.0, 0.0, 0.0, 1.0,
      0.0, 0.0, 1.0, 0.0};
  constexpr double imag[16] = {0.0};
  return make_two(control, target, real, imag);
}

}  // namespace

void build_simulation_gates(const circuit& circ, std::vector<sim_gate>& out) {
  out.clear();
  out.reserve(circ.gates().size() * 4);

  for (const auto& g : circ.gates()) {
    switch (g.kind) {
      case gate_kind::x:
        if (g.qubits.size() != 1) {
          throw std::logic_error("X gate requires one qubit");
        }
        out.push_back(make_x(checked_index(g.qubits.front())));
        break;
      case gate_kind::h:
        if (g.qubits.size() != 1) {
          throw std::logic_error("H gate requires one qubit");
        }
        out.push_back(make_h(checked_index(g.qubits.front())));
        break;
      case gate_kind::rz: {
        if (g.qubits.size() != 1) {
          throw std::logic_error("RZ gate requires one qubit");
        }
        const double angle = g.parametrized ? evaluate_expression(circ, g.expression) : g.angle;
        out.push_back(make_rz(checked_index(g.qubits.front()), angle));
        break;
      }
      case gate_kind::cnot:
        if (g.qubits.size() != 2) {
          throw std::logic_error("CNOT gate requires two qubits");
        }
        out.push_back(make_cnot(checked_index(g.qubits[0]), checked_index(g.qubits[1])));
        break;
      case gate_kind::pauli_rotation:
        if (!g.qubits.empty()) {
          sim_gate gate{};
          gate.op = sim_gate::kind::pauli;
          gate.pauli = g.pauli;
          gate.angle = evaluate_expression(circ, g.expression);
          gate.pauli_qubit_count = static_cast<std::uint32_t>(g.qubits.size());
          for (std::size_t i = 0; i < g.qubits.size() && i < gate.pauli_qubits.size(); ++i) {
            const std::size_t qubit = g.qubits[i];
            gate.pauli_qubits[i] = static_cast<std::uint32_t>(qubit);
            const bool x = (g.pauli.x_mask >> qubit) & 1ULL;
            const bool z = (g.pauli.z_mask >> qubit) & 1ULL;
            std::uint8_t basis = 0;
            if (x && z) {
              basis = 2;  // Y
            } else if (x) {
              basis = 1;  // X
            } else {
              basis = 3;  // Z
            }
            gate.pauli_basis[i] = basis;
          }
          out.push_back(gate);
        }
        break;
      default:
        throw std::runtime_error("unsupported gate kind in build_simulation_gates");
    }
  }

}

std::vector<sim_gate> build_simulation_gates(const circuit& circ) {
  std::vector<sim_gate> gates;
  build_simulation_gates(circ, gates);
  
  return gates;
}

}  // namespace vqe::backend
