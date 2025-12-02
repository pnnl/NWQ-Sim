#include "backend/gate_fusion.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace vqe::backend
{
  namespace
  {

    constexpr std::array<double, 4> kIdentityReal{{1.0, 0.0, 0.0, 1.0}};
    constexpr std::array<double, 4> kIdentityImag{{0.0, 0.0, 0.0, 0.0}};

    inline std::array<double, 4> extract_single_real(const sim_gate &gate)
    {
      std::array<double, 4> result{};
      std::copy_n(gate.real.begin(), 4, result.begin());
      return result;
    }

    inline std::array<double, 4> extract_single_imag(const sim_gate &gate)
    {
      std::array<double, 4> result{};
      std::copy_n(gate.imag.begin(), 4, result.begin());
      return result;
    }

    inline void kron(const std::array<double, 4> &a_real,
                     const std::array<double, 4> &a_imag,
                     const std::array<double, 4> &b_real,
                     const std::array<double, 4> &b_imag,
                     std::array<double, 16> &out_real,
                     std::array<double, 16> &out_imag)
    {
      for (std::size_t r = 0; r < 2; ++r)
      {
        for (std::size_t s = 0; s < 2; ++s)
        {
          const double ar = a_real[r * 2 + s];
          const double ai = a_imag[r * 2 + s];
          for (std::size_t v = 0; v < 2; ++v)
          {
            for (std::size_t w = 0; w < 2; ++w)
            {
              const std::size_t row = r * 2 + v;
              const std::size_t col = s * 2 + w;
              const std::size_t idx = row * 4 + col;
              const double br = b_real[v * 2 + w];
              const double bi = b_imag[v * 2 + w];
              out_real[idx] = ar * br - ai * bi;
              out_imag[idx] = ar * bi + ai * br;
            }
          }
        }
      }
    }

    inline void multiply_assign_single(const sim_gate &left, sim_gate &right)
    {
      double real[4] = {0.0};
      double imag[4] = {0.0};
      for (std::size_t m = 0; m < 2; ++m)
      {
        for (std::size_t n = 0; n < 2; ++n)
        {
          double r = 0.0;
          double i = 0.0;
          for (std::size_t k = 0; k < 2; ++k)
          {
            const double lr = left.real[m * 2 + k];
            const double li = left.imag[m * 2 + k];
            const double rr = right.real[k * 2 + n];
            const double ri = right.imag[k * 2 + n];
            r += lr * rr - li * ri;
            i += lr * ri + li * rr;
          }
          real[m * 2 + n] = r;
          imag[m * 2 + n] = i;
        }
      }
      std::copy(real, real + 4, right.real.begin());
      std::copy(imag, imag + 4, right.imag.begin());
    }

    inline void multiply_assign_two_left(const sim_gate &left, sim_gate &right)
    {
      double real[16] = {0.0};
      double imag[16] = {0.0};
      for (std::size_t m = 0; m < 4; ++m)
      {
        for (std::size_t n = 0; n < 4; ++n)
        {
          double r = 0.0;
          double i = 0.0;
          for (std::size_t k = 0; k < 4; ++k)
          {
            const double lr = left.real[m * 4 + k];
            const double li = left.imag[m * 4 + k];
            const double rr = right.real[k * 4 + n];
            const double ri = right.imag[k * 4 + n];
            r += lr * rr - li * ri;
            i += lr * ri + li * rr;
          }
          real[m * 4 + n] = r;
          imag[m * 4 + n] = i;
        }
      }
      std::copy(real, real + 16, right.real.begin());
      std::copy(imag, imag + 16, right.imag.begin());
    }

    inline void multiply_assign_two_right(sim_gate &left, const sim_gate &right)
    {
      double real[16] = {0.0};
      double imag[16] = {0.0};
      for (std::size_t m = 0; m < 4; ++m)
      {
        for (std::size_t n = 0; n < 4; ++n)
        {
          double r = 0.0;
          double i = 0.0;
          for (std::size_t k = 0; k < 4; ++k)
          {
            const double lr = left.real[m * 4 + k];
            const double li = left.imag[m * 4 + k];
            const double rr = right.real[k * 4 + n];
            const double ri = right.imag[k * 4 + n];
            r += lr * rr - li * ri;
            i += lr * ri + li * rr;
          }
          real[m * 4 + n] = r;
          imag[m * 4 + n] = i;
        }
      }
      std::copy(real, real + 16, left.real.begin());
      std::copy(imag, imag + 16, left.imag.begin());
    }

    inline void multiply_two_into(const sim_gate &left, const sim_gate &right, sim_gate &out)
    {
      out.op = sim_gate::kind::two;
      double real[16] = {0.0};
      double imag[16] = {0.0};
      for (std::size_t m = 0; m < 4; ++m)
      {
        for (std::size_t n = 0; n < 4; ++n)
        {
          double r = 0.0;
          double i = 0.0;
          for (std::size_t k = 0; k < 4; ++k)
          {
            const double lr = left.real[m * 4 + k];
            const double li = left.imag[m * 4 + k];
            const double rr = right.real[k * 4 + n];
            const double ri = right.imag[k * 4 + n];
            r += lr * rr - li * ri;
            i += lr * ri + li * rr;
          }
          real[m * 4 + n] = r;
          imag[m * 4 + n] = i;
        }
      }
      std::copy(real, real + 16, out.real.begin());
      std::copy(imag, imag + 16, out.imag.begin());
    }

    inline sim_gate make_swap_gate()
    {
      constexpr double real[16] = {
          1.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 1.0, 0.0,
          0.0, 1.0, 0.0, 0.0,
          0.0, 0.0, 0.0, 1.0};
      constexpr double imag[16] = {0.0};
      sim_gate gate{};
      gate.op = sim_gate::kind::two;
      gate.control = 0;
      gate.target = 1;
      std::copy(real, real + 16, gate.real.begin());
      std::copy(imag, imag + 16, gate.imag.begin());
      return gate;
    }

    inline sim_gate swap_control_target(const sim_gate &gate)
    {
      const sim_gate swap_gate = make_swap_gate();
      sim_gate tmp{};
      multiply_two_into(gate, swap_gate, tmp);
      sim_gate result{};
      multiply_two_into(swap_gate, tmp, result);
      result.control = gate.target;
      result.target = gate.control;
      return result;
    }

    inline sim_gate expand_to_target(const sim_gate &single, const sim_gate &two_gate)
    {
      sim_gate expanded{};
      expanded.op = sim_gate::kind::two;
      expanded.control = two_gate.control;
      expanded.target = two_gate.target;
      const auto sr = extract_single_real(single);
      const auto si = extract_single_imag(single);
      kron(kIdentityReal, kIdentityImag, sr, si, expanded.real, expanded.imag);
      return expanded;
    }

    inline sim_gate expand_to_control(const sim_gate &single, const sim_gate &two_gate)
    {
      sim_gate expanded{};
      expanded.op = sim_gate::kind::two;
      expanded.control = two_gate.control;
      expanded.target = two_gate.target;
      const auto sr = extract_single_real(single);
      const auto si = extract_single_imag(single);
      kron(sr, si, kIdentityReal, kIdentityImag, expanded.real, expanded.imag);
      return expanded;
    }

    inline void gate_fusion_1q(const std::vector<sim_gate> &input,
                               std::vector<sim_gate> &output,
                               std::size_t n_qubits)
    {
      std::vector<std::ptrdiff_t> table(n_qubits, -1);
      std::vector<bool> canfuse(n_qubits, false);
      output.clear();

      for (const auto &gate : input)
      {
        if (gate.is_single())
        {
          const std::size_t qubit = gate.target;
          if (!canfuse[qubit])
          {
            output.push_back(gate);
            canfuse[qubit] = true;
            table[qubit] = static_cast<std::ptrdiff_t>(output.size() - 1);
          }
          else
          {
            multiply_assign_single(gate, output[table[qubit]]);
          }
        }
        else
        {
          if (gate.control == gate.target)
          {
            throw std::logic_error("invalid two-qubit gate with identical qubits");
          }
          canfuse[gate.target] = false;
          canfuse[gate.control] = false;
          output.push_back(gate);
        }
      }
    }

    inline void gate_fusion_2q_absorb_1q_forward(const std::vector<sim_gate> &input,
                                                 std::vector<sim_gate> &output,
                                                 std::size_t n_qubits)
    {
      const std::size_t total = n_qubits * n_qubits;
      std::vector<std::ptrdiff_t> table(total, -1);
      std::vector<bool> canfuse(total, false);
      output.clear();

      for (const auto &gate : input)
      {
        if (gate.is_single())
        {
          const std::size_t qubit = gate.target;
          bool fused = false;
          for (std::size_t j = 0; j < n_qubits; ++j)
          {
            const std::size_t tgt_key = j * n_qubits + qubit;
            if (canfuse[tgt_key])
            {
              auto &two_gate = output[table[tgt_key]];
              const sim_gate expanded = expand_to_target(gate, two_gate);
              multiply_assign_two_left(expanded, two_gate);
              fused = true;
              break;
            }
            const std::size_t ctrl_key = qubit * n_qubits + j;
            if (canfuse[ctrl_key])
            {
              auto &two_gate = output[table[ctrl_key]];
              const sim_gate expanded = expand_to_control(gate, two_gate);
              multiply_assign_two_left(expanded, two_gate);
              fused = true;
              break;
            }
          }
          if (!fused)
          {
            for (std::size_t j = 0; j < n_qubits; ++j)
            {
              canfuse[j * n_qubits + qubit] = false;
              canfuse[qubit * n_qubits + j] = false;
            }
            output.push_back(gate);
          }
          continue;
        }

        sim_gate two_gate = gate;
        if (two_gate.control > two_gate.target)
        {
          two_gate = swap_control_target(two_gate);
        }

        for (std::size_t j = 0; j < n_qubits; ++j)
        {
          canfuse[two_gate.target * n_qubits + j] = false;
          canfuse[two_gate.control * n_qubits + j] = false;
          canfuse[j * n_qubits + two_gate.target] = false;
          canfuse[j * n_qubits + two_gate.control] = false;
        }

        output.push_back(two_gate);
        const std::size_t key = two_gate.control * n_qubits + two_gate.target;
        canfuse[key] = true;
        table[key] = static_cast<std::ptrdiff_t>(output.size() - 1);
      }
    }

    inline void gate_fusion_2q_absorb_1q_backward(const std::vector<sim_gate> &input,
                                                  std::vector<sim_gate> &output,
                                                  std::size_t n_qubits)
    {
      const std::size_t total = n_qubits * n_qubits;
      std::vector<std::ptrdiff_t> table(total, -1);
      std::vector<bool> canfuse(total, false);
      output.clear();

      for (std::size_t idx = input.size(); idx-- > 0;)
      {
        const auto &gate = input[idx];
        if (gate.is_single())
        {
          const std::size_t qubit = gate.target;
          bool fused = false;
          for (std::size_t j = 0; j < n_qubits; ++j)
          {
            const std::size_t tgt_key = j * n_qubits + qubit;
            if (canfuse[tgt_key])
            {
              auto &two_gate = output[table[tgt_key]];
              const sim_gate expanded = expand_to_target(gate, two_gate);
              multiply_assign_two_right(two_gate, expanded);
              fused = true;
              break;
            }
            const std::size_t ctrl_key = qubit * n_qubits + j;
            if (canfuse[ctrl_key])
            {
              auto &two_gate = output[table[ctrl_key]];
              const sim_gate expanded = expand_to_control(gate, two_gate);
              multiply_assign_two_right(two_gate, expanded);
              fused = true;
              break;
            }
          }
          if (!fused)
          {
            for (std::size_t j = 0; j < n_qubits; ++j)
            {
              canfuse[j * n_qubits + qubit] = false;
              canfuse[qubit * n_qubits + j] = false;
            }
            output.push_back(gate);
          }
          continue;
        }

        sim_gate two_gate = gate;
        if (two_gate.control > two_gate.target)
        {
          two_gate = swap_control_target(two_gate);
        }

        for (std::size_t j = 0; j < n_qubits; ++j)
        {
          canfuse[two_gate.target * n_qubits + j] = false;
          canfuse[two_gate.control * n_qubits + j] = false;
          canfuse[j * n_qubits + two_gate.target] = false;
          canfuse[j * n_qubits + two_gate.control] = false;
        }

        output.push_back(two_gate);
        const std::size_t key = two_gate.control * n_qubits + two_gate.target;
        canfuse[key] = true;
        table[key] = static_cast<std::ptrdiff_t>(output.size() - 1);
      }

      std::reverse(output.begin(), output.end());
    }

    inline void gate_fusion_2q(const std::vector<sim_gate> &input,
                               std::vector<sim_gate> &output,
                               std::size_t n_qubits)
    {
      const std::size_t total = n_qubits * n_qubits;
      std::vector<std::ptrdiff_t> table(total, -1);
      std::vector<bool> canfuse(total, false);
      output.clear();

      for (const auto &gate : input)
      {
        if (gate.is_single())
        {
          output.push_back(gate);
          continue;
        }

        sim_gate two_gate = gate;
        if (two_gate.control > two_gate.target)
        {
          two_gate = swap_control_target(two_gate);
        }

        const std::size_t ctrl = two_gate.control;
        const std::size_t tgt = two_gate.target;
        const std::size_t key = ctrl * n_qubits + tgt;

        if (!canfuse[key])
        {
          for (std::size_t j = 0; j < n_qubits; ++j)
          {
            canfuse[tgt * n_qubits + j] = false;
            canfuse[ctrl * n_qubits + j] = false;
            canfuse[j * n_qubits + tgt] = false;
            canfuse[j * n_qubits + ctrl] = false;
          }
          output.push_back(two_gate);
          canfuse[key] = true;
          table[key] = static_cast<std::ptrdiff_t>(output.size() - 1);
        }
        else
        {
          auto &existing = output[table[key]];
          multiply_assign_two_left(two_gate, existing);
        }
      }
    }

    void fuse_sequence(const std::vector<sim_gate> &gates,
                       std::size_t n_qubits,
                       std::vector<sim_gate> &tmp1,
                       std::vector<sim_gate> &tmp2,
                       std::vector<sim_gate> &tmp3,
                       std::vector<sim_gate> &out)
    {
      gate_fusion_1q(gates, tmp1, n_qubits);
      gate_fusion_2q_absorb_1q_forward(tmp1, tmp2, n_qubits);
      gate_fusion_2q_absorb_1q_backward(tmp2, tmp3, n_qubits);
      gate_fusion_2q(tmp3, out, n_qubits);
    }

    void fuse_simulation_gates_impl(const std::vector<sim_gate> &gates,
                                    std::size_t n_qubits,
                                    std::vector<sim_gate> &buffer,
                                    std::vector<sim_gate> &tmp1,
                                    std::vector<sim_gate> &tmp2,
                                    std::vector<sim_gate> &tmp3,
                                    std::vector<sim_gate> &chunk,
                                    std::vector<sim_gate> &fused)
    {
      buffer.clear();
      fused.clear();

      auto flush = [&]()
      {
        if (buffer.empty())
        {
          return;
        }
        fuse_sequence(buffer, n_qubits, tmp1, tmp2, tmp3, chunk);
        fused.insert(fused.end(), chunk.begin(), chunk.end());
        buffer.clear();
      };

      for (const auto &gate : gates)
      {
        if (gate.op == sim_gate::kind::pauli)
        {
          flush();
          fused.push_back(gate);
        }
        else
        {
          buffer.push_back(gate);
        }
      }

      flush();
    }

  } // namespace

  void fuse_simulation_gates(const std::vector<sim_gate> &gates,
                             std::size_t n_qubits,
                             std::vector<sim_gate> &buffer,
                             std::vector<sim_gate> &tmp1,
                             std::vector<sim_gate> &tmp2,
                             std::vector<sim_gate> &tmp3,
                             std::vector<sim_gate> &chunk,
                             std::vector<sim_gate> &fused)
  {
    fuse_simulation_gates_impl(gates, n_qubits, buffer, tmp1, tmp2, tmp3, chunk, fused);
  }

  std::vector<sim_gate> fuse_simulation_gates(const std::vector<sim_gate> &gates,
                                              std::size_t n_qubits)
  {
    std::vector<sim_gate> buffer;
    std::vector<sim_gate> tmp1;
    std::vector<sim_gate> tmp2;
    std::vector<sim_gate> tmp3;
    std::vector<sim_gate> chunk;
    std::vector<sim_gate> fused;
    fuse_simulation_gates_impl(gates, n_qubits, buffer, tmp1, tmp2, tmp3, chunk, fused);
    return fused;
  }

} // namespace vqe::backend
