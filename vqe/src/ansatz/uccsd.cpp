#include "ansatz/ansatz.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

#include <stdexcept>
#include "fermion_term.hpp"
#include "jw_transform.hpp"

namespace vqe
{

  namespace
  {

    constexpr double kCoefficientCutoff = 1e-10;
    const std::size_t kInvalidParam = std::numeric_limits<std::size_t>::max();

    std::size_t choose2(std::size_t value)
    {
      if (value < 2)
      {
        return 0;
      }
      return value * (value - 1) / 2;
    }

    fermion_op_kind to_fermion_kind(operator_kind op)
    {
      return op == operator_kind::creation ? fermion_op_kind::creation : fermion_op_kind::annihilation;
    }

    fermion_term build_fermion_term(const std::vector<fermion_operator> &ops,
                                    const molecular_environment &env)
    {
      fermion_term term;
      term.coefficient = 1.0;
      for (const auto &op : ops)
      {
        const std::size_t qubit = env.qubit_index(op.orbital_index, op.orbital, op.spin_state);
        term.operators.push_back({qubit, to_fermion_kind(op.op)});
        term.operators.back().index = qubit;
      }
      return term;
    }

    std::vector<pauli_term> hermitian_pauli_terms(const std::vector<fermion_operator> &ops,
                                                  const molecular_environment &env)
    {
      hamiltonian_data data;
      data.constant = 0.0;
      data.max_index = env.total_qubits() == 0 ? 0 : env.total_qubits() - 1;
      data.terms.push_back(build_fermion_term(ops, env));

      auto raw_terms = jordan_wigner_transform(data);
      std::vector<pauli_term> result;
      result.reserve(raw_terms.size());
      for (auto term : raw_terms)
      {
        const double value = 2.0 * term.coefficient.imag();
        if (std::abs(value) < kCoefficientCutoff)
        {
          continue;
        }
        term.coefficient = value;
        result.push_back(term);
      }
      return result;
    }

    std::string op_to_string(const fermion_operator &op, const molecular_environment &env)
    {
      std::ostringstream stream;
      stream << env.qubit_index(op.orbital_index, op.orbital, op.spin_state);
      if (op.op == operator_kind::creation)
      {
        stream << '^';
      }
      return stream.str();
    }

    std::string excitation_label(const std::vector<fermion_operator> &ops, const molecular_environment &env)
    {
      std::string label;
      bool first = true;
      for (const auto &op : ops)
      {
        const std::string token = op_to_string(op, env);
        if (!first)
        {
          label = " " + label;
        }
        else
        {
          first = false;
        }
        label = token + label;
      }
      return label;
    }

    parameter_expression build_expression(const excitation_info &excitation,
                                          const std::vector<std::size_t> &mapping,
                                          std::size_t parameter_offset)
    {
      parameter_expression expr;
      expr.terms.reserve(excitation.symmetry_terms.size());
      for (const auto &entry : excitation.symmetry_terms)
      {
        const std::size_t excitation_index = entry.first;
        if (excitation_index >= mapping.size())
        {
          throw std::out_of_range("symmetry term references invalid excitation index");
        }
        const std::size_t param_index = mapping[excitation_index];
        if (param_index == kInvalidParam)
        {
          throw std::runtime_error("symmetry references excitation without parameter");
        }
        expr.terms.push_back(parameter_term{param_index + parameter_offset, entry.second});
      }
      return expr;
    }

    bool is_non_trivial(const pauli_term &term)
    {
      return term.x_mask != 0 || term.z_mask != 0;
    }

  } // namespace

  uccsd_ansatz::uccsd_ansatz(const molecular_environment &env,
                             std::size_t trotter_steps,
                             std::size_t symmetry_level)
      : ansatz_base(env.total_qubits()),
        env_(env),
        trotter_steps_(trotter_steps),
        symmetry_level_(symmetry_level) {}

  void uccsd_ansatz::build()
  {
    excitations_.clear();
    excitation_to_parameter_.clear();
    excitation_to_parameter_label_.clear();
    unique_parameters_ = 0;

    const std::size_t n_occ = env_.occupied_orbitals();
    const std::size_t n_virt = env_.virtual_orbitals();
    const std::size_t n_singles = 2 * n_occ * n_virt;
    const std::size_t n_doubles = 2 * choose2(n_virt) * choose2(n_occ) + n_occ * n_virt * n_occ * n_virt;

    excitations_.reserve(n_singles + n_doubles);
    excitation_to_parameter_.reserve(n_singles + n_doubles);

    generate_single_excitations();
    generate_same_spin_doubles();
    generate_mixed_spin_doubles();

    pauli_generators_.clear();
    pauli_generators_.reserve(excitations_.size());
    for (const auto &excitation : excitations_)
    {
      pauli_generators_.push_back(hermitian_pauli_terms(excitation.operators, env_));
    }

    circuit_ = circuit(env_.total_qubits());
    const std::size_t total_parameters = unique_parameters_ * std::max<std::size_t>(trotter_steps_, 1);
    for (std::size_t idx = 0; idx < total_parameters; ++idx)
    {
      static_cast<void>(circuit_.add_parameter());
    }

    if (env_.total_qubits() > 0)
    {
      prepare_reference_state();
    }
    append_pauli_generators();
  }

  void uccsd_ansatz::generate_single_excitations()
  {
    const std::size_t n_occ = env_.occupied_orbitals();
    const std::size_t n_virt = env_.virtual_orbitals();
    for (std::size_t p = 0; p < n_occ; ++p)
    {
      const fermion_operator occ_up(p, orbital_kind::occupied, spin::up, operator_kind::annihilation);
      const fermion_operator occ_down(p, orbital_kind::occupied, spin::down, operator_kind::annihilation);
      for (std::size_t q = 0; q < n_virt; ++q)
      {
        const fermion_operator virt_up(q, orbital_kind::virtual_orbital, spin::up, operator_kind::creation);
        const fermion_operator virt_down(q, orbital_kind::virtual_orbital, spin::down, operator_kind::creation);
        if (symmetry_level_ >= 1)
        {
          const std::size_t term_index = excitations_.size();
          add_single_excitation(occ_up, virt_up, {{term_index, 1.0}}, true);
          add_single_excitation(occ_down, virt_down, {{term_index, 1.0}}, false);
        }
        else
        {
          add_single_excitation(occ_up, virt_up);
          add_single_excitation(occ_down, virt_down);
        }
      }
    }
  }

  void uccsd_ansatz::generate_same_spin_doubles()
  {
    const std::size_t n_occ = env_.occupied_orbitals();
    const std::size_t n_virt = env_.virtual_orbitals();
    for (std::size_t i = 0; i < n_occ; ++i)
    {
      const fermion_operator i_up(i, orbital_kind::occupied, spin::up, operator_kind::annihilation);
      const fermion_operator i_down(i, orbital_kind::occupied, spin::down, operator_kind::annihilation);
      for (std::size_t r = 0; r < n_virt; ++r)
      {
        const fermion_operator r_up(r, orbital_kind::virtual_orbital, spin::up, operator_kind::creation);
        const fermion_operator r_down(r, orbital_kind::virtual_orbital, spin::down, operator_kind::creation);
        for (std::size_t j = i + 1; j < n_occ; ++j)
        {
          const fermion_operator j_up(j, orbital_kind::occupied, spin::up, operator_kind::annihilation);
          const fermion_operator j_down(j, orbital_kind::occupied, spin::down, operator_kind::annihilation);
          for (std::size_t s = r + 1; s < n_virt; ++s)
          {
            const fermion_operator s_up(s, orbital_kind::virtual_orbital, spin::up, operator_kind::creation);
            const fermion_operator s_down(s, orbital_kind::virtual_orbital, spin::down, operator_kind::creation);
            if (symmetry_level_ >= 1)
            {
              const std::size_t term_index = excitations_.size();
              add_double_excitation(i_up, j_up, r_up, s_up, {{term_index, 1.0}}, true);
              add_double_excitation(i_down, j_down, r_down, s_down, {{term_index, 1.0}}, false);
            }
            else
            {
              add_double_excitation(i_up, j_up, r_up, s_up);
              add_double_excitation(i_down, j_down, r_down, s_down);
            }
          }
        }
      }
    }
  }

  void uccsd_ansatz::generate_mixed_spin_doubles()
  {
    const std::size_t n_occ = env_.occupied_orbitals();
    const std::size_t n_virt = env_.virtual_orbitals();

    for (std::size_t i = 0; i < n_occ; ++i)
    {
      const fermion_operator i_up(i, orbital_kind::occupied, spin::up, operator_kind::annihilation);
      for (std::size_t r = 0; r < n_virt; ++r)
      {
        const fermion_operator r_up(r, orbital_kind::virtual_orbital, spin::up, operator_kind::creation);
        for (std::size_t j = 0; j < n_occ; ++j)
        {
          const fermion_operator j_down(j, orbital_kind::occupied, spin::down, operator_kind::annihilation);
          if (!(symmetry_level_ < 2 || i == j || i < j))
          {
            continue;
          }
          for (std::size_t s = 0; s < n_virt; ++s)
          {
            if (!(symmetry_level_ < 2 || (i == j && r == s) || (i == j && r < s) || (i < j)))
            {
              continue;
            }
            const fermion_operator s_down(s, orbital_kind::virtual_orbital, spin::down, operator_kind::creation);
            if (symmetry_level_ < 2 || (i == j && r == s))
            {
              add_double_excitation(i_up, j_down, r_up, s_down);
            }
            else
            {
              const std::size_t term_index = excitations_.size();
              add_double_excitation(i_up, j_down, r_up, s_down, {{term_index, 1.0}}, true);

              const fermion_operator j_up(j, orbital_kind::occupied, spin::up, operator_kind::annihilation);
              const fermion_operator i_down(i, orbital_kind::occupied, spin::down, operator_kind::annihilation);
              const fermion_operator s_up(s, orbital_kind::virtual_orbital, spin::up, operator_kind::creation);
              const fermion_operator r_down(r, orbital_kind::virtual_orbital, spin::down, operator_kind::creation);
              add_double_excitation(j_up, i_down, s_up, r_down, {{term_index, 1.0}}, false);
            }
          }
        }
      }
    }
  }

  void uccsd_ansatz::add_single_excitation(const fermion_operator &from,
                                           const fermion_operator &to,
                                           const std::vector<std::pair<std::size_t, double>> &symmetry,
                                           bool owns_parameter)
  {
    excitation_info info;
    info.operators = {from, to};
    info.symmetry_terms = symmetry;
    if (info.symmetry_terms.empty())
    {
      info.symmetry_terms.emplace_back(excitations_.size(), 1.0);
    }
    info.label = excitation_label(info.operators, env_);
    excitations_.push_back(std::move(info));
    if (owns_parameter)
    {
      excitation_to_parameter_.push_back(unique_parameters_);
      excitation_to_parameter_label_[excitations_.back().label] = unique_parameters_;
      ++unique_parameters_;
    }
    else
    {
      excitation_to_parameter_.push_back(kInvalidParam);
    }
  }

  void uccsd_ansatz::add_double_excitation(const fermion_operator &i,
                                           const fermion_operator &j,
                                           const fermion_operator &r,
                                           const fermion_operator &s,
                                           const std::vector<std::pair<std::size_t, double>> &symmetry,
                                           bool owns_parameter)
  {
    excitation_info info;
    info.operators = {i, j, r, s};
    info.symmetry_terms = symmetry;
    if (info.symmetry_terms.empty())
    {
      info.symmetry_terms.emplace_back(excitations_.size(), 1.0);
    }
    info.label = excitation_label(info.operators, env_);
    excitations_.push_back(std::move(info));
    if (owns_parameter)
    {
      excitation_to_parameter_.push_back(unique_parameters_);
      excitation_to_parameter_label_[excitations_.back().label] = unique_parameters_;
      ++unique_parameters_;
    }
    else
    {
      excitation_to_parameter_.push_back(kInvalidParam);
    }
  }

  void uccsd_ansatz::add_single_excitation(const fermion_operator &from,
                                           const fermion_operator &to)
  {
    const std::size_t term_index = excitations_.size();
    add_single_excitation(from, to, {{term_index, 1.0}}, true);
  }

  void uccsd_ansatz::add_double_excitation(const fermion_operator &i,
                                           const fermion_operator &j,
                                           const fermion_operator &r,
                                           const fermion_operator &s)
  {
    const std::size_t term_index = excitations_.size();
    add_double_excitation(i, j, r, s, {{term_index, 1.0}}, true);
  }

  void uccsd_ansatz::prepare_reference_state()
  {
    const std::size_t n_occ = env_.occupied_orbitals();
    if (env_.xacc_indexing)
    {
      for (std::size_t i = 0; i < n_occ; ++i)
      {
        circuit_.add_x(i);
        circuit_.add_x(i + env_.n_spatial);
      }
    }
    else
    {
      for (std::size_t i = 0; i < 2 * n_occ; ++i)
      {
        circuit_.add_x(i);
      }
    }
  }

  void uccsd_ansatz::append_pauli_generators()
  {
    if (excitations_.empty())
    {
      return;
    }

    const std::size_t trotter = std::max<std::size_t>(trotter_steps_, 1);
    for (std::size_t step = 0; step < trotter; ++step)
    {
      const std::size_t offset = step * unique_parameters_;
      for (std::size_t idx = 0; idx < excitations_.size(); ++idx)
      {
        const auto &excitation = excitations_[idx];
        if (excitation.symmetry_terms.empty())
        {
          continue;
        }
        auto expr = build_expression(excitation, excitation_to_parameter_, offset);
        const auto &paulis = pauli_generators_[idx];
        bool used = false;
        for (const auto &term : paulis)
        {
          if (!is_non_trivial(term))
          {
            continue;
          }
          const double coeff = term.coefficient.real();
          if (std::abs(coeff) < kCoefficientCutoff)
          {
            continue;
          }
          circuit_.add_pauli_rotation(term, expr, 2.0 * coeff);
          used = true;
        }
        (void)used;
      }
    }
  }

} // namespace vqe
