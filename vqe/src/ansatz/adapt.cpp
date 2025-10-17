#include "ansatz/ansatz.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace vqe
{

  namespace
  {
    constexpr double kCoefficientCutoff = 1e-10;

    bool is_non_trivial(const pauli_term &term)
    {
      return term.x_mask != 0 || term.z_mask != 0;
    }

  } // namespace

  adapt_ansatz::adapt_ansatz(const molecular_environment &env, std::size_t symmetry_level)
      : ansatz_base(env.total_qubits()),
        env_(env),
        symmetry_level_(symmetry_level)
  {
    reset_circuit();
  }

  void adapt_ansatz::build_pool()
  {
    uccsd_ansatz generator(env_, 1, symmetry_level_);
    generator.build();
    const auto &excitations = generator.excitations();
    const auto &pauli_generators = generator.pauli_generators();

    // Get pool name from the generator
    pool_name_ = generator.pool_name();

    pool_excitations_.clear();
    pool_components_.clear();

    pool_excitations_.reserve(excitations.size());
    pool_components_.reserve(excitations.size());

    std::vector<std::vector<std::string>> group_labels;
    group_labels.reserve(excitations.size());

    const auto invalid_index = std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> excitation_to_group(excitations.size(), invalid_index);

    for (std::size_t idx = 0; idx < excitations.size(); ++idx)
    {
      const auto &info = excitations[idx];
      const auto &terms = pauli_generators[idx];

      const bool is_representative = info.symmetry_terms.empty() ||
                                     std::any_of(info.symmetry_terms.begin(), info.symmetry_terms.end(),
                                                 [idx](const auto &entry) { return entry.first == idx; });

      if (is_representative)
      {
        const std::size_t group_id = pool_components_.size();
        excitation_to_group[idx] = group_id;

        pool_components_.emplace_back();
        pool_components_.back().push_back(pool_operator_component{terms, 1.0});

        pool_excitations_.emplace_back();
        group_labels.emplace_back();
        group_labels.back().push_back(info.label);
        continue;
      }

      if (info.symmetry_terms.empty())
      {
        throw std::runtime_error("symmetry-linked excitation missing representative");
      }

      const std::size_t owner_idx = info.symmetry_terms.front().first;
      if (owner_idx >= excitation_to_group.size())
      {
        throw std::out_of_range("symmetry term references invalid excitation index");
      }
      const std::size_t group_id = excitation_to_group[owner_idx];
      if (group_id == invalid_index)
      {
        throw std::runtime_error("symmetry representative encountered after dependent excitation");
      }

      double scale = 0.0;
      for (const auto &entry : info.symmetry_terms)
      {
        if (entry.first == owner_idx)
        {
          scale += entry.second;
        }
      }
      if (scale == 0.0)
      {
        scale = 1.0;
      }

      pool_components_[group_id].push_back(pool_operator_component{terms, scale});
      group_labels[group_id].push_back(info.label);
      excitation_to_group[idx] = group_id;
    }

    for (std::size_t group_id = 0; group_id < pool_excitations_.size(); ++group_id)
    {
      auto &info = pool_excitations_[group_id];
      const auto &labels = group_labels[group_id];
      std::string combined;
      combined.reserve(labels.size() * 8);
      for (std::size_t i = 0; i < labels.size(); ++i)
      {
        if (i > 0)
        {
          combined += ", ";
        }
        combined += labels[i];
      }
      info.label = std::move(combined);
    }

    reset_circuit();
  }

  std::size_t adapt_ansatz::add_operator(std::size_t pool_index, double initial_parameter)
  {
    if (pool_index >= pool_components_.size())
    {
      throw std::out_of_range("pool index out of range");
    }
    const auto &components = pool_components_[pool_index];
    if (components.empty())
    {
      throw std::runtime_error("selected pool element has no symmetry components");
    }
    const std::size_t param_index = circuit_.add_parameter(initial_parameter);
    bool used = false;
    for (const auto &component : components)
    {
      if (component.terms.empty())
      {
        continue;
      }
      if (std::abs(component.parameter_scale) < 1e-12)
      {
        continue;
      }
      parameter_expression expr;
      expr.terms.push_back(parameter_term{param_index, component.parameter_scale});
      for (const auto &term : component.terms)
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
    }
    if (!used)
    {
      throw std::runtime_error("selected pool element produced no gates");
    }
    selected_indices_.push_back(pool_index);
    return param_index;
  }

  void adapt_ansatz::reset_circuit()
  {
    circuit_ = circuit(env_.total_qubits());
    selected_indices_.clear();
    prepare_reference_state();
  }

  void adapt_ansatz::prepare_reference_state()
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

} // namespace vqe
