#pragma once

#include <cstddef>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "pauli_term.hpp"
#include "core/circuit.hpp"
#include "core/environment.hpp"

namespace vqe
{

  struct excitation_info
  {
    std::vector<fermion_operator> operators;
    std::vector<std::pair<std::size_t, double>> symmetry_terms;
    std::string label;
  };

  struct pool_operator_component
  {
    std::vector<pauli_term> terms;
    double parameter_scale = 1.0;
  };

  class ansatz_base
  {
  public:
    explicit ansatz_base(std::size_t qubits) : circuit_(qubits) {}
    virtual ~ansatz_base() = default;

    [[nodiscard]] const circuit &get_circuit() const { return circuit_; }
    [[nodiscard]] circuit &mutable_circuit() { return circuit_; }

    [[nodiscard]] const std::vector<double> &parameters() const { return circuit_.parameters(); }
    [[nodiscard]] std::size_t parameter_count() const { return circuit_.parameters().size(); }

  protected:
    circuit circuit_;
  };

  class uccsd_ansatz : public ansatz_base
  {
  public:
    uccsd_ansatz(const molecular_environment &env,
                 std::size_t trotter_steps = 1,
                 std::size_t symmetry_level = 3);

    void build();

    [[nodiscard]] const molecular_environment &environment() const { return env_; }
    [[nodiscard]] std::size_t unique_parameter_count() const { return unique_parameters_; }
    [[nodiscard]] const std::vector<excitation_info> &excitations() const { return excitations_; }
    [[nodiscard]] const std::vector<std::vector<pauli_term>> &pauli_generators() const { return pauli_generators_; }
    [[nodiscard]] const std::unordered_map<std::string, std::size_t> &excitation_parameter_map() const
    {
      return excitation_to_parameter_label_;
    }

  private:
    void generate_single_excitations();
    void generate_same_spin_doubles();
    void generate_mixed_spin_doubles();
    void add_single_excitation(const fermion_operator &from,
                               const fermion_operator &to,
                               const std::vector<std::pair<std::size_t, double>> &symmetry,
                               bool owns_parameter);
    void add_double_excitation(const fermion_operator &i,
                               const fermion_operator &j,
                               const fermion_operator &r,
                               const fermion_operator &s,
                               const std::vector<std::pair<std::size_t, double>> &symmetry,
                               bool owns_parameter);

    void add_single_excitation(const fermion_operator &from,
                               const fermion_operator &to);
    void add_double_excitation(const fermion_operator &i,
                               const fermion_operator &j,
                               const fermion_operator &r,
                               const fermion_operator &s);

    void prepare_reference_state();
    void append_pauli_generators();

    const molecular_environment env_;
    std::size_t trotter_steps_ = 1;
    std::size_t symmetry_level_ = 0;

    std::vector<excitation_info> excitations_;
    std::vector<std::size_t> excitation_to_parameter_;
    std::unordered_map<std::string, std::size_t> excitation_to_parameter_label_;
    std::vector<std::vector<pauli_term>> pauli_generators_;
    std::size_t unique_parameters_ = 0;
  };

  class adapt_ansatz : public ansatz_base
  {
  public:
    explicit adapt_ansatz(const molecular_environment &env,
                          std::size_t symmetry_level = 3);

    void build_pool();
    std::size_t pool_size() const { return pool_components_.size(); }
    const std::vector<excitation_info> &pool_excitations() const { return pool_excitations_; }
    const std::vector<std::vector<pool_operator_component>> &pool_operator_components() const { return pool_components_; }
    const std::vector<std::size_t> &selected_indices() const { return selected_indices_; }

    std::size_t add_operator(std::size_t pool_index, double initial_parameter = 0.0);
    void reset_circuit();

  private:
    void prepare_reference_state();

    molecular_environment env_;
    std::size_t symmetry_level_ = 0;
    std::vector<excitation_info> pool_excitations_;
    std::vector<std::vector<pool_operator_component>> pool_components_;
    std::vector<std::size_t> selected_indices_;
  };

} // namespace vqe
