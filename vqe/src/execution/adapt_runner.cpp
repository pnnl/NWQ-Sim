#include "execution/adapt_runner.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <nlopt.hpp>
#include <chrono>
#ifdef VQE_ENABLE_MPI
#include <mpi.h>
#endif

#include "ansatz/ansatz.hpp"
#include "backend/statevector_cpu.hpp"
#if defined(VQE_ENABLE_CUDA) || defined(VQE_ENABLE_HIP)
#include "backend/statevector_gpu.hpp"
#endif
#include "hamiltonian_parser.hpp"
#include "jw_transform.hpp"
#include <sys/resource.h>

namespace vqe
{
#ifdef VQE_ENABLE_MPI
// Safe if MPI is already initialized elsewhere
struct MpiGuard 
{
    bool owner = false;
    MpiGuard()
	{
        int init = 0; MPI_Initialized(&init);
        if (!init) { MPI_Init(nullptr, nullptr); owner = true; }
    }
    ~MpiGuard() 
	{
        int finalized = 0; MPI_Finalized(&finalized);
        if (owner && !finalized) MPI_Finalize();
    }
};
#endif

  namespace
  {

    constexpr double kCoeffCutoff = 1e-10;

    struct parameter_summary
    {
      std::size_t count = 0;
      double max_abs = 0.0;
      double l2_norm = 0.0;
    };

    struct gate_tally
    {
      int single_qubit = 0;
      int two_qubit = 0;
    };

    std::string trim_copy(const std::string &value)
    {
      auto begin = value.begin();
      auto end = value.end();
      while (begin != end && std::isspace(static_cast<unsigned char>(*begin)))
      {
        ++begin;
      }
      while (end != begin && std::isspace(static_cast<unsigned char>(*(end - 1))) )
      {
        --end;
      }
      return std::string(begin, end);
    }

    std::string format_operator_label(const std::string &label)
    {
      std::vector<std::string> pieces;
      std::stringstream ss(label);
      std::string token;
      while (std::getline(ss, token, ','))
      {
        const auto trimmed = trim_copy(token);
        if (!trimmed.empty())
        {
          pieces.emplace_back("(" + trimmed + ")");
        }
      }
      if (pieces.empty())
      {
        if (label.empty())
        {
          return "<none>";
        }
        return "(" + label + ")";
      }
      std::ostringstream out;
      for (std::size_t idx = 0; idx < pieces.size(); ++idx)
      {
        if (idx > 0)
        {
          out << ", ";
        }
        out << pieces[idx];
      }
      return out.str();
    }

    int pauli_basis_single_qubit_cost(bool has_x, bool has_z)
    {
      if (has_x && has_z)
      {
        return 4;
      }
      if (has_x)
      {
        return 2;
      }
      if (has_z)
      {
        return 0;
      }
      return 0;
    }

    gate_tally pauli_rotation_tally(const gate &g)
    {
      gate_tally tally;
      for (std::size_t idx = 0; idx < g.qubits.size(); ++idx)
      {
        const std::size_t qubit = g.qubits[idx];
        const bool has_x = (g.pauli.x_mask >> qubit) & 1ULL;
        const bool has_z = (g.pauli.z_mask >> qubit) & 1ULL;
        tally.single_qubit += pauli_basis_single_qubit_cost(has_x, has_z);
      }
      if (!g.qubits.empty())
      {
        tally.single_qubit += 1;
      }
      if (g.qubits.size() > 1)
      {
        tally.two_qubit += static_cast<int>(2 * (g.qubits.size() - 1));
      }
      return tally;
    }

    gate_tally accumulate_gate_tally(const circuit &circ)
    {
      gate_tally tally;
      for (const auto &g : circ.gates())
      {
        switch (g.kind)
        {
        case gate_kind::x:
        case gate_kind::h:
        case gate_kind::rz:
          ++tally.single_qubit;
          break;
        case gate_kind::cnot:
          ++tally.two_qubit;
          break;
        case gate_kind::pauli_rotation:
        {
          const auto contrib = pauli_rotation_tally(g);
          tally.single_qubit += contrib.single_qubit;
          tally.two_qubit += contrib.two_qubit;
          break;
        }
        }
      }
      return tally;
    }

    parameter_summary summarize_parameters(const std::vector<double> &params)
    {
      parameter_summary summary;
      summary.count = params.size();
      double sum_sq = 0.0;
      for (double value : params)
      {
        const double abs_value = std::abs(value);
        if (abs_value > summary.max_abs)
        {
          summary.max_abs = abs_value;
        }
        sum_sq += value * value;
      }
      summary.l2_norm = std::sqrt(sum_sq);
      return summary;
    }

    std::string format_double(double value, int precision = 12)
    {
      std::ostringstream oss;
      oss << std::setprecision(precision) << value;
      return oss.str();
    }

    double current_rss_mebibytes()
    {
      struct rusage usage{};
      if (getrusage(RUSAGE_SELF, &usage) != 0)
      {
        return 0.0;
      }
#ifdef __APPLE__
      const double bytes = static_cast<double>(usage.ru_maxrss);
#else
      const double bytes = static_cast<double>(usage.ru_maxrss) * 1024.0;
#endif
      return bytes / (1024.0 * 1024.0);
    }

    // MZ: formatting time string for each ADAPT iteration
    std::string format_duration(double seconds)
    {
      if (seconds < 60.0)
      {
        // Less than 1 minute: show as seconds with 2 decimal places
        char buffer[16];
        std::snprintf(buffer, sizeof(buffer), "%.2fs", seconds);
        return std::string(buffer);
      }
      else
      {
        // 1 minute or more: show as "XXXm YYs"
        int minutes = static_cast<int>(seconds / 60.0);
        int secs = static_cast<int>(seconds) % 60;
        char buffer[16];
        std::snprintf(buffer, sizeof(buffer), "%dm %02ds", minutes, secs);
        return std::string(buffer);
      }
    }

    void save_adapt_parameters(std::ofstream &out,
                               std::size_t iteration,
                               const std::vector<double> &parameters,
                               const std::vector<std::string> &selected_labels,
                               const std::vector<std::size_t> &selected_indices,
                               double energy)
    {
      if (!out.is_open())
      {
        return;
      }
      out << std::setprecision(16);
      out << std::endl;
      out << "# Iteration " << iteration << std::endl;
      out << "# Energy: " << energy << std::endl;
      out << "# Number of parameters: " << parameters.size() << std::endl;
      out << "# Operator indices: ";
      for (std::size_t i = 0; i < selected_indices.size(); ++i)
      {
        if (i > 0) out << ", ";
        out << selected_indices[i];
      }
      out << std::endl;
      for (std::size_t i = 0; i < parameters.size(); ++i)
      {
        if (i < selected_labels.size())
        {
          out << "# " << selected_labels[i];
        }
        out << "::" << parameters[i] << std::endl;
      }
      out << std::endl;
      out.flush();
    }

    struct adapt_state
    {
      std::string pool_name;
      std::size_t symmetry_level = 0;
      std::size_t iteration = 0;
      double energy = 0.0;
      std::vector<std::size_t> operator_indices;
      std::vector<double> parameters;
      std::vector<std::string> operator_labels;
    };

    bool load_adapt_state(const std::string &filename, adapt_state &state, std::string &error)
    {
      std::ifstream file(filename);
      if (!file.is_open())
      {
        error = "Could not open file: " + filename;
        return false;
      }

      std::string line;
      adapt_state latest_state;
      bool found_pool = false;
      bool found_symmetry = false;
      bool in_iteration = false;

      while (std::getline(file, line))
      {
        // Skip empty lines
        if (line.empty()) {
          if (in_iteration) {
            // End of iteration block, save this as latest state
            state = latest_state;
            in_iteration = false;
          }
          continue;
        }

        // Parse header info
        if (line.find("# Pool: ") == 0)
        {
          latest_state.pool_name = line.substr(8);
          found_pool = true;
          continue;
        }

        if (line.find("# Symmetry level: ") == 0)
        {
          try {
            latest_state.symmetry_level = std::stoull(line.substr(18));
            found_symmetry = true;
          } catch (...) {
            error = "Invalid symmetry level in file";
            return false;
          }
          continue;
        }

        // Parse iteration block
        if (line.find("# Iteration ") == 0)
        {
          in_iteration = true;
          latest_state.operator_indices.clear();
          latest_state.parameters.clear();
          latest_state.operator_labels.clear();
          try {
            latest_state.iteration = std::stoull(line.substr(12));
          } catch (...) {
            error = "Invalid iteration number in file";
            return false;
          }
          continue;
        }

        if (line.find("# Energy: ") == 0)
        {
          try {
            latest_state.energy = std::stod(line.substr(10));
          } catch (...) {
            error = "Invalid energy value in file";
            return false;
          }
          continue;
        }

        if (line.find("# Operator indices: ") == 0)
        {
          std::string indices_str = line.substr(20);
          std::istringstream iss(indices_str);
          std::string token;
          while (std::getline(iss, token, ','))
          {
            // Trim whitespace
            token.erase(0, token.find_first_not_of(" \t"));
            token.erase(token.find_last_not_of(" \t") + 1);
            try {
              latest_state.operator_indices.push_back(std::stoull(token));
            } catch (...) {
              error = "Invalid operator index in file";
              return false;
            }
          }
          continue;
        }

        // Parse parameter lines (format: "# label::value")
        if (line.find("::") != std::string::npos)
        {
          std::size_t pos = line.find("::");
          std::string label_part = line.substr(0, pos);
          std::string value_part = line.substr(pos + 2);

          // Extract label (remove "# " prefix if present)
          std::string label;
          if (label_part.size() > 2 && label_part[0] == '#' && label_part[1] == ' ')
          {
            label = label_part.substr(2);
          }

          try {
            double value = std::stod(value_part);
            latest_state.parameters.push_back(value);
            if (!label.empty())
            {
              latest_state.operator_labels.push_back(label);
            }
          } catch (...) {
            error = "Invalid parameter value in file";
            return false;
          }
          continue;
        }
      }

      // Save final iteration if we were in one
      if (in_iteration)
      {
        state = latest_state;
      }

      if (!found_pool)
      {
        error = "File does not contain pool name";
        return false;
      }

      if (!found_symmetry)
      {
        error = "File does not contain symmetry level";
        return false;
      }

      if (state.iteration == 0)
      {
        error = "No valid iteration found in file";
        return false;
      }

      return true;
    }

    // Helper function to check if optimizer needs gradients based on algorithm name
    inline bool optimizer_needs_gradient(nlopt::algorithm algo)
    {
      const std::string name = nlopt_algorithm_name(static_cast<nlopt_algorithm>(algo));
      // Check if algorithm name contains "derivative" (not "no-derivative")
      // Naming: NLOPT_{G/L}{D/N}_* where D=derivative-based, N=no-derivative
      return (name.find("derivative-based") != std::string::npos);
    }

    template <typename Backend>
    struct objective_context
    {
      circuit *circ = nullptr;
      Backend *backend = nullptr;
      const std::vector<pauli_term> *terms = nullptr;
      std::size_t *eval_count = nullptr;
      bool compute_gradient = false;
      bool use_spsa_gradient = false;  // Use SPSA-style gradient estimation
      std::mt19937 *rng = nullptr;     // Random number generator for SPSA
      double gradient_step = 1e-5;     // Forward-difference step size
    };

    template <typename Backend>
    double objective_function_impl(const std::vector<double> &x,
                                   std::vector<double> &grad,
                                   void *data)
    {
      auto *ctx = static_cast<objective_context<Backend> *>(data);
      if (ctx == nullptr || ctx->circ == nullptr || ctx->backend == nullptr || ctx->terms == nullptr || ctx->eval_count == nullptr)
      {
        throw std::runtime_error("objective context not initialized");
      }

      auto &circuit_ref = *ctx->circ;
      for (std::size_t i = 0; i < x.size(); ++i)
      {
        circuit_ref.set_parameter(i, x[i]);
      }

      ctx->backend->reset();
      ctx->backend->apply(circuit_ref);
      const double energy = ctx->backend->expectation(*ctx->terms).real();
      ++(*ctx->eval_count);

      if (!std::isfinite(energy))
      {
        throw std::runtime_error("Non-finite energy value: " + std::to_string(energy));
      }

      // Compute gradients if needed (for derivative-based optimizers like LD_LBFGS)
      if (ctx->compute_gradient && !grad.empty())
      {
        if (ctx->use_spsa_gradient && ctx->rng != nullptr)
        {
          // SPSA gradient estimation: only 2 function evaluations regardless of dimensionality
          const double epsilon = 1e-4;  // SPSA perturbation size
          std::uniform_int_distribution<int> dist(0, 1);
          std::vector<double> delta(x.size());

          // Generate random perturbation direction (±1 for each parameter)
          for (std::size_t i = 0; i < x.size(); ++i)
          {
            delta[i] = (dist(*ctx->rng) == 0) ? -1.0 : 1.0;
          }

          // Evaluate at x + ε*delta
          for (std::size_t i = 0; i < x.size(); ++i)
          {
            circuit_ref.set_parameter(i, x[i] + epsilon * delta[i]);
          }
          ctx->backend->reset();
          ctx->backend->apply(circuit_ref);
          const double energy_plus = ctx->backend->expectation(*ctx->terms).real();
          ++(*ctx->eval_count);

          // Evaluate at x - ε*delta
          for (std::size_t i = 0; i < x.size(); ++i)
          {
            circuit_ref.set_parameter(i, x[i] - epsilon * delta[i]);
          }
          ctx->backend->reset();
          ctx->backend->apply(circuit_ref);
          const double energy_minus = ctx->backend->expectation(*ctx->terms).real();
          ++(*ctx->eval_count);

          // Compute gradient approximation: grad[i] ≈ (E+ - E-) / (2ε * delta[i])
          for (std::size_t i = 0; i < x.size(); ++i)
          {
            grad[i] = (energy_plus - energy_minus) / (2.0 * epsilon * delta[i]);
          }

          // Restore parameters
          for (std::size_t i = 0; i < x.size(); ++i)
          {
            circuit_ref.set_parameter(i, x[i]);
          }
        }
        else
        {
          // Standard finite difference gradient
          const double epsilon = ctx->gradient_step;  // finite difference step size
          if (epsilon <= 0.0)
          {
            throw std::runtime_error("Finite-difference step must be positive");
          }

          for (std::size_t i = 0; i < x.size(); ++i)
          {
            const double x_original = x[i];

            // Forward difference: grad[i] = [E(x + ε) - E(x)] / ε
            circuit_ref.set_parameter(i, x_original + epsilon);
            ctx->backend->reset();
            ctx->backend->apply(circuit_ref);
            const double energy_plus = ctx->backend->expectation(*ctx->terms).real();
            ++(*ctx->eval_count);

            if (!std::isfinite(energy_plus))
            {
              throw std::runtime_error("Non-finite energy_plus at param " + std::to_string(i));
            }

            grad[i] = (energy_plus - energy) / epsilon;
            circuit_ref.set_parameter(i, x_original);
          }
        }
      }
      else
      {
        grad.assign(x.size(), 0.0);
      }

      return energy;
    } // objective_function_impl() ends


    template <typename Backend>
    double compute_energy(Backend &backend,
                          const std::vector<pauli_term> &terms)
    {
      return backend.expectation(terms).real();
    }

    template <typename Backend>
    void apply_pool_operator(Backend &backend,
                             std::size_t num_qubits,
                             const std::vector<pauli_term> &pauli_terms,
                             double parameter)
    {
      if (std::abs(parameter) < 1e-12)
      {
        return;
      }
      circuit incremental(num_qubits);
      const std::size_t param_index = incremental.add_parameter(parameter);
      parameter_expression expr;
      expr.terms.push_back(parameter_term{param_index, 1.0});
      for (const auto &term : pauli_terms)
      {
        const double coeff = term.coefficient.real();
        if (std::abs(coeff) < kCoeffCutoff)
        {
          continue;
        }
        incremental.add_pauli_rotation(term, expr, 2.0 * coeff);
      }
      if (!incremental.gates().empty())
      {
        backend.apply(incremental);
      }
    }

  } // namespace

  template <typename Backend>
  adapt_result run_adapt_impl(const std::vector<pauli_term> &pauli_terms,
                              adapt_ansatz &ansatz,
                              const vqe_options &options,
                              const std::string &hamiltonian_path)
  {
#ifdef VQE_ENABLE_MPI
	MpiGuard guard;
	int world_rank = 0, world_size = 1;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
#endif

    Backend backend(ansatz.get_circuit().num_qubits());
    adapt_result result;
    backend.reset();
    backend.apply(ansatz.mutable_circuit());
    const double reference_energy = compute_energy(backend, pauli_terms);
    result.initial_energy = reference_energy;
    result.energy = reference_energy;
    const bool needs_grad = optimizer_needs_gradient(options.adapt_optimizer);

    // Load state from file if requested
    std::size_t start_iteration = 0;
    if (!options.adapt_load_state_file.empty())
    {
      adapt_state loaded_state;
      std::string load_error;
      if (!load_adapt_state(options.adapt_load_state_file, loaded_state, load_error))
      {
        throw std::runtime_error("Failed to load ADAPT state: " + load_error);
      }

      // Validate loaded state matches current configuration
      if (loaded_state.symmetry_level != ansatz.symmetry_level())
      {
        throw std::runtime_error("Symmetry level mismatch: file has " +
                                std::to_string(loaded_state.symmetry_level) +
                                ", current is " + std::to_string(ansatz.symmetry_level()));
      }

      // Rebuild circuit with loaded operators
      ansatz.reset_circuit();
      for (std::size_t idx : loaded_state.operator_indices)
      {
        if (idx >= ansatz.pool_size())
        {
          throw std::runtime_error("Invalid operator index " + std::to_string(idx) +
                                  " in state file (pool size is " + std::to_string(ansatz.pool_size()) + ")");
        }
        ansatz.add_operator(idx, 0.0);  // Add with zero initial parameter
      }

      // Set loaded parameters
      auto &circ = ansatz.mutable_circuit();
      if (loaded_state.parameters.size() != circ.parameters().size())
      {
        throw std::runtime_error("Parameter count mismatch: file has " +
                                std::to_string(loaded_state.parameters.size()) +
                                ", circuit has " + std::to_string(circ.parameters().size()));
      }
      for (std::size_t i = 0; i < loaded_state.parameters.size(); ++i)
      {
        circ.set_parameter(i, loaded_state.parameters[i]);
      }

      start_iteration = loaded_state.iteration;
      result.initial_energy = loaded_state.energy;
      result.energy = loaded_state.energy;

      if (options.verbose)
      {
        std::cout << "[adapt] Loaded state from: " << options.adapt_load_state_file << std::endl;
        std::cout << "[adapt] Resuming from iteration " << start_iteration
                  << " with energy " << format_double(loaded_state.energy, 12) << std::endl;
        std::cout << "[adapt] Loaded " << loaded_state.operator_indices.size()
                  << " operators and " << loaded_state.parameters.size() << " parameters" << std::endl;
      }
    }

    std::ofstream param_file;
    if (options.adapt_save_params)
    {
#ifdef VQE_ENABLE_MPI
      // Only rank 0 should write to file in MPI mode
      if (world_rank == 0)
#endif
      {
        std::string output_filename = hamiltonian_path + "-adapt_params.txt";
        // If resuming, append to existing file
        if (start_iteration > 0)
        {
          param_file.open(output_filename, std::ios::app);
        }
        else
        {
          param_file.open(output_filename);
          if (param_file.is_open())
          {
            param_file << "# ADAPT-VQE Parameters for " << hamiltonian_path << std::endl;
            param_file << "# Pool: " << ansatz.pool_name() << std::endl;
            param_file << "# Symmetry level: " << ansatz.symmetry_level() << std::endl;
            param_file << std::endl;
          }
        }
      }
    }

    if (options.verbose)
    {
      const auto stats = summarize_parameters(ansatz.mutable_circuit().parameters());
      std::cout << std::endl;
      std::cout << "[adapt] Initialization: pool_size=" << ansatz.pool_operator_components().size()
                << " # paulis=" << pauli_terms.size()
                << " gradient_step=" << format_double(options.adapt_gradient_step, 6)
                << " gradient_tol=" << format_double(options.adapt_gradient_tolerance, 6) << std::endl;
      std::cout << "[adapt] Reference energy: " << format_double(reference_energy, 12)
                << " parameters=" << stats.count
                << " max|θ|=" << format_double(stats.max_abs, 6)
                << " ||θ||₂=" << format_double(stats.l2_norm, 6) << std::endl;
      std::cout << "Running ADAPT-VQE..." << std::endl;
      std::cout << " Iter   Objective Value     # Evals  Grad Norm    Time   |  #1q Gates  #2q Gates  |  Selected Operator\n" << std::endl;
      std::cout << "-----------------------------------------------------------------------------------------------------------------------\n" << std::endl;
    }

    const auto &pool_components = ansatz.pool_operator_components();
    const auto pool_size = pool_components.size();
    const auto &pool_excitations = ansatz.pool_excitations();
    std::vector<bool> selected(pool_size, false);

    // Mark already selected operators if resuming
    if (start_iteration > 0)
    {
      for (std::size_t idx : ansatz.selected_indices())
      {
        if (idx < pool_size)
        {
          selected[idx] = true;
        }
      }
    }

    std::size_t total_energy_evals = 0; // MZ: I think what supposed to be printed is mistaken during the code re-organization
                                           // It was the number of optimizaton in each ADAPT round, not how many circuit evaluations
    bool converged = false;
    double previous_energy = (start_iteration > 0) ? result.energy : reference_energy;

    for (std::size_t iter = start_iteration; iter < start_iteration+options.adapt_max_iterations; ++iter)
    {
      auto iteration_start = std::chrono::steady_clock::now(); // MZ: timer start

      backend.reset();
      backend.apply(ansatz.mutable_circuit());
      double base_energy = compute_energy(backend, pauli_terms);
      ++total_energy_evals;

      if (options.adapt_log_memory)
      {
        std::cout << "[adapt] iteration " << iter << " RSS: " << current_rss_mebibytes() << " MiB" << std::endl;
      }

      double max_gradient = 0.0;
      std::size_t max_index = pool_size;
      Backend scratch(backend.num_qubits());

	  //=========================== MPI Parallelization =======================
	#ifdef VQE_ENABLE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
	  for (std::size_t idx = world_rank; idx < pool_size; idx += world_size) 
	  {
	#else
      for (std::size_t idx = 0; idx < pool_size; ++idx)
      {
	#endif
        if (selected[idx])
        {
          continue;
        }
        scratch.copy_state_from(backend);
        for (const auto &component : pool_components[idx])
        {
          apply_pool_operator(scratch, backend.num_qubits(), component.terms,
                              options.adapt_gradient_step * component.parameter_scale);
        }
        const double energy_plus = compute_energy(scratch, pauli_terms);
        ++total_energy_evals;

        const double gradient = (energy_plus - base_energy) / options.adapt_gradient_step;
        const double magnitude = std::abs(gradient);

        if (magnitude > max_gradient)
        {
          max_gradient = magnitude;
          max_index = idx;
        }
      }

	#ifdef VQE_ENABLE_MPI
	  // MPI_MAXLOC works with (double,int), so we map size_t to int temporarily
	  // If N might exceed INT_MAX, we need to use a custom reduction instead.
	  struct { double val; int idx; } in_pair, out_pair;
	  in_pair.val = max_gradient;
	  in_pair.idx = static_cast<int>(max_index);
	  MPI_Allreduce(&in_pair, &out_pair, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
	  max_gradient = out_pair.val;
	  max_index = static_cast<std::size_t>(out_pair.idx);
	#endif

	  //=========================================================================

      std::string selected_label = "<none>";
      if (max_index < pool_size && max_index < pool_excitations.size())
      {
        selected_label = pool_excitations[max_index].label;
      }

      if (max_index == pool_size || max_gradient < options.adapt_gradient_tolerance) // MZ: not sure when max_index == pool_size will happen
      {
        auto iteration_end = std::chrono::steady_clock::now();  // MZ: timer end
        std::chrono::duration<double> iteration_duration = iteration_end - iteration_start;
        double elapsed_seconds = iteration_duration.count();
        std::string time_str = format_duration(elapsed_seconds);

        converged = true;
        result.energy = base_energy;
        result.iterations = iter;
        if (options.verbose)
        {
          std::ostringstream oss;
          oss << std::setw(4) << iter+1 << "  "
              << std::setw(20) << std::fixed << std::setprecision(12) << base_energy << "  "
              << std::setw(7) << 0 << "  "
              << std::scientific << std::setprecision(3) << max_gradient
              << std::setw(9) << time_str << "  |  converged\n";
          std::cout << oss.str();
        }
        break;
      }

      ansatz.add_operator(max_index, 0.0);
      selected[max_index] = true;

      auto &circ = ansatz.mutable_circuit();
      auto current_params = circ.parameters();

      nlopt::opt opt(options.adapt_optimizer, current_params.size());
      const double lower = std::min(options.adapt_lower_bound, options.adapt_upper_bound);
      const double upper = std::max(options.adapt_lower_bound, options.adapt_upper_bound);
      opt.set_lower_bounds(std::vector<double>(current_params.size(), lower));
      opt.set_upper_bounds(std::vector<double>(current_params.size(), upper));
      if (options.adapt_max_evaluations > 0)
      {
        opt.set_maxeval(static_cast<int>(options.adapt_max_evaluations));
      }
      if (options.adapt_max_time > 0)
      {
        opt.set_maxtime(options.adapt_max_time);
      }
      if (options.adapt_relative_tolerance >= 0)
      {
        opt.set_ftol_rel(options.adapt_relative_tolerance);
      }
      if (options.adapt_absolute_tolerance >= 0)
      {
        opt.set_ftol_abs(options.adapt_absolute_tolerance);
      }
      if (std::isfinite(options.adapt_stop_value))
      {
        opt.set_stopval(options.adapt_stop_value);
      }

      for (const auto &param : options.adapt_algorithm_parameters)
      {
        opt.set_param(param.first.c_str(), param.second);
      }

      // Random number generator for SPSA (if enabled)
      static std::mt19937 rng(std::random_device{}());

      objective_context<Backend> ctx{&circ, &backend, &pauli_terms, &total_energy_evals,
                                      needs_grad, options.use_spsa_gradient, &rng,
                                      options.gradient_step};
      opt.set_min_objective(objective_function_impl<Backend>, &ctx);

      current_params = circ.parameters();
      if (options.adapt_initial_parameters.size() == current_params.size())
      {
        current_params = options.adapt_initial_parameters;
        for (std::size_t i = 0; i < current_params.size(); ++i)
        {
          circ.set_parameter(i, current_params[i]);
        }
      }
      double min_value = 0.0;
      int num_evals = 0;
      nlopt::result status;
      try
      {
        status = opt.optimize(current_params, min_value);
        num_evals = opt.get_numevals(); //MZ: number of optimization iterations is the iterative info to return
      }
      catch (const nlopt::roundoff_limited &ex)
      {
        // Roundoff errors - optimization may have converged as much as possible
        if (options.verbose)
        {
          std::cout << "[adapt][warning] NLopt roundoff limited: " << ex.what() << std::endl;
        }
        num_evals = opt.get_numevals();
        status = nlopt::ROUNDOFF_LIMITED;
      }
      catch (const nlopt::forced_stop &ex)
      {
        throw std::runtime_error(std::string("NLopt forced stop: ") + ex.what());
      }
      catch (const std::exception &ex)
      {
        throw std::runtime_error(std::string("NLopt optimization failed (iter=") + std::to_string(iter) +
                                 ", params=" + std::to_string(current_params.size()) +
                                 ", needs_grad=" + std::to_string(needs_grad) + "): " + ex.what());
      }

      for (std::size_t i = 0; i < current_params.size(); ++i)
      {
        circ.set_parameter(i, current_params[i]);
      }

      backend.reset();
      backend.apply(circ);
      const double expectation = backend.expectation(pauli_terms).real();
      double optimized_energy = expectation;
      ++total_energy_evals;

      result.energy = optimized_energy;
      if (status > 0)
      {
        result.converged = true;
      }

      auto iteration_end = std::chrono::steady_clock::now();  // MZ: timer end
      std::chrono::duration<double> iteration_duration = iteration_end - iteration_start;
      double elapsed_seconds = iteration_duration.count();
      std::string time_str = format_duration(elapsed_seconds);

      if (options.verbose)
      {
        const auto counts = accumulate_gate_tally(circ);
        const auto formatted_label = format_operator_label(selected_label);
        std::ostringstream row;
        row << std::setw(4) << iter+1 << "  ";
        row << std::fixed << std::setprecision(12) << std::setw(20) << optimized_energy << "  ";
        row << std::setw(7) << num_evals << "  ";
        row << std::scientific << std::setprecision(3) << std::setw(8) << max_gradient << std::setw(9) << time_str << "  |  ";
        row << std::defaultfloat;
        row << std::setw(9) << counts.single_qubit << "  "
            << std::setw(9) << counts.two_qubit << "  |  "
            << formatted_label;
        std::cout << row.str() << std::endl;
      }

      if (options.adapt_save_params)
      {
#ifdef VQE_ENABLE_MPI
        // Only rank 0 should write to file in MPI mode
        if (world_rank == 0)
#endif
        {
          std::vector<std::string> current_labels;
          for (auto idx : ansatz.selected_indices())
          {
            if (idx < pool_excitations.size())
            {
              current_labels.push_back(pool_excitations[idx].label);
            }
          }
          save_adapt_parameters(param_file, iter + 1, current_params, current_labels,
                              ansatz.selected_indices(), optimized_energy);
        }
      }

      if (options.adapt_energy_tolerance > 0.0 &&
          std::abs(previous_energy - optimized_energy) < options.adapt_energy_tolerance)
      {
        if (options.verbose)
        {
          std::cout << "[adapt][iter " << (iter + 1)
                    << "] energy change " << format_double(std::abs(previous_energy - optimized_energy), 6)
                    << " < tol=" << format_double(options.adapt_energy_tolerance, 6) << std::endl;

        }
        result.iterations = iter + 1;
        converged = true;
        break;
      }
      previous_energy = optimized_energy;

      if (iter + 1 == options.adapt_max_iterations)
      {
        result.iterations = iter + 1;
      }
    }

    if (!converged)
    {
      backend.reset();
      backend.apply(ansatz.mutable_circuit());
      result.energy = compute_energy(backend, pauli_terms);
      ++total_energy_evals;
    }

    result.parameters = ansatz.mutable_circuit().parameters();
    result.selected_indices = ansatz.selected_indices();
    result.selected_labels.reserve(result.selected_indices.size());
    for (auto idx : result.selected_indices)
    {
      if (idx >= pool_excitations.size())
      {
        throw std::out_of_range("selected operator index exceeds pool size");
      }
      result.selected_labels.push_back(pool_excitations[idx].label);
    }
    result.energy_evaluations = total_energy_evals;
    return result;
  }

  adapt_result run_adapt_vqe(const std::string &hamiltonian_path,
                             std::size_t n_particles,
                             const vqe_options &options)
  {
    const auto data = read_hamiltonian_file(hamiltonian_path);
    if (data.num_qubits() % 2 != 0)
    {
      throw std::runtime_error("Hamiltonian maps to an odd number of qubits; expected even for spin orbitals");
    }

    const auto pauli_terms = jordan_wigner_transform(data);
    const double constant_energy = data.constant.real();

    molecular_environment env;
    env.n_spatial = data.num_qubits() / 2;
    env.n_electrons = n_particles;
    env.xacc_indexing = options.use_xacc_indexing;
    env.constant_energy = constant_energy;

  adapt_ansatz ansatz(env, options.symmetry_level);
    ansatz.build_pool();

#if defined(VQE_ENABLE_CUDA) || defined(VQE_ENABLE_HIP)
    if (options.use_gpu)
    {
  return run_adapt_impl<backend::statevector_gpu>(pauli_terms, ansatz, options, hamiltonian_path);
    }
#else
    if (options.use_gpu)
    {
      throw std::runtime_error("vqe built without CUDA/HIP support; GPU backend unavailable");
    }
#endif

  return run_adapt_impl<backend::statevector_cpu>(pauli_terms, ansatz, options, hamiltonian_path);
  }

} // namespace vqe
