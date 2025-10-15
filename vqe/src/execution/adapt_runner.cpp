#include "execution/adapt_runner.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <nlopt.hpp>
#ifdef VQE_ENABLE_MPI
#include <mpi.h>
#endif

#include "ansatz/ansatz.hpp"
#include "backend/statevector_cpu.hpp"
#ifdef VQE_ENABLE_CUDA
#include "backend/statevector_gpu.hpp"
#endif
#include "hamiltonian_parser.hpp"
#include "jw_transform.hpp"
#include "nwq_util.hpp"
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
                               double energy)
    {
      if (!out.is_open())
      {
        return;
      }
      out << std::setprecision(16);
      out << "# Iteration " << iteration << std::endl;
      out << "# Energy: " << energy << std::endl;
      out << "# Number of parameters: " << parameters.size() << std::endl;
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

    template <typename Backend>
    struct objective_context
    {
      circuit *circ = nullptr;
      Backend *backend = nullptr;
      const std::vector<pauli_term> *terms = nullptr;
      std::size_t *eval_count = nullptr;
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

      ++(*ctx->eval_count);
      grad.assign(x.size(), 0.0);
      return ctx->backend->expectation(*ctx->terms).real();
    }

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

    std::ofstream param_file;
    if (options.adapt_save_interval > 0)
    {
      std::string output_filename = hamiltonian_path + "-adapt_params.txt";
      param_file.open(output_filename);
      if (param_file.is_open())
      {
        param_file << "# ADAPT-VQE Parameters for " << hamiltonian_path << std::endl;
        param_file << "# Save interval: every " << options.adapt_save_interval << " iterations" << std::endl;
        param_file << std::endl;
      }
    }

    if (options.verbose)
    {
      const auto stats = summarize_parameters(ansatz.mutable_circuit().parameters());
      NWQSim::safe_print("[adapt] Hamiltonian has %zu Pauli terms\n", pauli_terms.size());
      NWQSim::safe_print("[adapt] Initialization: pool_size=%zu gradient_step=%s gradient_tol=%s\n",
                         ansatz.pool_operator_components().size(),
                         format_double(options.adapt_gradient_step, 6).c_str(),
                         format_double(options.adapt_gradient_tolerance, 6).c_str());
      NWQSim::safe_print("[adapt] Reference energy: %s parameters=%zu max|θ|=%s ||θ||₂=%s\n",
                         format_double(reference_energy, 12).c_str(),
                         stats.count,
                         format_double(stats.max_abs, 6).c_str(),
                         format_double(stats.l2_norm, 6).c_str());
      NWQSim::safe_print("Running ADAPT-VQE...\n");
      NWQSim::safe_print(" Iter   Objective Value     # Evals  Grad Norm    Time     |  #1q Gates  #2q Gates  |  Selected Operator\n");
      NWQSim::safe_print("-----------------------------------------------------------------------------------------------------------------------\n");
    }

    const auto &pool_components = ansatz.pool_operator_components();
    const auto pool_size = pool_components.size();
    const auto &pool_excitations = ansatz.pool_excitations();
    std::vector<bool> selected(pool_size, false);

    std::size_t total_energy_evals = 0;
    bool converged = false;
    double previous_energy = reference_energy;

    for (std::size_t iter = 0; iter < options.adapt_max_iterations; ++iter)
    {
      auto iteration_start = std::chrono::steady_clock::now();

      backend.reset();
      backend.apply(ansatz.mutable_circuit());
      double base_energy = compute_energy(backend, pauli_terms);
      ++total_energy_evals;

      if (options.adapt_log_memory)
      {
        NWQSim::safe_print("[adapt] iteration %zu RSS: %.2f MiB\n", iter, current_rss_mebibytes());
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

        for (const auto &component : pool_components[idx])
        {
          apply_pool_operator(scratch, backend.num_qubits(), component.terms,
                              -2.0 * options.adapt_gradient_step * component.parameter_scale);
        }
        const double energy_minus = compute_energy(scratch, pauli_terms);
        ++total_energy_evals;

        const double gradient = (energy_plus - energy_minus) / (2.0 * options.adapt_gradient_step);
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

      if (max_index == pool_size || max_gradient < options.adapt_gradient_tolerance)
      {
        auto iteration_end = std::chrono::steady_clock::now();
        std::chrono::duration<double> iteration_duration = iteration_end - iteration_start;
        double elapsed_seconds = iteration_duration.count();
        std::string time_str = format_duration(elapsed_seconds);

        converged = true;
        result.energy = base_energy;
        result.iterations = iter;
        if (options.verbose)
        {
          NWQSim::safe_print("%4zu  %20.12f  %7zu  %.3e  %9s  |  converged\n",
                             iter,
                             base_energy,
                             total_energy_evals,
                             max_gradient,
                             time_str.c_str());
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

      objective_context<Backend> ctx{&circ, &backend, &pauli_terms, &total_energy_evals};
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
      nlopt::result status;
      try
      {
        status = opt.optimize(current_params, min_value);
      }
      catch (const std::exception &ex)
      {
        throw std::runtime_error(std::string("NLopt optimization failed: ") + ex.what());
      }

      for (std::size_t i = 0; i < current_params.size(); ++i)
      {
        circ.set_parameter(i, current_params[i]);
      }

      backend.reset();
      backend.apply(circ);
      const double expectation = backend.expectation(pauli_terms).real();
      const double optimized_energy = expectation;
      ++total_energy_evals;

      result.energy = optimized_energy;
      if (status > 0)
      {
        result.converged = true;
      }

      if (options.verbose)
      {
        auto iteration_end = std::chrono::steady_clock::now();
        std::chrono::duration<double> iteration_duration = iteration_end - iteration_start;
        double elapsed_seconds = iteration_duration.count();
        std::string time_str = format_duration(elapsed_seconds);

        const auto counts = accumulate_gate_tally(circ);
        const auto formatted_label = format_operator_label(selected_label);
        NWQSim::safe_print("%4zu  %20.12f  %7zu  %.3e  %9s  |  %9d  %9d  |  %s\n",
                           iter+1,
                           optimized_energy,
                           total_energy_evals,
                           max_gradient,
                           time_str.c_str(),
                           counts.single_qubit,
                           counts.two_qubit,
                           formatted_label.c_str());
      }

      if (options.adapt_save_interval > 0 && (iter + 1) % options.adapt_save_interval == 0)
      {
        std::vector<std::string> current_labels;
        for (auto idx : ansatz.selected_indices())
        {
          if (idx < pool_excitations.size())
          {
            current_labels.push_back(pool_excitations[idx].label);
          }
        }
        save_adapt_parameters(param_file, iter + 1, current_params, current_labels, optimized_energy);
      }

      if (options.adapt_energy_tolerance > 0.0 &&
          std::abs(previous_energy - optimized_energy) < options.adapt_energy_tolerance)
      {
        if (options.verbose)
        {
          NWQSim::safe_print("[adapt][iter %zu] energy change %s < tol=%s\n",
                             (iter + 1),
                             format_double(std::abs(previous_energy - optimized_energy), 6).c_str(),
                             format_double(options.adapt_energy_tolerance, 6).c_str());
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
    if (options.verbose && !result.selected_labels.empty())
    {
      NWQSim::safe_print("Final parameters:\n");
      for (std::size_t i = 0; i < result.selected_labels.size() && i < result.parameters.size(); ++i)
      {
        const auto formatted_label = format_operator_label(result.selected_labels[i]);
        NWQSim::safe_print("  %s :: %.16f\n", formatted_label.c_str(), result.parameters[i]);
      }
    }
    result.energy_evaluations = total_energy_evals;
    result.hamiltonian_terms = pauli_terms.size();
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

#ifdef VQE_ENABLE_CUDA
    if (options.use_gpu)
    {
  return run_adapt_impl<backend::statevector_gpu>(pauli_terms, ansatz, options, hamiltonian_path);
    }
#else
    if (options.use_gpu)
    {
      throw std::runtime_error("vqe built without CUDA support; GPU backend unavailable");
    }
#endif

  return run_adapt_impl<backend::statevector_cpu>(pauli_terms, ansatz, options, hamiltonian_path);
  }

} // namespace vqe
