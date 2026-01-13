#include "execution/vqe_runner.hpp"

#include <nlopt.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "ansatz/ansatz.hpp"
#include "backend/statevector_cpu.hpp"
#if defined(VQE_ENABLE_CUDA) || defined(VQE_ENABLE_HIP)
#include "backend/statevector_gpu.hpp"
#endif
#include "core/environment.hpp"
#include "hamiltonian_parser.hpp"
#include "jw_transform.hpp"

namespace vqe
{
  namespace
  {

    struct parameter_summary
    {
      std::size_t count = 0;
      double max_abs = 0.0;
      double l2_norm = 0.0;
    };

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

    struct iteration_logger
    {
      void configure(bool verbose_mode, double tolerance_value, std::size_t status_stride)
      {
        verbose = verbose_mode;
        tolerance = tolerance_value;
        status_interval = status_stride;
      }

      void initialize(double initial, const std::vector<double> &params)
      {
        initial_energy = initial;
        best_energy = initial;
        previous_best = initial;
        iteration = 0;
        last_print_eval = 0;
        if (!verbose)
        {
          return;
        }
        const auto stats = summarize_parameters(params);
        std::cout << "[vqe] Initial energy: " << format_double(initial, 12) << std::endl;
        std::cout << "[vqe] Initial parameter stats: count=" << stats.count
                  << " max|θ|=" << format_double(stats.max_abs, 6)
                  << " ||θ||₂=" << format_double(stats.l2_norm, 6) << std::endl;
      }

      void log_evaluation(std::size_t eval, double energy, const std::vector<double> &params)
      {
        if (!verbose)
        {
          return;
        }

        const bool improved = (energy + tolerance) < best_energy;
        if (!improved)
        {
          if (status_interval > 0 && (eval - last_print_eval) >= status_interval)
          {
            const auto stats = summarize_parameters(params);
            std::cout << "[vqe][status] eval=" << eval
                      << " energy=" << format_double(energy, 12)
                      << " best=" << format_double(best_energy, 12)
                      << " max|θ|=" << format_double(stats.max_abs, 6)
                      << " ||θ||₂=" << format_double(stats.l2_norm, 6) << std::endl;
            last_print_eval = eval;
          }
          return;
        }

        previous_best = best_energy;
        best_energy = energy;
        ++iteration;
        last_print_eval = eval;

        const auto stats = summarize_parameters(params);
        std::cout << "[vqe][iter " << iteration << "] eval=" << eval
                  << " energy=" << format_double(energy, 12)
                  << " ΔE=" << format_double(energy - previous_best, 12)
                  << " max|θ|=" << format_double(stats.max_abs, 6)
                  << " ||θ||₂=" << format_double(stats.l2_norm, 6) << std::endl;
      }

      void finalize(double final_energy)
      {
        best_energy = final_energy;
      }

      std::size_t iteration_count() const
      {
        return iteration;
      }

      double initial_value() const
      {
        return initial_energy;
      }

      double best_value() const
      {
        return best_energy;
      }

    private:
      bool verbose = false;
      double tolerance = 1e-8;
      std::size_t status_interval = 0;
      double initial_energy = std::numeric_limits<double>::quiet_NaN();
      double best_energy = std::numeric_limits<double>::infinity();
      double previous_best = std::numeric_limits<double>::infinity();
      std::size_t iteration = 0;
      std::size_t last_print_eval = 0;
    };

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
      uccsd_ansatz *ansatz = nullptr;
      Backend *backend = nullptr;
      const std::vector<pauli_term> *terms = nullptr;
      std::size_t *eval_count = nullptr;
      double *apply_time = nullptr;
      double *expectation_time = nullptr;
      iteration_logger *logger = nullptr;
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
      if (ctx == nullptr || ctx->ansatz == nullptr || ctx->backend == nullptr || ctx->terms == nullptr || ctx->eval_count == nullptr)
      {
        throw std::runtime_error("objective context not initialized");
      }

      auto &circuit = ctx->ansatz->mutable_circuit();
      for (std::size_t i = 0; i < x.size(); ++i)
      {
        circuit.set_parameter(i, x[i]);
      }

      const auto apply_start = std::chrono::steady_clock::now();
      ctx->backend->reset();
      ctx->backend->apply(ctx->ansatz->get_circuit());
      const auto apply_stop = std::chrono::steady_clock::now();

      const auto expect_start = std::chrono::steady_clock::now();
      const double expectation = ctx->backend->expectation(*ctx->terms).real();
      const auto expect_stop = std::chrono::steady_clock::now();
      ++(*ctx->eval_count);

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
            circuit.set_parameter(i, x[i] + epsilon * delta[i]);
          }
          ctx->backend->reset();
          ctx->backend->apply(ctx->ansatz->get_circuit());
          const double energy_plus = ctx->backend->expectation(*ctx->terms).real();
          ++(*ctx->eval_count);

          // Evaluate at x - ε*delta
          for (std::size_t i = 0; i < x.size(); ++i)
          {
            circuit.set_parameter(i, x[i] - epsilon * delta[i]);
          }
          ctx->backend->reset();
          ctx->backend->apply(ctx->ansatz->get_circuit());
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
            circuit.set_parameter(i, x[i]);
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
            circuit.set_parameter(i, x_original + epsilon);
            ctx->backend->reset();
            ctx->backend->apply(ctx->ansatz->get_circuit());
            const double energy_plus = ctx->backend->expectation(*ctx->terms).real();
            ++(*ctx->eval_count);

            grad[i] = (energy_plus - expectation) / epsilon;
            circuit.set_parameter(i, x_original);
          }
        }
      }
      else
      {
        grad.assign(x.size(), 0.0);
      }

      if (ctx->apply_time != nullptr)
      {
        *ctx->apply_time += std::chrono::duration<double>(apply_stop - apply_start).count();
      }
      if (ctx->expectation_time != nullptr)
      {
        *ctx->expectation_time += std::chrono::duration<double>(expect_stop - expect_start).count();
      }

      if (ctx->logger != nullptr)
      {
        ctx->logger->log_evaluation(*ctx->eval_count, expectation, x);
      }

      return expectation;
    }

    template <typename Backend>
    vqe_result run_default_vqe_backend(uccsd_ansatz &ansatz,
                                       const std::vector<pauli_term> &pauli_terms,
                                       const vqe_options &options)
    {
      const std::size_t param_count = ansatz.unique_parameter_count();
      if (param_count == 0)
      {
        throw std::runtime_error("UCCSDmin ansatz has zero parameters for this system");
      }

      Backend backend(ansatz.get_circuit().num_qubits());

      std::vector<double> parameters(param_count, 0.0);
      if (options.initial_parameters.size() == param_count)
      {
        parameters = options.initial_parameters;
        if (options.verbose)
        {
          std::cout << "[vqe] Using custom initial parameters (" << param_count << " values)" << std::endl;
        }
      }
      else if (options.verbose && !options.initial_parameters.empty())
      {
        std::cout << "[vqe] Warning: Initial parameters size mismatch (got "
                  << options.initial_parameters.size() << ", need " << param_count
                  << "). Using zeros." << std::endl;
      }

      if (options.verbose)
      {
        std::size_t display_count = std::min(param_count, std::size_t(10));
        std::cout << "[vqe] Initial parameters: [";
        for (std::size_t i = 0; i < display_count; ++i)
        {
          if (i > 0) std::cout << ", ";
          std::cout << parameters[i];
        }
        if (param_count > 10)
        {
          std::cout << ", ...]" << std::endl;
        }
        else
        {
          std::cout << "]" << std::endl;
        }
      }

      auto &mutable_circuit = ansatz.mutable_circuit();
      for (std::size_t i = 0; i < parameters.size(); ++i)
      {
        mutable_circuit.set_parameter(i, parameters[i]);
      }

      backend.reset();
      backend.apply(ansatz.get_circuit());
      const double initial_energy = backend.expectation(pauli_terms).real();

      iteration_logger logger;
      logger.configure(options.verbose, options.iteration_improvement_tolerance, options.status_interval);
      logger.initialize(initial_energy, parameters);

      nlopt::opt opt(options.optimizer, param_count);

      const double lower = std::min(options.lower_bound, options.upper_bound);
      const double upper = std::max(options.lower_bound, options.upper_bound);
      opt.set_lower_bounds(std::vector<double>(param_count, lower));
      opt.set_upper_bounds(std::vector<double>(param_count, upper));

      if (options.max_evaluations > 0)
      {
        opt.set_maxeval(static_cast<int>(options.max_evaluations));
      }
      if (options.max_time > 0)
      {
        opt.set_maxtime(options.max_time);
      }
      if (options.relative_tolerance >= 0)
      {
        opt.set_ftol_rel(options.relative_tolerance);
      }
      if (options.absolute_tolerance >= 0)
      {
        opt.set_ftol_abs(options.absolute_tolerance);
      }
      if (std::isfinite(options.stop_value))
      {
        opt.set_stopval(options.stop_value);
      }

      for (const auto &entry : options.algorithm_parameters)
      {
        opt.set_param(entry.first.c_str(), entry.second);
      }

      std::size_t eval_count = 0;
      double apply_time = 0.0;
      double expectation_time = 0.0;
      const bool needs_grad = optimizer_needs_gradient(options.optimizer);

      // Random number generator for SPSA (if enabled)
      static std::mt19937 rng(std::random_device{}());
      objective_context<Backend> context{&ansatz, &backend, &pauli_terms, &eval_count, &apply_time, &expectation_time, &logger,
                                          needs_grad, options.use_spsa_gradient, &rng,
                                          options.gradient_step};
      opt.set_min_objective(objective_function_impl<Backend>, &context);

      double min_value = 0.0;
      nlopt::result status = nlopt::FAILURE;

      try
      {
        status = opt.optimize(parameters, min_value);
      }
      catch (const std::exception &ex)
      {
        throw std::runtime_error(std::string("NLopt optimization failed: ") + ex.what());
      }

      logger.finalize(min_value);

      vqe_result result;
      result.energy = min_value;
      result.parameters = std::move(parameters);
      result.evaluations = eval_count;
      result.converged = (status > 0);
      result.iterations = logger.iteration_count();
      result.initial_energy = logger.initial_value();
      result.energy_delta = min_value - logger.initial_value();
      result.apply_time = apply_time;
      result.expectation_time = expectation_time;
      result.status = status;

      const char *profile_env = std::getenv("RUN_VQE_PROFILE");
      if (profile_env != nullptr && profile_env[0] != '\0')
      {
        std::cout << "[profile] total apply time (s): " << apply_time << std::endl;
        std::cout << "[profile] total expectation time (s): " << expectation_time << std::endl;
      }

      return result;
    }

  } // namespace

  vqe_result run_default_vqe(const std::string &hamiltonian_path,
                             std::size_t n_particles,
                             const vqe_options &options)
  {
    const auto data = read_hamiltonian_file(hamiltonian_path);
    if (data.num_qubits() % 2 != 0)
    {
      throw std::runtime_error("Hamiltonian maps to an odd number of qubits; expected even for spin orbitals");
    }

    const auto pauli_terms = jordan_wigner_transform(data);

    molecular_environment env;
    env.n_spatial = data.num_qubits() / 2;
    env.n_electrons = n_particles;
    env.xacc_indexing = options.use_xacc_indexing;
    env.constant_energy = data.constant.real();

    uccsd_ansatz ansatz(env, options.trotter_steps, options.symmetry_level);
    ansatz.build();

    return run_default_vqe_with_ansatz(ansatz, pauli_terms, options);
  }

} // namespace vqe

vqe::vqe_result vqe::run_default_vqe_with_ansatz(uccsd_ansatz &ansatz,
                                                 const std::vector<pauli_term> &pauli_terms,
                                                 const vqe_options &options)
{
#if defined(VQE_ENABLE_CUDA) || defined(VQE_ENABLE_HIP)
  if (options.use_gpu)
  {
    return run_default_vqe_backend<backend::statevector_gpu>(ansatz, pauli_terms, options);
  }
#else
  if (options.use_gpu)
  {
    throw std::runtime_error("vqe built without CUDA/HIP support; GPU backend unavailable");
  }
#endif

  return run_default_vqe_backend<backend::statevector_cpu>(ansatz, pauli_terms, options);
}
