#pragma once

#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include <nlopt.hpp>

namespace vqe
{

  enum class run_mode
  {
    standard,
    adapt
  };

  struct vqe_options
  {
    // General execution controls
    run_mode mode = run_mode::standard;
    bool verbose = false;
    bool use_gpu = false;
    bool use_xacc_indexing = true;
    std::optional<unsigned> random_seed{};

    // Shared physical configuration
    std::size_t symmetry_level = 0;
    std::size_t trotter_steps = 1;

    // Optimizer configuration for standard VQE
    nlopt::algorithm optimizer = nlopt::LN_COBYLA;
    double lower_bound = -2.0 * 3.14159265358979323846;
    double upper_bound = 2.0 * 3.14159265358979323846;
    std::size_t max_evaluations = 100;
    double relative_tolerance = -1.0;
    double absolute_tolerance = -1.0;
    double stop_value = -std::numeric_limits<double>::infinity();
    double max_time = -1.0;
    double iteration_improvement_tolerance = 1e-8;
    std::size_t status_interval = 25;
    std::unordered_map<std::string, double> algorithm_parameters;
    std::vector<double> initial_parameters;

    // ADAPT-VQE specific controls
    std::size_t adapt_max_iterations = 20;
    double adapt_gradient_step = 1e-2;
    double adapt_gradient_tolerance = 1e-3;
    double adapt_energy_tolerance = -1.0;
    bool adapt_log_memory = false;
    std::size_t adapt_save_interval = 0; // Save parameters every N iterations (0 = no saving)

    // Optimizer configuration for ADAPT inner VQE solves
    nlopt::algorithm adapt_optimizer = nlopt::LN_COBYLA;
    double adapt_lower_bound = -2.0 * 3.14159265358979323846;
    double adapt_upper_bound = 2.0 * 3.14159265358979323846;
    std::size_t adapt_max_evaluations = 100;
    double adapt_relative_tolerance = -1.0;
    double adapt_absolute_tolerance = -1.0;
    double adapt_stop_value = -std::numeric_limits<double>::infinity();
    double adapt_max_time = -1.0;
    double adapt_iteration_improvement_tolerance = 1e-8;
    std::size_t adapt_status_interval = 25;
    std::unordered_map<std::string, double> adapt_algorithm_parameters;
    std::vector<double> adapt_initial_parameters;
  };

} // namespace vqe
