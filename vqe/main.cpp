#include <algorithm>
#include <chrono>
#include <cctype>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <nlopt.h>
#include <nlopt.hpp>

#ifdef VQE_ENABLE_MPI
#include <mpi.h>
#endif

#include "ansatz/ansatz.hpp"
#include "execution/adapt_runner.hpp"
#include "execution/vqe_runner.hpp"
#include "hamiltonian_parser.hpp"
#include "jw_transform.hpp"

#include <private/nlohmann/json.hpp>

namespace
{
  using json = nlohmann::json;

  constexpr double kPi = 3.14159265358979323846;

  std::string nlopt_status_to_string(nlopt::result status)
  {
    switch (status)
    {
    case nlopt::FAILURE:
      return "FAILURE";
    case nlopt::INVALID_ARGS:
      return "INVALID_ARGS";
    case nlopt::OUT_OF_MEMORY:
      return "OUT_OF_MEMORY";
    case nlopt::ROUNDOFF_LIMITED:
      return "ROUNDOFF_LIMITED";
    case nlopt::FORCED_STOP:
      return "FORCED_STOP";
    case nlopt::SUCCESS:
      return "SUCCESS";
    case nlopt::STOPVAL_REACHED:
      return "STOPVAL_REACHED";
    case nlopt::FTOL_REACHED:
      return "FTOL_REACHED";
    case nlopt::XTOL_REACHED:
      return "XTOL_REACHED";
    case nlopt::MAXEVAL_REACHED:
      return "MAXEVAL_REACHED";
    case nlopt::MAXTIME_REACHED:
      return "MAXTIME_REACHED";
    default:
      return "UNKNOWN";
    }
  }

  std::string nlopt_algorithm_to_string(nlopt::algorithm algo)
  {
    try
    {
      nlopt::opt opt(algo, 1);
      return opt.get_algorithm_name();
    }
    catch (const std::exception &)
    {
      return "unknown";
    }
  }

  struct cli_config
  {
    bool show_help = false;
    bool list_backends = false;
    bool disable_fusion = false;
    std::optional<int> requested_threads;

    std::string hamiltonian_path;
    std::size_t n_particles = 0;
    bool have_particles = false;

    std::string backend = "CPU";
    std::string optimizer_config_path;
    std::string init_params_input;  // Initial parameters input (single value or file)
    bool save_params = false;  // Save optimized parameters to file

    vqe::vqe_options options;
    std::unordered_map<std::string, double> optimizer_params;
  };

  std::string to_upper(std::string value)
  {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c)
                   { return static_cast<char>(std::toupper(c)); });
    return value;
  }

  void print_help()
  {
    std::cout << "NWQ-VQE Options\n"
              << "INFORMATIONAL\n"
              << "  -h, --help            Show help menu.\n"
              << "  -l, --list-backends   List available backends and exit.\n"
              << "REQUIRED\n"
              << "  -f, --hamiltonian     Path to the input Hamiltonian file (sum of Fermionic operators).\n"
              << "  -p, -n, --nparticles  Number of electrons in molecule.\n"
              << "OPTIONAL (Hamiltonian, Ansatz and Backend)\n"
              << "  --ducc                Use DUCC indexing scheme (default).\n"
              << "  --xacc                Use XACC/Qiskit indexing scheme.\n"
              << "  --sym, --symm         UCCSD symmetry level (0->none, 1->spin, 2->orbital, 3->full) (default to 3).\n"
              << "  -b, --backend         Simulation backend (CPU|GPU). Defaults to CPU.\n"
              << "  --seed                Random seed for reproducibility.\n"
              << "OPTIONAL (Global Minimizer)\n"
              << "  -v, --verbose         Print additional progress information.\n"
              << "  -o, --optimizer       NLopt optimizer name (e.g. LN_COBYLA, LN_BOBYQA).\n"
              << "  --opt-config          Path to JSON file with optimizer parameter overrides.\n"
              << "  -lb, --lbound         Optimizer lower bound (default -π).\n"
              << "  -ub, --ubound         Optimizer upper bound (default π).\n"
              << "  --reltol              Relative tolerance termination criterion.\n"
              << "  --abstol              Absolute tolerance termination criterion.\n"
              << "  --stopval             Objective stop value.\n"
              << "  --maxeval             Maximum number of function evaluations (default 100).\n"
              << "  --maxtime             Maximum optimizer time (seconds).\n"
              << "  --spsa                Use SPSA gradient estimation (2 evals) instead of forward difference (N+1 evals).\n"
              << "  --init-params         [For VQE only] Initial parameters: single value (repeat for all) or file with comma-separated values.\n"
              << "  --save-params         [For VQE only] Save optimized parameters to {hamiltonian_path}-vqe_params.txt.\n"
              << "OPTIONAL (ADAPT-VQE)\n"
              << "  --adapt               Run ADAPT-VQE instead of standard VQE.\n"
              << "  -ag, --adapt-gradtol  Gradient norm tolerance (default 1e-3).\n"
              << "  -af, --adapt-fvaltol  Energy change tolerance (default disabled).\n"
              << "  -am, --adapt-maxeval  Maximum ADAPT iterations (default 50).\n"
              << "  -as, --adapt-saveint  Save parameters every N iterations (default 0 = no saving) to the same path as the Hamiltonian file.\n"
              << "SIMULATOR OPTIONS\n"
              << "  --num_threads         Specify number of threads (ignored in current backend).\n"
              << "  --disable_fusion      Disable gate fusion (ignored in current backend).\n";
  }

  void list_backends()
  {
    std::cout << "Available backends:\n  CPU";
#ifdef VQE_ENABLE_CUDA
    std::cout << "\n  GPU";
#endif
    std::cout << "\n";
  }

  void apply_optimizer_parameters(vqe::vqe_options &opts,
                                  const std::unordered_map<std::string, double> &params)
  {
    for (const auto &entry : params)
    {
      opts.algorithm_parameters[entry.first] = entry.second;
      opts.adapt_algorithm_parameters[entry.first] = entry.second;
    }
  }

  bool parse_opt_config(const std::string &path,
                        std::unordered_map<std::string, double> &out_map,
                        std::string &error)
  {
    std::ifstream in(path);
    if (!in.is_open())
    {
      error = "Failed to open optimizer config file: " + path;
      return false;
    }
    json doc;
    try
    {
      in >> doc;
    }
    catch (const std::exception &ex)
    {
      error = std::string("Failed to parse optimizer config: ") + ex.what();
      return false;
    }
    if (!doc.is_object())
    {
      error = "Optimizer config must be a JSON object mapping parameter names to numeric values";
      return false;
    }
    for (auto it = doc.begin(); it != doc.end(); ++it)
    {
      if (!it.value().is_number())
      {
        error = "Optimizer config entry '" + it.key() + "' is not numeric";
        return false;
      }
      out_map[it.key()] = it.value().get<double>();
    }
    return true;
  }

  bool parse_backend(const std::string &value, cli_config &config, std::string &error)
  {
    const std::string upper = to_upper(value);
    config.backend = upper;
    if (upper == "CPU")
    {
      config.options.use_gpu = false;
      return true;
    }
    if (upper == "GPU" || upper == "NVGPU" || upper == "NVGPU_MPI")
    {
      config.options.use_gpu = true;
      return true;
    }
    if (upper == "MPI")
    {
      error = "MPI backend is not supported in the new implementation";
      return false;
    }
    error = "Unknown backend: " + value;
    return false;
  }

  bool parse_optimizer_name(const std::string &name, vqe::vqe_options &opts, std::string &error)
  {
    const std::string upper = to_upper(name);
    nlopt_algorithm algo = nlopt_algorithm_from_string(upper.c_str());
    if (algo < 0 || algo >= NLOPT_NUM_ALGORITHMS)
    {
      error = "Unknown NLopt optimizer: " + name;
      return false;
    }
    opts.optimizer = static_cast<nlopt::algorithm>(algo);
    opts.adapt_optimizer = opts.optimizer;
    return true;
  }

  // Parse initial parameters: either a single value or a file with comma-separated values
  bool parse_initial_parameters(const std::string &input, std::size_t expected_size,
                                std::vector<double> &params, std::string &error)
  {
    params.clear();

    // Try to parse as a single double value
    try
    {
      double value = std::stod(input);
      // It's a single value - repeat for all parameters
      if (expected_size > 0)
      {
        params.assign(expected_size, value);
      }
      else
      {
        params.push_back(value);  // Will be repeated later when size is known
      }
      return true;
    }
    catch (...)
    {
      // Not a single value, try to read as a file
    }

    // Try to read from file
    std::ifstream file(input);
    if (!file.is_open())
    {
      error = "Cannot parse as number or open file: " + input;
      return false;
    }

    std::string line;
    while (std::getline(file, line))
    {
      // Split by comma
      std::stringstream ss(line);
      std::string token;
      while (std::getline(ss, token, ','))
      {
        // Trim whitespace
        token.erase(0, token.find_first_not_of(" \t\r\n"));
        token.erase(token.find_last_not_of(" \t\r\n") + 1);

        if (token.empty()) continue;

        try
        {
          params.push_back(std::stod(token));
        }
        catch (...)
        {
          error = "Invalid number in file: " + token;
          return false;
        }
      }
    }

    if (params.empty())
    {
      error = "No parameters found in file: " + input;
      return false;
    }

    // If file has fewer parameters than expected, fill with zeros
    if (expected_size > 0 && params.size() < expected_size)
    {
      std::cerr << "[warning] Parameter file has " << params.size() << " values, expected "
                << expected_size << ". Filling remaining with zeros." << std::endl;
      params.resize(expected_size, 0.0);
    }
    else if (expected_size > 0 && params.size() > expected_size)
    {
      std::cerr << "[warning] Parameter file has " << params.size() << " values, expected "
                << expected_size << ". Using first " << expected_size << " values." << std::endl;
      params.resize(expected_size);
    }

    return true;
  }

  void set_bounds(vqe::vqe_options &opts, double lower, double upper)
  {
    opts.lower_bound = lower;
    opts.upper_bound = upper;
    opts.adapt_lower_bound = lower;
    opts.adapt_upper_bound = upper;
  }

  void set_max_evaluations(vqe::vqe_options &opts, std::size_t value)
  {
    opts.max_evaluations = value;
    opts.adapt_max_evaluations = value;
  }

  void set_max_time(vqe::vqe_options &opts, double seconds)
  {
    opts.max_time = seconds;
    opts.adapt_max_time = seconds;
  }

  void set_relative_tolerance(vqe::vqe_options &opts, double tol)
  {
    opts.relative_tolerance = tol;
    opts.adapt_relative_tolerance = tol;
  }

  void set_absolute_tolerance(vqe::vqe_options &opts, double tol)
  {
    opts.absolute_tolerance = tol;
    opts.adapt_absolute_tolerance = tol;
  }

  void set_stop_value(vqe::vqe_options &opts, double value)
  {
    opts.stop_value = value;
    opts.adapt_stop_value = value;
  }

  bool parse_args(int argc, char **argv, cli_config &config, std::string &error)
  {
    config.options = vqe::vqe_options{};
    config.options.symmetry_level = 3; // MZ: changed to 3
    config.options.trotter_steps = 1;
    config.options.use_xacc_indexing = true;
    set_bounds(config.options, -1.0 * kPi, 1.0 * kPi);  // MZ: changed to [-pi, pi]
    set_relative_tolerance(config.options, -1.0);
    set_absolute_tolerance(config.options, -1.0);
    set_stop_value(config.options, -std::numeric_limits<double>::infinity());
    set_max_evaluations(config.options, 100);
    set_max_time(config.options, -1.0);
    config.options.optimizer = nlopt::LN_COBYLA;
    config.options.adapt_optimizer = nlopt::LN_COBYLA;
    config.options.adapt_gradient_tolerance = 1e-3;
    config.options.adapt_energy_tolerance = -1.0;
    config.options.adapt_max_iterations = 50;

    for (int i = 1; i < argc; ++i)
    {
      std::string arg = argv[i];
      if (arg == "-h" || arg == "--help")
      {
        config.show_help = true;
        return true;
      }
      if (arg == "-l" || arg == "--list-backends")
      {
        config.list_backends = true;
        continue;
      }
      if (arg == "-f" || arg == "--hamiltonian")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --hamiltonian";
          return false;
        }
        config.hamiltonian_path = argv[++i];
        continue;
      }
      if (arg == "-p" || arg == "-n" || arg == "--nparticles")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --nparticles";
          return false;
        }
        config.n_particles = static_cast<std::size_t>(std::stoull(argv[++i]));
        config.have_particles = true;
        continue;
      }
      if (arg == "--ducc")
      {
        config.options.use_xacc_indexing = false;
        continue;
      }
      if (arg == "--xacc")
      {
        config.options.use_xacc_indexing = true;
        continue;
      }
      if (arg == "--sym" || arg == "--symm")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --sym";
          return false;
        }
        std::size_t level = static_cast<std::size_t>(std::stoull(argv[++i]));
        config.options.symmetry_level = level;
        continue;
      }
      if (arg == "-b" || arg == "--backend")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --backend";
          return false;
        }
        if (!parse_backend(argv[++i], config, error))
        {
          return false;
        }
        continue;
      }
      if (arg == "--seed")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --seed";
          return false;
        }
        config.options.random_seed = static_cast<unsigned>(std::stoul(argv[++i]));
        continue;
      }
      if (arg == "-v" || arg == "--verbose")
      {
        config.options.verbose = true;
        continue;
      }
      if (arg == "--opt-config")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --opt-config";
          return false;
        }
        config.optimizer_config_path = argv[++i];
        continue;
      }
      if (arg == "-o" || arg == "--optimizer")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --optimizer";
          return false;
        }
        if (!parse_optimizer_name(argv[++i], config.options, error))
        {
          return false;
        }
        continue;
      }
      if (arg == "-lb" || arg == "--lbound")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --lbound";
          return false;
        }
        const double lower = std::stod(argv[++i]);
        set_bounds(config.options, lower, config.options.upper_bound);
        continue;
      }
      if (arg == "-ub" || arg == "--ubound")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --ubound";
          return false;
        }
        const double upper = std::stod(argv[++i]);
        set_bounds(config.options, config.options.lower_bound, upper);
        continue;
      }
      if (arg == "--reltol")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --reltol";
          return false;
        }
        set_relative_tolerance(config.options, std::stod(argv[++i]));
        continue;
      }
      if (arg == "--abstol")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --abstol";
          return false;
        }
        set_absolute_tolerance(config.options, std::stod(argv[++i]));
        continue;
      }
      if (arg == "--maxeval")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --maxeval";
          return false;
        }
        set_max_evaluations(config.options, static_cast<std::size_t>(std::stoull(argv[++i])));
        continue;
      }
      if (arg == "--maxtime")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --maxtime";
          return false;
        }
        set_max_time(config.options, std::stod(argv[++i]));
        continue;
      }
      if (arg == "--stopval")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --stopval";
          return false;
        }
        set_stop_value(config.options, std::stod(argv[++i]));
        continue;
      }
      if (arg == "--num_threads")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --num_threads";
          return false;
        }
        config.requested_threads = static_cast<int>(std::stoi(argv[++i]));
        continue;
      }
      if (arg == "--disable_fusion")
      {
        config.disable_fusion = true;
        continue;
      }
      if (arg == "--spsa")
      {
        config.options.use_spsa_gradient = true;
        continue;
      }
      if (arg == "--init-params")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --init-params";
          return false;
        }
        config.init_params_input = argv[++i];
        continue;
      }
      if (arg == "--save-params")
      {
        config.save_params = true;
        continue;
      }
      if (arg == "--adapt")
      {
        config.options.mode = vqe::run_mode::adapt;
        continue;
      }
      if (arg == "-ag" || arg == "--adapt-gradtol")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --adapt-gradtol";
          return false;
        }
        config.options.adapt_gradient_tolerance = std::stod(argv[++i]);
        continue;
      }
      if (arg == "-af" || arg == "--adapt-fvaltol")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --adapt-fvaltol";
          return false;
        }
        config.options.adapt_energy_tolerance = std::stod(argv[++i]);
        continue;
      }
      if (arg == "-am" || arg == "--adapt-maxeval")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --adapt-maxeval";
          return false;
        }
        config.options.adapt_max_iterations = static_cast<std::size_t>(std::stoull(argv[++i]));
        continue;
      }
      if (arg == "-as" || arg == "--adapt-saveint")
      {
        if (i + 1 >= argc)
        {
          error = "Missing value for --adapt-saveint";
          return false;
        }
        config.options.adapt_save_interval = static_cast<std::size_t>(std::stoull(argv[++i]));
        continue;
      }
      error = "Unrecognized option: " + arg;
      return false;
    }

    if (!config.optimizer_config_path.empty())
    {
      if (!parse_opt_config(config.optimizer_config_path, config.optimizer_params, error))
      {
        return false;
      }
    }

    apply_optimizer_parameters(config.options, config.optimizer_params);
    return true;
  }

  struct gate_stats
  {
    std::size_t total = 0;
    std::size_t single_qubit = 0;
    std::size_t two_qubit = 0;
    std::size_t pauli_rotations = 0;
    std::size_t x_gates = 0;
  };

  gate_stats analyze_gates(const vqe::circuit &circ)
  {
    gate_stats stats;
    stats.total = circ.gates().size();
    for (const auto &gate : circ.gates())
    {
      switch (gate.kind)
      {
      case vqe::gate_kind::x:
        ++stats.single_qubit;
        ++stats.x_gates;
        break;
      case vqe::gate_kind::h:
      case vqe::gate_kind::rz:
        ++stats.single_qubit;
        break;
      case vqe::gate_kind::cnot:
        ++stats.two_qubit;
        break;
      case vqe::gate_kind::pauli_rotation:
        ++stats.pauli_rotations;
        break;
      }
    }
    return stats;
  }

  std::vector<std::string> ordered_parameter_labels(const vqe::uccsd_ansatz &ansatz)
  {
    std::vector<std::string> labels;
    const auto &map = ansatz.excitation_parameter_map();
    for (const auto &excitation : ansatz.excitations())
    {
      if (map.find(excitation.label) != map.end() &&
          std::find(labels.begin(), labels.end(), excitation.label) == labels.end())
      {
        labels.push_back(excitation.label);
      }
    }
    return labels;
  }

  void print_parameter_summary(const vqe::uccsd_ansatz &ansatz,
                               const vqe::vqe_result &result,
                               std::size_t trotter_steps)
  {
    const auto &map = ansatz.excitation_parameter_map();
    auto labels = ordered_parameter_labels(ansatz);
    const std::size_t unique_params = ansatz.unique_parameter_count();
    const std::size_t per_step = unique_params == 0 ? result.parameters.size() : unique_params;
    if (labels.empty() && per_step == result.parameters.size())
    {
      labels.resize(per_step);
      for (std::size_t i = 0; i < per_step; ++i)
      {
        labels[i] = "param_" + std::to_string(i);
      }
    }
    std::cout << "Final parameters:" << std::endl;
    std::cout << std::setprecision(16);
    for (std::size_t step = 0; step < std::max<std::size_t>(trotter_steps, 1); ++step)
    {
      if (trotter_steps > 1)
      {
        std::cout << "  [Trotter step " << (step + 1) << "]" << std::endl;
      }
      for (const auto &label : labels)
      {
        const auto it = map.find(label);
        if (it == map.end())
        {
          continue;
        }
        const std::size_t base_index = it->second;
        const std::size_t parameter_index = base_index + step * per_step;
        if (parameter_index >= result.parameters.size())
        {
          continue;
        }
        std::cout << "    " << label << " :: " << result.parameters[parameter_index] << std::endl;
      }
    }
    std::cout << std::setprecision(6);
  }

  std::string format_duration(double seconds)
  {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4) << seconds << " seconds";
    return oss.str();
  }

  void emit_warnings(const cli_config &config)
  {
    if (config.disable_fusion)
    {
      std::cerr << "[warning] Gate fusion cannot be disabled in the new backend." << std::endl;
    }
    if (config.requested_threads.has_value())
    {
      std::cerr << "[warning] CPU backend does not currently support explicit thread configuration." << std::endl;
    }
  }








//-------------------------------------------------------------- run_vqe_mode() --------------------------------------------------------------//

  int run_vqe_mode(const cli_config &config)
  {
    const auto &opts = config.options;

    if (opts.verbose)
    {
      std::cout << "[vqe] Configuration summary:" << std::endl;
      std::cout << "  Hamiltonian path      : " << config.hamiltonian_path << std::endl;
      std::cout << "  Particles (electrons) : " << config.n_particles << std::endl;
      std::cout << "  Backend               : " << config.backend << std::endl;
      std::cout << "  Trotter steps         : " << opts.trotter_steps << std::endl;
      std::cout << "  Symmetry level        : " << opts.symmetry_level << std::endl;
      std::cout << "  Optimizer             : " << nlopt_algorithm_to_string(opts.optimizer);
      if (opts.use_spsa_gradient)
      {
        std::cout << " (SPSA gradient)";
      }
      std::cout << std::endl;
      std::cout << "  Parameter bounds      : [" << opts.lower_bound << ", " << opts.upper_bound << "]" << std::endl;
      if (opts.max_evaluations > 0)
      {
        std::cout << "  Max evaluations       : " << opts.max_evaluations << std::endl;
      }
      if (opts.max_time > 0)
      {
        std::cout << "  Max time (s)          : " << opts.max_time << std::endl;
      }
      if (opts.random_seed.has_value())
      {
        std::cout << "  Random seed           : " << *opts.random_seed << std::endl;
      }
    }

    std::cout << "Reading Hamiltonian..." << std::endl;
    const auto data = vqe::read_hamiltonian_file(config.hamiltonian_path);
    const auto pauli_terms = vqe::jordan_wigner_transform(data);
    if (opts.verbose)
    {
      std::cout << "Hamiltonian has " << pauli_terms.size() << " Pauli terms." << std::endl;
    }

    vqe::molecular_environment env;
    env.n_spatial = data.num_qubits() / 2;
    env.n_electrons = config.n_particles;
    env.xacc_indexing = opts.use_xacc_indexing;
    env.constant_energy = data.constant.real();

    if (env.n_spatial == 0)
    {
      throw std::runtime_error("Hamiltonian does not define any spatial orbitals");
    }

    vqe::uccsd_ansatz ansatz(env, opts.trotter_steps, opts.symmetry_level);
    ansatz.build();

    const std::size_t param_count = ansatz.get_circuit().parameters().size();

    if (opts.verbose)
    {
      std::cout << "Constructed ansatz with " << ansatz.excitations().size() << " generators and "
                << param_count << " parameters." << std::endl;
    }

    // Process initial parameters if specified
    vqe::vqe_options opts_with_params = opts;
    if (!config.init_params_input.empty())
    {
      std::string error;
      std::vector<double> init_params;

      if (!parse_initial_parameters(config.init_params_input, param_count, init_params, error))
      {
        std::cerr << "Error parsing initial parameters: " << error << std::endl;
        return EXIT_FAILURE;
      }

      opts_with_params.initial_parameters = init_params;
      if (opts.verbose)
      {
        std::cout << "Using custom initial parameters: " << init_params.size() << " values" << std::endl;
      }
    }

    auto start = std::chrono::steady_clock::now();
    auto result = vqe::run_default_vqe_with_ansatz(ansatz, pauli_terms, opts_with_params);
    auto stop = std::chrono::steady_clock::now();

    const double elapsed = std::chrono::duration<double>(stop - start).count();

    std::cout << "\n--------- Result Summary ---------" << std::endl;
    std::cout << "Method                 : VQE" << std::endl;
    std::cout << "Symmetry level         : " << opts.symmetry_level << std::endl;
    std::cout << "# Hamiltonian terms    : " << pauli_terms.size() << std::endl;
    std::cout << "# Evaluations          : " << result.evaluations << std::endl;
    std::cout << "Converged              : " << (result.converged ? "yes" : "no") << std::endl;
    std::cout << std::setprecision(16);
    std::cout << "Final objective value  : " << result.energy << std::endl;
    std::cout << "Initial objective value: " << result.initial_energy << std::endl;
    std::cout << "Objective delta        : " << result.energy_delta << std::endl;
    std::cout << std::setprecision(6);
    std::cout << "Evaluation time        : " << format_duration(elapsed) << std::endl;

    const auto stats = analyze_gates(ansatz.get_circuit());
    std::cout << "Operator stats         : " << ansatz.excitations().size() << " operators, "
              << ansatz.unique_parameter_count() << " unique parameters, "
              << ansatz.get_circuit().parameters().size() << " total parameters" << std::endl;
    std::cout << "Circuit stats          : " << stats.total << " gates (" << stats.single_qubit << " 1q, "
              << stats.two_qubit << " 2q, " << stats.pauli_rotations << " Pauli rotations)" << std::endl;

    print_parameter_summary(ansatz, result, opts.trotter_steps);

    if (opts.verbose)
    {
      std::cout << "\n[vqe] Optimization details:" << std::endl;
      std::cout << "  NLopt status          : " << nlopt_status_to_string(result.status) << std::endl;
      std::cout << "  Iterations recorded   : " << result.iterations << std::endl;
      std::cout << "  Evaluations           : " << result.evaluations << std::endl;
      if (result.evaluations > 0)
      {
        const double avg_apply = result.apply_time / static_cast<double>(result.evaluations);
        const double avg_expect = result.expectation_time / static_cast<double>(result.evaluations);
        std::cout << "  Total apply time (s)  : " << result.apply_time << " (avg " << avg_apply << ")" << std::endl;
        std::cout << "  Total expect time (s) : " << result.expectation_time << " (avg " << avg_expect << ")" << std::endl;
      }
      else
      {
        std::cout << "  Total apply time (s)  : " << result.apply_time << std::endl;
        std::cout << "  Total expect time (s) : " << result.expectation_time << std::endl;
      }
    }

    // Save optimized parameters to file if requested
    if (config.save_params)
    {
      std::string param_file_path = config.hamiltonian_path + "-vqe_params.txt";
      std::ofstream param_file(param_file_path);
      if (param_file.is_open())
      {
        for (std::size_t i = 0; i < result.parameters.size(); ++i)
        {
          if (i > 0) param_file << ", ";
          param_file << std::setprecision(16) << result.parameters[i];
        }
        param_file << std::endl;
        param_file.close();
        std::cout << "[vqe] Saved " << result.parameters.size()
                  << " optimized parameters to: " << param_file_path << std::endl;
      }
      else
      {
        std::cerr << "[vqe] Warning: Could not save parameters to: " << param_file_path << std::endl;
      }
    }

    return result.converged ? EXIT_SUCCESS : EXIT_FAILURE;
  }




//-------------------------------------------------------------- run_adapt_mode() --------------------------------------------------------------//

  int run_adapt_mode(const cli_config &config)
  {
    const auto &opts = config.options;

    if (opts.verbose)
    {
      std::cout << "[adapt] Configuration summary:" << std::endl;
      std::cout << "  Hamiltonian path      : " << config.hamiltonian_path << std::endl;
      std::cout << "  Particles (electrons) : " << config.n_particles << std::endl;
      std::cout << "  Backend               : " << config.backend << std::endl;
      std::cout << "  Symmetry level        : " << opts.symmetry_level << std::endl;
      std::cout << "  Max iterations        : " << opts.adapt_max_iterations << std::endl;
      std::cout << "  Gradient step         : " << opts.adapt_gradient_step << std::endl;
      std::cout << "  Gradient tolerance    : " << opts.adapt_gradient_tolerance << std::endl;
      if (opts.adapt_energy_tolerance > 0.0)
      {
        std::cout << "  Energy tolerance      : " << opts.adapt_energy_tolerance << std::endl;
      }
      std::cout << "  Optimizer             : " << nlopt_algorithm_to_string(opts.adapt_optimizer);
      if (opts.use_spsa_gradient)
      {
        std::cout << " (SPSA gradient)";
      }
      std::cout << std::endl;
      std::cout << "  Parameter bounds      : [" << opts.adapt_lower_bound << ", "
                << opts.adapt_upper_bound << "]" << std::endl;
      if (opts.adapt_max_evaluations > 0)
      {
        std::cout << "  Max evaluations       : " << opts.adapt_max_evaluations << std::endl;
      }
      if (opts.adapt_max_time > 0)
      {
        std::cout << "  Max time (s)          : " << opts.adapt_max_time << std::endl;
      }
      if (opts.random_seed.has_value())
      {
        std::cout << "  Random seed           : " << *opts.random_seed << std::endl;
      }
      std::cout << "Running ADAPT-VQE..." << std::endl;
    }

    auto start = std::chrono::steady_clock::now();
    auto result = vqe::run_adapt_vqe(config.hamiltonian_path,
                                     config.n_particles,
                                     opts);
    auto stop = std::chrono::steady_clock::now();
    const double elapsed = std::chrono::duration<double>(stop - start).count();

    std::cout << "\n--------- Result Summary ---------" << std::endl;
    std::cout << "Method                 : ADAPT-VQE" << std::endl;
    std::cout << "Symmetry level         : " << opts.symmetry_level << std::endl;
    std::cout << "# ADAPT iterations     : " << result.iterations << std::endl;
    std::cout << "Energy evaluations     : " << result.energy_evaluations << std::endl;
    std::cout << std::setprecision(16);
    std::cout << "Final objective value  : " << result.energy << std::endl;
    std::cout << "Initial objective value: " << result.initial_energy << std::endl;
    std::cout << "Objective delta        : " << (result.energy - result.initial_energy) << std::endl;
    std::cout << std::setprecision(6);
    std::cout << "Evaluation time        : " << format_duration(elapsed) << std::endl;
    std::cout << "Converged              : " << (result.converged ? "yes" : "no") << std::endl;

    std::cout << "\n[adapt] Optimization details:" << std::endl;
    std::cout << "  Iterations executed   : " << result.iterations << std::endl;
    std::cout << "  Energy evaluations    : " << result.energy_evaluations << std::endl;
    // std::cout << std::setprecision(16);
    // std::cout << "  Initial energy        : " << result.initial_energy << std::endl;
    // std::cout << "  Final energy          : " << result.energy << std::endl;
    // std::cout << "  ΔEnergy               : " << (result.energy - result.initial_energy) << std::endl;
    // std::cout << std::setprecision(6);
    std::cout << "  Selected operator cnt : " << result.selected_indices.size() << std::endl;
    std::cout << "  Parameter count       : " << result.parameters.size() << std::endl;


    return result.converged ? EXIT_SUCCESS : EXIT_FAILURE;
  }

} // namespace



//-------------------------------------------------------------- main() --------------------------------------------------------------//

int main(int argc, char **argv)
{
#ifdef VQE_ENABLE_MPI
  MPI_Init(&argc, &argv);
  //Disable printing for processes other than node-0:
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank != 0) 
  {
      // Redirect C stdio
      freopen("/dev/null", "w", stdout);
      freopen("/dev/null", "w", stderr);
      // Optional: also silence C++ streams
      std::cout.setstate(std::ios::failbit);
      std::cerr.setstate(std::ios::failbit);
  }
#endif

  cli_config config;
  std::string error;
  if (!parse_args(argc, argv, config, error))
  {
    std::cerr << "Error: " << error << "\n";
    print_help();
    return EXIT_FAILURE;
  }

  if (config.show_help)
  {
    print_help();
    return EXIT_SUCCESS;
  }

  if (config.list_backends)
  {
    list_backends();
    return EXIT_SUCCESS;
  }

  if (config.hamiltonian_path.empty())
  {
    std::cerr << "Error: Hamiltonian file not specified (--hamiltonian)." << std::endl;
    print_help();
    return EXIT_FAILURE;
  }

  if (!config.have_particles)
  {
    std::cerr << "Error: Number of particles not specified (--nparticles)." << std::endl;
    print_help();
    return EXIT_FAILURE;
  }

  emit_warnings(config);

  try
  {
    if (config.options.mode == vqe::run_mode::adapt)
    {
      return run_adapt_mode(config);
    }
    return run_vqe_mode(config);
  }
  catch (const std::exception &ex)
  {
    std::cerr << "error: " << ex.what() << std::endl;
    return EXIT_FAILURE;
  }
#ifdef VQE_ENABLE_MPI
  MPI_Finalize();
#endif

}
