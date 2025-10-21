#include "nwqsim_qflow.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <nlopt.hpp>

#include "ansatz/ansatz.hpp"
#include "execution/vqe_runner.hpp"
#include "fermion_term.hpp"
#include "hamiltonian_parser.hpp"
#include "jw_transform.hpp"

namespace
{
  constexpr double kPi = 3.14159265358979323846;

  vqe::hamiltonian_data build_data_from_ops(const std::vector<std::pair<std::string, std::complex<double>>> &ops)
  {
    vqe::hamiltonian_data data;
    for (const auto &entry : ops)
    {
      const std::string &config = entry.first;
      const auto &coeff = entry.second;
      if (config.empty())
      {
        data.constant += coeff;
        continue;
      }

      std::istringstream stream(config);
      std::string token;
      vqe::fermion_term term;
      term.coefficient = coeff;
      while (stream >> token)
      {
        if (token == "+")
        {
          continue;
        }
        bool creation = false;
        if (!token.empty() && token.back() == '^')
        {
          creation = true;
          token.pop_back();
        }
        if (token.empty())
        {
          continue;
        }
        const std::size_t index = static_cast<std::size_t>(std::stoull(token));
        data.max_index = std::max(data.max_index, index);
        term.operators.push_back({index, creation ? vqe::fermion_op_kind::creation : vqe::fermion_op_kind::annihilation});
      }
      if (term.operators.empty())
      {
        data.constant += coeff;
      }
      else
      {
        data.terms.push_back(std::move(term));
      }
    }
    return data;
  }

  bool is_gpu_backend(const std::string &backend)
  {
    std::string value = backend;
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c)
                   { return static_cast<char>(std::toupper(c)); });
    return value == "GPU" || value == "NVGPU" || value == "NVGPU_MPI";
  }

  std::vector<std::string> collect_parameter_labels(const vqe::uccsd_ansatz &ansatz)
  {
    std::vector<std::string> labels;
    const auto &map = ansatz.excitation_parameter_map();
    labels.reserve(map.size());
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

  std::vector<int> parse_parameter_label(const std::string &label)
  {
    std::vector<int> indices;
    std::istringstream stream(label);
    std::string token;
    while (stream >> token)
    {
      if (token == "+")
      {
        continue;
      }
      if (!token.empty() && token.back() == '^')
      {
        token.pop_back();
      }
      if (token.empty())
      {
        continue;
      }
      indices.push_back(static_cast<int>(std::stoi(token)));
    }
    std::reverse(indices.begin(), indices.end());
    return indices;
  }

} // namespace

std::vector<std::pair<std::string, std::complex<double>>> parseHamiltonianFile(const std::string &filename)
{
  const auto data = vqe::read_hamiltonian_file(filename);
  std::vector<std::pair<std::string, std::complex<double>>> result;
  if (std::abs(data.constant.real()) > 0.0 || std::abs(data.constant.imag()) > 0.0)
  {
    result.emplace_back("", data.constant);
  }
  for (const auto &term : data.terms)
  {
    std::ostringstream config;
    bool first = true;
    for (const auto &op : term.operators)
    {
      if (!first)
      {
        config << ' ';
      }
      first = false;
      config << op.index;
      if (op.kind == vqe::fermion_op_kind::creation)
      {
        config << '^';
      }
    }
    result.emplace_back(config.str(), term.coefficient);
  }
  return result;
}

std::pair<double, std::vector<std::pair<std::vector<int>, double>>> qflow_nwqsim(
    const std::vector<std::pair<std::string, std::complex<double>>> &hamiltonian_ops,
    int n_part,
    std::string backend,
    std::optional<vqe::vqe_options> options_override)
  // To use custom initial parameters, pass a vqe::vqe_options with initial_parameters set.
  // Use trial and error with verbose=true to find the correct parameter count:
  // Step 1: Run with verbose to see expected parameter count
  // vqe::vqe_options test_opts;
  // test_opts.verbose = true;
  // qflow_nwqsim(hamiltonian_ops, n_part, backend, test_opts);  // Check output for parameter count
  // Step 2: Create your initial parameters with the correct size
  // vqe::vqe_options custom_opts;
  // custom_opts.initial_parameters = {0.1, 0.2, 0.3, ...};  // Must match parameter count from step 1
  // custom_opts.verbose = true;
  // return qflow_nwqsim(hamiltonian_ops, n_part, backend, custom_opts);
{
  if (n_part < 0)
  {
    throw std::invalid_argument("Number of particles must be non-negative");
  }

  const auto data = build_data_from_ops(hamiltonian_ops);
  if (data.num_qubits() % 2 != 0)
  {
    throw std::runtime_error("Hamiltonian maps to an odd number of qubits; expected even for spin orbitals");
  }

  const auto pauli_terms = vqe::jordan_wigner_transform(data);

  vqe::molecular_environment env;
  env.n_spatial = data.num_qubits() / 2;
  env.n_electrons = static_cast<std::size_t>(n_part);
  env.constant_energy = data.constant.real();

  vqe::vqe_options options;
  if (options_override.has_value())
  {
    options = *options_override;
  }
  else
  {
    options = vqe::vqe_options{};
    options.verbose = false; // Set to true for more information
    options.trotter_steps = 1;
    options.symmetry_level = 3;
    options.lower_bound = -kPi; // -pi
    options.upper_bound = kPi; // pi
    options.max_evaluations = 100; // Max number of optimization iterations
    options.relative_tolerance = -1.0; // -1 for no relative tolerance
    options.absolute_tolerance = 1e-8; // -1 for no absolute tolerance
    options.stop_value = -std::numeric_limits<double>::infinity();
    options.max_time = -1.0; // -1 for no max time
    options.optimizer = nlopt::LD_LBFGS; // Use a derivative-based optimizer to enable gradient computation
  }

  options.mode = vqe::run_mode::standard;
  options.use_gpu = is_gpu_backend(backend);
  env.xacc_indexing = options.use_xacc_indexing;

  if (options.trotter_steps == 0)
  {
    throw std::invalid_argument("trotter_steps must be positive for qflow_nwqsim");
  }

  vqe::uccsd_ansatz ansatz(env, options.trotter_steps, options.symmetry_level);
  ansatz.build();

  auto result = vqe::run_default_vqe_with_ansatz(ansatz, pauli_terms, options);

  std::vector<std::pair<std::vector<int>, double>> parameter_output;
  const auto labels = collect_parameter_labels(ansatz);
  const auto &map = ansatz.excitation_parameter_map();
  for (const auto &label : labels)
  {
    const auto it = map.find(label);
    if (it == map.end())
    {
      continue;
    }
    const std::size_t index = it->second;
    if (index >= result.parameters.size())
    {
      continue;
    }
    parameter_output.emplace_back(parse_parameter_label(label), result.parameters[index]);
  }

  return {result.energy, parameter_output};
}

// MZ: "local gradient" (SPSA gradient) is not implemnted after re-organized anyway
// std::string get_termination_reason_local(int result)
// {
//   static const std::unordered_map<int, std::string> kReasons = {
//       {0, "Local gradient minimum is reached"},
//       {9, "Local Gradient Follower is not run"}};
//   const auto it = kReasons.find(result);
//   if (it != kReasons.end())
//   {
//     return it->second;
//   }
//   return "Unknown reason, code: " + std::to_string(result);
// }
