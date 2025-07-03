#include "include/nwqsim_qflow.hpp"
#include "include/utils.hpp"
#include "vqeBackendManager.hpp"
#include "utils.hpp"
#include "vqe_state.hpp"

#include "src/uccsdmin.cpp"
#include "src/singletgsd.cpp"

#include <sstream>
#include <chrono>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <regex>
#include <map>

// Callback function, requires signature (void*) (const std::vector<NWQSim::ValType>&, NWQSim::ValType, NWQSim::IdxType)
void silent_callback_function(const std::vector<NWQSim::ValType> &x, NWQSim::ValType fval, NWQSim::IdxType iteration) {}

std::string get_termination_reason_local(int result)
{
    static const std::map<int, std::string> reason_map_adapt = {
        {0, "Local gradient minimum is reached"},
        {9, "Local Gradient Follower is not run"}};
    auto it = reason_map_adapt.find(result);
    if (it != reason_map_adapt.end())
    {
        return it->second;
    }
    else
    {
        return "Unknown reason, code: " + std::to_string(result);
    }
}

double optimize_ansatz(const VQEBackendManager &manager,
                       const std::string &backend,
                       std::shared_ptr<NWQSim::VQE::Hamiltonian> hamil,
                       std::shared_ptr<NWQSim::VQE::Ansatz> ansatz,
                       NWQSim::VQE::OptimizerSettings &settings,
                       nlopt::algorithm &algo,
                       unsigned &seed,
                       int num_trials,
                       bool verbose,
                       double delta,
                       double eta)
{
    std::shared_ptr<NWQSim::VQE::VQEState> state = manager.create_vqe_solver(backend, ansatz, hamil, algo, silent_callback_function, seed, settings);

    std::vector<double> params;
    params.resize(ansatz->numParams());
    std::fill(params.begin(), params.end(), 0);

    double initial_ene, final_ene;
    long long num_iterations = 0;
    state->initialize();

    state->follow_fixed_gradient(params,
                                 initial_ene,
                                 final_ene,
                                 num_iterations,
                                 delta,
                                 eta,
                                 num_trials,
                                 verbose);

    return final_ene;
}

std::vector<std::pair<std::string, std::complex<double>>> parseHamiltonianFile(const std::string &filename)
{
    std::vector<std::pair<std::string, std::complex<double>>> result;
    std::ifstream file(filename);

    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string line;
    // Regex pattern to match (real, imag) followed by optional config string
    std::regex pattern(R"(\(\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*,\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\s*\)(.*)?)");

    while (std::getline(file, line))
    {
        // Skip empty lines and comment lines
        if (line.empty() || line[0] == 'c')
        {
            continue;
        }

        std::smatch match;
        if (std::regex_search(line, match, pattern))
        {
            if (match.size() >= 3)
            {
                // Parse the complex number
                double real_part = std::stod(match.str(1));
                double imag_part = std::stod(match.str(2));
                std::complex<double> coeff(real_part, imag_part);

                // Get the config string (everything after the closing parenthesis)
                std::string config = match.str(3);

                // Trim leading/trailing whitespace from config
                config.erase(0, config.find_first_not_of(" \t"));
                config.erase(config.find_last_not_of(" \t") + 1);

                result.emplace_back(config, coeff);
            }
        }
    }

    file.close();
    return result;
}

double qflow_nwqsim(const std::vector<std::pair<std::string, std::complex<double>>> &hamiltonian_ops, int n_part, std::string backend)
{
    VQEBackendManager manager;

    NWQSim::VQE::OptimizerSettings settings;
    settings.max_evals = 100;
    settings.lbound = -PI;
    settings.ubound = PI;

    nlopt::algorithm algo = (nlopt::algorithm)nlopt_algorithm_from_string("LN_COBYLA");
    double delta = 1e-4;
    double eta = 1e-3;
    unsigned seed = time(0);
    bool use_xacc = true;
    bool verbose = false;

    int n_trials = 1;

    std::shared_ptr<NWQSim::VQE::Hamiltonian> hamil = std::make_shared<NWQSim::VQE::Hamiltonian>(hamiltonian_ops, n_part, use_xacc, NWQSim::VQE::getJordanWignerTransform);

    std::shared_ptr<NWQSim::VQE::Ansatz> ansatz;

    ansatz = std::make_shared<NWQSim::VQE::UCCSDmin>(
        hamil->getEnv(),
        NWQSim::VQE::getJordanWignerTransform,
        1,
        0);

    ansatz->buildAnsatz();

    double final_energy = optimize_ansatz(manager, backend, hamil, ansatz, settings, algo, seed, n_trials, verbose, delta, eta);

    return final_energy;
}
