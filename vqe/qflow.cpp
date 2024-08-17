#include "vqeBackendManager.hpp"
#include "utils.hpp"
#include "vqe_state.hpp"
#include <unordered_map>
#include <sstream>
#define UNDERLINE "\033[4m"

#define CLOSEUNDERLINE "\033[0m"

int show_help() {
  std::cout << "NWQ-VQE QFlow Options" << std::endl;
  std::cout << UNDERLINE << "REQUIRED" << CLOSEUNDERLINE << std::endl;
  std::cout << "--hamiltonian, -f     Path to the input Hamiltonian file (formatted as a sum of Fermionic operators, see examples)" << std::endl;
  std::cout << "--nparticles, -n      Number of electrons in molecule" << std::endl;
  std::cout << UNDERLINE << "OPTIONAL" << CLOSEUNDERLINE << std::endl;
  std::cout << "--backend, -b         Simulation backend. Defaults to CPU" << std::endl;
  std::cout << "--list-backends, -l   List available backends and exit." << std::endl;
  std::cout << "--xacc                Use XACC indexing scheme, otherwise uses DUCC scheme." << std::endl;
  std::cout << "--seed                Random seed for initial point and empirical gradient estimation. Defaults to time(NULL)" << std::endl;
  std::cout << "--verbose             Print optimizer information on each iteration. Defaults to false" << std::endl;
  std::cout << UNDERLINE << "OPTIONAL (Global Minimizer)" << CLOSEUNDERLINE << std::endl;
  std::cout << "--optimizer           NLOpt optimizer name. Defaults to LN_COBYLA" << std::endl;
  std::cout << "--config              Path to config file for NLOpt optimizer parameters" << std::endl;
  std::cout << "--reltol              Relative tolerance termination criterion. Defaults to -1 (off)" << std::endl;
  std::cout << "--lbound              Optimizer lower bound. Defaults to -2" << std::endl;
  std::cout << "--ubound              Optimizer upper bound. Defaults to 2" << std::endl;
  std::cout << "--abstol              Relative tolerance termination criterion. Defaults to -1 (off)" << std::endl;
  std::cout << "--maxeval             Maximum number of function evaluations for optimizer. Defaults to 200" << std::endl;
  std::cout << "--maxtime             Maximum optimizer time (seconds). Defaults to -1.0 (off)" << std::endl;
  std::cout << "--stopval             Cutoff function value for optimizer. Defaults to -MAXFLOAT (off)" << std::endl;
  std::cout << UNDERLINE << "OPTIONAL (Local Gradient Follower)" << CLOSEUNDERLINE << std::endl;
  std::cout << "--local               Use local gradient follower pipeline." << std::endl;
  std::cout << "-g,--grad-samples     SPSA gradient samples." << std::endl;
  std::cout << "--delta               Perturbation magnitude for SPSA" << std::endl;
  std::cout << "--eta                 Gradient descent step size." << std::endl;
  return 1;
}

using json = nlohmann::json;
int parse_args(int argc, char** argv,
               VQEBackendManager& manager,
                std::string& hamilfile,
                std::string& backend,
                NWQSim::IdxType& n_particles,
                nlopt::algorithm& algo,
                NWQSim::VQE::OptimizerSettings& settings,
                int& n_trials,
                bool& use_xacc,
                bool& local,
                bool& verbose,
                unsigned& seed,
                double& delta,
                double& eta) {
  std::string config_file = "";
  std::string algorithm_name = "LN_COBYLA";
  hamilfile = "";
  backend = "CPU";
  n_particles = -1;
  n_trials = 1;
  settings.max_evals = 200;
  seed = time(0);
  delta = 1e-4;
  eta = 1e-3;
  use_xacc = false;
  local = false;
  verbose = false;
  settings.lbound = -2;
  settings.ubound = 2;
  for (size_t i = 1; i < argc; i++) {
    std::string argname = argv[i];
    if (argname == "-h" || argname == "--help") {
      return show_help();
    } if (argname == "-l" || argname == "--list-backends") {
      manager.print_available_backends();
      return 1;
    } else
    if (argname == "-b" || argname == "--backend") {
      backend = argv[++i];
      continue;
    } else
    if (argname == "-f" || argname == "--hamiltonian") {
      hamilfile = argv[++i];
      continue;
    } else 
    if (argname == "-n" || argname == "--nparticles") {
      n_particles = std::atoll(argv[++i]);
    } else 
    if (argname == "--local") {
      local = true;
    } else
    if (argname == "-g" || argname == "--grad-samples") {
      n_trials = std::atoll(argv[++i]);
    } else 
    if (argname == "--seed") {
      seed = (unsigned)std::atoi(argv[++i]);
    } else
    if (argname == "--xacc") {
      use_xacc = true;
    } else
    if (argname == "--verbose") {
      verbose = true;
    } else if (argname == "--delta") {
      delta = std::atof(argv[++i]);
    } else if (argname == "--eta") {
      eta = std::atof(argv[++i]);
    } else 
    if (argname == "--config") {
      config_file = argv[++i];
    } else 
    if (argname == "--optimizer") {
      algorithm_name = argv[++i];
    } else 
    if (argname == "--reltol") {
      settings.rel_tol = std::atof(argv[++i]);
    } else 
    if (argname == "--abstol") {
      settings.abs_tol = std::atof(argv[++i]);
    } else 
    if (argname == "--ubound") {
      settings.ubound = std::atof(argv[++i]);
    } else 
    if (argname == "--lbound") {
      settings.lbound = std::atof(argv[++i]);
    } else 
    if (argname == "--maxeval") {
      settings.max_evals = std::atoll(argv[++i]);
    } else if (argname == "--stopval") {
      settings.stop_val = std::atof(argv[++i]);
    } else if (argname == "--maxtime") {
      settings.max_time = std::atof(argv[++i]);
    } else {
      fprintf(stderr, "\033[91mERROR:\033[0m Unrecognized option %s, type -h or --help for a list of configurable parameters\n", argv[i]);
      return show_help();
    }
  }
  if (hamilfile == "") {
      fprintf(stderr, "\033[91mERROR:\033[0m Must pass a Hamiltonian file (--hamiltonian, -f)\n");

      return show_help();
  }
  if (n_particles == -1) {
      fprintf(stderr, "\033[91mERROR:\033[0m Must pass a particle count (--nparticles, -n)\n");
      return show_help();
  }
  algo = (nlopt::algorithm)nlopt_algorithm_from_string(algorithm_name.c_str());
  if (config_file != "") {
    std::ifstream f(config_file);
    json data = json::parse(f); 
    for (json::iterator it = data.begin(); it != data.end(); ++it) {
      settings.parameter_map[it.key()] = it.value().get<NWQSim::ValType>();
    }
  }
  return 0;
}

// Callback function, requires signature (void*) (const std::vector<NWQSim::ValType>&, NWQSim::ValType, NWQSim::IdxType)
void carriage_return_callback_function(const std::vector<NWQSim::ValType>& x, NWQSim::ValType fval, NWQSim::IdxType iteration) {
  printf("\33[2K\rEvaluation %lld, fval = %f", iteration, fval);fflush(stdout);
}

// Callback function, requires signature (void*) (const std::vector<NWQSim::ValType>&, NWQSim::ValType, NWQSim::IdxType)
void callback_function(const std::vector<NWQSim::ValType>& x, NWQSim::ValType fval, NWQSim::IdxType iteration) {
  printf("\33[2KEvaluation %lld, fval = %f\n", iteration, fval);fflush(stdout);
}

// Callback function, requires signature (void*) (const std::vector<NWQSim::ValType>&, NWQSim::ValType, NWQSim::IdxType)
void silent_callback_function(const std::vector<NWQSim::ValType>& x, NWQSim::ValType fval, NWQSim::IdxType iteration) {
  
}

void optimize_ansatz(const VQEBackendManager& manager,
                     const std::string& backend,
                     std::shared_ptr<NWQSim::VQE::Hamiltonian> hamil,
                     std::shared_ptr<NWQSim::VQE::Ansatz> ansatz,
                     NWQSim::VQE::OptimizerSettings& settings,
                     nlopt::algorithm& algo,
                     unsigned& seed,
                     int num_trials,
                     std::vector<double>& params,
                     bool local,
                     bool verbose,
                     double delta,
                     double eta) {
  NWQSim::VQE::Callback callback = verbose ? carriage_return_callback_function: silent_callback_function;
  std::shared_ptr<NWQSim::VQE::VQEState> state = manager.create_vqe_solver(backend, ansatz, hamil, algo, callback, seed, settings);
  if (!verbose) {
    NWQSim::Config::PRINT_SIM_TRACE = false;
  }  
  std::uniform_real_distribution<double> initdist(0, 2 * PI);
  std::mt19937_64 random_engine (seed);
  params.resize(ansatz->numParams());
  std::fill(params.begin(), params.end(), 0);
  double initial_ene, final_ene;
  long long num_iterations = 0;
  std::vector<std::pair<std::string, double> > param_tuple;
  if (local) {
  param_tuple = state->follow_fixed_gradient(params, 
                                              initial_ene, 
                                              final_ene, 
                                              num_iterations, 
                                              delta, 
                                              eta, 
                                              num_trials);

  } else {
    state->optimize(params, final_ene);
    param_tuple = ansatz->getFermionicOperatorParameters();
    num_iterations = state->get_iteration();
  }

  std::ostringstream strstream;
  manager.safe_print("\nFinished in %llu iterations. Initial Energy %f, Final Energy %f\nPrinting excitation amplitudes:\n", num_iterations, initial_ene, final_ene);
  for (auto& i: param_tuple) {
    strstream << i.first << ": " << i.second << std::endl;
  }
  manager.safe_print("%s", strstream.str().c_str());
}


int main(int argc, char** argv) {
  VQEBackendManager manager;
  std::string hamil_path, backend;
  NWQSim::IdxType n_part;
  NWQSim::VQE::OptimizerSettings settings;
  nlopt::algorithm algo;
  double delta;
  double eta;
  unsigned seed;
  bool use_xacc, local, verbose;
  int n_trials;
  if (parse_args(argc, argv, manager, hamil_path, backend, n_part, algo, settings, n_trials, use_xacc, local, verbose, seed, delta, eta)) {
    return 1;
  }
#ifdef MPI_ENABLED
  int i_proc;
  if (backend == "MPI" || backend == "NVGPU_MPI")
  {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &i_proc);
  }
#endif
  manager.safe_print("Reading Hamiltonian...\n");
  std::shared_ptr<NWQSim::VQE::Hamiltonian> hamil = std::make_shared<NWQSim::VQE::Hamiltonian>(hamil_path, n_part, use_xacc);
  manager.safe_print("Constructing UCCSD Ansatz...\n");

  std::shared_ptr<NWQSim::VQE::Ansatz> ansatz = std::make_shared<NWQSim::VQE::UCCSD>(
    hamil->getEnv(),
    NWQSim::VQE::getJordanWignerTransform,
    1
  );
  std::vector<double> params;
  manager.safe_print("Beginning VQE loop...\n");
  optimize_ansatz(manager, backend, hamil, ansatz, settings, algo, seed, n_trials, params, local, verbose, delta, eta);
#ifdef MPI_ENABLED
  if (backend == "MPI" || backend == "NVGPU_MPI")
  {
    MPI_Finalize();
  }
#endif
  return 0;
}