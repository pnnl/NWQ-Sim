#include "vqeBackendManager.hpp"
#include "utils.hpp"
#include <unordered_map>
#include <sstream>
#define UNDERLINE "\033[4m"

#define CLOSEUNDERLINE "\033[0m"

int show_help() {
  std::cout << "NWQ-VQE QFlow Options" << std::endl;
  std::cout << UNDERLINE << "REQUIRED" << CLOSEUNDERLINE << std::endl;
  std::cout << "--hamiltonian, -f     Path to the input Hamiltonian file (formatted as a sum of Fermionic operators, see examples)" << std::endl;
  std::cout << "--nparticles, -n      Number of electrons in molecule" << std::endl;
  std::cout << "--backend, -b         Simulation backend. Defaults to CPU" << std::endl;
  std::cout << "--list-backends, -l   List available backends and exit." << std::endl;
  std::cout << UNDERLINE << "OPTIONAL" << CLOSEUNDERLINE << std::endl;
  std::cout << "--seed                Random seed for initial point and empirical gradient estimation. Defaults to time(NULL)" << std::endl;
  std::cout << "--xacc                Use XACC indexing scheme, otherwise uses DUCC scheme." << std::endl;
  std::cout << "--delta               Magnitude of SPSA perturbation. Defaults to 1e-4" << std::endl;
  std::cout << "--eta                 Gradient descent step size. Defaults to 1e-3" << std::endl;
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
    } 
    if (argname == "-f" || argname == "--hamiltonian") {
      hamilfile = argv[++i];
      continue;
    } else 
    if (argname == "-n" || argname == "--nparticles") {
      n_particles = std::atoll(argv[++i]);
    } else 
    if (argname == "-g" || argname == "--grad-samples") {
      n_trials = std::atoll(argv[++i]);
    } else 
    if (argname == "--seed") {
      seed = (unsigned)std::atoi(argv[++i]);
    } else
    if (argname == "--xacc") {
      use_xacc = true;
    } else if (argname == "--delta") {
      delta = std::atof(argv[++i]);
    } else if (argname == "--eta") {
      eta = std::atof(argv[++i]);
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


void optimize_ansatz(const VQEBackendManager& manager,
                     const std::string& backend,
                     std::shared_ptr<NWQSim::VQE::Hamiltonian> hamil,
                     std::shared_ptr<NWQSim::VQE::Ansatz> ansatz,
                     NWQSim::VQE::OptimizerSettings& settings,
                     nlopt::algorithm& algo,
                     unsigned& seed,
                     int num_trials,
                     std::vector<double>& params,
                     double delta,
                     double eta) {
  std::shared_ptr<NWQSim::VQE::VQEState> state = manager.create_vqe_solver(backend, ansatz, hamil, algo, callback_function, seed, settings);  
  std::uniform_real_distribution<double> initdist(0, 2 * PI);
  std::mt19937_64 random_engine (seed);
  params.resize(ansatz->numParams());
  std::generate(params.begin(), params.end(), 
      [&random_engine, &initdist] () {return initdist(random_engine);});
  double initial_ene, final_ene;
  long long num_iterations = 0;
  std::vector<std::pair<std::string, double> > param_tuple = state->follow_fixed_gradient(params, initial_ene, final_ene, num_iterations, delta, eta, num_trials);
  std::ostringstream strstream;
  manager.safe_print("Finished in %llu iterations. Initial Energy %f, Final Energy %f\nPrinting excitation amplitudes:\n", num_iterations, initial_ene, final_ene);
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
  bool use_xacc;
  int n_trials;
  if (parse_args(argc, argv, manager, hamil_path, backend, n_part, algo, settings, n_trials, use_xacc, seed, delta, eta)) {
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
  optimize_ansatz(manager, backend, hamil, ansatz, settings, algo, seed, n_trials, params, delta, eta);
#ifdef MPI_ENABLED
  if (backend == "MPI" || backend == "NVGPU_MPI")
  {
    MPI_Finalize();
  }
#endif
  return 0;
}