#include "vqeBackendManager.hpp"
#include "utils.hpp"
#include <unordered_map>
#include <sstream>
#define UNDERLINE "\033[4m"

#define CLOSEUNDERLINE "\033[0m"

int show_help() {
  std::cout << "NWQ-VQE Options" << std::endl;
  std::cout << UNDERLINE << "REQUIRED" << CLOSEUNDERLINE << std::endl;
  std::cout << "--hamiltonian, -f     Path to the input Hamiltonian file (formatted as a sum of Fermionic operators, see examples)" << std::endl;
  std::cout << "--nparticles, -n      Number of electrons in molecule" << std::endl;
  std::cout << "--backend, -b         Simulation backend. Defaults to CPU" << std::endl;
  std::cout << "--list-backends, -l   List available backends and exit." << std::endl;
  std::cout << UNDERLINE << "OPTIONAL" << CLOSEUNDERLINE << std::endl;
  std::cout << "--seed                Random seed for initial point and empirical gradient estimation. Defaults to time(NULL)" << std::endl;
  std::cout << "--config              Path to config file for NLOpt optimizer parameters" << std::endl;
  std::cout << "--optimizer           NLOpt optimizer name. Defaults to LN_COBYLA" << std::endl;
  std::cout << "--reltol              Relative tolerance termination criterion. Defaults to -1 (off)" << std::endl;
  std::cout << "--abstol              Relative tolerance termination criterion. Defaults to -1 (off)" << std::endl;
  std::cout << "--maxeval             Maximum number of function evaluations for optimizer. Defaults to 200" << std::endl;
  std::cout << "--maxtime             Maximum optimizer time (seconds). Defaults to -1.0 (off)" << std::endl;
  std::cout << "--stopval             Cutoff function value for optimizer. Defaults to -MAXFLOAT (off)" << std::endl;
  std::cout << "--xacc                Use XACC indexing scheme, otherwise uses DUCC scheme." << std::endl;
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
                bool& use_xacc,
                unsigned& seed) {
  std::string config_file = "";
  std::string algorithm_name = "LN_COBYLA";
  hamilfile = "";
  backend = "CPU";
  n_particles = -1;
  settings.max_evals = 200;
  seed = time(NULL);
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
      backend = argv[++i];//-2.034241 -1.978760  -1.825736
      continue;
    } 
    if (argname == "-f" || argname == "--hamiltonian") {
      hamilfile = argv[++i];
      continue;
    } else 
    if (argname == "-n" || argname == "--nparticles") {
      n_particles = std::atoll(argv[++i]);
    } else 
    if (argname == "--seed") {
      seed = (unsigned)std::atoi(argv[++i]);
    } else if (argname == "--xacc") {
      use_xacc = true;
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


void optimize_ansatz(const VQEBackendManager& manager,
                     const std::string& backend,
                     std::shared_ptr<NWQSim::VQE::Hamiltonian> hamil,
                     std::shared_ptr<NWQSim::VQE::Ansatz> ansatz,
                     NWQSim::VQE::OptimizerSettings& settings,
                     nlopt::algorithm& algo,
                     unsigned& seed,
                     std::vector<double>& params,
                     double& fval) {
  std::shared_ptr<NWQSim::VQE::VQEState> state = manager.create_vqe_solver(backend, ansatz, hamil, algo, callback_function, seed, settings);  
  std::uniform_real_distribution<double> initdist(-2 * PI, 2 * PI);
  std::mt19937_64 random_engine (seed);
  params.resize(ansatz->numParams());
  std::generate(params.begin(), params.end(), 
      [&random_engine, &initdist] () {return initdist(random_engine);});
  state->optimize(params, fval);
}


int main(int argc, char** argv) {
  VQEBackendManager manager;
  std::string hamil_path, backend;
  NWQSim::IdxType n_part;
  NWQSim::VQE::OptimizerSettings settings;
  nlopt::algorithm algo;
  bool use_xacc;
  unsigned seed;
  std::cout << "Parsing args" << std::endl;
  if (parse_args(argc, argv, manager, hamil_path, backend, n_part, algo, settings, use_xacc, seed)) {
    return 1;
  }
#ifdef MPI_ENABLED
  std::cout << "Initializing MPI" << std::endl;
  int i_proc;
  if (backend == "MPI" || backend == "NVGPU_MPI")
  {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &i_proc);
  }
#endif
  manager.safe_print("Reading Hamiltonian...\n");
  std::shared_ptr<NWQSim::VQE::Hamiltonian> hamil = std::make_shared<NWQSim::VQE::Hamiltonian>(hamil_path, n_part, use_xacc);
  manager.safe_print("Constructed %lld Pauli Observables\n", hamil->num_ops());
  manager.safe_print("Constructing UCCSD Ansatz...\n");

  std::shared_ptr<NWQSim::VQE::Ansatz> ansatz = std::make_shared<NWQSim::VQE::UCCSD>(
    hamil->getEnv(),
    NWQSim::VQE::getJordanWignerTransform,
    1
  );

  manager.safe_print("%lld Gates with %lld parameters\n" ,ansatz->num_gates(), ansatz->numParams());
  std::vector<double> params;
  double fval;
  manager.safe_print("Beginning VQE loop...\n");
  optimize_ansatz(manager, backend, hamil, ansatz, settings, algo, seed, params, fval);
  
  std::string qasm_string = ansatz->toQASM3();
  std::ofstream outfile;
  outfile.open("../uccsd.qasm", std::fstream::out);
  outfile << qasm_string;
  outfile.close();
  std::vector<std::pair<std::string, double> > param_map = ansatz->getFermionicOperatorParameters();
  manager.safe_print("\nFinished VQE loop.\n\tFinal value: %e\n\tFinal parameters:\n", fval);
  for (auto& pair: param_map) {
    manager.safe_print("%s :: %e\n", pair.first.c_str(), pair.second);

  }
#ifdef MPI_ENABLED
  if (backend == "MPI" || backend == "NVGPU_MPI")
  {
    MPI_Finalize();
  }
#endif
  return 0;
}