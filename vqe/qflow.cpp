#include "include/utils.hpp"
#include "vqeBackendManager.hpp"
#include "utils.hpp"
#include "vqe_state.hpp"
#include <unordered_map>
#include "src/uccsdmin.cpp"
#include "src/singletgsd.cpp"
#include <sstream>
#include <chrono>
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
  std::cout << "--amplitudes          List of initial amplitudes." << std::endl;
  std::cout << "--xacc                Use XACC indexing scheme, otherwise uses DUCC scheme." << std::endl; // MZ: it was --ducc, but the code only detects -- ducc
  std::cout << "--seed                Random seed for initial point and empirical gradient estimation. Defaults to time(NULL)" << std::endl;
  std::cout << "--verbose             Print optimizer information on each iteration. Defaults to false" << std::endl;
  // std::cout << "--symm                Symmetry level (0->none, 1->single, 2->double non-mixed+single, 3->all). Defaults to 0." << std::endl;
  std::cout << "--symm                Symmetry level (0->none, 1->spin symmetry, 2->also orbital symmetry). Defaults to 0." << std::endl; // MZ: is this call orbital symmetry?
  std::cout << "--gsd                 Use singlet GSD ansatz for ADAPT-VQE. Default to false." << std::endl;
  std::cout << "--origin              Use old implementatin of UCCSD for VQE or ADAPT-VQE. Have duplicated operators and potential symmetry problem. Default to false." << std::endl;
  std::cout << UNDERLINE << "OPTIONAL (Global Minimizer)" << CLOSEUNDERLINE << std::endl;
  std::cout << "--optimizer           NLOpt optimizer name. Defaults to LN_COBYLA. Suggests to use LN_NEWUOA" << std::endl;
  std::cout << "--optimizer-config    Path to config file for NLOpt optimizer parameters" << std::endl;
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
                std::string& config_path,
                std::string& amplitudes,
                NWQSim::IdxType& n_particles,
                nlopt::algorithm& algo,
                NWQSim::VQE::OptimizerSettings& settings,
                int& n_trials,
                bool& use_xacc,
                bool& local,
                bool& verbose,
                uint64_t& symm_level,
                unsigned& seed,
                double& delta,
                double& eta,
                bool& origin_uccsd,
                bool& gsd) {
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
  symm_level = 0;
  settings.lbound = -2;
  settings.ubound = 2;
  origin_uccsd = false;
  gsd = false;
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
    if (argname == "--config") {
      config_file = argv[++i];
      continue;
    }  else 
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
    if (argname == "--symm") {
      symm_level = (unsigned)std::atoi(argv[++i]);;
    } else
    if (argname == "--xacc") {
      use_xacc = true;
    } else
    if (argname == "--amplitudes") {
      amplitudes = argv[++i];
    } else
    if (argname == "--verbose") {
      verbose = true;
    } else if (argname == "--delta") {
      delta = std::atof(argv[++i]);
    } else if (argname == "--eta") {
      eta = std::atof(argv[++i]);
    } else 
    if (argname == "--optimizer-config") {
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
    } else if (argname == "--origin") {
      origin_uccsd = true;
    } else if (argname == "--gsd") {
      gsd = true;
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

void print_header() {
    std::cout << "\n----- Iteration Summary -----\n" << std::left
              << std::setw(8) << " Iter."
              << std::setw(19) << "Objective Value"
              << std::setw(55) << "Parameters (first 5)"
              << std::endl;
    std::cout << std::string(80, '-') << std::endl;
}

// Callback function, requires signature (void*) (const std::vector<NWQSim::ValType>&, NWQSim::ValType, NWQSim::IdxType)
void callback_function_simple(const std::vector<NWQSim::ValType>& x, NWQSim::ValType fval, NWQSim::IdxType iteration) {
  if (iteration == 0) {
    print_header();
  }
  std::cout << std::left << " "
            << std::setw(7) << iteration
            << std::setw(19) << std::fixed << std::setprecision(12) << fval;
  
  std::cout << std::fixed << std::setprecision(6);
  for (size_t i = 0; i < std::min(x.size(), size_t(5)); ++i) {
      std::cout  << std::setw(11) << x[i];
  }
  std::cout << std::endl;
}

std::string get_termination_reason(nlopt::result result) {
    static const std::map<nlopt::result, std::string> reason_map = {
        {nlopt::SUCCESS, "Optimization converged successfully"},
        {nlopt::STOPVAL_REACHED, "Objective function value reached the specified stop value"},
        {nlopt::FTOL_REACHED, "Function tolerance reached"},
        {nlopt::XTOL_REACHED, "Variable tolerance reached"},
        {nlopt::MAXEVAL_REACHED, "Maximum number of evaluations reached"},
        {nlopt::MAXTIME_REACHED, "Maximum allowed time reached"},
        {nlopt::FAILURE, "Optimization failed"},
        {nlopt::INVALID_ARGS, "Invalid arguments"},
        {nlopt::OUT_OF_MEMORY, "Out of memory"},
        {nlopt::ROUNDOFF_LIMITED, "Roundoff limited"},
        {nlopt::FORCED_STOP, "Optimization was forcibly stopped by a callback function"}
    };
    auto it = reason_map.find(result);
    if (it != reason_map.end()) {
        return it->second;
    } else {
        return "Unknown reason, code: " + std::to_string(result);
    }
}

void optimize_ansatz(const VQEBackendManager& manager,
                     const std::string& backend,
                     const std::string& config,
                     const std::string& amplitudes,
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
  double fval;
  // NWQSim::VQE::Callback callback = verbose ? carriage_return_callback_function: silent_callback_function;
  NWQSim::VQE::Callback callback = verbose ? callback_function_simple: silent_callback_function;
  if (!verbose) {
    NWQSim::Config::PRINT_SIM_TRACE = false;
  }
  std::shared_ptr<NWQSim::VQE::VQEState> state = manager.create_vqe_solver(backend, config, ansatz, hamil, algo, callback, seed, settings); 
  params.resize(ansatz->numParams());
  std::cout << ansatz->numParams() << std::endl;
    std::fill(params.begin(), params.end(), 0);
  if (amplitudes != "") {
    std::cout << "Reading amplitudes from " << amplitudes << std::endl;
    NWQSim::VQE::read_amplitudes(amplitudes, params, ansatz->get_excitation_map());
  }
  std::cout << params << std::endl;
  double initial_ene, final_ene;
  long long num_iterations = 0;
  std::vector<std::pair<std::string, double> > param_tuple;
  state->initialize();

  auto start_time = std::chrono::high_resolution_clock::now(); // MZ: time the optimization
  if (local) {
  param_tuple = state->follow_fixed_gradient(params, 
                                              initial_ene, 
                                              final_ene, 
                                              num_iterations, 
                                              delta, 
                                              eta, 
                                              num_trials);

  } else {
    ansatz->setParams(params);
    if (settings.max_evals > 0) 
      state->optimize(params, final_ene);
    param_tuple = ansatz->getFermionicOperatorParameters();
    num_iterations = state->get_iteration();
  }
  auto end_time = std::chrono::high_resolution_clock::now(); // MZ: time the optimization
  double opt_duration = std::chrono::duration<double>(end_time - start_time).count(); // MZ: time the optimization
  state -> set_duration(opt_duration); // MZ: time the optimization

  std::ostringstream strstream;
  // manager.safe_print("\nFinished in %llu iterations. Initial Energy %f, Final Energy %f\nPrinting excitation amplitudes:\n", num_iterations, initial_ene, final_ene);
  // for (auto& i: param_tuple) {
  //   strstream << i.first << ": " << i.second << std::endl;
  // }
  // manager.safe_print("%s", strstream.str().c_str());
  try { // MZ: catch exception
      double total_seconds = state -> get_duration();
      // Calculate hours, minutes, and seconds
      int hours = static_cast<int>(total_seconds) / 3600;
      int minutes = (static_cast<int>(total_seconds) % 3600) / 60;
      double seconds = total_seconds - (hours * 3600 + minutes * 60);

      // Print out the Fermionic operators with their excitations                                        
      std::vector<std::pair<std::string, double> > param_map = ansatz->getFermionicOperatorParameters(); 
      manager.safe_print("\n--------- Result Summary ---------\n"); 
      manager.safe_print("Method                 : QFlow\n");
      manager.safe_print("Ansatz                 : %s\n", ansatz->getAnsatzName().c_str());  // MZ: don't want to scroll all the way up to see this
      manager.safe_print("# Ham. Pauli Strings   : %lld \n", hamil->num_ops());
      manager.safe_print("Circuit Stats          : %lld operators, %lld parameters, and %lld Gates\n" , ansatz->numOps(), ansatz->numParams(), ansatz->num_gates()); // MZ: don't want to scroll all the way up to see this
      std::string ter_rea = "Optimization terminated: "+get_termination_reason(state->get_optresult())+"\n";
      manager.safe_print(ter_rea.c_str());
      manager.safe_print("# iterations           : %d\n", num_iterations);
      manager.safe_print("Evaluation Time        : %d hrs %d mins %.4f secs\n", hours, minutes, seconds);
      manager.safe_print("Initial objective value: %.16f\n", initial_ene); 
      manager.safe_print("Final objective value  : %.16f\nFinal parameters:\n", final_ene); 
      for (auto& i: param_tuple) {
        strstream << i.first << ": " << i.second << std::endl;
      }
      manager.safe_print("%s", strstream.str().c_str());                                                                                            
    } 
    catch(std::exception &e) {
        std::cout << "Optimization failed: " << e.what() << std::endl;
    } // MZ: catch exception
}


int main(int argc, char** argv) {
  VQEBackendManager manager;
  std::string hamil_path, backend, config, amplitudes;
  NWQSim::IdxType n_part;
  NWQSim::VQE::OptimizerSettings settings;
  nlopt::algorithm algo;
  double delta;
  double eta;
  unsigned seed;
  uint64_t symm_level;
  bool use_xacc, local, verbose;
  int n_trials;
  bool origin_uccsd;
  bool gsd;

  if (parse_args(argc, argv, manager, hamil_path, backend, config, amplitudes,  n_part, algo, settings, 
                 n_trials, use_xacc, local, verbose, symm_level, 
                  seed, delta, eta, origin_uccsd, gsd)) {
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

  std::shared_ptr<NWQSim::VQE::Ansatz> ansatz;
  if (origin_uccsd) {
    ansatz = std::make_shared<NWQSim::VQE::UCCSD>(
      hamil->getEnv(),
      NWQSim::VQE::getJordanWignerTransform,
      1,
      symm_level
    );
  } else if (gsd) {
      ansatz  = std::make_shared<NWQSim::VQE::Singlet_GSD>(
      hamil->getEnv(),
      NWQSim::VQE::getJordanWignerTransform,
      1
    );
  } else {
    ansatz = std::make_shared<NWQSim::VQE::UCCSDmin>(
      hamil->getEnv(),
      NWQSim::VQE::getJordanWignerTransform,
      1,
      symm_level
    );
  }

  ansatz->buildAnsatz();
  std::vector<double> params;
  manager.safe_print("Beginning VQE loop...\n");
  optimize_ansatz(manager, backend, config, amplitudes, hamil, ansatz, settings, algo, seed, n_trials, params, local, verbose, delta, eta);
#ifdef MPI_ENABLED
  if (backend == "MPI" || backend == "NVGPU_MPI")
  {
    MPI_Finalize();
  }
#endif
  return 0;
}