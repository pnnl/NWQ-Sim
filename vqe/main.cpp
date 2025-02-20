#include "vqeBackendManager.hpp"
#include "utils.hpp"
#include <cmath>
#include <string>
#include "circuit/dynamic_ansatz.hpp"
#include "vqe_adapt.hpp"
#include "src/uccsdmin.cpp"
#include "src/singletgsd.cpp"
#include "src/ansatz_pool.cpp" // MZ: move the ansatz pool generation out of the src/utils.cpp
#include "nwq_util.hpp"
#include <chrono>

#define UNDERLINE "\033[4m"

#define CLOSEUNDERLINE "\033[0m"
using IdxType = NWQSim::IdxType;
using ValType = NWQSim::ValType;

struct VQEParams {
  /**
   * @brief  Structure to hold command line arguments for NWQ-VQE
   */
  // Problem Info
  std::string hamiltonian_path = "";
  IdxType nparticles = -1;
  bool xacc = true;

  // Simulator options
  std::string backend = "CPU";
  uint32_t seed;
  // Optimizer settings
  NWQSim::VQE::OptimizerSettings optimizer_settings;
  nlopt::algorithm algo;
  IdxType symm_level = 0;

  bool verbose = false; // MZ: if show optimization results in each iteration

  // ADAPT-VQE options
  bool adapt = false;
  bool qubit = false;
  IdxType adapt_maxeval = 100;
  ValType adapt_fvaltol = -1; // MZ: original ADAPT-VQE paper only used the gradient norm as the convergence criteria, 
                              // which is more reasonable as since fvaltol may give false convergence when the same operator is picked consecutively
  ValType adapt_gradtol = 1e-3;
  IdxType adapt_pool_size = -1;

  // Ansatz options
  NWQSim::VQE::PoolType pool = NWQSim::VQE::PoolType::Fermionic;
};
int show_help() {
  std::cout << "NWQ-VQE Options" << std::endl;
  std::cout << UNDERLINE << "INFORMATIONAL" << CLOSEUNDERLINE << std::endl;
  std::cout << "-h, --help            Show help menu." << std::endl;
  std::cout << "-l, --list-backends   List available backends and exit." << std::endl;
  std::cout << UNDERLINE << "REQUIRED" << CLOSEUNDERLINE << std::endl;
  std::cout << "-f, --hamiltonian     Path to the input Hamiltonian file (formatted as a sum of Fermionic operators, see examples)" << std::endl;
  std::cout << "-p, --nparticles      Number of electrons in molecule" << std::endl;
  std::cout << UNDERLINE << "OPTIONAL (Hamiltonian, Ansatz and Backend)" << CLOSEUNDERLINE << std::endl;
  std::cout << "--ducc                Use DUCC indexing scheme. Defaults to false (defaults to use XACC/Qiskit scheme)." << std::endl;
  // std::cout << "--symm                Symmetry level (0->none, 1->single, 2->double non-mixed+single, 3->all). Defaults to 0." << std::endl;
  std::cout << "--sym, --symm         UCCSD Symmetry level (0->none, 1->spin symmetry, 2->also orbital symmetry). Defaults to 0." << std::endl; // MZ: is this call orbital symmetry?
  std::cout << "--gsd                 Use singlet GSD ansatz for ADAPT-VQE. Default to false." << std::endl;
  std::cout << "-b, --backend         Simulation backend. Defaults to CPU" << std::endl;
  std::cout << "--seed                Random seed for initial point and empirical gradient estimation. Defaults to time(NULL)" << std::endl;
  // std::cout << "--config              Path to NWQ-Sim config file. Defaults to \"../default_config.json\"" << std::endl;
  std::cout << UNDERLINE << "OPTIONAL (Global Minimizer)" << CLOSEUNDERLINE << std::endl;
  std::cout << "-v, --verbose         Print optimization callback on each VQE iteration. Defaults to false" << std::endl;
  std::cout << "-o, --optimizer       NLOpt optimizer name. Defaults to LN_COBYLA. Some common optimizers are LN_BOBYQA, LN_NEWUOA, and LD_LBFGS." << std::endl;
  std::cout << "--opt-config          Path to config file for NLOpt optimizer parameters" << std::endl;
  std::cout << "-lb, --lbound         Optimizer lower bound. Defaults to -Pi" << std::endl;
  std::cout << "-ub, --ubound         Optimizer upper bound. Defaults to Pi" << std::endl;  std::cout << "--reltol              Relative tolerance termination criterion. Defaults to -1 (off)" << std::endl;
  std::cout << "--abstol              Relative tolerance termination criterion. Defaults to -1 (off)" << std::endl;
  std::cout << "--maxeval             Maximum number of function evaluations for optimizer (only for VQE). Defaults to 100" << std::endl;
  std::cout << "--maxtime             Maximum optimizer time (seconds). Defaults to -1.0 (off)" << std::endl;
  std::cout << "--stopval             Cutoff function value for optimizer. Defaults to -MAXFLOAT (off)" << std::endl;
  std::cout << UNDERLINE << "ADAPT-VQE OPTIONS" << CLOSEUNDERLINE << std::endl;
  std::cout << "--adapt               Use ADAPT-VQE for dynamic ansatz construction. Defaults to false" << std::endl;
  std::cout << "-ag, --adapt-gradtol  Cutoff absolute tolerance for operator gradient norm. Defaults to 1e-3" << std::endl;
  std::cout << "-af, --adapt-fvaltol  Cutoff absolute tolerance for function value. Defaults to -1 (off)" << std::endl;
  std::cout << "-am, --adapt-maxeval  Set a maximum iteration count for ADAPT-VQE. Defaults to 100" << std::endl;
  std::cout << "--qubit               Uses Qubit instead of Fermionic operators for ADAPT-VQE. Defaults to false" << std::endl;
  std::cout << "--adapt-pool          Sets the pool size for Qubit operators. Defaults to -1" << std::endl;
  std::cout << UNDERLINE << "LEGACY" << CLOSEUNDERLINE << std::endl;
  std::cout << "--xacc                Use XACC indexing scheme, otherwise uses DUCC scheme. (Deprecated, true by default)" << std::endl;
  std::cout << "--origin              Use old implementatin of UCCSD for VQE or ADAPT-VQE. Have duplicated operators and potential symmetry problem. Default to false." << std::endl;
  std::cout << "-n                    (Same as -p or --nparticules) Number of electrons in molecule" << std::endl; // MZ: this could be confusing as people usually use n for number of qubits (2*number of spartial orbitals)
  std::cout << UNDERLINE << "SIMULATOR OPTIONS" << CLOSEUNDERLINE << std::endl;
  std::cout << "--num_threads         Specify the number of OMP threads. Defaults to use all hardware threads" << std::endl;
  std::cout << "--disable_fusion      Disable gate fusion. Defaults to enabled" << std::endl;
  return 1;
}

using json = nlohmann::json;

/**
 * @brief  Parse command line arguments
 * @note   
 * @param  argc: Number of command line arguments passed
 * @param  argv: Pointer to command line arg C strings
 * @param  manager: Backend manager object
 * @param  params: Data structure to store command line arguments
 * @retval 
 */
int parse_args(int argc, char** argv,
               VQEBackendManager& manager,
               VQEParams& params) {
  std::string opt_config_file = "";
  std::string algorithm_name = "LN_COBYLA";
  NWQSim::VQE::OptimizerSettings& settings = params.optimizer_settings;
  params.seed = time(NULL);
  settings.lbound = -2*PI; //
  settings.ubound = 2*PI; //
  for (size_t i = 1; i < argc; i++) {
    std::string argname = argv[i];
    if (argname == "-h" || argname == "--help") {
      return show_help();
    } if (argname == "-l" || argname == "--list-backends") {
      manager.print_available_backends();
      return 1;
    } else
    if (argname == "-b" || argname == "--backend") {
      params.backend = argv[++i];//-2.034241 -1.978760  -1.825736
      continue;
    } else
    if (argname == "-f" || argname == "--hamiltonian") {
      params.hamiltonian_path = argv[++i];
      continue;
    } else 
    if (argname == "-p" || argname == "-n" || argname == "--nparticles") {
      params.nparticles = std::atoll(argv[++i]);
    } else 
    if (argname == "-v" || argname == "--verbose") {
      params.verbose = true;
    } else
    if (argname == "--seed") {
      params.seed = (unsigned)std::atoi(argv[++i]);
    }  else 
    if (argname == "--adapt-pool") {
      params.adapt_pool_size = (long long)std::atoi(argv[++i]);
    } else if (argname == "--xacc") {
      params.xacc = true;
    } else if (argname == "--ducc") {
      params.xacc = false;
    } else 
    if (argname == "--opt-config") {
      opt_config_file = argv[++i];
    } else  
    if (argname == "-o" || argname == "--optimizer") {
      algorithm_name = argv[++i];
    } else 
    if (argname == "--reltol") {
      params.optimizer_settings.rel_tol = std::atof(argv[++i]);
    } else 
    if (argname == "--abstol") {
      settings.abs_tol = std::atof(argv[++i]);
    }  else 
    if (argname == "--sym" || argname == "--symm") {
      params.symm_level = (unsigned)std::atoi(argv[++i]);;
    } else
    if (argname == "-ub" || argname == "--ubound") {
      settings.ubound = std::atof(argv[++i]);
    }  else 
    if (argname == "-lb" || argname == "--lbound") {
      settings.lbound = std::atof(argv[++i]);
    } else 
    if (argname == "-af" || argname == "--adapt-fvaltol") {
      params.adapt_fvaltol = std::atof(argv[++i]);
    } else 
    if (argname == "-ag" || argname == "--adapt-gradtol") {
      params.adapt_gradtol = std::atof(argv[++i]);
    }  else 
    if (argname == "-am" || argname == "--adapt-maxeval") {
      params.adapt_maxeval = std::atoll(argv[++i]);
    } else 
    if (argname == "--adapt") {
      params.adapt = true;
    }  else 
    if (argname == "--qubit") {
      // params.qubit = true;
      params.pool = NWQSim::VQE::PoolType::Pauli;
    } else 
    if (argname == "--maxeval") {
      settings.max_evals = std::atoll(argv[++i]);
    } else if (argname == "--stopval") {
      settings.stop_val = std::atof(argv[++i]);
    } else if (argname == "--maxtime") {
      settings.max_time = std::atof(argv[++i]);
    } else if (argname == "--gsd") {
      params.pool = NWQSim::VQE::PoolType::Singlet_GSD;
    } else if (argname == "--origin") {
      params.pool = NWQSim::VQE::PoolType::Fermionic_Origin;
    } else if (argname == "--num_threads") {
      NWQSim::Config::OMP_NUM_THREADS = std::atoi(argv[++i]);
    } else if (argname == "--disable_fusion") {
      NWQSim::Config::ENABLE_FUSION = false;
    } else {
      fprintf(stderr, "\033[91mERROR:\033[0m Unrecognized option %s, type -h or --help for a list of configurable parameters\n", argv[i]);
      return show_help();
    }
  }
  if (params.hamiltonian_path == "") {
      fprintf(stderr, "\033[91mERROR:\033[0m Must pass a Hamiltonian file (--hamiltonian, -f)\n");

      return show_help();
  }
  if (params.nparticles == -1) {
      fprintf(stderr, "\033[91mERROR:\033[0m Must pass a particle count (--nparticles, -n)\n");
      return show_help();
  }
  params.algo = (nlopt::algorithm)nlopt_algorithm_from_string(algorithm_name.c_str());
  if (opt_config_file != "") {
    std::ifstream f(opt_config_file);
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
  std::string paramstr = "[";
  for (auto i: x) {
    paramstr += std::to_string(i) + ", ";
  }
  printf("\33[2KEvaluation %lld, fval = %f, x=%s\n", iteration, fval, paramstr.c_str());fflush(stdout);
}
// Callback function, requires signature (void*) (const std::vector<NWQSim::ValType>&, NWQSim::ValType, NWQSim::IdxType)
void silent_callback_function(const std::vector<NWQSim::ValType>& x, NWQSim::ValType fval, NWQSim::IdxType iteration) {
  
}

//-------------------------------------------------------------------------------------------------
// MZ: sorry, original callback function is too much
void print_header() {
    std::cout << "\n----- Iteration Summary -----\n" << std::left
              << std::setw(8) << " Iter"
              << std::setw(29) << "Objective Value"
              << std::setw(55) << "Parameters (first 5)"
              << std::endl;
    std::cout << std::string(90, '-') << std::endl;
}

// Callback function, requires signature (void*) (const std::vector<NWQSim::ValType>&, NWQSim::ValType, NWQSim::IdxType)
void callback_function_simple(const std::vector<NWQSim::ValType>& x, NWQSim::ValType fval, NWQSim::IdxType iteration) {
  if (iteration == 0) {
    print_header();
  }
  std::cout << std::left << " "
            << std::setw(7) << iteration
            << std::setw(29) << std::fixed << std::setprecision(14) << fval;
  
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

std::string get_termination_reason_adapt(int result) {
    static const std::map<int, std::string> reason_map_adapt = {
        {0, "Abs. tol. for operator gradient norm is reached"},
        {1, "Abs. tol. for function value change is reached"},
        {2, "Max. number of iterations for ADAPT-VQE is reached"},
        {9, "ADAPT iteration is not run"}
    };
    auto it = reason_map_adapt.find(result);
    if (it != reason_map_adapt.end()) {
        return it->second;
    } else {
        return "Unknown reason, code: " + std::to_string(result);
    }
}
//-------------------------------------------------------------------------------------------------


/**
 * @brief  Optimized the UCCSD (or ADAPT-VQE) Ansatz and Report the Fermionic Excitations
 * @note   
 * @param  manager: API to handle calls to different backends
 * @param  params: Data structure with runtime-configurable options
 * @param  ansatz: Shared pointer to a parameterized quantum circuit
 * @param  hamil: Share pointer to a Hamiltonian observable
 * @param  x: Parameter vector (reference, output)
 * @param  fval: Energy value (reference, output)
 * @retval None
 */
std::shared_ptr<NWQSim::VQE::VQEState> optimize_ansatz(const VQEBackendManager& manager,
                     VQEParams params,
                     std::shared_ptr<NWQSim::VQE::Ansatz> ansatz,
                     std::shared_ptr<NWQSim::VQE::Hamiltonian> hamil,
                     std::vector<double>& x,
                     double& fval) {
  // Set the callback function (silent is default)
  // NWQSim::VQE::Callback callback = (params.adapt ? silent_callback_function : callback_function); // MZ: sorry, original callback function is too much
  // NWQSim::VQE::Callback callback = (params.adapt ? silent_callback_function : callback_function_simple); 
  NWQSim::VQE::Callback callback;
  if ( (params.adapt) || (!params.verbose) ) {
    callback = silent_callback_function;
  } else {
    callback = callback_function_simple;
  }
  
  std::shared_ptr<NWQSim::VQE::VQEState> state = manager.create_vqe_solver(params.backend,
                                                                           ansatz, 
                                                                           hamil, 
                                                                           params.algo, 
                                                                           callback, 
                                                                           params.seed, 
                                                                           params.optimizer_settings);  
  x.resize(ansatz->numParams());
  std::fill(x.begin(), x.end(), 0);

  
  if (params.adapt) {
    /***** ADAPT-VQE ******/
    // recast the ansatz pointer
    std::shared_ptr<NWQSim::VQE::DynamicAnsatz> dyn_ansatz = std::reinterpret_pointer_cast<NWQSim::VQE::DynamicAnsatz>(ansatz);
    // make the operator pool (either Fermionic or ADAPT)
    dyn_ansatz->make_op_pool(hamil->getTransformer(), params.seed, params.adapt_pool_size);
    // construct the AdaptVQE controller
    NWQSim::VQE::AdaptVQE adapt_instance(dyn_ansatz, state, hamil);
    // timing utilities
    auto start_time = std::chrono::high_resolution_clock::now();
    adapt_instance.make_commutators(); // Start making the commutators
    auto end_commutators = std::chrono::high_resolution_clock::now();
    // double commutator_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_commutators - start_time).count() / 1e9; // MZ: why nanoseconds?
    double commutator_time = std::chrono::duration<double>(end_commutators - start_time).count();
    // NWQSim::safe_print("Constructed ADAPT-VQE Commutators in %.2e seconds\n", commutator_time); // Report the commutator overhead
    state -> set_comm_duration(commutator_time); // MZ: time the commutator construction

    adapt_instance.optimize(x, fval, params.adapt_maxeval, params.adapt_gradtol, params.adapt_fvaltol); // MAIN OPTIMIZATION LOOP
    auto end_optimization = std::chrono::high_resolution_clock::now();
    // double optimization_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_optimization - end_commutators ).count() / 1e9; // MZ: why nanoseconds?
    double optimization_time = std::chrono::duration<double>(end_optimization - end_commutators ).count();
    // NWQSim::safe_print("Completed ADAPT-VQE Optimization in %.2e seconds\n", optimization_time); // Report the total time
    state -> set_duration(optimization_time); // MZ: time the optimization
    state -> set_numpauli(adapt_instance.get_numpauli()); // MZ: get the number of Pauli terms
    state -> set_numcomm(adapt_instance.get_numcomm()); // MZ: get the number of commuting cliques
//     double commutator_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_commutators - start_time).count() / 1e9;
//     NWQSim::safe_print("Constructed ADAPT-VQE Commutators in %.2e seconds\n", commutator_time); // Report the commutator overhead
//     adapt_instance.optimize(x, fval, params.adapt_maxeval, params.adapt_gradtol, params.adapt_fvaltol); // MAIN OPTIMIZATION LOOP
//     auto end_optimization = std::chrono::high_resolution_clock::now();
//     double optimization_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_optimization - end_commutators ).count() / 1e9;
//     NWQSim::safe_print("Completed ADAPT-VQE Optimization in %.2e seconds\n", optimization_time); // Report the total time
  } else {
    state->initialize(); // Initialize the state (AKA allocating measurement data structures and building the measurement circuit)

    auto start_time = std::chrono::high_resolution_clock::now(); // MZ: time the optimization
    state->optimize(x, fval); // MAIN OPTIMIZATION LOOP
    auto end_time = std::chrono::high_resolution_clock::now(); // MZ: time the optimization
    double opt_duration = std::chrono::duration<double>(end_time - start_time).count(); // MZ: time the optimization
    state -> set_duration(opt_duration); // MZ: time the optimization
  }
  return state; // MZ: I need information
}


int main(int argc, char** argv) {
  VQEBackendManager manager;
  VQEParams params;

  if (parse_args(argc, argv, manager, params)) {
    return 1;
  }
  #ifdef MPI_ENABLED
    int i_proc;
    if (params.backend == "MPI" || params.backend == "NVGPU_MPI")
    {
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &i_proc);
  }
  #endif

  // Get the Hamiltonian from the external file
  NWQSim::safe_print("Reading Hamiltonian...\n");
  std::shared_ptr<NWQSim::VQE::Hamiltonian> hamil = std::make_shared<NWQSim::VQE::Hamiltonian>(params.hamiltonian_path, 
                                                                                               params.nparticles,
                                                                                               params.xacc);
  NWQSim::safe_print("Constructed %lld Pauli observables\n", hamil->num_ops());
  NWQSim::safe_print("Constructing the ansatz...\n");
  
  // Build the parameterized ansatz
  std::shared_ptr<NWQSim::VQE::Ansatz> ansatz;
  if (params.adapt)
  {
    ansatz = std::make_shared<NWQSim::VQE::DynamicAnsatz>(hamil->getEnv(), params.pool);
  } else {
    // Static ansatz
    if (params.pool == NWQSim::VQE::PoolType::Fermionic_Origin) {
        ansatz  = std::make_shared<NWQSim::VQE::UCCSD>(
        hamil->getEnv(),
        NWQSim::VQE::getJordanWignerTransform,
        1,
        params.symm_level
      );
    } else if (params.pool == NWQSim::VQE::PoolType::Singlet_GSD) {
        ansatz  = std::make_shared<NWQSim::VQE::Singlet_GSD>(
        hamil->getEnv(),
        NWQSim::VQE::getJordanWignerTransform,
        1
      );
    } else {
      ansatz  = std::make_shared<NWQSim::VQE::UCCSDmin>(
        hamil->getEnv(),
        NWQSim::VQE::getJordanWignerTransform,
        1,
        params.symm_level
      );
    }
  }
  ansatz->buildAnsatz();

  NWQSim::safe_print("Starting with %lld Gates and %lld parameters\n" ,ansatz->num_gates(), ansatz->numParams());
  std::vector<double> x;
  double fval;
  if (params.adapt) {
    NWQSim::safe_print("Beginning ADAPT-VQE loop...\n");
  } else {
    NWQSim::safe_print("Beginning VQE loop...\n");
  }

//   // Print out the Fermionic operators with their excitations
//   std::vector<std::pair<std::string, double> > param_map = ansatz->getFermionicOperatorParameters();
//   NWQSim::safe_print("\nFinished VQE loop.\n\tFinal value: %e\n\tFinal parameters:\n", fval);
//   for (auto& pair: param_map) {
//     NWQSim::safe_print("%s :: %e\n", pair.first.c_str(), pair.second);

  // Print out the Fermionic operators with their excitations                                        
  // std::vector<std::pair<std::string, double> > param_map = ansatz->getFermionicOperatorParameters();  // MZ: comment out for better printing
  // NWQSim::safe_print("\nFinished VQE loop.\n\tFinal value: %e\n\tFinal parameters:\n", fval);         // MZ: comment out for better printing
  // for (auto& pair: param_map) {                                                                       // MZ: comment out for better printing
  //   NWQSim::safe_print("%s :: %e\n", pair.first.c_str(), pair.second);                                // MZ: comment out for better printing
  // }                                                                                                   // MZ: comment out for better printing


  std::shared_ptr<NWQSim::VQE::VQEState> opt_info = optimize_ansatz(manager, params, ansatz,  hamil, x, fval); // MZ: add opt_result

  double total_seconds = opt_info -> get_duration();
  int hours = static_cast<int>(total_seconds) / 3600;
  int minutes = (static_cast<int>(total_seconds) % 3600) / 60;
  double seconds = total_seconds - (hours * 3600 + minutes * 60);

  // Print out the Fermionic operators with their excitations                                        
  std::vector<std::pair<std::string, double> > param_map = ansatz->getFermionicOperatorParameters(); 
  // MZ: A better summary at the end so I don't scroll all the way up, especially with gradient-free optimizater    
  NWQSim::safe_print("\n--------- Result Summary ---------\n");
  if (params.adapt) {
    NWQSim::safe_print("Method                 : ADAPT-VQE\n");
  } else {
    NWQSim::safe_print("Method                 : VQE, Symmetry Level = %d\n", params.symm_level);
  }
  NWQSim::safe_print("Ansatz                 : %s\n", ansatz->getAnsatzName().c_str());  // MZ: don't want to scroll all the way up to see this
  NWQSim::safe_print("# Ham. Pauli Strings   : %lld \n", hamil->num_ops()); // MZ: don't want to scroll all the way up to see this
  if (params.adapt) {
    // MZ: time
    double comm_total_secs = opt_info -> get_comm_duration();
    int comm_hours = static_cast<int>(comm_total_secs) / 3600;
    int comm_minutes = (static_cast<int>(comm_total_secs) % 3600) / 60;
    double comm_seconds = comm_total_secs - (comm_hours * 3600 + comm_minutes * 60);
    NWQSim::safe_print("# commutators   (ADAPT): %d\n", opt_info->get_numcomm());
    NWQSim::safe_print("# Pauli Strings (ADAPT): %d\n", opt_info->get_numpauli());
    NWQSim::safe_print("Commutator Time (ADAPT): %d hrs %d mins %.4f secs\n", comm_hours, comm_minutes, comm_seconds);
  }
    NWQSim::safe_print("Operator Stats         : %lld operators, %lld parameters, and %lld Gates\n" , ansatz->numOps(), ansatz->numParams(), ansatz->num_gates()); // MZ: don't want to scroll all the way up to see this
    NWQSim::CircuitMetrics final_metrics = ansatz -> circuit_metrics();
    NWQSim::safe_print("Circuit Stats          : %lld depth, %lld 1q gates, %lld 2q gates, %.3f gate density\n", final_metrics.depth, final_metrics.one_q_gates, final_metrics.two_q_gates, final_metrics.gate_density);
    std::string ter_rea;
    if (params.adapt) {
      ter_rea = "Optimization terminated: "+get_termination_reason_adapt(opt_info->get_adaptresult())+"\n";
    } else {
      ter_rea = "Optimization terminated: "+get_termination_reason(opt_info->get_optresult())+"\n";
    }
    NWQSim::safe_print(ter_rea.c_str());
    if (params.adapt) {
      NWQSim::safe_print("# ADAPT rounds         : %d\n", opt_info->get_adaptrounds()); //MZ: ADAPT round
    } else {
      NWQSim::safe_print("# function eval.       : %d\n", opt_info->get_numevals());
    }
    NWQSim::safe_print("Evaluation Time        : %d hrs %d mins %.4f secs\n", hours, minutes, seconds);
    NWQSim::safe_print("Final objective value  : %.16f\nFinal parameters:\n", fval); 
    for (auto& pair: param_map) {                                                                       
      NWQSim::safe_print("  %s :: %.16f\n", pair.first.c_str(), pair.second);  
    }
    // // MZ: print the QASM3 circuit
    // if (params.adapt) {
    //   NWQSim::safe_print( (ansatz->toQASM3()).c_str() );
    // }

#ifdef MPI_ENABLED
  if (params.backend == "MPI" || params.backend == "NVGPU_MPI")
  {
    MPI_Finalize();
  }
#endif
  return 0;
}