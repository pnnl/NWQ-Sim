#include "argparse.hpp"
#include "svsim_vqe/sv_cpu_vqe.hpp"
#include "private/nlohmann/json.hpp"
#include "utils.hpp"



using json = nlohmann::json;
int parse_args(int argc, char** argv,
                std::string& hamilfile,
                NWQSim::IdxType& n_particles,
                nlopt::algorithm& algo,
                NWQSim::VQE::OptimizerSettings& settings,
                unsigned& seed) {
  argparse::ArgumentParser parser("NWQ-VQE");
  parser.add_argument("--hamiltonian", "-f").help("Path to Fermionic Hamiltonian file");
  parser.add_argument("--nparticles", "-n").help("Number of particles to model").scan<'i', int>();
  parser.add_argument("--config", "-c").help("Optimizer configuration file").default_value(std::string(""));
  parser.add_argument("--maxeval").help("Maximum number of function evaluations").scan<'i', long long>().default_value(settings.max_evals);
  parser.add_argument("--seed").help("Random seed").scan<'u', unsigned>().default_value((unsigned)time(NULL));
  parser.add_argument("--reltol").help("Relative function value tolerance cutoff").scan<'g', double>().default_value(settings.rel_tol);
  parser.add_argument("--abstol").help("Absolute function value tolerance cutoff").scan<'g', double>().default_value(settings.abs_tol);
  parser.add_argument("--maxtime").help("Optimizer timeout (seconds)").scan<'g', double>().default_value(settings.max_time);
  parser.add_argument("--stopval").help("Optimizer value cutoff").scan<'g', double>().default_value(settings.stop_val);
  parser.add_argument("--optimizer").help(
    "Optimization algorithm name (NLOpt key, see [NLOpt Algorithms](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/))"
    ).default_value(std::string("LN_COBYLA"));
  try {
    parser.parse_args(argc, argv);
  }
  catch (const std::exception& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << parser;
    return 1;
  }
  hamilfile = parser.get<std::string>("--hamiltonian");
  n_particles = parser.get<int>("--nparticles");
  seed = parser.get<unsigned>("--seed");
  std::string config_file = parser.get<std::string>("--config");
  std::string algorithm_name = parser.get<std::string>("--optimizer");
  settings.rel_tol = parser.get<double>("--reltol");
  settings.abs_tol = parser.get<double>("--abstol");
  settings.max_evals = parser.get<long long>("--maxeval");
  settings.stop_val = parser.get<double>("--stopval");
  settings.max_time = parser.get<double>("--maxtime");
  algo = (nlopt::algorithm)nlopt_algorithm_from_string(algorithm_name.c_str());
  std::cout << nlopt::algorithm_name(algo) << std::endl;
  if (config_file != "") {
    std::ifstream f(config_file);
    json data = json::parse(f); 
    for (json::iterator it = data.begin(); it != data.end(); ++it) {
      settings.parameter_map[it.key()] = it.value().get<NWQSim::ValType>();
    }
  }

}


// Callback function, requires signature (void*) (const std::vector<NWQSim::ValType>&, NWQSim::ValType, NWQSim::IdxType)
void carriage_return_callback_function(const std::vector<NWQSim::ValType>& x, NWQSim::ValType fval, NWQSim::IdxType iteration) {
  printf("\33[2K\rEvaluation %lld, fval = %f", iteration, fval);fflush(stdout);
}

// Callback function, requires signature (void*) (const std::vector<NWQSim::ValType>&, NWQSim::ValType, NWQSim::IdxType)
void callback_function(const std::vector<NWQSim::ValType>& x, NWQSim::ValType fval, NWQSim::IdxType iteration) {
  printf("\33[2K\rEvaluation %lld, fval = %f", iteration, fval);fflush(stdout);
}
int main(int argc, char** argv) {
  std::string hamil_path;
  NWQSim::IdxType n_part;
  NWQSim::VQE::OptimizerSettings settings;
  nlopt::algorithm algo;
  unsigned seed;
  parse_args(argc, argv, hamil_path, n_part, algo, settings, seed);
  NWQSim::VQE::Hamiltonian hamil(hamil_path, n_part);

  std::shared_ptr<NWQSim::VQE::Ansatz> ansatz = std::make_shared<NWQSim::VQE::UCCSD>(
    hamil.getEnv(),
    NWQSim::VQE::getJordanWignerTransform,
    1
  );
  NWQSim::VQE::SV_CPU_VQE state(ansatz, hamil, algo, carriage_return_callback_function, seed, settings);
  
  std::uniform_real_distribution<double> initdist(0, 2 * PI);
  std::mt19937_64 random_engine (seed);
  std::vector<double> params(ansatz->numParams());
  std::generate(params.begin(), params.end(), 
      [&random_engine, &initdist] () {return initdist(random_engine);});

  double fval;
  state.optimize(params, fval);
  std::cout << std::endl << "Final Energy: " << fval << std::endl;
  std::cout << "Final Parameters: " << params << std::endl;
  return 0;
}