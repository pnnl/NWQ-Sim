#include "vqeBackendManager.hpp"
#include "utils.hpp"
#include <unordered_map>
#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

namespace py = pybind11;

// Callback function, requires signature (void*) (const std::vector<NWQSim::ValType>&, NWQSim::ValType, NWQSim::IdxType)
void null_callback_function(const std::vector<NWQSim::ValType>& x, NWQSim::ValType fval, NWQSim::IdxType iteration) {
}

// Callback function, requires signature (void*) (const std::vector<NWQSim::ValType>&, NWQSim::ValType, NWQSim::IdxType)
void callback_function(const std::vector<NWQSim::ValType>& x, NWQSim::ValType fval, NWQSim::IdxType iteration) {
  printf("\33[2KEvaluation %lld, fval = %f\n", iteration, fval);fflush(stdout);
}


std::vector<std::pair<std::string, double> > optimize_ansatz(const VQEBackendManager& manager,
                     const std::string& backend,
                     std::shared_ptr<NWQSim::VQE::Hamiltonian> hamil,
                     std::shared_ptr<NWQSim::VQE::Ansatz> ansatz,
                     NWQSim::VQE::OptimizerSettings& settings,
                     nlopt::algorithm& algo,
                     unsigned& seed,
                     double eta,
                     double delta,
                     std::vector<double>& params,
                     double& fval) {
  py::print("Started function");
  std::shared_ptr<NWQSim::VQE::VQEState> state = manager.create_vqe_solver(backend, ansatz, hamil, algo, null_callback_function, seed, settings);  
  std::uniform_real_distribution<double> initdist(0, 2 * PI);
  py::print("Calling optimization routine");

  std::vector<std::pair<std::string, double> > param_tuple = state->follow_fixed_gradient(params, fval, eta, delta);
  py::print("Exited optimization routine");
  return param_tuple;
}

std::vector<std::pair<std::string, double> > 
// std::string
optimize_effective_hamiltonian(
  const std::vector<std::pair<std::string, std::complex<double>>>& fermionic_operators,
  NWQSim::VQE::IdxType n_particles,
  const std::vector<double>& init_params,
  std::string backend = "CPU",
  uint32_t seed = 0,
  double delta = 1e-4,
  double eta = 1e-4,
  int n_trotter = 1) {
  seed = seed;
  std::shared_ptr<NWQSim::VQE::Hamiltonian> hamil = std::make_shared<NWQSim::VQE::Hamiltonian>(fermionic_operators, n_particles, NWQSim::VQE::getJordanWignerTransform);
  std::shared_ptr<NWQSim::VQE::Ansatz> ansatz = std::make_shared<NWQSim::VQE::UCCSD>(
                                                    hamil->getEnv(),
                                                    NWQSim::VQE::getJordanWignerTransform,
                                                    n_trotter
                                                  );
  if (init_params.size() != ansatz->numParams()) {
    throw std::runtime_error("Not enough initial parameters provided to ansatz, please pass " + std::to_string(ansatz->numParams()) + " parameters\n");
  }
  std::vector<double> local_params (init_params);
  double fval;
  VQEBackendManager manager;

  NWQSim::VQE::OptimizerSettings settings;
  nlopt::algorithm algo = nlopt::algorithm::LN_COBYLA; // default value for ctor, not actually used
  std::vector<std::pair<std::string, double> > result = optimize_ansatz(manager, backend, hamil, ansatz, settings, algo, seed, eta, delta, local_params, fval);
    // throw std::runtime_error("Done\n");
  // std::vector<std::pair<std::string, double> > result = {{"test", 0.0}};
  return result;
}
PYBIND11_MODULE(nwqflow, m) {
    m.doc() = "QFlow backend based on NWQ-Sim"; // optional module docstring

    m.def("optimize_effective_hamiltonian", 
    &optimize_effective_hamiltonian, 
    "Perform single-direction gradient descent using an SPSA-estimated gradient and return the locally-optimal parameters along with the associated Fermionic excitation.");
}
/*
int main(int argc, char** argv) {
  VQEBackendManager manager;
  std::string hamil_path, backend;
  NWQSim::IdxType n_part;
  NWQSim::VQE::OptimizerSettings settings;
  nlopt::algorithm algo;
  unsigned seed;
  if (parse_args(argc, argv, manager, hamil_path, backend, n_part, algo, settings, seed)) {
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
  std::shared_ptr<NWQSim::VQE::Hamiltonian> hamil = std::make_shared<NWQSim::VQE::Hamiltonian>(hamil_path, n_part);
  manager.safe_print("Constructing UCCSD Ansatz...\n");

  std::shared_ptr<NWQSim::VQE::Ansatz> ansatz = std::make_shared<NWQSim::VQE::UCCSD>(
    hamil->getEnv(),
    NWQSim::VQE::getJordanWignerTransform,
    1
  );
  std::vector<double> params;
  double fval;
  manager.safe_print("Beginning VQE loop...\n");
  optimize_ansatz(manager, backend, hamil, ansatz, settings, algo, seed, params, fval);
  std::ostringstream paramstream;
  paramstream << params;
  manager.safe_print("\nFinished VQE loop.\n\tFinal value: %e\n\tFinal parameters: %s\n", fval, paramstream.str().c_str());
#ifdef MPI_ENABLED
  if (backend == "MPI" || backend == "NVGPU_MPI")
  {
    MPI_Finalize();
  }
#endif
  return 0;
}*/