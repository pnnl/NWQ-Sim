#include "vqeBackendManager.hpp"
#include "utils.hpp"
#include <unordered_map>
#include <sstream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <cmath>

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
                     int n_evals,
                     int& seed,
                     double eta,
                     double delta,
                     std::vector<double>& params,
                     double& fval) {
  // py::print("Started function");
  std::shared_ptr<NWQSim::VQE::VQEState> state = manager.create_vqe_solver(backend, ansatz, hamil, algo, null_callback_function, seed, settings);  
  std::uniform_real_distribution<double> initdist(0, 2 * PI);
  // py::print("Calling optimization routine");

  std::vector<std::pair<std::string, double> > param_tuple = state->follow_fixed_gradient(params, fval, eta, delta, n_evals);
  // py::print("Exited optimization routine");
  return param_tuple;
}

std::vector<std::pair<std::string, double> > 
// std::string
optimize_effective_hamiltonian(
  const std::vector<std::pair<std::string, std::complex<double>>>& fermionic_operators,
  NWQSim::VQE::IdxType n_particles,
  const std::vector<double>& init_params,
  std::string backend = "CPU",
  bool use_xacc = false,
  int seed = -1,
  double delta = 1e-4,
  double eta = 1e-4,
  int n_trotter = 1,
  int n_samples = 1) {
  seed = seed;
  std::shared_ptr<NWQSim::VQE::Hamiltonian> hamil = std::make_shared<NWQSim::VQE::Hamiltonian>(fermionic_operators, n_particles, use_xacc, NWQSim::VQE::getJordanWignerTransform);
  std::shared_ptr<NWQSim::VQE::Ansatz> ansatz = std::make_shared<NWQSim::VQE::UCCSD>(
                                                    hamil->getEnv(),
                                                    NWQSim::VQE::getJordanWignerTransform,
                                                    n_trotter
                                                  );
  if (init_params.size() != ansatz->numParams()) {
    throw std::runtime_error("Not enough initial parameters provided to ansatz, please pass " + std::to_string(ansatz->numParams()) + " parameters\n");
  }
  if (seed == -1) {
    seed = (int)time(NULL);
  }
  std::vector<double> local_params (init_params);
  double fval;
  VQEBackendManager manager;

  NWQSim::VQE::OptimizerSettings settings;
  nlopt::algorithm algo = nlopt::algorithm::LN_COBYLA; // default value for ctor, not actually used
  std::vector<std::pair<std::string, double> > result = optimize_ansatz(manager, backend, hamil, ansatz, settings, algo, n_samples, seed, eta, delta, local_params, fval);
    // throw std::runtime_error("Done\n");
  // std::vector<std::pair<std::string, double> > result = {{"test", 0.0}};
  return result;
}

PYBIND11_MODULE(nwqflow, m) {
    m.doc() = "QFlow backend based on NWQ-Sim"; // optional module docstring

    m.def("optimize_effective_hamiltonian", 
    &optimize_effective_hamiltonian, 
    "Perform single-direction gradient descentreturn the locally-optimal parameters\n"
    "\tArguments:\n"
    "\t\toperators (Iterable[(str, complex)]): List of xacc-formatted operator strings with coefficients\n"
    "\t\tnum_particles (int): Number of electrons (assumed to be equal number of alpha/beta)\n"
    "\t\tx0 (Iterable[float]): Initial parameter values\n"
    "\t\tbackend (str): NWQ-Sim backend for simulation. Defaults to \"CPU\"\n"
    "\t\txacc (bool): Use XACC operator indexing, otherwise use DUCC. Defaults to False\n"
    "\t\tseed (int): Random seed for optimizer and SPSA perturbation. Defaults to time(NULL)\n"
    "\t\tdelta (float): Magnitude of SPSA perturbation. Defaults to 1e-3\n"
    "\t\teta (float): Gradient descent stepsize. Defaults to 1e-3\n"
    "\t\tnum_trotter (int): Number of Trotter steps (linearly increases number of parameters). Defaults to 1\n"
    "\t\tnum_samples (int): Number of gradient samples for SPSA average. Defaults to 1\n",
    py::arg("operators"), 
    py::arg("num_particles"), 
    py::arg("x0"), 
    py::arg("backend") = "CPU", 
    py::arg("xacc") = false, 
    py::arg("seed") = -1, 
    py::arg("delta") = 1e-3, 
    py::arg("eta") = 1e-3, 
    py::arg("num_trotter") = 1, 
    py::arg("tnum_samples") = 1);


    m.def("get_param_count",
    [] (int num_spatial_orbitals, int num_particles) {
      int n_occ = num_particles / 2;
      int n_virt = num_spatial_orbitals - n_occ;
      // technically not the full fermi op count, but missing contribution accounts for the symmetry terms
      return n_occ * n_virt + n_occ * n_occ * n_virt * n_virt;
    },
    "Get the parameter count for given orbital/particle counts\n"
    "\tArguments:\n"
    "\t\tnum_spatial_orbitals (int): Number of spatial orbitals (must be greater than the number of particles / 2)\n"
    "\t\tnum_particles (int): Number of electrons (assumed to be equal number of alpha/beta)\n",
    py::arg("num_spatial_orbitals"), 
    py::arg("num_particles")
    );
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