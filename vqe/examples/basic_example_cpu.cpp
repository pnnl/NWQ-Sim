#include "transform/transform.hpp"
#include "observable/fermionic_operator.hpp"
#include "observable/pauli_operator.hpp"
#include "circuit/ansatz.hpp"
#include "circuit/measurement.hpp"
#include <fstream>
#include "utils.hpp"
#include "state.hpp"
#include "svsim_vqe/sv_cpu_vqe.hpp"
#include "gradient/sa_gradient.hpp"


using namespace NWQSim::VQE;


// Callback function, requires signature (void*) (const std::vector<NWQSim::ValType>&, NWQSim::ValType, NWQSim::IdxType)
double last = 0.0;
void callback_function(const std::vector<NWQSim::ValType>& x, NWQSim::ValType fval, NWQSim::IdxType iteration) {
  double relchange = abs(last) > 1e-6 ? 0.0 : (fval - last) / last;
  std::cout << "Iteration " << iteration << ", fval = " << fval << ", relchange " << relchange << std::endl;
  last = fval;
}


// H2O Example
int main(int argc, char** argv) {
  NWQSim::IdxType n_particles = 4; // Set the number of particles
  // Note: path relative to presumed build directory
  std::string hamiltonian_path = "../vqe/example_hamiltonians/H4_4_0.9_xacc.hamil"; //  Effective Hamiltonian file path

  std::shared_ptr<Hamiltonian> hamil = std::make_shared<Hamiltonian>(hamiltonian_path, n_particles, true); // Build the Hamiltonian object (used for energy calculation)
  Transformer jw_transform = getJordanWignerTransform; // Choose a transformation function

  // Build the ansatz circuit using the Hamiltonian Molecular environment and JW mapping
  //      (shared_ptr used to match baseline NWQ-Sim functionality)
  std::shared_ptr<Ansatz> ansatz = std::make_shared<UCCSD>(hamil->getEnv(), jw_transform, 1);
  ansatz->buildAnsatz();
  std::string config_path = "../default_config.json";
  // Build the Quantum State object
  NWQSim::VQE::SV_CPU_VQE state(ansatz, // reference to ansatz
                                hamil,  // reference to Hamiltonian
                                nlopt::algorithm::LN_COBYLA, // NLOpt algorithm for optimization
                                callback_function, // Callback function for each energy evaluation
                                config_path, // path to config file
                                0 // Random seed (passed to the SPSA gradient estimator for random perturbations)
                                );
  // All zero starting point. Note that the function modifies `parameters` inplace
  std::vector<double> parameters(ansatz->numParams(), 0.0);

  // Return destination for the function value
  double fval;

  // Start the VQE optimization
  state.optimize(parameters, fval);
  std::cout << "Final Parameters: " << parameters << std::endl;
  std::cout << "Final Energy: " << fval << std::endl;

  return 0;
}