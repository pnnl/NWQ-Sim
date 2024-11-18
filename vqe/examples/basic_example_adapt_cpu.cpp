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
#include "circuit/dynamic_ansatz.hpp"
#include "vqe_adapt.hpp"


using namespace NWQSim::VQE;


// Callback function, requires signature (void*) (const std::vector<NWQSim::ValType>&, NWQSim::ValType, NWQSim::IdxType)
double last = 0.0;
void callback_function(const std::vector<NWQSim::ValType>& x, NWQSim::ValType fval, NWQSim::IdxType iteration) {
  
}

using Pauli = NWQSim::VQE::PauliOperator;
// H2O Example
int main(int argc, char** argv) {
  NWQSim::IdxType n_particles = 6; // Set the number of particles
  // Note: path relative to presumed build directory
  std::string hamiltonian_path = "../vqe/example_hamiltonians/BeH2_6_2.4_xacc.hamil"; //  Effective Hamiltonian file path
  std::vector<std::vector<Pauli> > op = {{Pauli("IIIIIIIIIIII")}};
  std::shared_ptr<Hamiltonian> hamil = std::make_shared<Hamiltonian>(op, n_particles, true); // Build the Hamiltonian object (used for energy calculation)
  Transformer jw_transform = getJordanWignerTransform; // Choose a transformation function

  // Build the ansatz circuit using the Hamiltonian Molecular environment and JW mapping
  //      (shared_ptr used to match baseline NWQ-Sim functionality)
  std::shared_ptr<Ansatz> ansatz = std::make_shared<NWQSim::VQE::DynamicAnsatz>(hamil->getEnv());
  ansatz->buildAnsatz();
  // Build the Quantum State object
  OptimizerSettings osets;
  osets.max_evals = 1;
  std::shared_ptr<NWQSim::VQE::VQEState> state = std::make_shared<NWQSim::VQE::SV_CPU_VQE>(ansatz, // reference to ansatz
                                hamil,  // reference to Hamiltonian
                                nlopt::algorithm::LN_COBYLA, // NLOpt algorithm for optimization
                                callback_function, // Callback function for each energy evaluation,
                                0, // Random seed (passed to the SPSA gradient estimator for random perturbations)
                                osets);
  // Random initial parameters (sets to all 0 by default if not provided). Note that the function modifies `parameters` inplace
  std::vector<double> parameters(ansatz->numParams(), 1.0);
  std::uniform_real_distribution<double> dist(0.0, 2 * PI);
  std::mt19937_64 rand_device(342);
  std::generate(parameters.begin(), parameters.end(), [&] () {return dist(rand_device);});

  // Return destination for the function value
  double fval;

  // Start the VQE optimization   
  std::shared_ptr<NWQSim::VQE::DynamicAnsatz> dyn_ansatz = std::reinterpret_pointer_cast<NWQSim::VQE::DynamicAnsatz>(ansatz);
    dyn_ansatz->make_op_pool(hamil->getTransformer(), 0, 40);

    NWQSim::VQE::AdaptVQE adapt_instance(dyn_ansatz, state, hamil);
    adapt_instance.make_commutators();
    adapt_instance.optimize(parameters, fval, 100);
  std::cout << "Final Energy: " << fval << std::endl;
  auto param_map = ansatz->getFermionicOperatorParameters();
  for (auto pair: param_map) {
    std::cout << pair.first << " " << pair.second << std::endl;
  }
  return 0;
}