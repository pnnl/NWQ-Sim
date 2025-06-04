
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
  std::string hamiltonian_path = "../vqe/example_hamiltonians/LiH_3_xacc.hamil"; //  Effective Hamiltonian file path

  std::shared_ptr<Hamiltonian> hamil = std::make_shared<Hamiltonian>(hamiltonian_path, n_particles, false); // Build the Hamiltonian object (used for energy calculation) using the DUCC mapping
  Transformer jw_transform = getJordanWignerTransform; // Choose a transformation function

  // Build the ansatz circuit using the Hamiltonian Molecular environment and JW mapping
  //      (shared_ptr used to match baseline NWQ-Sim functionality)
  std::shared_ptr<Ansatz> ansatz = std::make_shared<UCCSD>(hamil->getEnv(), jw_transform, 1);
  ansatz->buildAnsatz();
  
  // Now we get a bit fancier. We can pass an `OptimizerSettings` object to specify termination criteria and optimizer parameters
  OptimizerSettings settings;
  settings.abs_tol = 1e-5;  // absolution function value tolerance cutoff
  settings.rel_tol = 1e-3;  // relative function value tolerance cutoff
  settings.max_evals = 100; // max number of function calls (circuit simulations)
  settings.max_time = 1000; // timeout (in seconds)
  settings.stop_val = -76.387; // ground state for our problem
  // algorithm-specific parameters, see NLOpt docs for specifics. `inner_maxeval` specifies the maximum number
  //    of function evalutions for the inner loop of the Method of Moving Averages (MMA) and Conservative Convex Separable Approximation (CCSA)
  //    gradient-based algorithms
  settings.parameter_map = {
    {"inner_maxeval", 20}
  }; 
  
  
  // Build the Quantum State object
  NWQSim::VQE::SV_CPU_VQE state(ansatz, // reference to ansatz
                                hamil,  // reference to Hamiltonian
                                nlopt::algorithm::LD_MMA, // NLOpt algorithm for optimization
                                callback_function, // Callback function for each energy evaluation
                                0, // Random seed (passed to the SPSA gradient estimator for random perturbations)
                                settings
                                );
  // Random initial parameters (sets to all 0 by default if not provided). Note that the function modifies `parameters` inplace
  std::vector<double> parameters(ansatz->numParams(), 1.0);
  std::uniform_real_distribution<double> dist(0.0, 2 * PI);
  std::mt19937_64 rand_device(342);
  std::generate(parameters.begin(), parameters.end(), [&] () {return dist(rand_device);});

  // Return destination for the function value
  double fval;

  // Start the VQE optimization
  state.optimize(parameters, fval);
  std::cout << "Final Parameters: " << parameters << std::endl;
  std::cout << "Final Energy: " << fval << std::endl;

  return 0;
}