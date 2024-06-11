#include "transform/transform.hpp"
#include "observable/fermionic_operator.hpp"
#include "observable/pauli_operator.hpp"
#include "circuit/ansatz.hpp"
#include "circuit/measurement.hpp"
#include <fstream>
#include "utils.hpp"
#include "state.hpp"
#include "svsim_vqe/sv_cuda_mpi_vqe.cuh"
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
  int i_proc;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &i_proc);
  NWQSim::IdxType n_particles = 10; // Set the number of particles
  // Note: path relative to presumed build directory
  std::string hamiltonian_path = "../vqe/example_hamiltonians/h2O.hamil"; //  Effective Hamiltonian file path

  Hamiltonian hamil(hamiltonian_path, n_particles); // Build the Hamiltonian object (used for energy calculation)
  Transformer jw_transform = getJordanWignerTransform; // Choose a transformation function

  // Build the ansatz circuit using the Hamiltonian Molecular environment and JW mapping
  //      (shared_ptr used to match baseline NWQ-Sim functionality)
  std::shared_ptr<Ansatz> ansatz = std::make_shared<UCCSD>(hamil.getEnv(), jw_transform, 1);
  // Build the Quantum State object
  NWQSim::VQE::SV_CUDA_MPI_VQE state(ansatz, // reference to ansatz
                                hamil,  // reference to Hamiltonian
                                nlopt::algorithm::LN_COBYLA, // NLOpt algorithm for optimization
                                callback_function, // Callback function for each energy evaluation
                                0 // Random seed (passed to the SPSA gradient estimator for random perturbations)
                                );
  // Random initial parameters (sets to all 0 by default if not provided). Note that the function modifies `parameters` inplace
  std::vector<double> parameters(ansatz->numParams(), 1.0);
  if (i_proc == 0) {
    std::uniform_real_distribution<double> dist(0.0, 2 * PI);
    std::mt19937_64 rand_device(342);
    std::generate(parameters.begin(), parameters.end(), [&] () {return dist(rand_device);});
  }

  // Return destination for the function value
  double fval;

  // Start the VQE optimization
  state.optimize(parameters, fval);
  if (i_proc == 0) {
    std::cout << "Final Parameters: " << parameters << std::endl;
    std::cout << "Final Energy: " << fval << std::endl;
  }
  MPI_Finalize();
  return 0;
}