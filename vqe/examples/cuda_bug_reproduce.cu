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


/**
    * Steps to reproduce bug on Perlmutter
    * 1. salloc --nodes 1 --qos interactive -t 60 --constraint gpu --account=m4243
    * 2. cd NWQ-Sim/build
    * 3. source ../environment/setup_perlmutter.sh
    * 4. cmake ..;make -j8
    * 5. srun -N 1 -n <n_gpus> -c 1 --gpus-per-task=1 --gpu-bind=single:1 vqe/examples/cuda_bug_reproduce
    In all cases, the expected output includes:
    `Identity expectation value (should be 1.0): 1`
    
    For n_gpus=1, this is the case. However, we get the following output for n_gpus=2:
    `Identity expectation value (should be 1.0): 0.976758`
    and n_gpus=4:
    `Identity expectation value (should be 1.0): 0.87145`
    

*/
void example () {
  int i_proc;

  MPI_Comm_rank(MPI_COMM_WORLD, &i_proc);
  NWQSim::IdxType n_particles = 14; // Set the number of particles
  NWQSim::IdxType n_spatial = 10; // Set the number of active spaces
  // Test case: From N_2 molecule environment, just testing the identity operator. Expectation should be 1.0 regardless of circuit 
  NWQSim::VQE::MolecularEnvironment env (n_spatial, n_particles, true);
  std::vector<std::vector<PauliOperator> > pauli_operators = {{PauliOperator("IIIIIIIIIIIIIIIIIIII")}};
  std::shared_ptr<Hamiltonian> hamil = std::make_shared<Hamiltonian>(env, pauli_operators); // Build the Hamiltonian object (used for energy calculation)
  Transformer jw_transform = getJordanWignerTransform; // Choose a transformation function

  // Build the ansatz circuit using the Hamiltonian Molecular environment and JW mapping
  //      (shared_ptr used to match baseline NWQ-Sim functionality)
  std::shared_ptr<Ansatz> ansatz = std::make_shared<UCCSD>(hamil->getEnv(), jw_transform, 1);
  OptimizerSettings settings;
    settings.max_evals = 1; // only run one eval for test
  // Build the Quantum State object
  NWQSim::VQE::SV_CUDA_MPI_VQE state(ansatz, // reference to ansatz
                                hamil,  // reference to Hamiltonian
                                nlopt::algorithm::LN_COBYLA, // NLOpt algorithm for optimization
                                callback_function, // Callback function for each energy evaluation
                                0, // Random seed (passed to the SPSA gradient estimator for random perturbations)
                                settings);
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
    std::cout << "Identity expectation value (should be 1.0): " << fval << std::endl;
  }
}

// H2O Example
int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  example();
  MPI_Finalize();
  return 0;
}