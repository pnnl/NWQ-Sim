#ifndef VQE_CUDA_MPI_STATE
#define VQE_CUDA_MPI_STATE
#include "svsim/sv_cuda_mpi.cuh"
#include "vqe_state.hpp"
#include "observable/pauli_operator.hpp"
#include "utils.hpp"
#include "circuit/ansatz.hpp"
#include "circuit/measurement.hpp"
#include "gradient/sa_gradient.hpp"
#include "observable/hamiltonian.hpp"
#include "nlopt.hpp"
#include <memory>
#include <cmath>

namespace NWQSim
{
  namespace VQE {
    class SV_CUDA_VQE: public VQEState, public SV_CUDA_MPI {
      public:
        SV_CUDA_VQE(std::shared_ptr<Ansatz> a, 
                   const Hamiltonian& h, 
                   nlopt::algorithm optimizer_algorithm,
                   Callback _callback,
                   IdxType seed = 0,
                   OptimizerSettings opt_settings = OptimizerSettings()): 
                                      SV_CUDA_MPI(a->num_qubits()),
                                      VQEState(a, h, optimizer_algorithm, _callback, seed, opt_settings) {};
      virtual void call_simulator(std::shared_ptr<Ansatz> ansatz) override {        
        
        reset_state();
        sim(ansatz);
      };
   
  };
}; // namespace NWQSim
};
#endif