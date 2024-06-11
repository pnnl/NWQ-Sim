#ifndef VQE_MPI_STATE
#define VQE_MPI_STATE
#include "svsim/sv_mpi.hpp"
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
  enum STATUS {
    CALL_SIMULATOR,
    WAIT,
    EXIT_LOOP
  };
  namespace VQE {
    class SV_MPI_VQE: public VQEState, public SV_MPI {
      public:
        SV_MPI_VQE(std::shared_ptr<Ansatz> a, 
                   const Hamiltonian& h, 
                   nlopt::algorithm optimizer_algorithm,
                   Callback _callback,
                   IdxType seed = 0,
                   OptimizerSettings opt_settings = OptimizerSettings()): 
                                      SV_MPI(a->num_qubits()),
                                      VQEState(a, h, optimizer_algorithm, _callback, seed, opt_settings) {};
      virtual void call_simulator(std::shared_ptr<Ansatz> _ansatz) override {  
        std::vector<ValType> xparams;  
        if (i_proc == 0) {
          const std::vector<ValType>* ansatz_params = _ansatz->getParams();
          stat = CALL_SIMULATOR;
          for(IdxType i = 1; i < n_cpus; i++) {
            MPI_Send(&stat, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(ansatz_params, _ansatz->numParams(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
          }
        } else {
          xparams.resize(_ansatz->numParams());
          MPI_Recv(xparams.data(), _ansatz->numParams(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          _ansatz->setParams(xparams);
        }
        reset_state();
        sim(_ansatz);
        if (i_proc != 0) {
          stat = WAIT;
        }
      };
      virtual void process_loop(std::shared_ptr<Ansatz> _ansatz) {
        assert(i_proc != 0);
        while(stat != EXIT_LOOP) {
          MPI_Recv(&stat, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          if (stat == CALL_SIMULATOR) {
            call_simulator(_ansatz);
          }
        }

      }
      virtual void optimize(std::vector<ValType>& parameters, ValType& final_ene) {
          iteration = 0;
          Config::PRINT_SIM_TRACE = false;
          if (parameters.size() == 0) {
            parameters = std::vector<ValType>(ansatz->numParams(), 0.0);
          }
          if (i_proc == 0) {
            nlopt::result optimization_result = optimizer.optimize(parameters, final_ene);
            stat = EXIT_LOOP;
            for(IdxType i = 1; i < n_cpus; i++) {
              MPI_Send(&stat, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }

          } else {
            process_loop(ansatz);
          }
      }
      
      protected:
        std::vector<ValType> current_params;
        STATUS stat;
    };
  };
} // namespace NWQSim

#endif