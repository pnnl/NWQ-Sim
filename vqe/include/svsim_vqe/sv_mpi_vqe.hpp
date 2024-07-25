#ifndef VQE_MPI_STATE
#define VQE_MPI_STATE
#include "state.hpp"
#include "svsim/sv_mpi.hpp"
#include "vqe_state.hpp"
#include "observable/pauli_operator.hpp"
#include "utils.hpp"
#include "circuit/ansatz.hpp"
#include "circuit/measurement.hpp"
#include "gradient/sa_gradient.hpp"
#include "observable/hamiltonian.hpp"
#include "nlopt.hpp"
#include <iostream>
#include <memory>
#include <cmath>

namespace NWQSim
{
  namespace VQE {
    class SV_MPI_VQE: public VQEState, public SV_MPI {
      public:
        SV_MPI_VQE(std::shared_ptr<Ansatz> a, 
                   std::shared_ptr<Hamiltonian> h, 
                   nlopt::algorithm optimizer_algorithm,
                   Callback _callback,
                   IdxType seed = 0,
                   OptimizerSettings opt_settings = OptimizerSettings()): 
                                      SV_MPI(a->num_qubits()),
                                      VQEState(a, h, optimizer_algorithm, _callback, seed, opt_settings) {
        expvals.resize(1);
        initialize();
      };

      virtual void fill_obslist(IdxType index) override {
        ObservableList*& obs = obsvec[index];
        obs = new ObservableList;
        obs->coeffs = coeffs[index].data();
        obs->zmasks = zmasks[index].data();
        obs->exp_output = 0;
        obs->numterms = xmasks[index].size();
        ansatz->EXPECT(obs); 
      };
      virtual void call_simulator() override {  
        if (iteration > 0){
          Config::PRINT_SIM_TRACE = false;
        }
        std::vector<ValType> xparams;  
        if (iteration > 0){
          Config::PRINT_SIM_TRACE = false;
        }
        if (i_proc == 0) {
          std::vector<double>* ansatz_params = ansatz->getParams();
          stat = CALL_SIMULATOR;
          for(IdxType i = 1; i < n_cpus; i++) {
            MPI_Send(&stat, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
            MPI_Send(ansatz_params->data(), ansatz->numParams(), MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
          }
        } else {
          xparams.resize(ansatz->numParams());
          MPI_Recv(xparams.data(), ansatz->numParams(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          ansatz->setParams(xparams);
        }
        BARR_MPI;
        std::fill(expvals.begin(), expvals.end(), 0);
        reset_state();
        sim(ansatz);
        std::vector<ValType> exptemp(expvals.begin(), expvals.end());
        MPI_Reduce(exptemp.data(), expvals.data(), expvals.size(),  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (i_proc != 0) {
          stat = WAIT;
        }
      };
      virtual void process_loop() {
        assert(i_proc != 0);
        while(stat != EXIT_LOOP) {
          MPI_Recv(&stat, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          if (stat == CALL_SIMULATOR) {
            call_simulator();
            iteration++;
          }
        }

      }

      virtual void optimize(std::vector<ValType>& parameters, ValType& final_ene) override {
          iteration = 0;
          if (parameters.size() == 0) {
            parameters = std::vector<ValType>(ansatz->numParams(), 0.0);
          }
          if (i_proc == 0) {
            nlopt::result optimization_result = optimizer.optimize(parameters, final_ene);
            // energy(parameters);
            stat = EXIT_LOOP;
            for(IdxType i = 1; i < n_cpus; i++) {
              MPI_Send(&stat, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
            }

          } else {
            process_loop();
          }
      }
      
      protected:
        std::vector<ValType> current_params;
        STATUS stat;
    };
  };
} // namespace NWQSim

#endif