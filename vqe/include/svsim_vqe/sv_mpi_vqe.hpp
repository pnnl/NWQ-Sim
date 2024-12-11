#ifndef VQE_MPI_STATE
#define VQE_MPI_STATE
#include "mpi.h"
#include "state.hpp"
#include "svsim/sv_mpi.hpp"
#include "vqe_state.hpp"
#include "utils.hpp"
#include "circuit/ansatz.hpp"
#include "observable/hamiltonian.hpp"
#include "nlopt.hpp"
#include <iostream>
#include <memory>
#include <cmath>

namespace NWQSim
{
  namespace VQE {
    /**
     * @brief  MPI-Specific backend for VQE routine
     * @note   
     * @retval None
     */
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
        size_t index = 0;
        size_t curr_ptr = 0;
        int rank = 0;
        MPI_Comm_rank(comm_global, &rank);
        process_rank = rank;
      };

      /**
       * @brief  Allocate an ObservableList object for a commuting group
       * @note   Need to free objects before termination, see destructor
       * @param  index: Index of the QWC group
       * @retval None
       */
      virtual void fill_obslist(IdxType index) override {
        ObservableList*& obs = obsvec[index];
        obs = new ObservableList;
        obs->coeffs = coeffs[index].data();
        obs->zmasks = zmasks[index].data();
        obs->numterms = zmasks[index].size();
        measurement->EXPECT(obs); 
      };
      /**
       * @brief  Call to NWQ-Sim
       * @note   
       * @param  _measurement: 
       * @retval None
       */
      virtual void call_simulator(std::shared_ptr<Ansatz> _measurement, bool reset) override { 
        Config::PRINT_SIM_TRACE = false; 
        if (iteration > 0){
          Config::PRINT_SIM_TRACE = false;
        }
        if (reset) {
          MPI_Barrier(comm_global);
          std::vector<ValType> xparams;
          if (i_proc == 0) {
            std::vector<double>* ansatz_params = ansatz->getParams();
            stat = CALL_SIMULATOR;
            for(IdxType i = 1; i < n_cpus; i++) {
              MPI_Send(ansatz_params->data(), ansatz->numParams(), MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
            }
          } else {
            xparams.resize(ansatz->numParams());
            MPI_Recv(xparams.data(),
                    ansatz->numParams(), 
                          MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            ansatz->setParams(xparams);
          }
          BARR_MPI;
          reset_state();
          BARR_MPI;
          sim(ansatz);
        }

        BARR_MPI;
        sim(_measurement);
        BARR_MPI;
      };

      virtual void get_exp_values(const std::vector<ObservableList*>& observables, std::vector<IdxType> sizes, std::vector<ValType>& output) override {
        std::vector<ValType> temp(output.size(), 0);
        for (size_t i = 0; i < observables.size(); i++) {
          for (size_t j = 0; j < sizes[i]; j++) {
            temp[i] += observables[i][j].exp_output;
          }
          // output.at(i) = observables.at(i)->exp_output;
        }
        MPI_Allreduce(temp.data(), output.data(), temp.size(), MPI_DOUBLE, MPI_SUM, comm_global);
      };
      ~SV_MPI_VQE() {

        for (auto i: obsvec) {
            delete i;
        }
      }
      virtual void call_simulator() override { 
        /**
         * @brief Prepares the trial state and computes the Hamiltonian expectation value
         * 
         * */
        if (iteration > 0){
          Config::PRINT_SIM_TRACE = false;
        }
        std::vector<ValType> xparams;
        if (i_proc == 0) {
          std::vector<double>* ansatz_params = ansatz->getParams();
          stat = CALL_SIMULATOR;
          for(IdxType i = 1; i < n_cpus; i++) {
            MPI_Send(&stat, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
            MPI_Send(ansatz_params->data(), ansatz->numParams(), MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
          }
        } else {
          xparams.resize(ansatz->numParams());
          MPI_Recv(xparams.data(),
                   ansatz->numParams(), 
                        MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          ansatz->setParams(xparams);
        }
        BARR_MPI;
        reset_state();
        BARR_MPI;
        // for (auto& i: obsvec) {
          // iexp_output = 0;
        // }
        sim(ansatz);
        BARR_MPI;
        sim(measurement);
        if (i_proc != 0) {
          stat = WAIT;
        }
        expectation_value = 0;
        for (auto i: obsvec) {
          expectation_value += i->exp_output;
        }
        double temp = expectation_value;
        expectation_value = 0;
        MPI_Allreduce(&temp, &expectation_value, 1, MPI_DOUBLE, MPI_SUM, comm_global);
        BARR_MPI;
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
          stat = WAIT;

          BARR_MPI;
          nlopt::opt optimizer = nlopt::opt(optimizer_algorithm, ansatz->numParams());
          optimizer.set_min_objective(nl_opt_function, (void*)this);
          // std::vector<double> lower_bounds(ansatz->numParams(), -2 * PI); // MZ: Why?
          // std::vector<double> upper_bounds(ansatz->numParams(), 2 * PI); // MZ: Why?
          // std::vector<double> lower_bounds(ansatz->numParams(), optimizer_settings.lbound); //MZ: my fix
          // std::vector<double> upper_bounds(ansatz->numParams(), optimizer_settings.ubound);  //MZ: my fix
          // optimizer.set_lower_bounds(lower_bounds);
          // optimizer.set_upper_bounds(upper_bounds);
          optimizer.set_lower_bounds(optimizer_settings.lbound); // MZ: use the overload since bounds for all parameters are the same (https://nlopt.readthedocs.io/en/latest/NLopt_C-plus-plus_Reference/#bound-constraints)
          optimizer.set_upper_bounds(optimizer_settings.ubound); // MZ: same as above
          // Set the termination criteria
          optimizer.set_maxeval(optimizer_settings.max_evals);
          optimizer.set_maxtime(optimizer_settings.max_time);
          optimizer.set_ftol_abs(optimizer_settings.abs_tol);
          optimizer.set_ftol_rel(optimizer_settings.rel_tol);
          optimizer.set_stopval(optimizer_settings.stop_val);
          // Set any specified optimizer parameters
          for (auto& kv_pair: optimizer_settings.parameter_map) {
              optimizer.set_param(kv_pair.first.c_str(), kv_pair.second);
          }
          iteration = 0;
          if (parameters.size() == 0) {
            parameters = std::vector<ValType>(ansatz->numParams(), 0.0);
          }
          MPI_Bcast(parameters.data(), parameters.size(), MPI_DOUBLE, 0, comm_global);
          if (i_proc == 0) {
            // nlopt::result optimization_result = optimizer.optimize(parameters, final_ene); //MZ
            optimizer.optimize(parameters, final_ene); //MZ
            opt_result = optimizer.last_optimize_result(); // MZ: this is the correct way to get the result, otherwise always give 0
            num_evals = optimizer.get_numevals(); // MZ: get numebr of function evaluations
            // energy(parameters);
            stat = EXIT_LOOP;
            for(IdxType i = 1; i < n_cpus; i++) {
              MPI_Send(&stat, 1, MPI_INT, i, 3, MPI_COMM_WORLD);
            }
          } else {
            process_loop();
          }
          parameters = ansatz->getParamRef();
          MPI_Bcast(&final_ene, 1, MPI_DOUBLE, 0, comm_global);
          BARR_MPI;
      }
      
      protected:
        std::vector<ValType> current_params;
        STATUS stat;
    };
  };
} // namespace NWQSim

#endif