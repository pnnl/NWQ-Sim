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
    class SV_CUDA_MPI_VQE: public VQEState, public SV_CUDA_MPI {
      public:
        SV_CUDA_MPI_VQE(std::shared_ptr<Ansatz> a, 
                   std::shared_ptr<Hamiltonian> h, 
                   nlopt::algorithm optimizer_algorithm,
                   Callback _callback,
                   IdxType seed = 0,
                   OptimizerSettings opt_settings = OptimizerSettings()): 
                                      SV_CUDA_MPI(a->num_qubits()),
                                      VQEState(a, h, optimizer_algorithm, _callback, seed, opt_settings) {
        iteration = 0;
        // Pauli term data sizes
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        n_cpus = size;    
        int rank = 0;
        MPI_Comm_rank(comm_global, &rank);
        process_rank = rank;
      };
      virtual void allocate_observables(IdxType size) override {
        SAFE_ALOC_GPU(obsvals_dev, size * sizeof(ObservableList));
        obsvals.resize(size);
        obsvec.resize(size);
      };
virtual void fill_obslist(IdxType index) override {
        /**
         * @brief Add an entry to the ObservableList data structure
         */
        ObservableList obs;
        obs.numterms = zmasks[index].size();
        IdxType isize = obs.numterms * sizeof(IdxType);
        IdxType vsize = obs.numterms * sizeof(ValType);
        SAFE_ALOC_GPU(obs.zmasks, isize);
        obs.exp_output = 0;
        SAFE_ALOC_GPU(obs.coeffs, vsize);
        // allocate memory for the zmasks and coefficients
        cudaSafeCall(cudaMemcpy(obs.zmasks, zmasks[index].data(), isize,
                                cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(obs.coeffs, coeffs[index].data(), vsize,
                                cudaMemcpyHostToDevice));
          // std::cout << "Filling list: " << index << " " <<(obsvals_dev + index) << std::endl;
        ObservableList* obs_device = (ObservableList*)((void*)obsvals_dev + index * sizeof(ObservableList));
        obsvals[index] = obs;
        obsvec[index] = obs_device;
        measurement->EXPECT(obs_device); 
      };
      virtual void set_exp_gate(std::shared_ptr<Ansatz> circuit, ObservableList* o, std::vector<IdxType>& _zmasks, std::vector<ValType>& _coeffs) override {
        ObservableList obs;
        
        obs.numterms = _zmasks.size();
        IdxType isize = obs.numterms * sizeof(IdxType);
        IdxType vsize = obs.numterms * sizeof(ValType);
        SAFE_ALOC_GPU(obs.zmasks, isize);
        obs.exp_output = 0;
        SAFE_ALOC_GPU(obs.coeffs, vsize);
        cudaSafeCall(cudaMemcpy(obs.zmasks, _zmasks.data(), isize,
                                cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(obs.coeffs, _coeffs.data(), _coeffs.size() * sizeof(ValType),
                                cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(o, &obs, sizeof(ObservableList),
                                cudaMemcpyHostToDevice));
        circuit->EXPECT(o);
      };
      
      virtual void call_simulator() override {  
        cudaSafeCall(cudaMemcpy(obsvals_dev, obsvals.data(), obsvals.size() * sizeof(ObservableList),
                              cudaMemcpyHostToDevice));
        std::vector<ValType> xparams;   
        if (iteration > 0){
          Config::PRINT_SIM_TRACE = false;
        }
        if (i_proc == 0) {
          std::vector<double>* ansatz_params = ansatz->getParams();
          stat = CALL_SIMULATOR;
          for(IdxType i = 1; i < n_gpus; i++) {
            MPI_Send(&stat, 1, MPI_INT, i, iteration, MPI_COMM_WORLD);
            MPI_Send(ansatz_params->data(), ansatz->numParams(), MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
          }
        } else {
          xparams.resize(ansatz->numParams());
          MPI_Recv(xparams.data(), ansatz->numParams(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          ansatz->setParams(xparams);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        reset_state();
        sim(ansatz);
        sim(measurement);
        cudaDeviceSynchronize();

        cudaSafeCall(cudaMemcpy(obsvals.data(), obsvals_dev, obsvals.size() * sizeof(ObservableList),
                                    cudaMemcpyDeviceToHost));
        cudaDeviceSynchronize();
        // for (size_t i = 0; i < dim; i++) {
        //   std::cout << "(" << sv_real_cpu[i] * sv_real_cpu[i] + sv_imag_cpu[i] * sv_imag_cpu[i] << "), ";
        // }
        // std::cout << std::endl;
        double temp = 0;
        size_t index = 0;
        for (auto o: obsvals) {
          temp += o.exp_output;
          index++;
        }
        MPI_Allreduce(&temp, &expectation_value, 1, MPI_DOUBLE, MPI_SUM, comm_global);
        if (i_proc != 0) {
          stat = WAIT;
        }
        MPI_Barrier(MPI_COMM_WORLD);     
        
      };
      virtual void call_simulator(std::shared_ptr<Ansatz> _measurement, bool reset) override {  
        std::vector<ValType> xparams;   
        if (iteration > 0){
          Config::PRINT_SIM_TRACE = false;
        }
        if (reset) {
          if (i_proc == 0) {
            std::vector<double>* ansatz_params = ansatz->getParams();
            stat = CALL_SIMULATOR;
            for(IdxType i = 1; i < n_gpus; i++) {
              MPI_Send(ansatz_params->data(), ansatz->numParams(), MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
            }
          } else {
            xparams.resize(ansatz->numParams());
            MPI_Recv(xparams.data(), ansatz->numParams(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            ansatz->setParams(xparams);
          }
          MPI_Barrier(MPI_COMM_WORLD);
          reset_state();
          sim(ansatz);
        }
        sim(_measurement);
        cudaDeviceSynchronize();
        
      };
      virtual void process_loop() {
        assert(i_proc != 0);
        while(stat != EXIT_LOOP) {
          MPI_Recv(&stat, 1, MPI_INT, 0, iteration, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          if (stat == CALL_SIMULATOR) {
            call_simulator();
            iteration++;
          }
        }

      }
      virtual void optimize(std::vector<ValType>& parameters, ValType& final_ene) override {
          nlopt::opt optimizer = nlopt::opt(optimizer_algorithm, ansatz->numParams());
          optimizer.set_min_objective(nl_opt_function, (void*)this);
          // std::vector<double> lower_bounds(ansatz->numParams(), -2 * PI); // MZ: this is 100% error
          // std::vector<double> upper_bounds(ansatz->numParams(), 2 * PI); // MZ: this is 100% error
          // std::vector<double> lower_bounds(ansatz->numParams(), optimizer_settings.lbound); // MZ: need someone test it
          // std::vector<double> upper_bounds(ansatz->numParams(), optimizer_settings.ubound); // MZ: need someone test it
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
          if (i_proc == 0) {
            // nlopt::result optimization_result = optimizer.optimize(parameters, final_ene); // MZ
            optimizer.optimize(parameters, final_ene); // MZ
            opt_result = optimizer.last_optimize_result(); // MZ: the correct way to get the result of the optimization, otherwise gives 0
            num_evals = optimizer.get_numevals(); // MZ: get numebr of function evaluations
            
            // std::cout << energy(parameters) << std::endl;iteration = 1;
            // std::cout << zmasks[1376] << std::endl;
            stat = EXIT_LOOP;
            for(IdxType i = 1; i < n_cpus; i++) {
              MPI_Send(&stat, 1, MPI_INT, i, iteration, MPI_COMM_WORLD);
            }

          } else {
            stat = WAIT;
            process_loop();
          }
      }
      ~SV_CUDA_MPI_VQE()
        {
          cudaSafeCall(cudaMemcpy(obsvals.data(), obsvals_dev, obsvals.size() * sizeof(ObservableList),
                                    cudaMemcpyDeviceToHost));
          for (auto& obs: obsvals) {
              SAFE_FREE_GPU(obs.zmasks);
              SAFE_FREE_GPU(obs.coeffs);
          }
          SAFE_FREE_GPU(obsvals_dev);
            
            // SAFE_FREE_GPU(randoms_gpu);
            // SAFE_FREE_GPU(gates_gpu);
            // SAFE_FREE_HOST_CUDA(randoms);
            // SAFE_FREE_HOST_CUDA(results);
        }
         virtual void get_exp_values(const std::vector<ObservableList*>& observables, std::vector<IdxType> sizes, std::vector<ValType>& output) override {
          std::vector<double> temp(output.size(), 0);
          for (size_t i = 0; i < observables.size(); i++) {
            std::vector<ObservableList> host_temp(sizes[i]);
            cudaSafeCall(cudaMemcpy(host_temp.data(), observables[i], sizes[i] * sizeof(ObservableList),
                                      cudaMemcpyDeviceToHost));
            
            for (ObservableList obs: host_temp) {
              temp[i] += obs.exp_output;
            }
            // output.at(i) = observables.at(i)->exp_output;
          }
          MPI_Allreduce(temp.data(), output.data(), temp.size(), MPI_DOUBLE, MPI_SUM, comm_global);
      };

      virtual void allocate_observables(ObservableList*& observables, IdxType size) override {
        SAFE_ALOC_GPU(observables, size * sizeof(ObservableList));
      };
    protected:
        IdxType n_cpus;
        STATUS stat; 
        ObservableList* obsvals_dev;
        std::vector<ObservableList> obsvals;
   
  };
}; // namespace NWQSim
};
#endif