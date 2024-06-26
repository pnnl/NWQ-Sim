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
                   OptimizerSettings opt_settings = OptimizerSettings(),
                   MPI_Comm _comm = MPI_COMM_WORLD): 
                                      SV_CUDA_MPI(a->num_qubits(), _comm),
                                      VQEState(a, h, optimizer_algorithm, _callback, seed, opt_settings) {
        
        // Pauli term data sizes
        IdxType nterms = xmasks.size();
        obs.numterms = nterms;
        IdxType isize = nterms * sizeof(IdxType);
        IdxType vsize = nterms * sizeof(ValType);

        // Allocate space on the GPU for the Pauli expectation information
        //  For now, duplicate info on each GPU (i.e. don't use nvshmem)
        SAFE_ALOC_GPU(obs.xmasks, isize);
        SAFE_ALOC_GPU(obs.zmasks, isize);
        SAFE_ALOC_GPU(obs.x_index_sizes, isize);
        SAFE_ALOC_GPU(obs.x_indices, x_indices.size() * sizeof(IdxType));
        SAFE_ALOC_GPU(obs.exp_output, vsize);

        // Copy Pauli masks and indices from host to device
        cudaSafeCall(cudaMemcpy(obs.xmasks, xmasks.data(), isize,
                                    cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(obs.zmasks, zmasks.data(), isize,
                                    cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(obs.x_index_sizes, x_index_sizes.data(), isize,
                                    cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(obs.x_indices, x_indices.data(), x_indices.size() * sizeof(IdxType),
                                    cudaMemcpyHostToDevice));
        ObservableList* obs_device;
        SAFE_ALOC_GPU(obs_device, sizeof(ObservableList));
        cudaSafeCall(cudaMemcpy(obs_device, &obs, sizeof(ObservableList),
                                    cudaMemcpyHostToDevice));
        ansatz->EXPECT(obs_device);
        int size;
        MPI_Comm_size(comm_global, &size);
        n_cpus = size;

      };
      virtual void call_simulator() override {  
        std::vector<ValType> xparams;   
        if (iteration > 0){
          Config::PRINT_SIM_TRACE = false;
        }
        if (i_proc == 0) {
          std::vector<double>* ansatz_params = ansatz->getParams();
          stat = CALL_SIMULATOR;
          for(IdxType i = 1; i < n_gpus; i++) {
            MPI_Send(&stat, 1, MPI_INT, i, iteration, comm_global);
            MPI_Send(ansatz_params->data(), ansatz->numParams(), MPI_DOUBLE, i, 1, comm_global);
          }
        } else {
          xparams.resize(ansatz->numParams());
          MPI_Recv(xparams.data(), ansatz->numParams(), MPI_DOUBLE, 0, 1, comm_global, MPI_STATUS_IGNORE);
          ansatz->setParams(xparams);
        }
        MPI_Barrier(comm_global);
        reset_state();
        sim(ansatz);
        cudaDeviceSynchronize();
        cudaSafeCall(cudaMemcpy(expvals.data(), obs.exp_output, expvals.size() * sizeof(ValType),
                                    cudaMemcpyDeviceToHost));
        MPI_Barrier(comm_global);
        if (i_proc != 0) {
          MPI_Reduce(expvals.data(), expvals.data(), expvals.size(), MPI_DOUBLE, MPI_SUM, 0, comm_global);
          stat = WAIT;
        } else {
          std::vector<ValType> temp(expvals.size());
          MPI_Reduce(expvals.data(), temp.data(), expvals.size(), MPI_DOUBLE, MPI_SUM, 0, comm_global);
          std::copy(temp.begin(), temp.end(), expvals.begin());
        }
        MPI_Barrier(comm_global);              
        
      };
      virtual std::vector<std::pair<std::string, ValType>> follow_fixed_gradient(const std::vector<ValType>& x0, 
                                                                                 ValType& final_ene, 
                                                                                 ValType delta, 
                                                                                 ValType eta, 
                                                                                 IdxType n_grad_est) override {
        Config::PRINT_SIM_TRACE = false;
        if (i_proc != 0) {
          process_loop();
        }
        else {
          std::vector<ValType> gradient (x0.size(),1.0);
          std::vector<ValType> params(x0);
          std::vector<ValType> minima_params(x0);
          ValType ene_prev = MAXFLOAT;
          ValType ene_curr = energy(params);
          // gradient
          // get the single-direction starting vector
          g_est.estimate([&] (const std::vector<double>& xval) { return energy(xval);}, params, gradient, delta, n_grad_est);
          IdxType step = 0;
          // for (auto& i: gradient) {
            // std::cout << i <<  " ";
          // }
          // std::cout << std::endl;
          // auto s1 = std::chrono::high_resolution_clock::now();
          // follow the starting vector until we hit a global minimum
          do {
            for (size_t i = 0; i < params.size(); i++) {
              params[i] -= eta * gradient[i];
            }
            // auto s1 =  std::chrono::high_resolution_clock::now();
            // ene_curr = 0;
            ene_curr = energy(params);
            // std::cout << step << " " << ene_curr << " " << ene_prev << std::endl;
            if (ene_curr >= ene_prev) {
              for (size_t i = 0; i < params.size(); i++) {
                params[i] += eta * gradient[i];
              }
              break;
            } else {
              ene_prev = ene_curr;
            }
            step++;
            // auto s2 =  std::chrono::high_resolution_clock::now();
            // std::cout << (s2-s1).count()/1e9 << std::endl;
          } while(true);
          // std::cout << "Ended loop\n" << std::endl;
          stat = EXIT_LOOP;
            for(IdxType i = 1; i < n_cpus; i++) {
              MPI_Send(&stat, 1, MPI_INT, i, 3, comm_global);
            }
        }
        std::vector<std::pair<std::string, ValType>> result = ansatz->getFermionicOperatorParameters();
        return result;
      }
      virtual void process_loop() {
        assert(i_proc != 0);
        while(stat != EXIT_LOOP) {
          MPI_Recv(&stat, 1, MPI_INT, 0, iteration, comm_global, MPI_STATUS_IGNORE);
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
            
            // std::cout << energy(parameters) << std::endl;iteration = 1;
            // std::cout << zmasks[1376] << std::endl;
            stat = EXIT_LOOP;
            for(IdxType i = 1; i < n_cpus; i++) {
              MPI_Send(&stat, 1, MPI_INT, i, iteration, comm_global);
            }

          } else {
            process_loop();
          }
      }
      ~SV_CUDA_MPI_VQE()
        {
            SAFE_FREE_GPU(obs.xmasks);
            SAFE_FREE_GPU(obs.zmasks);
            SAFE_FREE_GPU(obs.x_index_sizes);
            SAFE_FREE_GPU(obs.x_indices);
            SAFE_FREE_GPU(obs.exp_output);
            SAFE_FREE_GPU(ansatz->gates->back().data);
            // SAFE_FREE_GPU(randoms_gpu);
            // SAFE_FREE_GPU(gates_gpu);
            // SAFE_FREE_HOST_CUDA(randoms);
            // SAFE_FREE_HOST_CUDA(results);
        }
    protected:
        IdxType n_cpus;
        STATUS stat; 
   
  };
}; // namespace NWQSim
};
#endif