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
        
        // Pauli term data sizes
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        n_cpus = size;        
        expvals.resize(1);
        SAFE_ALOC_GPU(expvals_dev, expvals.size()* sizeof(ValType));
        initialize();

      };
      virtual void fill_obslist(IdxType index) override {
        ObservableList& obs = obsvec[index];
        
        obs.numterms = zmasks[index].size();
        IdxType isize = obs.numterms * sizeof(IdxType);
        IdxType vsize = obs.numterms * sizeof(ValType);
        SAFE_ALOC_GPU(obs.zmasks, isize);
        obs.exp_output = 0;
        SAFE_ALOC_GPU(obs.coeffs, vsize);
        cudaSafeCall(cudaMemcpy(obs.zmasks, zmasks[index].data(), isize,
                                    cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(obs.coeffs, coeffs[index].data(), coeffs[index].size() * sizeof(ValType),
                                    cudaMemcpyHostToDevice));
        ObservableList* obs_device;
        SAFE_ALOC_GPU(obs_device, sizeof(ObservableList));
        cudaSafeCall(cudaMemcpy(obs_device, &obs, sizeof(ObservableList),
                                    cudaMemcpyHostToDevice));
        ansatz->EXPECT(obs_device); 
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
            MPI_Send(&stat, 1, MPI_INT, i, iteration, MPI_COMM_WORLD);
            MPI_Send(ansatz_params->data(), ansatz->numParams(), MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
          }
        } else {
          xparams.resize(ansatz->numParams());
          MPI_Recv(xparams.data(), ansatz->numParams(), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          ansatz->setParams(xparams);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        cudaSafeCall(cudaMemset(expvals_dev, 0, expvals.size() * sizeof(ValType)));
        reset_state();
        sim(ansatz);
        cudaDeviceSynchronize();
        cudaSafeCall(cudaMemcpy(expvals.data(), expvals_dev, expvals.size() * sizeof(ValType),
                                    cudaMemcpyDeviceToHost));
        MPI_Barrier(MPI_COMM_WORLD);
        std::vector<ValType> temp(expvals.begin(), expvals.end());
          MPI_Reduce(temp.data(), expvals.data(), expvals.size(), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (i_proc != 0) {
          stat = WAIT;
        }else {
        
        MPI_Barrier(MPI_COMM_WORLD);              
        
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
              MPI_Send(&stat, 1, MPI_INT, i, iteration, MPI_COMM_WORLD);
            }

          } else {
            process_loop();
          }
      }
      ~SV_CUDA_MPI_VQE()
        {
          for (auto& obs: obsvec) {
              SAFE_FREE_GPU(obs.zmasks);
              SAFE_FREE_GPU(obs.coeffs);
          }
            SAFE_FREE_GPU(expvals_dev);
          for (auto g: *ansatz->gates) {
              if (g.op_name == OP::EXPECT) {
                  SAFE_FREE_GPU(g.data);
              }
          }
            
            // SAFE_FREE_GPU(randoms_gpu);
            // SAFE_FREE_GPU(gates_gpu);
            // SAFE_FREE_HOST_CUDA(randoms);
            // SAFE_FREE_HOST_CUDA(results);
        }
    protected:
        IdxType n_cpus;
        STATUS stat; 
        ValType* expvals_dev;
   
  };
}; // namespace NWQSim
};
#endif