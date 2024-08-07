#ifndef VQE_CUDA_STATE
#define VQE_CUDA_STATE
#include "state.hpp"
#include "svsim/sv_cuda.cuh"
#include "vqe_state.hpp"
#include "observable/pauli_operator.hpp"
#include "utils.hpp"
#include "circuit/ansatz.hpp"
#include "circuit/measurement.hpp"
#include "gradient/sa_gradient.hpp"
#include "observable/hamiltonian.hpp"
#include "nlopt.hpp"
#include "private/cuda_util.cuh"
#include <memory>
#include <cmath>

namespace NWQSim
{
  namespace VQE {
    class SV_CUDA_VQE: public VQEState, public SV_CUDA {
      public:
        SV_CUDA_VQE(std::shared_ptr<Ansatz> a, 
                   std::shared_ptr<Hamiltonian> h, 
                   nlopt::algorithm optimizer_algorithm,
                   Callback _callback,
                   const std::string& configpath,
                   IdxType seed = 0,
                   OptimizerSettings opt_settings = OptimizerSettings()): 
                                      SV_CUDA(a->num_qubits(), configpath),
                                      VQEState(a, h, optimizer_algorithm, _callback, seed, opt_settings) {
        // resize the expectation value data structure (to-do, just make a double)
        initialize(); // initialize ansatz, measurement circuit
        // allocate space for host, device memory for observable lists
        SAFE_ALOC_GPU(obsvals_dev, zmasks.size() * sizeof(ObservableList));
        obsvals_host = new ObservableList[zmasks.size()];

      };
      virtual void fill_obslist(IdxType index) override {
        ObservableList obs;
        
        obs.numterms = xmasks[index].size();
        IdxType isize = obs.numterms * sizeof(IdxType);
        IdxType vsize = obs.numterms * sizeof(ValType);
        SAFE_ALOC_GPU(obs.zmasks, isize);
        obs.exp_output = 0;
        SAFE_ALOC_GPU(obs.coeffs, vsize);
        cudaSafeCall(cudaMemcpy(obs.zmasks, zmasks[index].data(), isize,
                                cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(obs.coeffs, coeffs[index].data(), coeffs[index].size() * sizeof(ValType),
                                cudaMemcpyHostToDevice));
        ObservableList* obs_device = obsvals_dev + index;
        obsvec[index] = obs_device;
        cudaSafeCall(cudaMemcpy(obs_device, &obs, sizeof(ObservableList),
                                cudaMemcpyHostToDevice));
        measurement->EXPECT(obs_device); 
      };
      ~SV_CUDA_VQE()
        {
            cudaSafeCall(cudaMemcpy(obsvals_host, obsvals_dev, zmasks.size() * sizeof(ObservableList),
                                    cudaMemcpyDeviceToHost));
            for (size_t i = 0; i < zmasks.size(); i++) {
              ObservableList o = obsvals_host[i];
              SAFE_FREE_GPU(o.zmasks);
              SAFE_FREE_GPU(o.coeffs);
            }
            SAFE_FREE_GPU(obsvals_dev);
            SAFE_FREE_GPU(ansatz->gates->back().data);
        }
      virtual void call_simulator() override {        
        reset_state();
        sim(ansatz);
        sim(measurement);
        cudaDeviceSynchronize();
        cudaSafeCall(cudaMemcpy(obsvals_host, obsvals_dev, obsvals_host.size() * sizeof(ObservableList),
                                    cudaMemcpyDeviceToHost));
        expectation_value = 0;
        for (size_t i = 0; i < zmasks.size(); i++) {
          expectation_value += obsvals_host[i].exp_output;
        }
      };
      virtual void set_exp_gate(std::shared_ptr<Ansatz> circuit, ObservableList* o, std::vector<IdxType>& zmasks, std::vector<ValType>& coeffs) override {
        ObservableList obs;
        
        obs.numterms = xmasks[index].size();
        IdxType isize = obs.numterms * sizeof(IdxType);
        IdxType vsize = obs.numterms * sizeof(ValType);
        SAFE_ALOC_GPU(obs.zmasks, isize);
        obs.exp_output = 0;
        SAFE_ALOC_GPU(obs.coeffs, vsize);
        cudaSafeCall(cudaMemcpy(obs.zmasks, zmasks[index].data(), isize,
                                cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(obs.coeffs, coeffs[index].data(), coeffs[index].size() * sizeof(ValType),
                                cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(o, &obs, sizeof(ObservableList),
                                cudaMemcpyHostToDevice));
        measurement->EXPECT(obs_device);
      };
      virtual void delete_observables(ObservableList* observables) override {
        SAFE_FREE_GPU(observables);
      };

      virtual void get_exp_values(const std::vector<ObservableList*>& observables, std::vector<IdxType> sizes, std::vector<ValType>& output) override {
        for (size_t i = 0; i < observables.size(); i++) {
          std::vector<ObservableList> host_temp(sizes[i]);
          cudaSafeCall(cudaMemcpy(host_temp.data(), observables[i], sizes[i] * sizeof(ObservableList),
                                    cudaMemcpyDeviceToHost));
          
          for (ObservableList obs: host_temp) {
            output[i] += obs.exp_output;
          }
          // output.at(i) = observables.at(i)->exp_output;
        }
      };

      virtual void allocate_observables(ObservableList*& observables, IdxType size) override {
        SAFE_ALOC_GPU(observables, size * sizeof(ObservableList));
      };
    protected:
        ObservableList* obsvals_dev;
        ObservableList* obsvals_host;

};
  };
} // namespace NWQSim

#endif