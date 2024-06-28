#ifndef VQE_CUDA_STATE
#define VQE_CUDA_STATE
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
                   IdxType seed = 0,
                   OptimizerSettings opt_settings = OptimizerSettings()): 
                                      SV_CUDA(a->num_qubits()),
                                      VQEState(a, h, optimizer_algorithm, _callback, seed, opt_settings) {
        IdxType nterms = xmasks.size();
        obs.numterms = nterms;
        IdxType isize = nterms * sizeof(IdxType);
        IdxType vsize = nterms * sizeof(ValType);
        expvals.resize(1);
        SAFE_ALOC_GPU(obs.xmasks, isize);
        SAFE_ALOC_GPU(obs.zmasks, isize);
        SAFE_ALOC_GPU(obs.x_index_sizes, isize);
        SAFE_ALOC_GPU(obs.x_indices, x_indices.size() * sizeof(IdxType));
        SAFE_ALOC_GPU(obs.exp_output, expvals.size()* sizeof(ValType));
        SAFE_ALOC_GPU(obs.coeffs, vsize);
        cudaSafeCall(cudaMemcpy(obs.xmasks, xmasks.data(), isize,
                                    cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(obs.zmasks, zmasks.data(), isize,
                                    cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(obs.x_index_sizes, x_index_sizes.data(), isize,
                                    cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(obs.x_indices, x_indices.data(), x_indices.size() * sizeof(IdxType),
                                    cudaMemcpyHostToDevice));
        cudaSafeCall(cudaMemcpy(obs.coeffs, coeffs.data(), coeffs.size() * sizeof(ValType),
                                    cudaMemcpyHostToDevice));
        ObservableList* obs_device;
        SAFE_ALOC_GPU(obs_device, sizeof(ObservableList));
        cudaSafeCall(cudaMemcpy(obs_device, &obs, sizeof(ObservableList),
                                    cudaMemcpyHostToDevice));
        ansatz->EXPECT(obs_device);
      };

      ~SV_CUDA_VQE()
        {
            SAFE_FREE_GPU(obs.xmasks);
            SAFE_FREE_GPU(obs.zmasks);
            SAFE_FREE_GPU(obs.x_index_sizes);
            SAFE_FREE_GPU(obs.x_indices);
            SAFE_FREE_GPU(obs.coeffs);
            SAFE_FREE_GPU(obs.exp_output);
            SAFE_FREE_GPU(ansatz->gates->back().data);
        }
      virtual void call_simulator() override {        
        reset_state();
        cudaSafeCall(cudaMemset(obs.exp_output, 0, expvals.size() * sizeof(ValType)));
        sim(ansatz);
        cudaDeviceSynchronize();
        cudaSafeCall(cudaMemcpy(expvals.data(), obs.exp_output, expvals.size() * sizeof(ValType),
                                    cudaMemcpyDeviceToHost));
      };


  };
} // namespace NWQSim

#endif