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
        SAFE_ALOC_GPU(obs.xmasks, isize);
        SAFE_ALOC_GPU(obs.zmasks, isize);
        SAFE_ALOC_GPU(obs.x_index_sizes, isize);
        SAFE_ALOC_GPU(obs.x_indices, x_indices.size() * sizeof(IdxType));
        SAFE_ALOC_GPU(obs.exp_output, vsize);
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
      };

      ~SV_CUDA_VQE()
        {
            SAFE_FREE_GPU(obs.xmasks);
            SAFE_FREE_GPU(obs.zmasks);
            SAFE_FREE_GPU(obs.x_index_sizes);
            SAFE_FREE_GPU(obs.x_indices);
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

      virtual ValType getPauliExpectation(const PauliOperator& op) override {
          IdxType qubit = 0;
          IdxType xmask = 0;
          IdxType zmask = 0;
          IdxType y_phase = 0;
          IdxType max_x = 0;
          for (auto pauli_op: *op.getOps()) {
            switch (pauli_op)
            {
              case PauliOp::X:
                xmask = xmask | (1 << qubit);
                max_x = qubit;
                break;
              case PauliOp::Y:
                xmask = xmask | (1 << qubit);
                zmask = zmask | (1 << qubit);
                max_x = qubit;
                y_phase += 1;
                break;
              case PauliOp::Z:
                zmask = zmask | (1ll << qubit);
                break;
              default:
                break;
            }
            qubit++;
          }
          ValType expectation = 0.0;
          if (xmask == 0) {
            for (IdxType i = 0; i < dim; i++) {
              ValType local_exp = sv_real[i] * sv_real[i] - sv_imag[i] * sv_imag[i];
              if (count_ones(zmask & i) & 1) {
                local_exp *= -1;
              }
              expectation += local_exp;
            }
            return expectation;
          }
          ValType sign = (y_phase / 2) % 2 ? -1: 1;
          ValType phase = y_phase % 2;
          size_t mask_u = ~((1lu << (max_x + 1)) - 1);
          size_t mask_l = (1lu << (max_x)) - 1;
          for (IdxType i = 0; i < half_dim; i++) {
            IdxType idx0 = ((i << 1) & mask_u) | (i & mask_l);
            IdxType idx1 = xmask ^ idx0;
            ValType v0, v1;
            if (phase) {
              v0 = -sv_imag[idx1] * sv_real[idx0] + sv_imag[idx0] * sv_real[idx1];
              v1 = -sv_imag[idx0] * sv_real[idx1] + sv_imag[idx1] * sv_real[idx0];
            } else {
              v0 = sv_real[idx1] * sv_real[idx0] + sv_imag[idx1] * sv_imag[idx0];
              v1 = sv_real[idx1] * sv_real[idx0] + sv_imag[idx0] * sv_imag[idx1];
            }
            v0 *= sign;
            v1 *= sign;
            ValType thisval = v0;
            if ((count_ones(idx0 & zmask) & 1) != 0) {
              thisval *= -1;
            }
            if ((count_ones(idx1 & zmask) & 1) != 0) {
              thisval -= v1;
            } else {
              thisval += v1;
            }
            expectation += thisval;
          }
          return expectation;
        }       

    };
  };
} // namespace NWQSim

#endif