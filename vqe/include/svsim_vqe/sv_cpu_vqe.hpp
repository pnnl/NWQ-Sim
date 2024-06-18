#ifndef VQE_CPU_STATE
#define VQE_CPU_STATE
#include "svsim/sv_cpu.hpp"
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
    class SV_CPU_VQE: public VQEState, public SV_CPU {
      public:
        SV_CPU_VQE(std::shared_ptr<Ansatz> a, 
                   std::shared_ptr<Hamiltonian> h, 
                   nlopt::algorithm optimizer_algorithm,
                   Callback _callback,
                   IdxType seed = 0,
                   OptimizerSettings opt_settings = OptimizerSettings()): 
                                      SV_CPU(a->num_qubits()),
                                      VQEState(a, h, optimizer_algorithm, _callback, seed, opt_settings) {
        obs.xmasks = xmasks.data();
        obs.zmasks = zmasks.data();
        obs.numterms = xmasks.size();
        obs.exp_output = expvals.data();
        obs.x_indices = x_indices.data();
        obs.x_index_sizes = x_index_sizes.data();
        ansatz->EXPECT(&obs);
      };
      virtual void call_simulator() override {        
        reset_state();
        sim(ansatz);
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
          double val = 0.0;
          for (IdxType i = 0; i < half_dim; i++) {
            IdxType idx0 = ((i << 1) & mask_u) | (i & mask_l);
            IdxType idx1 = xmask ^ idx0;
            ValType v0, v1;
            val += sv_imag[idx1] * sv_imag[idx1]  + sv_real[idx1] * sv_real[idx1];
            val += sv_imag[idx0] * sv_imag[idx0]  + sv_real[idx0] * sv_real[idx0];
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