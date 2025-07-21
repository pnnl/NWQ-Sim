#ifndef VQE_CPU_STATE
#define VQE_CPU_STATE
#include "state.hpp"
#include "svsim/sv_cpu.hpp"
#include "vqe_state.hpp"
#include "observable/pauli_operator.hpp"
#include "utils.hpp"
#include "circuit/ansatz.hpp"
#include "circuit/measurement.hpp"
#include "gradient/sa_gradient.hpp"
#include "observable/hamiltonian.hpp"
#include "nlopt/nlopt.hpp"
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
      };
      virtual void call_simulator() override {      
        reset_state();
        expectation_value = 0;
        sim(ansatz);
        sim(measurement);
        for (auto i: obsvec) {
          expectation_value += i->exp_output;
        }
      };

      virtual void call_simulator(std::shared_ptr<Ansatz> _measurement, bool reset) override {    
        if (reset) {
          reset_state();
          sim(ansatz);
        }
        sim(_measurement);
      };
      virtual void fill_obslist(IdxType index) override {
        ObservableList*& obs = obsvec[index];
        obs = new ObservableList;
        obs->coeffs = coeffs[index].data();
        obs->zmasks = zmasks[index].data();
        obs->numterms = zmasks[index].size();
        measurement->EXPECT(obs);
      };
      virtual ValType getPauliExpectation(const PauliOperator& op) override {
          IdxType qubit = op.get_dim();
          IdxType xmask = op.get_xmask();
          IdxType zmask = op.get_zmask();
          IdxType y_phase = 0;
          IdxType max_x = 0;
          for (IdxType i = 0; i < qubit; i++) {
            bool xbit = (xmask >> i) & 1;
            bool zbit = (zmask >> i) & 1;
            y_phase += (xbit && zbit);
            if (xbit)
              max_x = std::max(i, max_x);
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

      ~SV_CPU_VQE() {

        for (auto i: obsvec) {
            delete i;
        }
      }
    };
  };
} // namespace NWQSim

#endif