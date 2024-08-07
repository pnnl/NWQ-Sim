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
                   const std::string& configpath,
                   IdxType seed = 0,
                   OptimizerSettings opt_settings = OptimizerSettings()): 
                                      SV_CPU(a->num_qubits(), configpath),
                                      VQEState(a, h, optimizer_algorithm, _callback, seed, opt_settings) {
        /*1.75702 6.23383 0.24862 3.75417 3.40728 0.359146 3.96801 2.66137 5.29208 5.69446 2.65872 4.11879 5.84607 2.16384 1.59979 2.79395 4.96491 4.36946 3.50401 1.80598 
        * 51(0.124291, -2.32453e-16), 
53(0.0854489, 1.04083e-16), 
54(-0.0145972, 1.02349e-16), 
57(0.164585, -6.93889e-17), 
58(0.0932238, 3.1225e-17), 
60(0.199378, -2.77556e-17), 
83(0.178024, 2.77556e-17), 
85(-0.02072, -8.06646e-17), 
86(0.208928, 3.67761e-16), 
89(0.00308507, -6.91721e-17), 
90(-0.072424, -1.52656e-16), 
92(0.267061, 1.38778e-16), 
99(-0.200104, 0), 
101(0.0221777, 7.02563e-17), 
102(0.253261, 4.996e-16), 
105(-0.0218919, 1.85615e-16), 
106(-0.152765, -3.19189e-16), 
108(-0.119388, -6.245e-17), 
147(-0.00111163, 2.08275e-16), 
149(0.125115, -1.31839e-16), 
150(0.126866, -1.11022e-16), 
153(0.139077, -3.46945e-16), 
154(0.0320659, 1.12757e-16), 
156(-0.0506285, 1.73472e-17), 
163(0.0617692, 1.97758e-16), 
165(0.0704726, 9.71445e-17), 
166(0.00355785, 7.56773e-17), 
169(-0.241969, 4.16334e-17), 
170(-0.0461118, -3.79904e-16), 
172(-0.335234, -2.77556e-16), 
195(0.199378, 2.63678e-16), 
197(-0.080827, -4.51028e-16), 
198(0.356338, 9.71445e-16), 
201(0.137854, -1.59595e-16), 
202(0.125496, 2.35922e-16), 
204(0.402323, 9.71445e-16),  */
      };
      virtual void call_simulator() override {      
        reset_state();
        expectation_value = 0;
        sim(ansatz);
        sim(measurement);
        for (auto i: obsvec) {
          std::cout << i->exp_output << std::endl;
          expectation_value += i->exp_output;
        }
        std::cout << expectation_value << std::endl;
      };

      virtual void call_simulator(std::shared_ptr<Ansatz> _measurement) override {    
        reset_state();
        sim(ansatz);
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