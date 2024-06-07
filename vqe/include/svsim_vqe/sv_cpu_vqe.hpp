#ifndef VQE_STATE
#define VQE_STATE
#include "svsim/sv_cpu.hpp"
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
    inline
    ValType normsq(ValType real, ValType imag) {
      return (real * real + imag * imag);
    }

    
    inline
    IdxType count_ones(IdxType val) {
      IdxType count = 0;
      IdxType mask = 1;
      for (size_t i = 0; i < sizeof(IdxType) * 8 - 1; i++) {
        count += (val & mask) > 0;
        mask = mask << 1;
      }
      return count;
    }

    typedef std::function<void(const std::vector<ValType>& x, ValType ene, IdxType iter)> Callback;

    double nl_opt_function(const std::vector<double>& x, std::vector<double>& gradient, void* val);
    class SV_CPU_VQE: public SV_CPU {
      public:
      
        SV_CPU_VQE(std::shared_ptr<Ansatz> a, 
                   const Hamiltonian& h, 
                   nlopt::algorithm optimizer_algorithm,
                   Callback _callback,
                  IdxType seed = 0,
                  OptimizerSettings opt_settings = OptimizerSettings()): 
                                      SV_CPU(a->num_qubits()),
                                      hamil(h),
                                      ansatz(a),
                                      callback(_callback),
                                      g_est(seed),
                                      optimizer_settings(opt_settings)
                                      {
          optimizer = nlopt::opt(optimizer_algorithm, ansatz->numParams());
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

          // Check if the chosen algorithm requires derivatives
          compute_gradient = std::string(optimizer.get_algorithm_name()).find("no-derivative") == std::string::npos;
          optimizer.set_min_objective(nl_opt_function, (void*)this);
          std::vector<double> lower_bounds(ansatz->numParams(), 0);
          std::vector<double> upper_bounds(ansatz->numParams(), 2 * PI);
          optimizer.set_lower_bounds(lower_bounds);
          optimizer.set_upper_bounds(upper_bounds);
          
        };
      // function for the NLOpt plugin
      double cost_function(const std::vector<double> x, std::vector<double>& gradient) {
        if (compute_gradient) {
          gradient.resize(x.size());
          g_est.estimate([&] (const std::vector<double>& xval) { return energy(xval);}, x, gradient, 1e-4);
        }
        double ene = energy(x);
        if (callback != NULL) {
          callback(x, ene, iteration);
        }
        iteration++;
        return ene;
      }
        void optimize(std::vector<ValType>& parameters, ValType& final_ene) {
          iteration = 0;
          Config::PRINT_SIM_TRACE = false;
          if (parameters.size() == 0) {
            parameters = std::vector<ValType>(ansatz->numParams(), 0.0);
          }
          nlopt::result optimization_result = optimizer.optimize(parameters, final_ene);
        }

      ValType energy(const std::vector<double>& x) {
        ansatz->setParams(x);
        reset_state();
        sim(ansatz);

        ExpectationMap emap;
        auto& pauli_operators = hamil.getPauliOperators();
        for (auto& pauli_list: pauli_operators) {
          for (const PauliOperator& pauli: pauli_list) {
            ValType expectation = getPauliExpectation(pauli);
            emap[pauli] = expectation;
          }
        }
        ValType energy = hamil.expectation(emap);
        return energy;
      }
      protected:
        std::shared_ptr<Ansatz> ansatz;
        const Hamiltonian& hamil;
        SPSA g_est;
        nlopt::opt optimizer;
        bool compute_gradient;
        Callback callback;
        OptimizerSettings optimizer_settings;
        IdxType iteration;


      


        ValType getPauliExpectation(const PauliOperator& op) {
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
    double nl_opt_function(const std::vector<double>& x, std::vector<double>& gradient, void* val) {
        return ((SV_CPU_VQE*)val)->cost_function(x, gradient);
    };
  };
} // namespace NWQSim

#endif