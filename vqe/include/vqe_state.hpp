#ifndef VQE_STATE
#define VQE_STATE
#include "state.hpp"
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
    class VQEState {
      public:
      
        VQEState(std::shared_ptr<Ansatz> a, 
                   const Hamiltonian& h, 
                   nlopt::algorithm optimizer_algorithm,
                   Callback _callback,
                  IdxType seed = 0,
                  OptimizerSettings opt_settings = OptimizerSettings()): 
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
          auto& pauli_operators = hamil.getPauliOperators();
          for (auto& pauli_list: pauli_operators) {
            for (const PauliOperator& pauli: pauli_list) {
              const std::vector<IdxType>& xinds = pauli.get_xindices();
              xmasks.push_back(pauli.get_xmask());
              zmasks.push_back(pauli.get_zmask());
              x_index_sizes.push_back(xinds.size());
              for (auto i: xinds) {
                x_indices.push_back(i);
              }
            }
          }
          expvals.resize(hamil.num_ops());
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
      virtual void optimize(std::vector<ValType>& parameters, ValType& final_ene) {
          iteration = 0;
          if (parameters.size() == 0) {
            parameters = std::vector<ValType>(ansatz->numParams(), 0.0);
          }
          nlopt::result optimization_result = optimizer.optimize(parameters, final_ene);
      }
      virtual void call_simulator() {};
      virtual ValType energy(const std::vector<double>& x) {
        ansatz->setParams(x);

        if (iteration > 0){
          Config::PRINT_SIM_TRACE = false;
        } else {
          Config::PRINT_SIM_TRACE = true;
        }
        call_simulator();

        ExpectationMap emap;
        auto& pauli_operators = hamil.getPauliOperators();
        IdxType index = 0;
        for (auto& pauli_list: pauli_operators) {
          for (const PauliOperator& pauli: pauli_list) {
            emap[pauli] = expvals[index++];
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
        ObservableList obs;
        std::vector<IdxType> x_index_sizes;
        std::vector<IdxType> xmasks;
        std::vector<IdxType> zmasks;
        std::vector<IdxType> x_indices;
        std::vector<ValType> expvals;

      


      virtual ValType getPauliExpectation(const PauliOperator& op) {
        throw std::runtime_error("Pauli expectation not implemented for this backend");
      };

    };
    double nl_opt_function(const std::vector<double>& x, std::vector<double>& gradient, void* val) {
        return ((VQEState*)val)->cost_function(x, gradient);
    };
  };
} // namespace NWQSim

#endif