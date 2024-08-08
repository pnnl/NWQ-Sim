#ifndef VQE_STATE
#define VQE_STATE
#include "state.hpp"
#include "observable/pauli_operator.hpp"
#include "utils.hpp"
#include "circuit/ansatz.hpp"
#include "circuit/measurement.hpp"
#include "gradient/sa_gradient.hpp"
#include "observable/hamiltonian.hpp"
#include "circuit_pass/fusion.hpp"
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

    


    typedef std::function<void(const std::vector<ValType>& x, ValType ene, IdxType iter)> Callback;

    double nl_opt_function(const std::vector<double>& x, std::vector<double>& gradient, void* val);
    class VQEState {
      public:
      
        VQEState(std::shared_ptr<Ansatz> a, 
                   std::shared_ptr<Hamiltonian> h, 
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
        }
        void initialize() {
            // 0 3.273514e-01
// Evaluation 0, fval = 21.949690

// Finished VQE loop.
        // Final value: 2.194969e+01
         //   0.519016 (0.00324322 + 0i)IIIIZIIIIIXX
         //-0.263524 (-0.00210224 + 0i)YZZZZYZIIIII
            // -0.264174 (-5.01215e-05 + 0i)YZZZZYIIIZII
            // 0.00287325 (0.00444125 + 0i)IIIXXIIIIXZX
            // 0.158064 (3.65849e-18 + 0i)IYYIIIIIIIIIIIIIIIII
            // 0.0981487 (3.65855e-18 + 0i)IIIIIIIIIIIYYIIIIIII
            // 0.829777 (0.094532 + 0i)IIIIIIIIIIIZIZIIIIII
            // 0.0108362 (0.0108748 + 0i)XZXIIIIIIIXZXIIIIIII
          const std::vector<std::vector<PauliOperator> >& pauli_operators = hamil->getPauliOperators();  
          IdxType index = 0;    
          obsvec.resize(pauli_operators.size());
          xmasks.resize(pauli_operators.size());
          zmasks.resize(pauli_operators.size());
          x_index_sizes.resize(pauli_operators.size());
          x_indices.resize(pauli_operators.size());
          coeffs.resize(pauli_operators.size());
          std::vector<IdxType> mapping (ansatz->num_qubits());
          std::iota(mapping.begin(), mapping.end(), 0);
          for (auto& pauli_list: pauli_operators) {
            xmasks[index].reserve(pauli_list.size());
            zmasks[index].reserve(pauli_list.size());
            coeffs[index].reserve(pauli_list.size());
            x_index_sizes[index].reserve(pauli_list.size());
            IdxType composite_xmask = 0;
            IdxType composite_zmask = 0;
            IdxType ncommute = pauli_list.size();
            for (const PauliOperator& pauli: pauli_list) {
              std::vector<IdxType> xinds;
              pauli.get_xindices(xinds);
              coeffs[index].push_back(pauli.getCoeff().real());
              composite_xmask |= pauli.get_xmask();
              xmasks[index].push_back(0);
              coeffs[index].back() *= (pauli.count_y() % 2) ? -1.0 : 1.0;
              x_index_sizes[index].push_back(0);
              zmasks[index].push_back(pauli.get_zmask() | pauli.get_xmask());
              composite_zmask |= pauli.get_zmask();
            }

            PauliOperator common(composite_xmask, composite_zmask, ansatz->num_qubits());
            
            Measurement circ1(common);
            ansatz->compose(circ1, mapping);
            
            fill_obslist(index);
            
            Measurement circ2(common, true);
            ansatz->compose(circ2, mapping);
            index++;
          }
                                          
          // Check if the chosen algorithm requires derivatives
          compute_gradient = std::string(optimizer.get_algorithm_name()).find("no-derivative") == std::string::npos;
          optimizer.set_min_objective(nl_opt_function, (void*)this);
          std::vector<double> lower_bounds(ansatz->numParams(), -2 * PI);
          std::vector<double> upper_bounds(ansatz->numParams(), 2 * PI);
          optimizer.set_lower_bounds(lower_bounds);
          optimizer.set_upper_bounds(upper_bounds);
          
        };
      virtual void fill_obslist(IdxType index) {};
      // function for the NLOpt plugin
      double cost_function(const std::vector<double> x, std::vector<double>& gradient) {
        if (iteration > 0){
          Config::PRINT_SIM_TRACE = false;
        }
        if (compute_gradient) {
          gradient.resize(x.size());
          g_est.estimate([&] (const std::vector<double>& xval) { return energy(xval);}, x, gradient, 1e-4);
        }
        if (iteration > 0){
          Config::PRINT_SIM_TRACE = false;
        }
        double ene = energy(x);
        if (callback != NULL) {
          callback(x, ene, iteration);
        }
        iteration++;
        return ene;
      }
      virtual std::vector<std::pair<std::string, ValType>> follow_fixed_gradient(const std::vector<ValType>& x0, ValType& initial_ene, ValType& final_ene, IdxType& num_iterations, ValType delta, ValType eta, IdxType n_grad_est) {
        Config::PRINT_SIM_TRACE = false;
        std::vector<ValType> gradient (x0.size(),1.0);
        std::vector<ValType> params(x0);
        std::vector<ValType> minima_params(x0);
        ValType ene_prev = MAXFLOAT;
        ValType ene_curr = energy(params);
        initial_ene = ene_curr;
        // gradient
        // get the single-direction starting vector
        g_est.estimate([&] (const std::vector<double>& xval) { return energy(xval);}, params, gradient, delta, n_grad_est);
        iteration = 0;
        do {
          for (size_t i = 0; i < params.size(); i++) {
            params[i] -= eta * gradient[i];
          }
          // auto s1 =  std::chrono::high_resolution_clock::now();
          // ene_curr = 0;
          ene_curr = energy(params);
          // for (auto d: gradient) {
          //   std::cout << d << " ";
          // }
          // std::cout << std::endl;
          
          // std::cout << step << " " << ene_curr << " " << ene_prev << std::endl;
          if (ene_curr >= ene_prev) {
            for (size_t i = 0; i < params.size(); i++) {
              params[i] += eta * gradient[i];
            }
            break;
          } else {
            ene_prev = ene_curr;
          }
          iteration++;
          // auto s2 =  std::chrono::high_resolution_clock::now();
          // std::cout << (s2-s1).count()/1e9 << std::endl;
        } while(true);
        num_iterations = iteration;
        // std::cout << "Ended loop\n" << std::endl;
        std::vector<std::string> fermi_strings = ansatz->getFermionicOperatorStrings();
        final_ene = ene_curr;
        std::vector<std::pair<std::string, ValType>> result = ansatz->getFermionicOperatorParameters();
        return result;
      }
      virtual void optimize(std::vector<ValType>& parameters, ValType& final_ene) {
          iteration = 0;
          if (parameters.size() == 0) {
            parameters = std::vector<ValType>(ansatz->numParams(), 0.0);
          }
          // final_ene = energy(parameters);
          nlopt::result optimization_result = optimizer.optimize(parameters, final_ene);
      }
      virtual void call_simulator() {};
      virtual ValType energy(const std::vector<double>& x) {
        ansatz->setParams(x);

        call_simulator();

        // ExpectationMap emap;

      
        // const std::vector<std::vector<PauliOperator> >& pauli_operators = hamil->getPauliOperators();    
        // auto& pauli_operators = hamil->getPauliOperators();
        // double expval = 0.0;
        // for (auto& clique: pauli_operators) {
        //   for (auto& pauli: clique) {
        //     expval += getPauliExpectation(pauli) * pauli.getCoeff().real();
        //   }
        // }
        IdxType index = 0;
        ValType expectation = hamil->getEnv().constant + expvals.front();
        // ValType ene = 0.0;
        // ValType ene = hamil.expectation(emap);
        return expectation;
      }
      IdxType get_iteration() const {return iteration;};
      protected:
        std::shared_ptr<Ansatz> ansatz;
        std::shared_ptr<Hamiltonian> hamil;
        SPSA g_est;
        nlopt::opt optimizer;
        bool compute_gradient;
        Callback callback;
        OptimizerSettings optimizer_settings;
        IdxType iteration;
        std::vector<ObservableList> obsvec;
        std::vector<std::vector<IdxType> >  x_index_sizes;
        std::vector<std::vector<IdxType> > xmasks;
        std::vector<std::vector<IdxType> > zmasks;
        std::vector<std::vector<IdxType> > x_indices;
        std::vector<std::vector<ValType> > coeffs;
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