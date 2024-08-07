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
#include <vector>

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
                   nlopt::algorithm _optimizer_algorithm,
                   Callback _callback,
                  IdxType seed = 0,
                  OptimizerSettings opt_settings = OptimizerSettings()): 
                                      hamil(h),
                                      ansatz(a),
                                      callback(_callback),
                                      g_est(seed),
                                      optimizer_algorithm(_optimizer_algorithm),
                                      optimizer_settings(opt_settings)
                                      {
          
        }
        ~VQEState(){
        for (auto i: obsvec) {
            delete i;
        }
      }
        void initialize() {
          const std::vector<std::vector<PauliOperator> >& pauli_operators = hamil->getPauliOperators(); 
          IdxType index = 0;    
          measurement.reset(new Ansatz(ansatz->num_qubits()));
          obsvec.clear();
          zmasks.clear();
          coeffs.clear();
          obsvec.resize(pauli_operators.size());
          zmasks.resize(pauli_operators.size());
          coeffs.resize(pauli_operators.size());
          std::vector<IdxType> mapping (ansatz->num_qubits());
          std::iota(mapping.begin(), mapping.end(), 0);
          for (auto& pauli_list: pauli_operators) {
            PauliOperator common = make_common_op(pauli_list, 
                                                  zmasks[index], 
                                                  coeffs[index]);
            Measurement circ1(common);
            measurement->compose(circ1, mapping);
            fill_obslist(index);
            Measurement circ2(common, true);
            measurement->compose(circ2, mapping);
            index++;
          }
                                          
          // Check if the chosen algorithm requires derivatives
          compute_gradient = std::string(nlopt::algorithm_name(optimizer_algorithm)).find("no-derivative") == std::string::npos;
          
          
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
      virtual void get_current_gradient(std::vector<double>& gradient, ValType delta, ValType eta, IdxType n_grad_est) {
        g_est.estimate([&] (const std::vector<double>& xval) { return energy(xval);}, ansatz->getParamRef(), gradient, delta, n_grad_est);
      }
      void set_ansatz(std::shared_ptr<Ansatz> new_a) {
        ansatz = new_a;
      }
      void swap_hamil(std::shared_ptr<Hamiltonian>& _hamil) {
        hamil.swap(_hamil);
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
          nlopt::opt optimizer = nlopt::opt(optimizer_algorithm, ansatz->numParams());
          optimizer.set_min_objective(nl_opt_function, (void*)this);
          std::vector<double> lower_bounds(ansatz->numParams(), -2 * PI);
          std::vector<double> upper_bounds(ansatz->numParams(), 2 * PI);
          optimizer.set_lower_bounds(lower_bounds);
          optimizer.set_upper_bounds(upper_bounds);
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
          iteration = 0;
          nlopt::result optimization_result = optimizer.optimize(parameters, final_ene);
      }
      virtual void call_simulator() {};
      virtual void call_simulator(std::shared_ptr<Ansatz> _measurement) {};
      virtual void set_exp_gate(std::shared_ptr<Ansatz> circuit, ObservableList* o, std::vector<IdxType>& zmasks, std::vector<ValType>& coeffs) {
        o->zmasks = zmasks.data();
        o->coeffs = coeffs.data();
        o->numterms = coeffs.size();
        circuit->EXPECT(o);
      };
      virtual void get_exp_values(const std::vector<ObservableList*>& observables, std::vector<IdxType> sizes, std::vector<ValType>& output) {
        for (size_t i = 0; i < observables.size(); i++) {
          for (size_t j = 0; j < sizes[i]; j++) {
            output[i] += observables[i][j].exp_output;
          }
          // output.at(i) = observables.at(i)->exp_output;
        }
      };
      virtual void allocate_observables(ObservableList*& observables, IdxType size) {
        observables = new ObservableList[size];
      };
      virtual void delete_observables(ObservableList* observables) {
        delete[] observables;
      };
      virtual ValType energy(const std::vector<double>& x) {
        ansatz->setParams(x);

        call_simulator();
        ValType expectation = hamil->getEnv().constant + expectation_value;
        return expectation;
      }
      
      std::shared_ptr<Hamiltonian> get_hamiltonian() const { return hamil; }
      IdxType get_iteration() const {return iteration;};
      protected:
        std::shared_ptr<Ansatz> ansatz;
        std::shared_ptr<Ansatz> measurement;
        std::shared_ptr<Hamiltonian> hamil;
        SPSA g_est;
        bool compute_gradient;
        Callback callback;
        OptimizerSettings optimizer_settings;
        IdxType iteration;
        std::vector<ObservableList*> obsvec;
        std::vector<std::vector<IdxType> >  x_index_sizes;
        std::vector<std::vector<IdxType> > xmasks;
        std::vector<std::vector<IdxType> > zmasks;
        std::vector<std::vector<IdxType> > x_indices;
        std::vector<std::vector<ValType> > coeffs;
        double expectation_value;
        nlopt::algorithm optimizer_algorithm;

      


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