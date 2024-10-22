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
#include <vector>

namespace NWQSim
{
  namespace VQE {
    inline
    ValType normsq(ValType real, ValType imag) {
      return (real * real + imag * imag);
    }

    

    // callback function for printing updates to terminal
    typedef std::function<void(const std::vector<ValType>& x, ValType ene, IdxType iter)> Callback;

    // function forward declaration to pass to NLOpt
    double nl_opt_function(const std::vector<double>& x, std::vector<double>& gradient, void* val);
    class VQEState {
      /**
       * @brief  VQE Backend Base class
       * @note   Needs to be inherited by backend-specific implementations
       * @param  a: Ansatz to optimize
       * @param  h: Hamiltonian whose expectation value we want to minimize
       * @param  _optimizer_algorithm: NLOpt optimization algorithm to use
       * @param  _callback: Callback function for updates
       * @param  seed: Random seed for SPSA, stochastic optimization
       * @param  opt_settings: Settings for NLOpt
       * @retval None
       */
      public:
      
        // Ctor
        VQEState( std::shared_ptr<Ansatz> a, 
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
          process_rank = 0;
        }

        // Dtor
        ~VQEState(){}

        /**
         * @brief  Initialize the VQE backend
         * @note   Responsible for constructing the measurement circuit and allocating storage space for the VQE 
         * @retval None
         */
        void initialize() {
          // Get Hamiltonian operators, assumes already grouped into QWC cliques
          const std::vector<std::vector<PauliOperator> >& pauli_operators = hamil->getPauliOperators(); 

          IdxType index = 0;    

          // Allocate memory for measurement circuit, ObservableList instances, and observable vectors
          measurement.reset(new Ansatz(ansatz->num_qubits()));
          allocate_observables(pauli_operators.size()); // Backend specific memory allocation routine

          zmasks.resize(pauli_operators.size());
          coeffs.resize(pauli_operators.size());

          std::vector<IdxType> mapping (ansatz->num_qubits()); // make a dummy layout
          std::iota(mapping.begin(), mapping.end(), 0);
          // For each commuting group
          for (auto& pauli_list: pauli_operators) {
            // Make the common Pauli operator by taking a bitwise "or" over stabilizer states
            PauliOperator common = make_common_op(pauli_list, 
                                                  zmasks[index], 
                                                  coeffs[index]);
            // make a measurement circuit for the common operator and append
            Measurement circ1(common);
            measurement->compose(circ1, mapping);
            fill_obslist(index); // backend-specific observable construction/initialization routines
            Measurement circ2(common, true);
            measurement->compose(circ2, mapping); // append the measurement inverse to undo diagonalization
            index++;
          }
                            
          // Check if the chosen algorithm requires derivatives
          compute_gradient = std::string(nlopt::algorithm_name(optimizer_algorithm)).find("no-derivative") == std::string::npos;     
        };
      virtual void fill_obslist(IdxType index) {};
      // function for the NLOpt plugin
      double cost_function(const std::vector<double>& x, std::vector<double>& gradient) {
        // only print NWQ-Sim stats on the first iteration 
        if (iteration > 0){
          Config::PRINT_SIM_TRACE = false;
        }
        // if gradient-based optimizer, use SPSA to approximate gradient
        if (compute_gradient) {
          gradient.resize(x.size());
          g_est.estimate([&] (const std::vector<double>& xval) { return energy(xval);}, x, gradient, 1e-4, 1);
        }
        // compute the Hamiltonian expectation value (Calls simulator)
        double ene = energy(x);
        // callback message to terminal
        if (callback != NULL) {
          callback(x, ene, iteration);
        }
        iteration++;
        return ene;
      }

      /**
       * @brief  QFlow local minimum finder
       * @note   Follows a single gradient vector until a local minimum is reached, used in QFlow iterative schema 
       * @param  x0: Initial parameter vector (input)
       * @param  initial_ene: Initial energy (output)
       * @param  final_ene: Final energy of local minimum (output)
       * @param  num_iterations: Number of iterations before convergence (output)
       * @param  delta: Perturbation magnitude for SPSA
       * @param  eta: Gradient descent step size
       * @param  n_grad_est: Number of gradient estimates for SPSA average (higher=more accurate but slower)
       * @retval 
       */
      virtual std::vector<std::pair<std::string, ValType>> follow_fixed_gradient(const std::vector<ValType>& x0, ValType& initial_ene, ValType& final_ene, IdxType& num_iterations, ValType delta, ValType eta, IdxType n_grad_est) {
        Config::PRINT_SIM_TRACE = false;
        // initialize data structures
        std::vector<ValType> gradient (x0.size(),1.0);
        std::vector<ValType> params(x0);
        std::vector<ValType> minima_params(x0);
        ValType ene_prev = MAXFLOAT;
        ValType ene_curr = energy(params); // starting energy
        initial_ene = ene_curr;
        // get the single-direction starting vector
        g_est.estimate([&] (const std::vector<double>& xval) { return energy(xval);}, params, gradient, delta, n_grad_est);
        iteration = 0;
        // until we reach a local minimum, descend along the starting gradient 
        do {
          // descent step
          for (size_t i = 0; i < params.size(); i++) {
            params[i] -= eta * gradient[i];
          }
          ene_curr = energy(params);
          if (ene_curr >= ene_prev) {
            // local minimum reached, undo the last update and break
            for (size_t i = 0; i < params.size(); i++) {
              params[i] += eta * gradient[i];
            }
            break;
          } else {
            ene_prev = ene_curr;
          }
          iteration++;
        } while(true);
        num_iterations = iteration;
        std::vector<std::string> fermi_strings = ansatz->getFermionicOperatorStrings();
        final_ene = ene_curr;
        std::vector<std::pair<std::string, ValType>> result = ansatz->getFermionicOperatorParameters();
        return result;
      }

      /**
       * @brief  Main Optimization Loop
       * @note  May be overloaded by specific backends (e.g. MPI, NVGPU_MPI)
       * @param  parameters: Parameter vector reference (output)
       * @param  final_ene: Energy output
       * @retval None
       */
      virtual void optimize(std::vector<ValType>& parameters, ValType& final_ene) {
          // Initialize optimizer 
          nlopt::opt optimizer = nlopt::opt(optimizer_algorithm, ansatz->numParams());
          // set the objective function
          optimizer.set_min_objective(nl_opt_function, (void*)this);
          std::vector<double> lower_bounds(ansatz->numParams(), optimizer_settings.lbound);
          std::vector<double> upper_bounds(ansatz->numParams(), optimizer_settings.ubound);  
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
          // Call the optimizer
          // final_ene = energy(parameters);
          // nlopt::result optimization_result = optimizer.optimize(parameters, final_ene); // MZ
          optimizer.optimize(parameters, final_ene); // MZ
          opt_result = optimizer.last_optimize_result(); // MZ: this is the correct way to get the result, otherwise always give 0
          num_evals = optimizer.get_numevals(); // MZ: get numebr of function evaluations
          
      }

      // Function declarations (overloaded by backends) 
      virtual void call_simulator() {};
      virtual void call_simulator(std::shared_ptr<Ansatz> _measurement, bool reset) {};
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
      virtual void allocate_observables(IdxType size) {
        obsvec.resize(size);
      };
      virtual void delete_observables(ObservableList* observables, IdxType size) {
        delete[] observables;
      };
      virtual ValType energy(const std::vector<double>& x) {
        ansatz->setParams(x);

        call_simulator();
        ValType expectation = hamil->getEnv().constant + expectation_value;
        return expectation;
      }
      IdxType get_process_rank() const { return process_rank; }
      
      std::shared_ptr<Hamiltonian> get_hamiltonian() const { return hamil; }
      IdxType get_iteration() const {return iteration;};
      nlopt::result get_optresult() const {return opt_result;}; // MZ
      unsigned int get_numevals() const { return num_evals;}; // MZ
      double get_duration() const { return opt_duration;}; // MZ
      void set_duration(double value) { opt_duration = value; }; //MZ
      protected:
        std::shared_ptr<Ansatz> ansatz;                    // state preparation circuit
        std::shared_ptr<Ansatz> measurement;               // circuit to measure expectation values
        std::shared_ptr<Hamiltonian> hamil;                // target system Hamiltonian
        SPSA g_est;                                        // stochastic gradient estimator
        IdxType process_rank;                              // process rank (used by MPI/NVGPU_MPI backends)
        bool compute_gradient;                             // flag to compute gradient (depends on optimizer selected)
        Callback callback;                                 // callback function for terminal updates
        OptimizerSettings optimizer_settings;              // NLOpt optimizer settings (bounds, termination critera)
        IdxType iteration;                                 // current iteration counter
        std::vector<ObservableList*> obsvec;               // vector of pointers to ObservableList structures for clique expectation value calculation
        std::vector<std::vector<IdxType> > zmasks;         // vector of diagonalized operator zmasks
        std::vector<std::vector<ValType> > coeffs;         // vector of diagonalized operator coefficients
        double expectation_value;                          // last computed expectation value
        nlopt::algorithm optimizer_algorithm;              // NLOpt optimization algorithm for circuit updates 
        nlopt::result opt_result;                          //MZ: optimzation success or fail and the reason
        unsigned int num_evals;                            //MZ: total number of evaluations
        double opt_duration;                               //MZ: time the optimization

      

      // deprecated function, only implemented by CPU backend for debugging
      virtual ValType getPauliExpectation(const PauliOperator& op) {
        throw std::runtime_error("Pauli expectation not implemented for this backend");
      };

    };
    // Target function for NLOpt object, has to match function prototype
    double nl_opt_function(const std::vector<double>& x, std::vector<double>& gradient, void* val) {
        return ((VQEState*)val)->cost_function(x, gradient);
    };
  };
} // namespace NWQSim

#endif