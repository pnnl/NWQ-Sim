#ifndef VQE_STATE
#define VQE_STATE
#include <random>
#include "state.hpp"
#include "observable/pauli_operator.hpp"
#include "utils.hpp"
#include "circuit/ansatz.hpp"
#include "circuit/measurement.hpp"
#include "gradient/sa_gradient.hpp"
#include "observable/hamiltonian.hpp"
#include "nlopt/nlopt.hpp"
#include <memory>
#include <cmath>
#include <vector>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <algorithm>

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
                                      g_est(seed ? seed : static_cast<IdxType>(std::random_device{}())),
                                      // g_est(),
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
      virtual std::vector<std::pair<std::string, ValType>> follow_fixed_gradient(const std::vector<ValType>& x0, ValType& initial_ene, ValType& final_ene, IdxType& num_iterations, 
                                                                                ValType delta, ValType eta, IdxType n_grad_est, bool verbose) {
        Config::PRINT_SIM_TRACE = false;
        local_result = 9; // MZ: initial flag
        // initialize data structures
        std::vector<ValType> gradient (x0.size(),1.0);
        std::vector<ValType> params(x0);
        std::vector<ValType> minima_params(x0);
        ValType ene_prev = MAXFLOAT;
        ValType ene_curr = energy(params); // starting energy
        initial_ene = ene_curr;
        // get the single-direction starting vector
        
        // Reseed SPSA so each call uses a fresh random direction unless caller passed a nonzero seed
        // g_est.reseed(static_cast<uint64_t>(std::random_device{}()));
        g_est.estimate([&] (const std::vector<double>& xval) { return energy(xval);}, params, gradient, delta, n_grad_est);
        iteration = 0;
        // until we reach a local minimum, descend along the starting gradient 
        do {
          // descent step
          for (size_t i = 0; i < params.size(); i++) {
            params[i] -= eta * gradient[i];
          }
          ene_curr = energy(params);
          // MZ: some printout
          if ( (process_rank == 0) && verbose) {
            if (iteration == 0) {
              std::cout << "\n----------- Iteration Summary -----------\n" << std::left
                        << std::setw(8) << " Iter."
                        << std::setw(20) << "Objective Value"
                        << std::setw(12) << "Grad. Norm"
                        << std::endl;
              std::cout << std::string(41, '-') << std::endl;
            }
            double grad_norm = std::sqrt(std::accumulate(gradient.begin(), gradient.end(), 0.0, [] (ValType a, ValType b) {
              return a + b * b;
            }));
            std::cout << std::left << " "
                      << std::setw(7) << iteration
                      << std::setw(20) << std::fixed << std::setprecision(12) << ene_curr
                      << std::setw(12) << std::fixed << std::setprecision(8) << grad_norm
                      << std::endl;
          }
          //
          if (ene_curr >= ene_prev) {
            // local minimum reached, undo the last update and break
            for (size_t i = 0; i < params.size(); i++) {
              params[i] += eta * gradient[i];
            }
            break;
            local_result = 0; // MZ: reach local graident minimum
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
          // std::vector<double> lower_bounds(ansatz->numParams(), optimizer_settings.lbound);
          // std::vector<double> upper_bounds(ansatz->numParams(), optimizer_settings.ubound); 
          // optimizer.set_lower_bounds(lower_bounds);
          // optimizer.set_upper_bounds(upper_bounds);
          optimizer.set_lower_bounds(optimizer_settings.lbound); // MZ: use the overload since bounds for all parameters are the same (https://nlopt.readthedocs.io/en/latest/NLopt_C-plus-plus_Reference/#bound-constraints)
          optimizer.set_upper_bounds(optimizer_settings.ubound); // MZ: same as above
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

          // std::cout << "Lower bound" << std::endl;
          // for (auto& i: optimizer.get_lower_bounds()) {
          //   std::cout << i << std::endl;
          // }
          // std::cout << "Upper bound" << std::endl;
          // for (auto& i: optimizer.get_upper_bounds()) {
          //   std::cout << i << std::endl;
          // }
      }


      //------------------------------------------------------------------------




      /**
       * @brief  True Gradient Descent using Finite Differences
       * @note   Similar to follow_fixed_gradient but uses true gradients instead of SPSA
       * @param  x0: Initial parameter vector (input)
       * @param  initial_ene: Initial energy (output)
       * @param  final_ene: Final energy of local minimum (output)
       * @param  num_iterations: Number of iterations before convergence (output)
       * @param  delta: Step size for finite difference gradient computation
       * @param  eta: Gradient descent step size
       * @param  max_iterations: Maximum number of iterations (default: 1000)
       * @param  verbose: Print iteration details 
       * @retval Parameter map with final values
       */
      virtual std::vector<std::pair<std::string, ValType>> follow_true_gradient(const std::vector<ValType>& x0, ValType& initial_ene, ValType& final_ene, 
                                                                                  IdxType& num_iterations, ValType delta, ValType eta, 
                                                                                  IdxType max_iterations = 100, bool verbose = false) {
        Config::PRINT_SIM_TRACE = false;
        local_result = 9; // MZ: initial flag

        // Initialize finite difference gradient estimator
        FiniteDifference fd_estimator;

        // Initialize data structures
        std::vector<ValType> gradient(x0.size());
        std::vector<ValType> params(x0);
        ValType ene_curr = energy(params); // starting energy
        ValType ene_prev = ene_curr; // initialize to current energy
        initial_ene = ene_curr;

        // Adaptive learning rate parameters
        ValType eta_current = eta; // current learning rate
        ValType eta_min = eta * 1e-5; // minimum learning rate
        ValType eta_max = eta * 10.0; // maximum learning rate
        ValType eta_increase_factor = 1.1; // factor to increase eta when making progress
        ValType eta_decrease_factor = 0.5; // factor to decrease eta when overshooting
        int consecutive_improvements = 0; // track consecutive energy improvements

        iteration = 0;

        // Store initial gradient norm for relative convergence check
        ValType initial_grad_norm = 0.0;

        // Gradient descent loop
        do {
          // Compute true gradient using finite differences
          fd_estimator.estimate([&] (const std::vector<double>& xval) { return energy(xval);}, params, gradient, delta, 1);

          // Compute gradient norm for convergence check
          double grad_norm = std::sqrt(std::accumulate(gradient.begin(), gradient.end(), 0.0, [] (ValType a, ValType b) {
            return a + b * b;
          }));

          // Store initial gradient norm for relative threshold
          if (iteration == 0) {
            initial_grad_norm = grad_norm;
          }

          // Print iteration info if verbose
          if ((process_rank == 0) && verbose) {
            if (iteration == 0) {
              std::cout << "\n----------- Finite Difference Gradient Descent -----------\n" << std::left
                        << std::setw(8) << " Iter."
                        << std::setw(20) << "Objective Value"
                        << std::setw(12) << "Grad. Norm"
                        << std::setw(12) << "Learn. Rate"
                        << std::endl;
              std::cout << std::string(53, '-') << std::endl;
            }
            std::cout << std::left << " "
                      << std::setw(7) << iteration
                      << std::setw(20) << std::fixed << std::setprecision(12) << ene_curr
                      << std::setw(12) << std::fixed << std::setprecision(8) << grad_norm
                      << std::setw(12) << std::fixed << std::setprecision(8) << eta_current
                      << std::endl;
          }

          // Check for convergence using relative gradient norm threshold
          ValType relative_grad_threshold = 1e-6 * std::max(initial_grad_norm, 1e-8);
          if (grad_norm < relative_grad_threshold) {
            if (verbose && process_rank == 0) {
              std::cout << "Converged: gradient norm below relative threshold (" << relative_grad_threshold << ")" << std::endl;
            }
            local_result = 0;
            break;
          }

          // Gradient descent step with current learning rate
          for (size_t i = 0; i < params.size(); i++) {
            params[i] -= eta_current * gradient[i];
          }

          // Store previous energy before computing new one
          ene_prev = ene_curr;
          ene_curr = energy(params);

          // Adaptive learning rate adjustment
          if (ene_curr < ene_prev) {
            // Energy decreased - good step
            consecutive_improvements++;

            // If we've had several consecutive improvements, increase learning rate
            if (consecutive_improvements >= 3 && eta_current < eta_max) {
              eta_current = std::min(eta_current * eta_increase_factor, eta_max);
              if (verbose && process_rank == 0) {
                std::cout << "Increasing learning rate to " << eta_current << std::endl;
              }
            }
          } else {
            // Energy increased or stayed same - bad step, reduce learning rate
            consecutive_improvements = 0;

            // Undo the step
            for (size_t i = 0; i < params.size(); i++) {
              params[i] += eta_current * gradient[i];
            }

            // Reduce learning rate
            eta_current = std::max(eta_current * eta_decrease_factor, eta_min);

            if (verbose && process_rank == 0) {
              std::cout << "Energy increased, reducing learning rate to " << eta_current << std::endl;
            }

            // Check if learning rate is too small (convergence)
            if (eta_current <= eta_min) {
              ene_curr = ene_prev; // restore previous energy
              if (verbose && process_rank == 0) {
                std::cout << "Learning rate too small, converged to local minimum" << std::endl;
              }
              local_result = 0; // reach local gradient minimum
              break;
            }

            // Try the step again with reduced learning rate
            for (size_t i = 0; i < params.size(); i++) {
              params[i] -= eta_current * gradient[i];
            }
            ene_curr = energy(params);

            // If still no improvement, accept current position and continue
            if (ene_curr >= ene_prev) {
              for (size_t i = 0; i < params.size(); i++) {
                params[i] += eta_current * gradient[i]; // undo step
              }
              ene_curr = ene_prev; // restore previous energy
              if (verbose && process_rank == 0) {
                std::cout << "Local minimum reached: energy stopped decreasing" << std::endl;
              }
              local_result = 0; // reach local gradient minimum
              break;
            }
          }

          iteration++;

        } while(iteration < max_iterations);

        num_iterations = iteration;
        final_ene = ene_curr;

        // Set final parameters in ansatz
        ansatz->setParams(params);

        if (verbose && process_rank == 0) {
          std::cout << "True gradient descent completed after " << iteration << " iterations" << std::endl;
          std::cout << "Final energy: " << final_ene << std::endl;
          std::cout << "Final learning rate: " << eta_current << " (initial: " << eta << ")" << std::endl;
        }

        return ansatz->getFermionicOperatorParameters();
      }

      //------------------------------------------------------------------------



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
      int get_numevals() const { return num_evals;}; // MZ
      double get_duration() const { return opt_duration;}; // MZ
      void set_duration(double value) { opt_duration = value; }; //MZ

      double get_comm_duration() const { return comm_duration;}; // MZ, for adapt-vqe
      void set_comm_duration(double value) { comm_duration = value; }; //MZ, for adapt-vqe

      size_t get_numpauli() const { return num_pauli_terms_total;}; // MZ, for adapt-vqe
      void set_numpauli(size_t value) { num_pauli_terms_total = value; }; //MZ, for adapt-vqe

      size_t  get_numcomm() const { return num_comm_cliques;}; // MZ, for adapt-vqe
      void set_numcomm(size_t value) { num_comm_cliques = value; }; //MZ, for adapt-vqe

      int get_adaptrounds() const {return num_adapt_rounds;}; // MZ: for adapt-vqe
      void set_adaptrounds(int rounds) { num_adapt_rounds = rounds;}; // MZ: for adapt-vqe

      int get_adaptresult() const {return adapt_result;}; // MZ: for get adapt result flag
      void set_adaptresult(int res) { adapt_result = res;}; // MZ: for set adapt result flag

      int get_localresult() const {return local_result;}; // MZ: for get Local Gradient Follower flag
      void set_localresult(int res) { local_result = res;}; // MZ: for set Local Gradient Follower flag
      protected:
        std::shared_ptr<Ansatz> ansatz;                    // state preparation circuit
        std::shared_ptr<Ansatz> measurement;               // circuit to measure expectation values
        std::shared_ptr<Hamiltonian> hamil;                // target system Hamiltonian
        SPSA g_est;                                        // stochastic gradient estimator
        // FiniteDifference g_est;                            // finite difference gradient estimator
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
        nlopt::result opt_result;                          // MZ: optimzation success or fail and the reason
        int num_evals;                                     // MZ: total number of evaluations
        double opt_duration;                               // MZ: time the optimization
        double comm_duration;                              // MZ: time construction of ADAPT-VQE Commutators
        size_t num_pauli_terms_total;                      // MZ: Total number of Pauli terms in the commutator, for ADAPT-VQE 
        size_t num_comm_cliques;                           // MZ: Total number of commuting cliques, for ADAPT-VQE 
        int num_adapt_rounds;                              // MZ: total number of ADAPT rounds
        int adapt_result;                                  // MZ: save adapt results: 
                                                           //    0 -> Reach gradient norm tolerance
                                                           //    1 -> Reach function tolerance
                                                           //    2 -> Reach maximum # ADAPT iterations
                                                           //    9 -> ADAPT iteration is not run successfully
        int local_result;                                  // MZ: save local gradient flow results: 
                                                           //    0 -> Reach local gradient minimum
                                                           //    9 -> Local Gradient Follower is not run
      

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