# VQE States
The `VQEState` class defined in vqe/include/vqe_state.hpp implements the VQE iteration loop. It takes the problem, ansatz, and assorted parameters as input and returns an optimized circuit (though probably not globally optimal). 

The base `VQEState` class is a general API, however it requires subclasses to implement specific backends. Backend architectures currently supported are:
1. CPU (vqe/include/svsim_vqe/sv_cpu_vqe.hpp)
2. NVIDIA GPU (vqe/include/svsim_vqe/sv_cuda_vqe.hpp)
3. CPU MPI (vqe/include/svsim_vqe/sv_mpi_vqe.hpp)
4. NVIDIA GPU MPI (via NVSHMEM) (vqe/include/svsim_vqe/sv_cuda_mpi_vqe.hpp)

Each of these has backend-specific quirks, which we address in turn  

## Base API
We start with the base API, then move into the backend-specific implementations.

### Properties and Data Structures
The following are `protected` properties inherited by subclasses, and are used for general VQE functions.

The `ansatz` circuit is the general state preparation circuit, described in circuits.md. `measurement` is a composition of measurement circuits and EXPECT gates (see measurement.md and expectation_values.md) to compute the expectation values for the target Hamiltonian `hamil`.
```c++
std::shared_ptr<Ansatz> ansatz;                    // state preparation circuit
std::shared_ptr<Ansatz> measurement;               // circuit to measure expectation values
std::shared_ptr<Hamiltonian> hamil;                // target system Hamiltonian
```

The `obsvec`, `zmasks`, and `coeffs` structures are used for expectation estimation. `zmasks` and `coeffs` stores the Z stabilizer arrays and sign-adjusted coefficients for each diagonalized Pauli operator. `obsvec` contains an array of pointers to individual `ObservablesList` structs, each of which computes the expectation value for a QWC Pauli group.
```c++
std::vector<ObservableList*> obsvec;               // vector of pointers to ObservableList structures for clique expectation value calculation
std::vector<std::vector<IdxType> > zmasks;         // vector of diagonalized operator zmasks
std::vector<std::vector<ValType> > coeffs;         // vector of diagonalized operator coefficients
``` 

Gradient-based optimization is handled by an SPSA estimator (see optimization.md), with the `compute_gradient` flag used to indicate whether the chosen optimizer was gradient-based (e.g. LD_MMA) or gradient-free (e.g. LN_COBYLA).
```c++
SPSA g_est;                                        // stochastic gradient estimator
bool compute_gradient;                             // flag to compute gradient (depends on optimizer selected)
```
`process_rank` is used to ensure that only the root node prints messages in MPI settings. Non-MPI subclasses initialize `process_rank=0`. `callback` is a function to print optimization status messages back to the user, with the prototype:

`typedef std::function<void(const std::vector<ValType>& x, ValType ene, IdxType iter)> Callback;` 
```c++
IdxType process_rank;                              // process rank (used by MPI/NVGPU_MPI backends)
Callback callback;                                 // callback function for terminal updates
```

The remaining properties are used for optimization. 
```c++
OptimizerSettings optimizer_settings;              // NLOpt optimizer settings (bounds, termination critera)
IdxType iteration;                                 // current iteration counter

double expectation_value;                          // last computed expectation value
nlopt::algorithm optimizer_algorithm;              // NLOpt optimization algorithm for circuit updates 
```

The `OptimizerSettings` struct (defined in vqe/include/utils.hpp) contains criteria for optimizer operation and termination, shown below. The `parameter_map` data structure holds algorithm-specific options, see `vqe/mma_config.json` for an example pertaining to the LD_MMA optimizer.
```c++
struct OptimizerSettings {
  /**
  * @brief  Data structure for NLOpt optimizer settings
  */
  ValType rel_tol; // relative tolerance cutoff
  ValType abs_tol; // absolute tolerance cutoff
  ValType stop_val; //
  IdxType max_evals; // Max number of function evaluations
  ValType max_time; // Optimizer timeout (seconds)
  ValType lbound; // Lower bound for optimization
  ValType ubound; // Upper bound for optimization
  std::unordered_map<std::string, ValType> parameter_map; // map for setting optimizer-specific parameters
  // Defaults (disables all of the settings, except for the max_eval ceiling)
  OptimizerSettings(): rel_tol(-1), 
                        abs_tol(-1),
                        stop_val(-MAXFLOAT),
                        max_evals(50),
                        max_time(-1),
                        lbound(-PI),
                        ubound(PI) {}
};
```
### Constructor
The `VQEState` constructor is defined as follows:
```c++
VQEState( std::shared_ptr<Ansatz> a, 
                  std::shared_ptr<Hamiltonian> h, 
                  nlopt::algorithm _optimizer_algorithm,
                  Callback _callback,
                  IdxType seed = 0,
                  OptimizerSettings opt_settings = OptimizerSettings())
```
Note that since we store a pointer to the `Ansatz` and `Hamiltonian` objects, we can modify them externally without re-contructing the `VQEState`. This is particularly useful in ADAPT-VQE.

To fully initialize the VQEState, the `initialize()` function is called after the base class constructor to build the measurement circuit and build the data structures required for energy estimation.

Note that `initialize` calls the following functions which backends are required to implement:
1. `void allocate_observables(IdxType size)`: Make space for the base `ObservableList` instances and any subclass-specific data structures. Necessary to properly allocate device memory for NVGPU-based backends
2. `void fill_obslist(IdxType index)`: Initialize a single `ObservableList` entry and add the EXPECT gate to the measurement circuit. This NEEDS to occur in a specific order:
```c++
// Make the common Pauli operator by taking a bitwise "or" over stabilizer states
PauliOperator common = make_common_op(pauli_list, 
                                      zmasks[index], 
                                      coeffs[index]);
Measurement circ1(common);
measurement->compose(circ1, mapping);
fill_obslist(index); // backend-specific observable construction/initialization routines
Measurement circ2(common, true);
measurement->compose(circ2, mapping); // append the measurement inverse to undo diagonalization
```
where the call is sandwiched between a Pauli diagonalization circuit and its inverse. Here `pauli_list` is a QWC group obtained from the `Hamiltonian`.

### VQE Loop
The VQE algorithm begins with a call to the `optimize` function:
```c++
...
std::vector<double> parameters(ansatz.numParams(), 0);
double fval;
VQEState state(ansatz, hamiltonian, algorithm, callback, seed, settings);
state.optimize(parameters, fval);
```
which sets the parameters for an `nlopt::opt optimizer` object from `optimizer_settings` and then enters the NLOpt-library loop:
```c++
// Initialize optimizer 
nlopt::opt optimizer = nlopt::opt(optimizer_algorithm, ansatz->numParams());
// set the objective function
optimizer.set_min_objective(nl_opt_function, (void*)this);
  ...
nlopt::result optimization_result = optimizer.optimize(parameters, fval);
```
The optimizer uses the target function
```c++
// Target function for NLOpt object, has to match function prototype
double nl_opt_function(const std::vector<double>& x, std::vector<double>& gradient, void* val) {
    return ((VQEState*)val)->cost_function(x, gradient);
};
```
which calls the `cost_function` method of the base class. `cost_function` estimates the gradient (if necessary), calls `energy(x)`, then prints calls the `callback` function.

`energy` first assigns the parameters `x` to the `Ansatz`, then invokes the backend-specific `call_simulator()` function. 

`call_simulator()` first simulates `ansatz` to prepare the trial state, then `measurement` to diagonalize and measure the Pauli strings.
```c++
virtual ValType energy(const std::vector<double>& x) {
  ansatz->setParams(x);

  call_simulator();
  ValType expectation = hamil->getEnv().constant + expectation_value;
  return expectation;
}
```

This loop continues until an NLOpt termination criterion is met, such as the iteration limit or some convergence criterion.

### Backend-Specific Prototypes
The following functions are merely prototypes or CPU defaults defined in the base class, and are expected to be overriden by subclasses:
```c++
virtual void call_simulator() {};
```
Calls the `SVSim` instance associated with the subclass.

```c++
virtual void call_simulator(std::shared_ptr<Ansatz> _measurement, bool reset) {};
```
Calls the `SVSim` instance associated with the subclass, but use the provided measurement circuit instead of the class property `measurement`. Used by ADAPT-VQE for commutator gradient calculations.

```c++
virtual void set_exp_gate(std::shared_ptr<Ansatz> circuit, ObservableList* o, std::vector<IdxType>& zmasks, std::vector<ValType>& coeffs)
```
Append an EXPECT gate to the circuit and initialize the `ObservableList` object with the stabilizer and coefficient data. Used by ADAPT-VQE for commutator gradient calculations.

```c++
virtual void allocate_observables(ObservableList*& observables, IdxType size)
```
Allocate space for `size` entries (on host and/or device) for the provided  `observables` pointer.


```c++
virtual void allocate_observables(IdxType size)
```
Allocate space for `size` entries (on host and/or device) for the class `obs_vec` and other observable-related data structures (may be subclass specific).


```c++
virtual void delete_observables(ObservableList* observables, IdxType size)
```
Delete the data associated with `observables`, could be host or device storage depending on backend.
