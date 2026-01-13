# NWQ-VQE

Production-ready Variational Quantum Eigensolver (VQE) module built on the NWQ-Sim platform.

## Features

- **High Performance**: GPU-accelerated quantum circuit simulation with MPI parallelization for ADAPT-VQE
- **Multiple Algorithms**: UCCSD and ADAPT-VQE implementations with adaptive operator selection
- **Checkpoint and Resumption**: Save and restore VQE/ADAPT-VQE state for long-running optimizations
- **Flexible Optimization**: Integration with NLopt optimizers including SPSA gradient estimation and numerical gradients
- **Parameter Management**: Custom initial parameters from files or values, parameter saving to disk
- **Python API**: Interface with error handling and convenience functions
- **Advanced Logging**: Iteration details, memory tracking, convergence monitoring, and performance profiling
- **Documentation**: API documentation with examples

## Building

This module will not automatically build together with NWQ-Sim. Follow the project-wide instructions in `doc/user_manual.md` to configure and compile the tree. If you cloned without `--recursive`, fetch the bundled NLopt and pybind11 sources before running CMake:

```bash
git submodule update --init --recursive
```

Importantly, `-DNWQSIM_ENABLE_VQE=ON` should be added in the `make` command to enable VQE function. For example,

```bash
cmake .. -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC \
 -DNWQSIM_ENABLE_VQE=ON \        # enable VQE
 -DVQE_ENABLE_MPI=ON \           # enable MPI
 -DCMAKE_CUDA_HOST_COMPILER=CC \ # enable NVDIA GPU backend
 -DCMAKE_BUILD_TYPE=Release
```

The build system will take care of compiling NLopt during the normal configure/build step.

### Python Installation

The VQE module is packaged as an optimized Python package. Install with:

```bash
pip install .
```

This automatically handles all dependencies including `scikit-build-core`, `pybind11`, CMake, and compilers. The installation creates the `nwqsim_vqe` package with proper library bundling and rpath configuration.

**Optional Dependencies:**
- For development: `pip install .[dev]`
- For SciPy optimization: `pip install scipy numpy` (user responsibility)

**Installation Options:**
- Production: `pip install .` (recommended)
- Development: Configure with CMake and add `build/vqe` to `PYTHONPATH`
- C++ only: Set `-DVQE_BUILD_PYTHON=OFF` during CMake configuration

Running `cmake --build <build-dir>` produces:

- `build/vqe/nwq_vqe` – command-line VQE driver.
- `build/vqe/libvqe.{so,dylib}` – shared library exposing the solver primitives.
- `build/vqe/libnwqsim.{so,dylib}` – QFlow-facing wrapper compatible with third-party integrations.

## Command-Line Usage

Launch VQE via the `nwq_vqe` executable. The two mandatory arguments are the Hamiltonian file and electron count.

```
"INFORMATIONAL"
"  -h, --help            Show help menu."
"  -l, --list-backends   List available backends and exit."
"REQUIRED"
"  -f, --hamiltonian     Path to the input Hamiltonian file (sum of Fermionic operators)."
"  -p, -n, --nparticles  Number of electrons in molecule."
"OPTIONAL (Hamiltonian, Ansatz and Backend)"
"  --xacc                Enable XACC/Qiskit indexing scheme (default)."
"  --ducc                Enable DUCC indexing scheme."
"  --sym, --symm         UCCSD symmetry level (0->none, 1->spin, 2->orbital, 3->full). Default to 3."
"  -b, --backend         Simulation backend (CPU, NVGPU, AMDGPU, MPI (not supported)). Defaults to CPU."
"  --seed                Random seed for reproducibility."
"OPTIONAL (Global Minimizer)"
"  -v, --verbose         Print additional progress information."
"  -o, --optimizer       NLopt optimizer (e.g. LN_COBYLA, LN_BOBYQA, LN_NEWUOA, LD_LBFGS)."
"  --opt-config          Path to JSON file with optimizer parameter overrides."
"  -lb, --lbound         Optimizer lower bound (default -π)."
"  -ub, --ubound         Optimizer upper bound (default π)."
"  --reltol              Relative tolerance termination criterion."
"  --abstol              Absolute tolerance termination criterion (default 1e-6)."
"  --grad-step           Forward-difference gradient step for derivative-based optimizers (VQE and ADAPT inner solves) (default 1e-5)."
"  --stopval             Objective stop value."
"  --maxeval             Maximum number of function evaluations in VQE optimization (default 100)."
"  --maxtime             Maximum VQE optimizer time (seconds)."
"  --spsa                Enable SPSA gradient estimation (2 evals) instead of forward difference (N+1 evals)."
"  --init-params         [For VQE only] Initial parameters: single value (repeat for all) or file with comma-separated values."
"  --save-params         [For VQE only] Save optimized parameters to {hamiltonian_path}-vqe_params.txt."
"OPTIONAL (ADAPT-VQE)"
"  --adapt               Enable ADAPT-VQE instead of standard VQE."
"  -ag, --adapt-gradtol  Operator Gradient norm tolerance (Default 1e-3)."
"  --adapt-grad-step     Central-difference step for ADAPT operator gradients (default 1e-4)."
"  -af, --adapt-fvaltol  Energy change tolerance (default disabled)."
"  -am, --adapt-maxeval  Maximum ADAPT iterations (default 50)."
"  -ak, --adapt-batch-k  Maximum operators appended per ADAPT iteration (default 1)."
"  -at, --adapt-tau      Gradient threshold fraction for batched selection (default 1.0)."
"  -as, --adapt-save     Save parameters every iteration to {hamiltonian_path}-adapt_params.txt."
"  -al, --adapt-load     Load ADAPT-VQE state from file to resume optimization."
"SIMULATOR OPTIONS"
"  --num_threads         Specify number of threads (ignored in current backend)."
"  --disable_fusion      Disable gate fusion (ignored in current backend).";
```

Example – plain VQE run on the supplied `H4_4_0.9_xacc.hamil` benchmark:

```bash
./build/vqe/nwq_vqe -f vqe/example_hamiltonians/H4_4_0.9_xacc.hamil -p 4 --maxeval 1000
```

Example – VQE with custom initial parameters and parameter saving:

```bash
./build/vqe/nwq_vqe -f vqe/example_hamiltonians/H4_4_0.9_xacc.hamil -p 4 --init-params init_params.txt --maxeval 1000 --save-params
```

Example – Fermionic ADAPT-VQE with the same problem:

```bash
./build/vqe/nwq_vqe -f vqe/example_hamiltonians/H4_4_0.9_xacc.hamil -p 4 --adapt --adapt-maxeval 50
```

Example – ADAPT-VQE with checkpointing (the parameter values are saved under the same parent folder as the Hamiltonian file):

```bash
./build/vqe/nwq_vqe -f vqe/example_hamiltonians/H4_4_0.9_xacc.hamil -p 4 --adapt --adapt-save --adapt-log-memory
# To resume from checkpoint:
./build/vqe/nwq_vqe -f vqe/example_hamiltonians/H4_4_0.9_xacc.hamil -p 4 --adapt --adapt-load qe/example_hamiltonians/H4_4_0.9_xacc.hamil-adapt_params.txt
```

Example – ADAPT-VQE with SPSA gradient and GPU acceleration:

```bash
./build/vqe/nwq_vqe -f vqe/example_hamiltonians/H4_4_0.9_xacc.hamil -p 4 --adapt -b NVGPU --spsa --adapt-maxeval 50
```

Example - Batched-ADAPT-VQE: include at most `-ak` number of operators in each outter ADAPT iteration if their operator gradient norms are larger than $\tau |g_{max}|$ where $ \tau > 0$ is defined by `-at`.

```bash
./build/vqe/nwq_vqe -f vqe/example_hamiltonians/H4_4_0.9_xacc.hamil -p 4 --adapt -ak 10 -at 0.5
```



## Python API

The `nwqsim_vqe` package provides a Python interface with error handling and type hints.

### Quick Start

```python
import nwqsim_vqe as vqe

# Check system capabilities
print("Available optimizers:", vqe.list_optimizers())
print("GPU support:", vqe.check_gpu_support())
print("Version:", vqe.__version__)

# 1) Create molecular system (convenience function)
system = vqe.system(n_spatial_orbitals=2, n_electrons=2)

# 2) Load and transform Hamiltonian  
hamiltonian = vqe.load_hamiltonian("h2_sto3g.txt")
pauli_terms = vqe.transform(hamiltonian)

# 3) Build UCCSD ansatz
ansatz = vqe.UCCSD(system, trotter_steps=1, symmetry_level=3)
ansatz.build()

# 4) Configure optimization with custom initial parameters
options = vqe.Options()
options.optimizer = "L-BFGS-B"
options.max_evaluations = 1000
options.use_gpu = True
options.relative_tolerance = 1e-8
options.initial_parameters = [0.1, -0.2, 0.05]  # Custom starting parameters (must match ansatz.parameter_count)

# 5) Run VQE
result = vqe.run(ansatz, pauli_terms, options)
print(f"Ground state energy: {result.energy:.8f}")
print(f"Converged: {result.success} ({result.evaluations} evaluations)")
print(f"Final parameters: {result.parameters}")
```

### SciPy Integration

Use the `EnergyEvaluator` directly with SciPy optimizers for maximum flexibility:

```python
import numpy as np
from scipy.optimize import minimize

# Create energy evaluator
evaluator = vqe.EnergyEvaluator(ansatz, pauli_terms, use_gpu=True)

# Set up initial parameters
x0 = np.zeros(evaluator.parameter_count)

# Define objective function
def objective(params):
    return evaluator.energy(params.tolist())

# Use any SciPy optimizer
result = minimize(objective, x0, method="L-BFGS-B", 
                 options={"maxiter": 1000})

print(f"Optimal energy: {result.fun:.8f}")
print(f"Converged: {result.success}")
print(f"Function evaluations: {result.nfev}")
```

### Usage Examples

#### Energy Evaluation
```python
# Create reusable energy evaluator for optimization loops
evaluator = vqe.EnergyEvaluator(ansatz, pauli_terms, use_gpu=True)

# Energy evaluation  
energy = evaluator.energy([0.1, -0.2, 0.3])

# Or make it callable
energy_func = evaluator  # EnergyEvaluator is callable
energy = energy_func([0.1, -0.2, 0.3])
```

#### Data Structures
```python
# All classes have meaningful __repr__ and __str__
print(system)  # System(orbitals=2, electrons=2, qubits=4)
print(ansatz)  # UCCSD(qubits=4, parameters=2)
print(result)  # VQEResult(energy=-1.85727503, evaluations=247, converged=True)

# Pauli terms with readable output  
for term in pauli_terms[:3]:
    print(term)  # 0.174073 IIXX, -0.174073 IIXY, etc.
```

#### Configuration Classes
```python
options = vqe.Options()
try:
    options.max_evaluations = -100  # Raises ValueError
except ValueError as e:
    print(f"Validation error: {e}")

# Configuration display
print(options)  # Options(optimizer='L-BFGS-B', max_eval=1000, gpu=False)
```

### ADAPT-VQE Configuration

ADAPT-VQE provides adaptive circuit construction with configuration options:

```python
import nwqsim_vqe as vqe

# Load system using convenience functions
system = vqe.system(n_spatial_orbitals=4, n_electrons=4)
hamiltonian = vqe.load_hamiltonian("h4_square.txt")

# Configure ADAPT with the unified options object
adapt_options = vqe.Options()
adapt_options.mode = "adapt"
adapt_options.use_xacc_indexing = True
adapt_options.adapt_max_iterations = 50
adapt_options.adapt_gradient_tolerance = 1e-4
adapt_options.adapt_energy_tolerance = 1e-8

# Configure inner VQE optimization for each ADAPT step
adapt_options.adapt_optimizer = "L-BFGS-B"
adapt_options.adapt_max_evaluations = 500
adapt_options.adapt_max_time = 600.0    # Time limit per ADAPT iteration in seconds
adapt_options.use_gpu = True

# Run ADAPT-VQE
result = vqe.adapt(
  hamiltonian_file="h4_square.txt",
  n_electrons=4,
  options=adapt_options
)

# Result analysis
print(f"Final energy: {result.energy:.8f}")
print(f"ADAPT iterations: {result.iterations}")
print(f"Total evaluations: {result.energy_evaluations}")
print(f"Selected operators: {result.selected_labels}")
print(f"Converged: {result.success}")
```


### Error Handling

The API provides error handling:

```python
try:
    # Invalid system configuration  
    system = vqe.system(n_spatial_orbitals=0, n_electrons=2)
except ValueError as e:
    print(f"Configuration error: {e}")

try:
    # File not found with context
    hamiltonian = vqe.load_hamiltonian("nonexistent.txt") 
except FileNotFoundError as e:
    print(f"File error: {e}")

try:
    # Parameter check
    options = vqe.Options()
    options.relative_tolerance = -1.0  # Invalid
except ValueError as e:
    print(f"Validation error: {e}")
```

## API Reference Summary

### Core Classes
- **`System`**: Molecular system specification
- **`UCCSD`**: UCCSD ansatz builder with parameter management  
- **`Options`**: Unified configuration options
- **`VQEResult`** / **`AdaptResult`**: Results with metadata
- **`EnergyEvaluator`**: Energy evaluation
- **`HamiltonianData`** / **`PauliTerm`**: Data structures with string representations

### Key Functions
- **`system()`**: Convenience constructor
- **`load_hamiltonian()`**: File loading with error context
- **`run_vqe()` / `run_vqe_with_ansatz()`**: Main VQE entry points
- **`adapt()`**: ADAPT-VQE with adaptive circuit construction
- **`EnergyEvaluator`**: Interface for external optimizers
- **`list_optimizers()` / `check_gpu_support()`**: System introspection

### Enumerations (Python-style)
- **`Spin.UP`** / **`Spin.DOWN`**: Electron spin states
- **`OrbitalKind.OCCUPIED`** / **`OrbitalKind.VIRTUAL`**: Orbital types
- **`OperatorKind.CREATION`** / **`OperatorKind.ANNIHILATION`**: Fermion operators

### Key Features
- **Error handling**
- **Input validation** for all configuration parameters  
- **String representations** (`__repr__`, `__str__`) for all classes
- **Property-based access** with getters/setters instead of raw attributes
- **Parallel computation support**
- **Clean external optimizer interface** via `EnergyEvaluator`
- **GPU auto-fallback** with warnings instead of crashes
- **Type hints and documentation** throughout the API



## C++ Library Usage

The shared library `libnwqsim` exports the QFlow ABI defined in `include/nwqsim_qflow.hpp`. It provides:

- `parseHamiltonianFile` – load a Fermionic Hamiltonian into the coefficient/operator form.
- `qflow_nwqsim` – run a UCCSD VQE and obtain the ground-state energy plus the parameter list. Pass an optional `vqe::vqe_options` instance (from `vqe_options.hpp`) to override defaults such as optimizer, bounds, GPU usage, or initial parameters.

Link against `libvqe` for native access to the solver primitives (environment, ansatz builders, state-vector backends, and the `execution` runners).


Example with custom options:

```cpp
#include "nwqsim_qflow.hpp"
#include "vqe_options.hpp"

vqe::vqe_options opts;
opts.verbose = false; // Set to true for more information
opts.trotter_steps = 1;
opts.symmetry_level = 3;
opts.lower_bound = -kPi; // -pi
opts.upper_bound = kPi; // pi
opts.max_evaluations = 100; // Max number of optimization iterations
opts.relative_tolerance = -1.0; // -1 for no relative tolerance
opts.absolute_tolerance = 1e-8; // -1 for no absolute tolerance
opts.max_time = -1.0; // -1 for no max time
opts.optimizer = nlopt::LD_LBFGS; // Use a derivative-based optimizer to enable gradient computation
opts.initial_parameters = {0.1, 0.2, 0.3, ...};  // Must match parameter count for the used ansatz (UCCSD)

auto [energy, parameters] = qflow_nwqsim(ham_ops, n_electrons, "NVGPU", opts);
```




## Advanced CLI Features

### Parameter Management

Initial parameters can be provided as a starting point for VQE optimization. This is useful when you have a good initial guess from prior optimization or another method.

**From file (CLI):**
```bash
./build/vqe/nwq_vqe -f hamiltonian.txt -p 4 --init-params params.txt
```

The file should contain comma-separated values, one per parameter. If the file has fewer parameters than required, missing values are filled with zeros. If it has more, only the first N values are used.

**From single value (CLI):**
```bash
./build/vqe/nwq_vqe -f hamiltonian.txt -p 4 --init-params 0.5
```

This repeats the value for all parameters in the ansatz.

**Saving optimized parameters (CLI):**
```bash
./build/vqe/nwq_vqe -f hamiltonian.txt -p 4 --maxeval 1000 --save-params
```

Creates `hamiltonian.txt-vqe_params.txt` with the final optimized parameters.

### Gradient Methods

The module supports multiple gradient estimation methods for optimization:

**Numerical gradients (default for derivative-based optimizers):**
- Uses forward finite difference with step size 1e-5
- Requires N+1 function evaluations per gradient
- More accurate for smooth functions

**SPSA gradient (Simultaneous Perturbation Stochastic Approximation):**
```bash
./build/vqe/nwq_vqe -f hamiltonian.txt -p 4 --spsa
```

- Only 2 function evaluations per gradient regardless of dimension
- Beneficial for large parameter spaces
- Useful when evaluation cost is very high

### ADAPT-VQE Checkpointing

Long-running ADAPT-VQE optimizations can be saved and resumed:

**Saving parameters every iteration:**
```bash
./build/vqe/nwq_vqe -f hamiltonian.txt -p 4 --adapt --adapt-save
```

Creates `hamiltonian.txt-adapt_params.txt` with full ADAPT state including:
- Selected operators and their indices
- Optimized parameters at each iteration
- Energy values
- Metadata (pool type, symmetry level, iteration number)

**Resuming from checkpoint:**
```bash
./build/vqe/nwq_vqe -f hamiltonian.txt -p 4 --adapt --adapt-load hamiltonian.txt-adapt_params.txt
```

The optimizer validates that symmetry settings match and resumes from the saved iteration.

### Memory and Performance Monitoring

**Memory tracking:**
```bash
./build/vqe/nwq_vqe -f hamiltonian.txt -p 4 --adapt --adapt-log-memory
```

Logs RSS memory usage per ADAPT iteration to track memory growth.

**Verbose output:**
```bash
./build/vqe/nwq_vqe -f hamiltonian.txt -p 4 --verbose
```

Provides detailed iteration information including:
- Parameter statistics (max absolute value, L2 norm)
- Energy changes between iterations
- Selected operators (for ADAPT)
- Gate count statistics
- Timing information

### Parallel ADAPT-VQE

Operator gradient computation in ADAPT-VQE can be parallelized using MPI:

```bash
mpirun -np 8 ./build/vqe/nwq_vqe -f hamiltonian.txt -p 4 --adapt --adapt-maxeval 50
```

Each MPI rank evaluates a subset of operators, reducing wall-clock time for large operator pools.

## Configuration Reference

### CLI Defaults and Behavior

**Symmetry Level:** Default is now 3 (full symmetry), forcing equivalentance of parameters corresponding to symmetric excitation operators. Values range from 0 (no symmetry) to 3 (full).

**Parameter Bounds:** Default bounds are -π to π. Through the range [-0.9, 0.9] is suitable for most UCCSD ansatze.

**ADAPT-VQE Iterations:** Default maximum is 50 iterations.

**Tolerances:**
- Gradient tolerance (ADAPT): 1e-3 (convergence when gradient norm is below this)
- Energy tolerance (ADAPT): Disabled by default. When enabled, ADAPT stops when energy change between iterations drops below this threshold.

**Time Limits:** Specified in seconds for VQE optimizations in full UCCSD VQE function and ADAPT-VQE function. A value less than or equal to 0 means unlimited. Useful for resource-constrained environments.

### VQE-specific vqe_options

- `initial_parameters`: Custom starting parameters for optimization
- `max_evaluations`: Function evaluation limit
- `max_time`: Wall-clock time limit
- `optimizer`: NLopt algorithm (LD_LBFGS, LN_COBYLA, etc.)
- `relative_tolerance` / `absolute_tolerance`: Convergence criteria

### ADAPT-specific vqe_options

- `adapt_initial_parameters`: Starting parameters for inner VQE in first ADAPT iteration
- `adapt_max_evaluations`: Function evaluations per inner VQE solve
- `adapt_max_time`: Time limit per inner VQE solve
- `adapt_gradient_step`: Step size for central difference gradient (default 1e-4)
- `adapt_log_memory`: Enable memory tracking
- `adapt_save_params`: Auto-save state every iteration
- `adapt_load_state_file`: Resume from checkpoint
- `use_spsa_gradient`: Enable SPSA gradient method
