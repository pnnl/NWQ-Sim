# NWQ-VQE

Production-ready Variational Quantum Eigensolver (VQE) module built on the NWQ-Sim platform.

## Features

- **High Performance**: GPU-accelerated quantum circuit simulation
- **Multiple Algorithms**: UCCSD and ADAPT-VQE implementations  
- **Python API**: Interface with error handling
- **Flexible Optimization**: Integration with NLopt and SciPy optimizers

- **Documentation**: API documentation with examples

## Building

This module builds together with NWQ-Sim. Follow the project-wide instructions in `doc/user_manual.md` to configure and compile the tree. If you cloned without `--recursive`, fetch the bundled NLopt source before running CMake:

```bash
git submodule update --init --recursive vqe/nlopt
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
INFORMATIONAL
  -h, --help            Show usage information and exit.
  -l, --list-backends   List compiled back-ends and exit.

REQUIRED
  -f, --hamiltonian     Path to Hamiltonian file (sum of Fermionic operators).
  -p, -n, --nparticles  Number of electrons in the system.

HAMILTONIAN / ANSATZ / BACKEND OPTIONS
  --ducc                Use DUCC orbital indexing (default).
  --xacc                Use XACC/Qiskit orbital indexing.
  --sym, --symm         UCCSD symmetry level (0 none, 1 spin, 2 orbital, 3 full).
  -b, --backend         Simulation backend (CPU | GPU). Defaults to CPU.
  --seed                Seed for initial parameters and stochastic gradients.

OPTIMISER OPTIONS
  -v, --verbose         Print additional optimisation progress.
  -o, --optimizer       NLopt algorithm name (default LN_COBYLA).
  --opt-config          JSON file with NLopt parameter overrides.
  -lb, --lbound         Lower bound for all parameters (default -2π).
  -ub, --ubound         Upper bound for all parameters (default  2π).
  --reltol              Relative objective tolerance (default disabled).
  --abstol              Absolute objective tolerance (default disabled).
  --maxeval             Maximum objective evaluations (default 100).
  --maxtime             Wall-clock limit for the optimiser (seconds).
  --stopval             Objective value threshold for early exit.

SIMULATION OPTIONS
  --num_threads         Requested CPU thread count (advisory).
  --disable_fusion      Disable gate fusion (advisory).

ADAPT-VQE OPTIONS
  --adapt               Enable ADAPT-VQE loop.
  -ag, --adapt-gradtol  Gradient-norm termination threshold (default 1e-3).
  -af, --adapt-fvaltol  Energy-change termination threshold (default disabled).
  -am, --adapt-maxeval  Maximum ADAPT iterations (default 20).
```

Example – plain VQE run on the supplied `H4_4_0.9_xacc.hamil` benchmark:

```bash
./build/vqe/nwq_vqe -f vqe/example_hamiltonians/H4_4_0.9_xacc.hamil -p 4 --maxeval 1000
```

Example – Fermionic ADAPT-VQE with the same problem:

```bash
./build/vqe/nwq_vqe -f vqe/example_hamiltonians/H4_4_0.9_xacc.hamil -p 4 --adapt --adapt-maxeval 50
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

# 4) Configure optimization
options = vqe.Options()
options.optimizer = "L-BFGS-B"
options.max_evaluations = 1000
options.use_gpu = True
options.relative_tolerance = 1e-8

# 5) Run VQE
result = vqe.run(ansatz, pauli_terms, options)
print(f"Ground state energy: {result.energy:.8f}")
print(f"Converged: {result.success} ({result.evaluations} evaluations)")
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
adapt_options.adapt_max_iterations = 20
adapt_options.adapt_gradient_tolerance = 1e-4
adapt_options.adapt_energy_tolerance = 1e-8

# Configure inner VQE optimization for each ADAPT step
adapt_options.adapt_optimizer = "L-BFGS-B"
adapt_options.adapt_max_evaluations = 500
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
print(f"Circuit depth: {result.circuit_depth}")
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
- ✅ **Error handling**
- ✅ **Input validation** for all configuration parameters  
- ✅ **String representations** (`__repr__`, `__str__`) for all classes
- ✅ **Property-based access** with getters/setters instead of raw attributes
- ✅ **Parallel computation support**
- ✅ **Clean external optimizer interface** via `EnergyEvaluator`
- ✅ **GPU auto-fallback** with warnings instead of crashes
- ✅ **Type hints and documentation** throughout the API



## C++ Library Usage

The shared library `libnwqsim` exports the QFlow ABI defined in `include/nwqsim_qflow.hpp`. It provides:

- `parseHamiltonianFile` – load a Fermionic Hamiltonian into the coefficient/operator form.
- `qflow_nwqsim` – run a UCCSD-min VQE and obtain the ground-state energy plus the parameter list. Pass an optional `vqe::vqe_options` instance (from `vqe_options.hpp`) to override defaults such as optimizer, bounds, or GPU usage.
- `get_termination_reason_local` – translate ADAPT termination codes to human-readable descriptions.

Link against `libvqe` for native access to the solver primitives (environment, ansatz builders, state-vector backends, and the `execution` runners).

Example with custom options:

```cpp
#include "nwqsim_qflow.hpp"
#include "vqe_options.hpp"

vqe::vqe_options opts;
opts.max_evaluations = 500;
opts.optimizer = nlopt::LN_BOBYQA;

auto [energy, parameters] = qflow_nwqsim(ham_ops, n_electrons, "GPU", opts);
```
