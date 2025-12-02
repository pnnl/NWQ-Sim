# Python VQE API Reference

This document provides a complete reference for the `nwqsim_vqe` Python package.

## Installation

```bash
cd vqe/
pip install -e .
```

## Quick Start

```python
import nwqsim_vqe

# Create molecular system
system = nwqsim_vqe.system(n_spatial_orbitals=2, n_electrons=2)

# Build ansatz
ansatz = nwqsim_vqe.UCCSD(system, trotter_steps=1)
ansatz.build()

# Load Hamiltonian and transform
hamiltonian = nwqsim_vqe.load("h2.txt")
pauli_terms = nwqsim_vqe.transform(hamiltonian)

# Configure and run VQE
options = nwqsim_vqe.Options()
result = nwqsim_vqe.run(ansatz, pauli_terms, options)
```

---

## Core Functions

### `run(ansatz, hamiltonian, options=Options())`

Run VQE optimization with pre-built ansatz.

**Parameters:**
- `ansatz` : `UCCSD` - Prepared UCCSD ansatz
- `hamiltonian` : `list[PauliTerm]` - Qubit Hamiltonian as Pauli terms  
- `options` : `Options`, optional - VQE configuration options

**Returns:**
- `VQEResult` - Optimization results

**Example:**
```python
result = nwqsim_vqe.run(ansatz, pauli_terms)
print(f"Energy: {result.energy}")
```

### `adapt(path, n_particles, options=Options())`

Run ADAPT-VQE with adaptive circuit construction.

**Parameters:**
- `path` : `str` - Path to Hamiltonian file
- `n_particles` : `int` - Number of particles in the system
- `use_xacc` : `bool`, optional - Use XACC orbital indexing (default: True)
- `options` : `Options`, optional - configuration (mode is forced to `"adapt"`)

**Returns:**
- `AdaptResult` - ADAPT-VQE results with circuit evolution

**Example:**
```python
result = nwqsim_vqe.adapt("h2.txt", n_particles=2)
print(f"Final energy: {result.energy}")
```

### `run_from_file(path, n_particles, use_xacc=True, options=Options())`

Run VQE directly from Hamiltonian file.

**Parameters:**
- `path` : `str` - Path to Hamiltonian file
- `n_particles` : `int` - Number of particles
- `use_xacc` : `bool`, optional - Use XACC orbital indexing (default: True)
- `options` : `Options`, optional - VQE configuration options

**Returns:**
- `VQEResult` - Optimization results

### `energy(ansatz, hamiltonian, parameters, use_gpu=False)`

Evaluate energy for given ansatz and parameters.

**Parameters:**
- `ansatz` : `UCCSD` - Prepared UCCSD ansatz
- `hamiltonian` : `list[PauliTerm]` - Qubit Hamiltonian as Pauli terms
- `parameters` : `list[float]` - Variational parameters
- `use_gpu` : `bool`, optional - Use GPU acceleration (default: False)

**Returns:**
- `float` - Energy value

---

## Core Classes

### `System(n_spatial_orbitals, n_electrons)`

Molecular system specification.

**Parameters:**
- `n_spatial_orbitals` : `int` - Number of spatial orbitals
- `n_electrons` : `int` - Number of electrons

**Properties:**
- `n_spatial` : `int` - Number of spatial orbitals
- `n_electrons` : `int` - Number of electrons  
- `n_qubits` : `int` - Number of qubits (2 * n_spatial)

**Example:**
```python
system = nwqsim_vqe.System(2, 2)
print(system)  # System(orbitals=2, electrons=2, qubits=4)
```

### `UCCSD(system, trotter_steps=1, symmetry_level=3)`

UCCSD ansatz for VQE calculations.

**Parameters:**
- `system` : `System` - Molecular system specification
- `trotter_steps` : `int`, optional - Number of Trotter steps (default: 1)
- `symmetry_level` : `int`, optional - Symmetry level (default: 3)

**Methods:**
- `build()` - Build the UCCSD circuit
- `get_circuit()` - Get the quantum circuit
- `get_parameter_count()` - Get number of variational parameters

**Example:**
```python
ansatz = nwqsim_vqe.UCCSD(system, trotter_steps=1)
ansatz.build()
print(f"Parameters: {ansatz.get_parameter_count()}")
```

### `Options()`

Unified configuration for both VQE and ADAPT-VQE workflows.

**Core fields:**
- `mode` : `str` - Solver mode (`"vqe"` or `"adapt"`)
- `use_gpu` : `bool` - Enable GPU acceleration (when compiled with CUDA)
- `use_xacc_indexing` : `bool` - Select XACC/Qiskit orbital ordering
- `random_seed` : `Optional[int]` - Seed for reproducible randomness
- `optimizer` : `str` - NLopt optimizer name (e.g., `"LN_COBYLA"`)
- `max_evaluations`, `relative_tolerance`, `absolute_tolerance`, `stop_value` - Convergence controls
- `lower_bound`, `upper_bound` - Global parameter bounds
- `trotter_steps`, `symmetry_level` - UCCSD ansatz configuration

**ADAPT-specific fields (used when `mode="adapt"`):**
- `adapt_max_iterations` - Maximum ADAPT iterations
- `adapt_gradient_step` - Finite-difference step for gradients
- `adapt_gradient_tolerance`, `adapt_energy_tolerance` - Convergence thresholds
- `adapt_optimizer` - NLopt optimizer for inner VQE solves
- `adapt_max_evaluations`, `adapt_relative_tolerance`, `adapt_absolute_tolerance`, `adapt_stop_value`, `adapt_max_time`
- `adapt_algorithm_parameters`, `adapt_initial_parameters` - Per-iteration optimizer configuration
- `adapt_log_memory` - Enable memory usage logging each iteration

**Example:**
```python
options = nwqsim_vqe.Options()
options.mode = "adapt"
options.use_gpu = True
options.adapt_optimizer = "L-BFGS-B"
options.adapt_max_iterations = 25
options.adapt_gradient_tolerance = 1e-4
```

---

## Result Classes

### `VQEResult`

VQE optimization results.

**Properties:**
- `energy` : `float` - Ground state energy
- `parameters` : `list[float]` - Optimal variational parameters
- `evaluations` : `int` - Number of function evaluations
- `success` : `bool` - Whether optimization converged
- `message` : `str` - Optimization status message

### `AdaptResult`

ADAPT-VQE results with circuit evolution.

**Properties:**
- `energy` : `float` - Final ground state energy
- `parameters` : `list[float]` - Final variational parameters
- `iterations` : `int` - Number of ADAPT iterations
- `circuit_depth` : `int` - Final circuit depth
- `success` : `bool` - Whether ADAPT converged

---

## Utility Functions

### `load(path)`

Load Hamiltonian from file.

**Parameters:**
- `path` : `str` - Path to Hamiltonian file

**Returns:**
- `HamiltonianData` - Loaded Hamiltonian data

### `transform(hamiltonian)`

Transform fermionic Hamiltonian to qubit Hamiltonian using Jordan-Wigner transformation.

**Parameters:**
- `hamiltonian` : `HamiltonianData` - Fermionic Hamiltonian

**Returns:**
- `list[PauliTerm]` - Qubit Hamiltonian as Pauli terms

### `pauli_to_string(pauli_term)`

Convert Pauli term to string representation.

**Parameters:**
- `pauli_term` : `PauliTerm` - Pauli term

**Returns:**
- `str` - String representation

### `list_optimizers()`

Get list of available optimization algorithms.

**Returns:**
- `list[str]` - Available optimizer names

### `check_gpu_support()`

Check if GPU acceleration is available.

**Returns:**
- `bool` - True if GPU support is available

---

## Convenience Functions

### `system(n_spatial_orbitals, n_electrons)`

Create a molecular system.

**Parameters:**
- `n_spatial_orbitals` : `int` - Number of spatial orbitals
- `n_electrons` : `int` - Number of electrons

**Returns:**
- `System` - Molecular system object

**Raises:**
- `ValueError` - If parameters are invalid

### `load_hamiltonian(filename)`

Load Hamiltonian with error handling.

**Parameters:**
- `filename` : `str` - Path to Hamiltonian file

**Returns:**
- `HamiltonianData` - Loaded Hamiltonian

**Raises:**
- `FileNotFoundError` - If file does not exist
- `ValueError` - If file format is invalid

---

## Data Structures

### `HamiltonianData`

Container for Hamiltonian data.

**Properties:**
- `terms` : `list[FermionTerm]` - Fermionic terms
- `constant` : `float` - Constant energy offset

### `PauliTerm`

Pauli term in qubit Hamiltonian.

**Properties:**
- `coefficient` : `complex` - Complex coefficient
- `pauli_string` : `str` - Pauli string (e.g., "XXYY")

### `FermionTerm`

Fermionic term in molecular Hamiltonian.

**Properties:**
- `coefficient` : `complex` - Complex coefficient
- `operators` : `list[FermionOperator]` - List of fermionic operators

---

## Enums

### `Spin`

Electron spin states.

**Values:**
- `Spin.UP` - Spin up
- `Spin.DOWN` - Spin down

### `OrbitalKind`

Orbital types.

**Values:**
- `OrbitalKind.OCCUPIED` - Occupied orbital
- `OrbitalKind.VIRTUAL` - Virtual orbital

### `FermionOpKind`

Fermionic operator types.

**Values:**
- `FermionOpKind.CREATION` - Creation operator
- `FermionOpKind.ANNIHILATION` - Annihilation operator

---

## Error Handling

The API provides error handling for common issues:

```python
try:
    # Invalid system specification
    system = nwqsim_vqe.system(n_spatial_orbitals=0, n_electrons=2)
except ValueError as e:
    print(f"Input error: {e}")

try:
    # Missing file
    hamiltonian = nwqsim_vqe.load_hamiltonian("nonexistent.txt")
except FileNotFoundError as e:
    print(f"File error: {e}")

try:
    # Invalid parameter
    options = nwqsim_vqe.Options()
    options.relative_tolerance = -1.0
except ValueError as e:
    print(f"Validation error: {e}")
```---

## Examples

### Basic VQE Calculation

```python
import nwqsim_vqe

# 1) Create molecular system
system = nwqsim_vqe.system(n_spatial_orbitals=2, n_electrons=2)

# 2) Load and transform Hamiltonian
hamiltonian = nwqsim_vqe.load("h2.txt")
pauli_terms = nwqsim_vqe.transform(hamiltonian)

# 3) Build UCCSD ansatz  
ansatz = nwqsim_vqe.UCCSD(system, trotter_steps=1, symmetry_level=3)
ansatz.build()

# 4) Configure optimization
options = nwqsim_vqe.Options()
options.optimizer = "L-BFGS-B"
options.max_evaluations = 1000
options.use_gpu = True

# 5) Run VQE
result = nwqsim_vqe.run(ansatz, pauli_terms, options)
print(f"Ground state energy: {result.energy:.8f}")
print(f"Converged: {result.success} ({result.evaluations} evaluations)")
```

### ADAPT-VQE Calculation

```python
import nwqsim_vqe

# Configure ADAPT-VQE
adapt_options = nwqsim_vqe.Options()
adapt_options.mode = "adapt"
adapt_options.adapt_max_iterations = 50
adapt_options.adapt_gradient_tolerance = 1e-3
adapt_options.use_xacc_indexing = True

# Run ADAPT-VQE
result = nwqsim_vqe.adapt(
    "h2.txt", 
    n_electrons=2, 
    options=adapt_options
)

print(f"Final energy: {result.energy:.8f}")
print(f"Circuit depth: {result.circuit_depth}")
print(f"ADAPT iterations: {result.iterations}")
```

### Energy Evaluation

```python
import nwqsim_vqe

# Create evaluator for external optimizers
evaluator = nwqsim_vqe.EnergyEvaluator(ansatz, pauli_terms, use_gpu=True)

def objective_function(parameters):
    return evaluator.evaluate(parameters)

# Use with SciPy
from scipy.optimize import minimize
result = minimize(objective_function, initial_params, method='L-BFGS-B')
```

---

## Version Information

```python
import nwqsim_vqe

print(f"Version: {nwqsim_vqe.__version__}")
print(f"Build info: {nwqsim_vqe.__build_info__}")
```