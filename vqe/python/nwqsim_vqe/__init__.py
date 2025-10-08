"""
NWQ-Sim VQE Module
==================

A high-performance Variational Quantum Eigensolver (VQE) implementation with support 
for UCCSD ansätze, ADAPT-VQE, and GPU acceleration.

Key Features:
- UCCSD and ADAPT-VQE algorithms
- GPU acceleration via CUDA
- Integration with SciPy optimizers
- Efficient Jordan-Wigner transformation
- Thread-safe energy evaluation

Basic Usage:
    >>> import nwqsim_vqe as vqe
    >>> # Load Hamiltonian and run VQE
    >>> result = vqe.run_vqe("h2.txt", n_electrons=2)
    >>> print(f"Ground state energy: {result.energy}")
"""

from __future__ import annotations

from typing import Optional, Sequence, Union, List, Dict, Any
import warnings

from . import _core

# Re-export optimized C++ bindings
from ._core import (  # noqa: F401
    # Configuration classes
    Options,
    
    # Result classes  
    AdaptResult,
    VQEResult,
    
    # Core data structures
    BasicFermionOp,
    FermionOperator, 
    FermionTerm,
    HamiltonianData,
    System,
    PauliTerm,
    UCCSD,
    
    # Enums
    FermionOpKind,
    OperatorKind,
    OrbitalKind,
    Spin,
    
    # Energy evaluation
    EnergyEvaluator,
    
    # Core functions
    energy,
    transform,
    pauli_to_string,
    load,
    adapt,
    run_from_file,
    run,
    
    # Utilities
    check_gpu_support,
    list_optimizers,
)

# Backwards compatibility aliases (with deprecation warnings)
def FermionOperatorKind(*args, **kwargs):
    warnings.warn("FermionOperatorKind is deprecated, use FermionOpKind instead", 
                  DeprecationWarning, stacklevel=2)
    return FermionOpKind(*args, **kwargs)

__all__ = [
    # Main API
    "run",
    "adapt", 
    "run_from_file",
    
    # Configuration
    "Options",
    
    # Results
    "VQEResult", 
    "AdaptResult",
    
    # Core classes
    "UCCSD",
    "System",
    "EnergyEvaluator",
    
    # Data structures
    "HamiltonianData",
    "PauliTerm",
    "FermionTerm",
    "FermionOperator",
    "BasicFermionOp",
    
    # Enums
    "Spin",
    "OrbitalKind", 
    "OperatorKind",
    "FermionOpKind",
    
    # Utilities  
    "load",
    "transform", 
    "energy",
    "pauli_to_string",
    "check_gpu_support",
    "list_optimizers",
    
    # Convenience functions
    "system",
    "load_hamiltonian",
    
    # Backwards compatibility
    "FermionOperatorKind",
]





def system(
    n_spatial_orbitals: int,
    n_electrons: int, 
    *,
    xacc_indexing: bool = True,
    nuclear_repulsion: float = 0.0
) -> System:
    """
    Create a molecular system specification.
    
    Create a molecular system with
    validation and sensible defaults.
    
    Parameters
    ----------
    n_spatial_orbitals : int
        Number of spatial molecular orbitals.
    n_electrons : int
        Total number of electrons.
    xacc_indexing : bool, default=True
        Use XACC orbital ordering convention.
    nuclear_repulsion : float, default=0.0
        Nuclear repulsion energy offset.
        
    Returns
    -------
    System
        Configured molecular environment.
        
    Raises
    ------
    ValueError
        If orbital or electron counts are invalid.
        
    Examples
    --------
    >>> # H2 molecule
    >>> env = system(n_spatial_orbitals=2, n_electrons=2)
    >>> print(f"System requires {env.total_qubits()} qubits")
    """
    
    if n_spatial_orbitals <= 0:
        raise ValueError("Number of spatial orbitals must be positive")
    if n_electrons <= 0:
        raise ValueError("Number of electrons must be positive")  
    if n_electrons > 2 * n_spatial_orbitals:
        raise ValueError("Too many electrons for given number of orbitals")
    
    return System(
        n_spatial_orbitals, 
        n_electrons, 
        xacc_indexing, 
        nuclear_repulsion
    )


def load_hamiltonian(filename: str) -> HamiltonianData:
    """
    Load fermionic Hamiltonian from file with error handling.
    
    Wrapper around read_hamiltonian_file with improved error messages.
    
    Parameters
    ----------
    filename : str
        Path to Hamiltonian file.
        
    Returns
    -------
    HamiltonianData
        Loaded Hamiltonian data.
        
    Raises  
    ------
    FileNotFoundError
        If file does not exist.
    ValueError  
        If file format is invalid.
    """
    
    import os
    
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"Hamiltonian file not found: {filename}")
    
    try:
        return load(filename)
    except Exception as e:
        raise ValueError(f"Failed to parse Hamiltonian file '{filename}': {e}") from e


# Module metadata
__version__ = getattr(_core, "__version__", "unknown")
__build_info__ = getattr(_core, "__build_info__", {})
