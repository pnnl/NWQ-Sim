# VQE Functionality Reference

This document outlines the functionalities supported by the current VQE implementation. It is intended to be used as a reference for code refactoring to ensure all features are correctly ported.

## 1. VQE Modes

The VQE implementation supports three main execution modes:

1.  **Standard VQE:** A conventional VQE loop using a fixed, predefined ansatz.
    *   **Executable:** `main` (compiled from `main.cpp`)
    *   **Activation:** This is the default mode when the `--adapt` flag is NOT used.

2.  **ADAPT-VQE:** An adaptive VQE algorithm that dynamically grows the ansatz from an operator pool.
    *   **Executable:** `main` (compiled from `main.cpp`)
    *   **Activation:** Enabled with the `--adapt` flag.

3.  **QFlow:** A specific workflow that uses a "local gradient follower" to optimize the ansatz parameters. This appears to be a custom gradient descent implementation.
    *   **Executable:** `qflow` (compiled from `qflow.cpp`)
    *   **Activation:** Enabled with flags like `--local`, `--localspsa`, or `--localfd`.

## 2. Supported Ansatze

The following parameterized circuits (ansatze) are available:

1.  **UCCSD (Unitary Coupled-Cluster Singles and Doubles):**
    *   `UCCSDmin`: A memory-efficient implementation of UCCSD. This is the default for standard VQE.
    *   `UCCSD`: An older implementation, available via the `--origin` flag.
    *   `Singlet_GSD`: A variation using a singlet generalized singles and doubles ansatz, available via the `--gsd` flag.

2.  **DynamicAnsatz:**
    *   Used for ADAPT-VQE. The structure of this ansatz is built dynamically during the optimization.

## 3. Operator Pools for ADAPT-VQE

When using the ADAPT-VQE mode, the following operator pools can be used to generate the ansatz:

1.  **Fermionic Pool:** The default pool, consisting of fermionic operators.
2.  **Pauli Pool (Qubit Pool):** A pool of Pauli operators (qubit operators), enabled with the `--qubit` flag.

## 4. Classical Optimizers

The VQE code uses two main types of classical optimizers:

1.  **NLopt Library:** A library of non-linear optimization algorithms.
    *   **Usage:** This is the default for standard VQE and ADAPT-VQE.
    *   **Configuration:** The specific algorithm can be chosen with the `-o` or `--optimizer` flag (e.g., `LN_COBYLA`, `LN_BOBYQA`, `LD_LBFGS`). Other parameters like tolerance and max evaluations are also configurable.

2.  **Local Gradient Follower:** A custom implementation of gradient descent.
    *   **Usage:** Used in the "QFlow" mode.
    *   **Gradient Estimation:**
        *   **Finite Differences (FD):** Enabled with `--localfd`.
        *   **SPSA (Simultaneous Perturbation Stochastic Approximation):** Enabled with `--localspsa`.

## 5. Hamiltonian Format

*   The Hamiltonian is specified via an input file (`-f` or `--hamiltonian`).
*   The file format is a sum of fermionic operators, with each line representing a term.

## 6. Key Configuration Flags

*   `-f, --hamiltonian`: **Required.** Path to the Hamiltonian file.
*   `-p, --nparticles`: **Required.** Number of electrons.
*   `-b, --backend`: Specifies the simulation backend (e.g., `CPU`, `MPI`).
*   `--adapt`: Enables ADAPT-VQE mode.
*   `--local`, `--localspsa`, `--localfd`: Enables the QFlow mode with a specific gradient estimation method.
*   `-o, --optimizer`: Selects the NLopt optimization algorithm.
*   `--symm`: Sets the symmetry level for the UCCSD ansatz.
*   `--qubit`: Use a qubit operator pool for ADAPT-VQE.
*   `--origin`: Use the older UCCSD implementation.
*   `--gsd`: Use the Singlet GSD ansatz.
