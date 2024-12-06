# NWQ-VQE
These documents describe the structure and functionality of the NWQ-VQE simulator. 


## Contents
- [Operators](components/operators.md): Utilities for creating and manipulating Hamiltonians, Pauli operators, and Fermionic operators. Covers:
  - vqe/include/hamiltonian.hpp
  - vqe/include/pauli_operator.hpp
  - vqe/include/observable/fermionic_operator.hpp
  - vqe/include/environment.hpp
- [Fermion to Qubit Transformers](components/transformers.md): Transformer function prototype and the Jordan-Wigner mapper
  - vqe/include/transform/transform.cpp
- [Ansatzë](components/circuits.md): Ansatz base and UCCSD classes
  - vqe/include/ansatz.hpp
  - vqe/src/uccsd.cpp
  - vqe/src/uccsdmin.cpp
  - vqe/src/singletgsd.cpp
- [VQE States](components/vqe_states.md): Covers the base VQEState interface as well as the backend-specific implementations
  - vqe/include/vqe_state.hpp
  - vqe/include/svsim_vqe/sv_cpu_vqe.hpp
  - vqe/include/svsim_vqe/sv_cuda_vqe.hpp
  - vqe/include/svsim_vqe/sv_mpi_vqe.hpp
  - vqe/include/svsim_vqe/sv_cuda_mpi_vqe.hpp
- [Measurement](components/measurement.md): Construction of QWC groups via Sorted Insertion, measurement circuit construction, `ObservableList` initialization
  - vqe/include/circuit/measurement.hpp
  - Related functions in vqe/include/utils.hpp, vqe/src/utils.cpp
  - Related functions in vqe_state.hpp and assorted backends
- [Expectation Value Calculation](components/expectation_values.md): Expectation value calculation functions, `ObservableList` usage, 
  - Modifications to include/gates.hpp
  - EXPECT_GATE implementations in include/svsim/* backends
  - Related functions in vqe/include/utils.hpp, vqe/src/utils.cpp
  - Related functions in vqe_state.hpp and assorted backends
- [Function Optimization](components/optimization.md): Interface with NLOpt and associated options
  - Related code in vqe/include/vqe_state.hpp
- [adapt.md](components/adapt.md): Fermionic-ADAPT and Qubit-ADAPT optimizers. Covers:
  - [vqe_adapt.hpp](../vqe/include/vqe_adapt.hpp)
  - related functions in [utils.cpp](../vqe/src/utils.cpp)
- [Utilities](components/utils.md): General utilities defined in [utils.cpp](../vqe/src/utils.cpp) and [utils.hpp](../vqe/include/utils.hpp)
