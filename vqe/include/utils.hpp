#ifndef VQE_UTILS
#define VQE_UTILS
#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <math.h>
#include <sstream>


// Templated print function for std::vector
template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& target) {
  out << "[";
  size_t len = target.size();
  if (len > 0) {
      for (size_t i = 0; i < len - 1; i++) {
          out << target[i] << ", ";
      }
      out << target[len-1];
  }
  out << "]";
  return out;
}

// Status enum for MPI processes
enum STATUS {
    CALL_SIMULATOR,
    WAIT,
    EXIT_LOOP
};
template <typename T, typename S>
std::stringstream& operator<<(std::stringstream& out, const std::pair<T, S>& target) {
  out << "(" << target.first << ", " << target.second << ")";
  return out;
}
template <typename T, typename S>
std::ostream& operator<<(std::ostream& out, const std::pair<T, S>& target) {
  out << "(" << target.first << ", " << target.second << ")";
  return out;
}
template <typename T>
std::stringstream& operator<<(std::stringstream& out, const std::vector<T>& target) {
  out << "[";
  size_t len = target.size();
  if (len > 0) {
      for (size_t i = 0; i < len - 1; i++) {
          out << target[i] << ", ";
      }
      out << target[len-1];
  }
  out << "]";
  return out;
}

namespace NWQSim{
  namespace VQE{
  using IdxType = long long;
  using ValType = double;
  struct OptimizerSettings {
    // Structure to pass settings to the optimizer
    ValType rel_tol; // relative tolerance cutoff
    ValType abs_tol; // absolute tolerance cutoff
    ValType stop_val; //
    IdxType max_evals; // Max number of function evaluations
    ValType max_time; // Optimizer timeout (seconds)
    std::unordered_map<std::string, ValType> parameter_map; // map for setting optimizer-specific parameters
    // Defaults (disables all of the settings, except for the max_eval ceiling)
    OptimizerSettings(): rel_tol(-1), 
                         abs_tol(-1),
                         stop_val(-MAXFLOAT),
                         max_evals(50),
                         max_time(-1) {}
  };

  inline
  IdxType getQubitIndex(IdxType orbital_index, IdxType spin, IdxType orb_type, IdxType n_occ, IdxType n_virt) {
    // Flattened indexing scheme
    // IdxType index = (orbital_index) \
          + (orb_type * spin * n_virt + (!orb_type) * spin * n_occ) \
          + orb_type * 2 * n_occ;
    // Qiskit indexing scheme
    IdxType index = (orbital_index) \
            + (orb_type * n_occ) \
            + spin * (n_occ + n_virt);
    return index;
  }


  inline
  void getFermiInfoFromQubit(IdxType qubit_idx, IdxType& orbital_index, IdxType& spin, IdxType& orb_type, IdxType n_occ, IdxType n_virt) {
    // Flattened indexing scheme (reversed). Extracts orbital/operator/spin properties from the qubit index
    // orbital_index = qubit_idx;
    // orb_type = (qubit_idx >= (2 * n_occ));
    // orbital_index -= orb_type * (2 * n_occ);
    // if (orb_type) {
    //   spin = orbital_index >= n_virt;
    //   orbital_index -= spin * n_virt;
    // } else {
    //   spin = orbital_index >= n_occ;
    //   orbital_index -= spin * n_occ;
    // }
    orbital_index = qubit_idx;
    spin = (orbital_index >= (n_virt + n_occ));
    orbital_index -= spin * (n_virt + n_occ);
    orb_type = orbital_index >= n_occ;
    orbital_index -= orb_type * n_occ;
    // Below is the reverse indexing for comparison with Qiskit
    // std::cout << qubit_idx << " " << orbital_index << " " << spin << " " << orb_type << " " << n_occ << " " << n_virt << std::endl;
  }

  // Convert an integer to an  `n_qubits`-digit binary string
  std::string to_binary_string(IdxType val, IdxType n_qubits);

};};

#endif