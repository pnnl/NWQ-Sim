#ifndef VQE_UTILS
#define VQE_UTILS
#include "nwq_util.hpp"
#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <math.h>
#include <sstream>
#include <list>
// #include "nwq_util.hpp"





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
// Templated print function for std::pair
template <typename T, typename S>
std::ostream& operator<<(std::ostream& out, const std::pair<T, S>& target) {
  std::cout << "(" << target.first << ", " << target.second << ")" << std::endl;
  return out;
}
template <typename T>
inline
T choose2(T n) {
    if (n >= 2) {
        return n * (n-1) / 2;
    }
    return 0;
}
// Status enum for MPI processes
enum STATUS {
    CALL_SIMULATOR,
    WAIT,
    EXIT_LOOP
};

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

  class FermionOperator;

  class PauliOperator;
  struct MolecularEnvironment;
  using IdxType = long long;
  using ValType = double;
  enum class Commute {
    GC,  // general commutativity (aka "FC")
    QWC, // qubit-wise commutativity
    TRC  // topology-restricted commutativity
  };


  struct OptimizerSettings {
   /**
    * @brief  Data structure for NLOpt optimizer settings
    */
    ValType rel_tol; // relative tolerance cutoff
    ValType abs_tol; // absolute tolerance cutoff
    ValType stop_val; //
    IdxType max_evals; // Max number of function evaluations
    ValType max_time; // Optimizer timeout (seconds)
    ValType lbound; // Lower bound
    ValType ubound; // Upper bound
    std::unordered_map<std::string, ValType> parameter_map; // map for setting optimizer-specific parameters
    // Defaults (disables all of the settings, except for the max_eval ceiling)
    OptimizerSettings(): rel_tol(-1), 
                         abs_tol(-1),
                         stop_val(-MAXFLOAT),
                         max_evals(100), // MZ: changed from 50 to 100, to match the default in main.cpp.
                         max_time(-1),
                         lbound(-PI+0.00001), // MZ: avoid bad numerical issues
                         ubound(PI-0.00001) {}
  };


 /**
  * @brief  Get the qubit index of a spin orbital
  * @note   
  * @param  orbital_index: The numerical index of an orbital (starting from 0) within either occupied or virtual
  * @param  spin: Spin type (Spin::Up or Spin::Down)
  * @param  orb_type: Orbital type (virtual or occupied)
  * @param  n_occ: Number of occupied orbitals in system
  * @param  n_virt: Number of virtual orbitals in system
  * @param  xacc_scheme: Flag for XACC ordering (true) or canonical (false)
  */
  inline
  IdxType getQubitIndex(IdxType orbital_index, IdxType spin, IdxType orb_type, IdxType n_occ, IdxType n_virt, bool xacc_scheme) {
    // Flattened indexing scheme
    IdxType index;
    if (!xacc_scheme) {
      // DUCC scheme
      index = (orbital_index) \
          + (orb_type * spin * n_virt + (!orb_type) * spin * n_occ) \
          + orb_type * 2 * n_occ;  // MZ: not for PySCF, temporarily commented out
      // index = 2*(orbital_index+ (orb_type * n_occ)) + spin; // MZ: DEBUG: for PySCF and original ADAPT-VQE code, TODO: keep it or recover it
    } else {
      // Qiskit/XACC scheme
      index = (orbital_index) \
            + (orb_type * n_occ) \
            + spin * (n_occ + n_virt);
    }
    return index;
  }


/**
 * @brief  Count the ones in a long long bitmask
 * @note   Used for Pauli parity operations
 * @param  val: Bitmask to count
 * @retval 
 */
inline
IdxType count_ones(IdxType val) {
  IdxType count = 0;
  IdxType mask = 1;
  for (size_t i = 0; i < sizeof(IdxType) * 8 - 1; i++) {
    count += (val & mask) > 0;
    mask = mask << 1;
  }
  return count;
}
 /**
  * @brief  Extract orbital information from a qubit index
  * @note   Inverse of `getQubitIndex`
  * @param  qubit_idx: Qubit index (input)
  * @param  orbital_index: Numerical index of orbital (output)
  * @param  spin: Spin type (Spin::Up or Spin::Down) (output)
  * @param  orb_type: Type of orbital (Occupied or Virtual) (output)
  * @param  n_occ: # Occupied orbitals
  * @param  n_virt: # Virtual orbitals
  * @param  xacc_scheme: Flag on whether to use XACC ordering (true) or canonical ordering (false)
  * @retval None
  */
  inline
  void getFermiInfoFromQubit(IdxType qubit_idx, IdxType& orbital_index, IdxType& spin, IdxType& orb_type, IdxType n_occ, IdxType n_virt, bool xacc_scheme) {
    // Flattened indexing scheme (reversed). Extracts orbital/operator/spin properties from the qubit index
    if (!xacc_scheme) {
      // DUCC reverse indexing
      orbital_index = qubit_idx;
      orb_type = (qubit_idx >= (2 * n_occ));
      orbital_index -= orb_type * (2 * n_occ);
      if (orb_type) {
        spin = orbital_index >= n_virt;
        orbital_index -= spin * n_virt;
      } else {
        spin = orbital_index >= n_occ;
        orbital_index -= spin * n_occ;
      }
    } else {
      // Qiskit/XACC reverse indexing
      orbital_index = qubit_idx;
      spin = (orbital_index >= (n_virt + n_occ));
      orbital_index -= spin * (n_virt + n_occ);
      orb_type = orbital_index >= n_occ;
      orbital_index -= orb_type * n_occ;
    }
  }


/**
  * @brief MZ: Counting the number of double excitations for singlet GSD
  * @note MZ: yeah, I know it's stupid
  * @param n_spatial 
  * @return int 
  */
  inline
  int counting_doubles(int n_spatial) {
    int doublt_term1_counter = 0;
    int doublt_term2_counter = 0;
    int doublt_term4_counter = 0;
    int doublt_term6_counter = 0;

    // Singlets and Triplets
    int rs = -1;
    for (IdxType r = 0; r < n_spatial; r++) {
      for (IdxType s = r; s < n_spatial; s++) {
        int pq = -1;
        for (IdxType p = 0; p < n_spatial; p++) {
          for (IdxType q = p; q < n_spatial; q++) {
            pq += 1;
            if (rs > pq) continue;
            if ( (p == r) && (q == s) ) continue;
            if (p!=q) {
              if (r == s) { // only singlet with two terms, not sure Group 3 or Group 4 cases
                doublt_term2_counter += 1;
                continue;
              }
              if ( ((r!=s)&&(q!=r)) || ((p==r)&&(q!=s)&(r!=s)) || ((p==s)&&(q!=r)&(r!=s)) || ((q==s)&&(p!=r)) || ((q==r)&&(p!=s))) {
                doublt_term4_counter += 1;
                doublt_term6_counter += 1;
              }
            }
            // Group 3 to 5
            if (p == q) {
              // Group 3
              if ((q != r)&&(r!=s)) {
                // only singlet
                doublt_term2_counter += 1;
                continue;
              }
              // Group 4
              if ((q == r)&&(r!=s)) {
                // only singlet
                doublt_term2_counter += 1;
                continue;
              }
              // Group 5
              if ((q != r)&&(r==s)) {
                // only singlet
                doublt_term1_counter += 1;
                continue;
              }
            } // p == q
          } // q
        } // p
      } // s
    } // r
    return doublt_term1_counter+2*doublt_term2_counter+4*doublt_term4_counter+6*doublt_term6_counter;
  }





  // Forward declaration for PauliOperator class
  class PauliOperator;

  // Make a common operator for a QWC class
  PauliOperator make_common_op(const std::vector<PauliOperator>& pauli_list, 
                               std::vector<IdxType>& zmasks,
                               std::vector<ValType>& coeffs);
  // Heuristic to partition operators into QWC-compatible cliques
  void sorted_insertion(const std::vector<PauliOperator>& paulilist, std::list<std::vector<IdxType> >& cliques, bool overlap);
  // Convert an integer to an  `n_qubits`-digit binary string
  std::string to_binary_string(IdxType val, IdxType n_qubits);
  std::string to_fermionic_string(const std::vector<FermionOperator>& product, const MolecularEnvironment& env);
  void read_amplitudes(std::string fpath, std::vector<ValType>& params, const std::unordered_map<std::string, IdxType>& idx_map);

};};

#endif