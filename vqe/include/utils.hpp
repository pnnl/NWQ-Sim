#ifndef VQE_UTILS
#define VQE_UTILS
#include <iostream>
#include <vector>
namespace NWQSim{
  using IdxType = long long;
  using ValType = double;
  namespace VQE{
  inline
  IdxType getQubitIndex(IdxType orbital_index, IdxType spin, IdxType orb_type, IdxType n_occ, IdxType n_virt) {
    // Flattened indexing scheme
    IdxType index = (orbital_index) \
            + (orb_type * spin * n_virt + (!orb_type) * spin * n_occ) \
            + orb_type * 2 * n_occ;
    // Qiskit indexing scheme
    // IdxType index = (orbital_index) \
            + (orb_type * n_occ) \
            + spin * (n_occ + n_virt);
    return index;

  }
  inline
  void getFermiInfoFromQubit(IdxType qubit_idx, IdxType& orbital_index, IdxType& spin, IdxType& orb_type, IdxType n_occ, IdxType n_virt) {
    // Flattened indexing scheme (reversed)
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
    // orbital_index = qubit_idx;
    // spin = (orbital_index >= (n_virt + n_occ));
    // orbital_index -= spin * (n_virt + n_occ);
    // orb_type = orbital_index >= n_occ;
    // orbital_index -= orb_type * n_occ;
    // std::cout << qubit_idx << " " << orbital_index << " " << spin << " " << orb_type << " " << n_occ << " " << n_virt << std::endl;
  }
  std::string to_binary_string(NWQSim::IdxType val, NWQSim::IdxType n_qubits);
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
};};

#endif