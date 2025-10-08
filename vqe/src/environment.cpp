#include "core/environment.hpp"

namespace vqe {

std::size_t molecular_environment::qubit_index(std::size_t orbital_index,
                                               orbital_kind kind,
                                               spin spin_state) const {
  const std::size_t n_occ = occupied_orbitals();
  const std::size_t n_virt = virtual_orbitals();
  if (!xacc_indexing) {
    const std::size_t spin_val = (spin_state == spin::down) ? 1 : 0;
    const std::size_t orb_val = (kind == orbital_kind::virtual_orbital) ? 1 : 0;
    return orbital_index + (orb_val * spin_val * n_virt + (1 - orb_val) * spin_val * n_occ) + orb_val * 2 * n_occ;
  }

  const std::size_t spin_val = (spin_state == spin::down) ? 1 : 0;
  const std::size_t orb_val = (kind == orbital_kind::virtual_orbital) ? 1 : 0;
  return orbital_index + orb_val * n_occ + spin_val * (n_occ + n_virt);
}

}  // namespace vqe
