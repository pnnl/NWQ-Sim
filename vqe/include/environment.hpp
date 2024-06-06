#ifndef MOLECULE
#define MOLECULE
#include <vector>
#include <string>
#include "utils.hpp"


namespace NWQSim {
  using IdxType = long long;
  using ValType = double;
  namespace VQE {
    struct MolecularEnvironment {
      IdxType n_spatial;
      IdxType n_part; 
      IdxType n_occ;
      IdxType n_virt;
      ValType constant;
      MolecularEnvironment() {};
      MolecularEnvironment(IdxType _nspatial,
                           IdxType _npart,
                           ValType _const = 0.0): n_spatial(_nspatial),
                                            n_part (_npart),
                                            n_occ (n_part / 2),
                                            n_virt (n_spatial - n_occ),
                                            constant(_const) {}
      IdxType addr(IdxType orbital, bool is_virtual, bool spin_down) const {
        return getQubitIndex(orbital, spin_down, is_virtual, n_occ, n_virt);
      }
    };
  };// namespace vqe
};// namespace nwqsim
#endif