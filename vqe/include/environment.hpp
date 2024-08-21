#ifndef MOLECULE
#define MOLECULE
#include "utils.hpp"


namespace NWQSim {
  using IdxType = long long;
  using ValType = double;
  namespace VQE {

    struct MolecularEnvironment {
      /**
       * Data structure to store environment information (orbital/particle config, energy offsets)
      */
      IdxType n_spatial; // # spatial orbitals
      IdxType n_part;    // # electrons
      IdxType n_occ;     // # occupied orbitals
      IdxType n_virt;    // # virtual orbitals
      bool xacc_scheme;  // Whether to use XACC (true) or canonical (false) qubit indexing
      ValType constant;  // constant offset
      
      // Ctors
      MolecularEnvironment() {};

      MolecularEnvironment(IdxType _nspatial,
                           IdxType _npart,
                           bool _xacc_scheme,
                           ValType _const = 0.0): n_spatial(_nspatial),
                                            n_part (_npart),
                                            n_occ (n_part / 2),
                                            n_virt (n_spatial - n_occ),
                                            xacc_scheme(_xacc_scheme),
                                            constant(_const) {}

      IdxType addr(IdxType orbital, bool is_virtual, bool spin_down) const {
        // Return the qubit index for a given orbital index, virtual/occupied status, and spin value
        return getQubitIndex(orbital, spin_down, is_virtual, n_occ, n_virt, xacc_scheme);
      }
    };

  };// namespace vqe
};// namespace nwqsim
#endif