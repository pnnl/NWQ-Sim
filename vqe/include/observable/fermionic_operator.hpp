#ifndef __FERMI_OP__
#define __FERMI_OP__
#include <vector>
#include <complex>
#include "utils.hpp"
namespace NWQSim {
  /* Basic data type for indices */
  using IdxType = long long int;
  /* Basic data type for value */
  using ValType = double;
  namespace VQE {
    enum FermionOpType {
      Creation,
      Annihilation
    };
    enum Spin {
      Up,
      Down
    };
    enum OrbitalType {
      Occupied,
      Virtual
    };
    class FermionOperator {
      // Basic IR before transformation operations (JW, BK, etc.)
      protected:
      IdxType orbital_index;
      OrbitalType orb_type;
      Spin spin;
      FermionOpType type;
      std::complex<ValType> coeff;
      public:
      FermionOperator(IdxType _idx,
                      OrbitalType _orb_type,
                      Spin _spin,
                      FermionOpType _optype,
                      std::complex<ValType> _coeff = 1.0): orbital_index(_idx),
                                       orb_type(_orb_type),
                                       spin(_spin),
                                       type(_optype),
                                       coeff(_coeff) { }
      FermionOperator(const FermionOperator& other): orbital_index(other.orbital_index),
                                                     orb_type(other.orb_type),
                                                     spin(other.spin),
                                                     type(other.type),
                                                     coeff(other.coeff) {}
      FermionOperator& operator=(const FermionOperator& other) {
        orbital_index = other.orbital_index;
        orb_type = other.orb_type;
        spin = other.spin;
        type = other.type;
        coeff = other.coeff;
        return *this;
      }
      IdxType qubitIndex(IdxType n_occ, IdxType n_virt) const {
        // Flattened indexing scheme
        return getQubitIndex(orbital_index, spin, orb_type, n_occ, n_virt);
      }
      std::complex<ValType> getCoeff() const { return coeff; }
      FermionOpType getType() const { return type; }
      std::string toString (IdxType n_occ, IdxType n_virt) const {
        std::stringstream ss;
        ss << ((type==Creation) ? "^": "") << qubitIndex(n_occ, n_virt);
        return ss.str();
      };
    };
  };// namespace vqe
};// namespace nwqsim

#endif