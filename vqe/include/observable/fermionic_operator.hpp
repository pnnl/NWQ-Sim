#ifndef __FERMI_OP__
#define __FERMI_OP__
#include <vector>
#include <complex>
#include "utils.hpp"
#include "environment.hpp"
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
      bool xacc_scheme;
      std::complex<ValType> coeff;
      public:
      FermionOperator(IdxType _idx,
                      OrbitalType _orb_type,
                      Spin _spin,
                      FermionOpType _optype,
                      bool _xacc_scheme,
                      std::complex<ValType> _coeff = 1.0): orbital_index(_idx),
                                       orb_type(_orb_type),
                                       spin(_spin),
                                       type(_optype),
                                       xacc_scheme(_xacc_scheme),
                                       coeff(_coeff) { }
      FermionOperator(const FermionOperator& other): orbital_index(other.orbital_index),
                                                     orb_type(other.orb_type),
                                                     spin(other.spin),
                                                     type(other.type),
                                                     xacc_scheme(other.xacc_scheme),
                                                     coeff(other.coeff) {}
      FermionOperator& operator=(const FermionOperator& other) {
        orbital_index = other.orbital_index;
        orb_type = other.orb_type;
        spin = other.spin;
        type = other.type;
        coeff = other.coeff;
        xacc_scheme = other.xacc_scheme;
        return *this;
      }
      bool operator==(const FermionOperator& other)const  {
        return (orbital_index == other.orbital_index && spin == other.spin && orb_type == other.orb_type && type == other.type);
      }
      bool operator!=(const FermionOperator& other)const  {
        return !(orbital_index != other.orbital_index || spin != other.spin || orb_type != other.orb_type || type != other.type);
      }
      IdxType qubitIndex(IdxType n_occ, IdxType n_virt) const {
        // Flattened indexing scheme
        return getQubitIndex(orbital_index, spin, orb_type, n_occ, n_virt, xacc_scheme);
      }
      std::complex<ValType> getCoeff() const { return coeff; }
      FermionOpType getType() const { return type; }
      Spin getSpin() const { return spin; }
      OrbitalType getOrbitalType() const { return orb_type; }
      std::string toString (IdxType n_occ, IdxType n_virt) const {
        std::stringstream ss;
        ss  << qubitIndex(n_occ, n_virt) << ((type==Creation) ? "^": "");
        return ss.str();
      };
    };
    void read_fermion_operators(std::string input_path,
                                std::vector<std::vector<FermionOperator>>& fermi_operators,
                                MolecularEnvironment& env,
                                size_t n_particles,
                                bool xacc_scheme);
  };// namespace vqe
};// namespace nwqsim

#endif