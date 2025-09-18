#ifndef PAULI
#define PAULI
#include <vector>
#include <string>
#include <cassert>
#include <sstream>
#include <unordered_map>
#include <iostream>
#include <complex>
#include <algorithm>
#include <memory>
#include <numeric>
#include <list>
#include "utils.hpp"
namespace NWQSim {
  /* Basic data type for indices */
  using IdxType = long long int;
  /* Basic data type for value */
  using ValType = double;
  namespace VQE {
    const std::complex<ValType> Imag = {0.0, 1.0};
    
    // also determines phase, all anti-commuting operators add a pi/2 phase
    const bool commutationRelations[4][4] = {
      {1, 1, 1, 1},
      {1, 1, 0, 0},
      {1, 0, 1, 0},
      {1, 0, 0, 1}
    };
    // Sign relations for pauli products
    const int signRelations[4][4] = {
      {1,  1,   1,   1},   // I commutes with everything
      {1,  1,   -1,  1}, // XI=X, XX=I, XZ=-iY, XZ=iZ
      {1,  1,  1,   -1}, // ZI=Z, ZX=iY, ZZ=I, ZY=-iX
      {1, -1 ,  1,   1} // YI=Y, YX=-iZ, YY=I, YZ=iX
    };
 
    const char *const PAULI_OP_NAMES[] = {
      "I",
      "X",
      "Z",
      "Y"
    };

    class PauliOperator {
      private:
        IdxType dim;
        bool non_trivial;
        std::complex<ValType> coeff;
        IdxType xmask;
        IdxType zmask;
        std::vector<IdxType> x_indices;
        // std::shared_ptr<std::vector<PauliOp> > ops; // Store in big-endian order v, i.e. qubit 0 is in position 0
      public:
        PauliOperator() {};
        PauliOperator(IdxType _xmask,
                      IdxType _zmask,
                      IdxType _dim,
                      std::complex<ValType> _coeff = 1.0): dim(_dim), coeff(_coeff), xmask(_xmask), zmask(_zmask) {
          non_trivial = (_xmask | _zmask) != 0;
        } 
        
        PauliOperator(std::string _opstring,
                      std::complex<ValType> _coeff = 1.0): dim(_opstring.length()), coeff(_coeff) {
          // ops = std::make_shared<std::vector<PauliOp> >(dim);
          size_t index = 0;
          non_trivial = false;
          xmask = 0;
          zmask = 0;
          for (char i: _opstring) {
            assert(index < dim);
            switch (i)
            {
            case 'X':
              xmask |= (1 << (dim - (index) - 1));
              break;
            case 'Y':
              xmask |= (1 << (dim - (index) - 1));
              zmask |= (1 << (dim - (index) - 1));
              break;
            case 'Z':
              zmask |= (1 << (dim - (index) - 1));
              break;
            case 'I':
              break;
            
            default:
              assert(false);
              break;
            }
            index++;
          }
          non_trivial = (xmask | zmask) != 0;
        } 

        PauliOperator& operator=(const PauliOperator& other) {
          dim = other.dim;
          coeff = other.coeff;
          non_trivial = other.non_trivial;
          xmask = other.xmask;
          zmask = other.zmask;
          x_indices = other.x_indices;

          // assert(coeff. > 0.0);
          return *this;
        }
        PauliOperator(const PauliOperator& other) {
          dim = other.dim;
          // ops = std::make_shared<std::vector<PauliOp> >(
          //   other.ops.get()->begin(), other.ops.get()->end()); 
          coeff = other.coeff;
          non_trivial = other.non_trivial;
          xmask = other.xmask;
          zmask = other.zmask;
          x_indices = other.x_indices;
          // assert(coeff.real() > 0.0 || coeff.imag() > 0.0);
        }
        PauliOperator conj() const {
          return PauliOperator(xmask, zmask, dim, std::complex<ValType>(coeff.real(), -coeff.imag()));
        }
        PauliOperator operator*(const PauliOperator& other) const {
          std::complex<ValType> new_coeff = coeff * other.coeff;
          IdxType dim1 = dim;
          IdxType dim2 = other.dim;
          if (dim1 != dim2) {
            fprintf(stderr, "Pauli strings of different size: %s %s\n", pauliToString().c_str(), other.pauliToString().c_str());
            assert(false);
          }
          IdxType newxmask = xmask ^ other.xmask;
          IdxType newzmask = zmask ^ other.zmask;
          IdxType acmask = (xmask & other.zmask) ^ (other.xmask & zmask);
          for (IdxType t = 0; t < dim; t++) {
            if (acmask & (1 << t)) {
              int xbit1 = (xmask & (1 << t)) >> t;
              int xbit2 = (other.xmask & (1 << t))  >> t;
              int zbit1 = (zmask & (1 << t))  >> (t);
              int zbit2 = (other.zmask & (1 << t))  >> (t);
              new_coeff *= std::complex<ValType>(0, signRelations[xbit1 + 2 * zbit1][xbit2 + 2 * zbit2]);
            }
          }
          return PauliOperator(newxmask, newzmask, dim, new_coeff);
        }
        PauliOperator& operator*=(std::complex<ValType> scalar){
          coeff *= scalar;
          return *this;
        }
        PauliOperator operator*(std::complex<ValType> scalar) const {
          PauliOperator newop (*this);
          newop.coeff *= scalar;
          return newop;
        }
        void get_xindices(std::vector<IdxType>& indices) const {
          for (int i = 0; i < dim; i++) {
            if (xmask & (1 << i)) {
              indices.push_back(i);
            }
          }
        }
        IdxType get_dim() const {
          return dim;
        }
        IdxType get_xmask() const {
          return xmask;
        }
        IdxType get_zmask() const {
          return zmask;
        }
        // Dump the Pauli operator to string
        std::string pauliToString(bool print_coeff = true) const {
          std::stringstream ss;
          if (coeff != 1.0 || !print_coeff) {
            ss << "(" << coeff.real() << " ";
            if (coeff.imag() >= 0) {
              ss << "+ " << coeff.imag() << "i)";
            } else {
              ss << "- " << fabs(coeff.imag()) << "i)";
            }
            // ss
          }
          // Reverse to ensure correct order (little endian)
          for (int i = dim-1; i >= 0; i--) {
            int xbit = (xmask & (1 << i)) >> i;
            int zbit = (zmask & (1 << i)) >> i;
            ss << PAULI_OP_NAMES[xbit + 2 * zbit];
          }
          return ss.str();
        }
        bool parity(const PauliOperator& other, ValType& sign) const {
          IdxType mindim = std::min(other.dim, dim);
          IdxType n_ones = count_ones((xmask & other.zmask) | (zmask & other.xmask));
          sign = (n_ones / 2) % 2 ? -1.0 : 1.0;
          return (n_ones % 2) == 0;
          // return parity_xor;
        }
        bool parity(IdxType other_zmask, ValType& sign) const {
          IdxType n_ones = count_ones((xmask & other_zmask) | zmask );
          sign = (n_ones / 2) % 2 ? -1.0 : 1.0;
          return (n_ones % 2) == 0;
          // return parity_xor;
        }
        bool parity(const PauliOperator& other) const {
          IdxType mindim = std::min(other.dim, dim);
          IdxType n_ones = count_ones((xmask & other.zmask) | (zmask & other.xmask));
          return (n_ones % 2) == 0;
          // return parity_xor;
        }
        bool isNonTrivial() const { return non_trivial;}
        bool QWC(const PauliOperator& other) const {          
          // return ((xmask & other.zmask) | (zmask & other.xmask)) == 0;
          return ((xmask & other.zmask) ^ (zmask & other.xmask)) == 0; // MZ: noticed by Sean Garner
        }
        bool GC(const PauliOperator& other) const {
          return parity(other);
        }
        IdxType count_y() const {
          return count_ones(xmask & zmask);
        }
        bool TRC(const PauliOperator& other, 
                 IdxType group_mask, 
                 const std::vector<std::vector<IdxType> >& distances,
                 IdxType tolerance) const {
          IdxType sympprod = (xmask & other.zmask) | (zmask & other.xmask);
          IdxType total_anticomm = sympprod | group_mask;
          IdxType n_anticomm_total = 0;
          IdxType n_anticomm = 0;
          std::list<IdxType> anticomm;
          IdxType D = 0;
          for (IdxType i = 0; i < dim; i++) {
            if (total_anticomm & (1 << i)) {
              n_anticomm_total++;
              for (auto other: anticomm) {
                D = std::max(D, distances[i][other]);
              }
              anticomm.push_back(i);

              if (sympprod & (1 << i))
                n_anticomm++;
            }
          }
          bool gc = (n_anticomm % 2) == 0;
          bool within_tol = (D * n_anticomm_total * n_anticomm_total) < tolerance;
          return within_tol && gc;
        }
        bool commutes(const PauliOperator& other, 
                      Commute relation) const {
          switch(relation) {
            case Commute::QWC:
              return QWC(other);
            case Commute::GC:
              return GC(other);
            default:
              return false;
          }
        }
        std::complex<ValType> getCoeff() const {return coeff;}
        void setCoeff (std::complex<ValType> new_coeff) {coeff = new_coeff;}
        // const std::shared_ptr<std::vector<PauliOp> > getOps () const {return ops;};
        // const PauliOp& at (IdxType _ix) const {return ops->at(_ix);}
        bool operator==(const PauliOperator& other) const {
          
          return (xmask == other.xmask) && (zmask == other.zmask);
        }
        ValType sample_expectation(std::unordered_map<IdxType, IdxType> counts) const {
          ValType expect = 0.0;
          IdxType shots = 0;
          for (auto& i: counts) {
            IdxType count = i.second;
            ValType sign;
            bool par = parity(i.first, sign);
            shots += count;
            expect += sign * (par ? -count: count);
          }
          return expect / shots;
        }
        IdxType numQubits() const {return dim;};
    };
    struct PauliHash {
      std::size_t operator()(const PauliOperator& op) const {
        return std::hash<int>{} (op.get_xmask()) * std::hash<int>{} (op.get_zmask());
      }
    };
    std::ostream& operator<<(std::ostream& out, const PauliOperator& op);

    using PauliMap = std::unordered_map<PauliOperator, std::complex<ValType>, PauliHash>;
  };// namespace vqe
};// namespace nwqsim
#endif