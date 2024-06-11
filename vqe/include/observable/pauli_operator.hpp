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
namespace NWQSim {
  /* Basic data type for indices */
  using IdxType = long long int;
  /* Basic data type for value */
  using ValType = double;
  namespace VQE {
    const std::complex<ValType> Imag = {0.0, 1.0};
    enum PauliOp {
      I,
      X,
      Y,
      Z
    };
    // also determines phase, all anti-commuting operators add a pi/2 phase
    const bool commutationRelations[4][4] = {
      {1, 1, 1, 1},
      {1, 1, 0, 0},
      {1, 0, 1, 0},
      {1, 0, 0, 1}
    };
    // also determines phase, all anti-commuting operators add a pi/2 phase
    const PauliOp products[4][4] = {
      {PauliOp::I, PauliOp::X, PauliOp::Y, PauliOp::Z},
      {PauliOp::X, PauliOp::I, PauliOp::Z, PauliOp::Y},
      {PauliOp::Y, PauliOp::Z, PauliOp::I, PauliOp::X},
      {PauliOp::Z, PauliOp::Y, PauliOp::X, PauliOp::I}
    };
    // Sign relations for pauli products
    const int signRelations[4][4] = {
      {1,  1,   1,   1},   // I commutes with everything
      {1,  1,   1,  -1}, // XI=X, XX=I, XY=iZ, XZ=-iY
      {1, -1 ,  1,   1}, // YI=Y, YX=-iZ, YY=I, YZ=iX
      {1,  1,  -1,   1} // ZI=Z, ZX=iY, ZY=-iX, ZZ=I
    };
    
    const char *const PAULI_OP_NAMES[] = {
      "I",
      "X",
      "Y",
      "Z"
    };

    class PauliOperator {
      private:
        IdxType dim;
        bool non_trivial;
        std::complex<ValType> coeff;
        IdxType xmask;
        IdxType zmask;
        ValType sign;
        std::vector<IdxType> x_indices;
        bool phase;
        std::shared_ptr<std::vector<PauliOp> > ops; // Store in big-endian order v, i.e. qubit 0 is in position 0
      public:
        PauliOperator() {};
        PauliOperator(const std::vector<PauliOp>& _ops,
                      std::complex<ValType> _coeff = 1.0): dim(_ops.size()), coeff(_coeff) {
          ops = std::make_shared<std::vector<PauliOp> >(_ops);
          setExpectationVariables();
        } 
        void setExpectationVariables() {
          xmask = 0;
          zmask = 0;
          phase = 0;
          IdxType index = 0;
          for (auto op: *ops) {
            switch (op)
            {
            case X:
              xmask |= 1 << index;
              x_indices.push_back(index);
              break;
            case Y:
              xmask |= 1 << index;
              zmask |= 1 << index;
              x_indices.push_back(index);
              phase += 1;
              break;
            case Z:
              zmask |= 1 << index;
              break;
            
            default:
              break;
            }
            index++;
          }
          // if even number of imaginary contributions, sign is negative
          sign = (phase >> 1) & 1 ? 1 : -1;
          // either i or 0
          phase = phase & 1;
          non_trivial = (xmask | zmask) > 0;
        }
        PauliOperator(std::string _opstring,
                      std::complex<ValType> _coeff = 1.0): dim(_opstring.length()), coeff(_coeff) {
          ops = std::make_shared<std::vector<PauliOp> >(dim);
          size_t index = 0;
          non_trivial = false;
          for (char i: _opstring) {
            assert(index < dim);
            switch (i)
            {
            case 'X':
              ops->at(dim - (index++) - 1) = PauliOp::X;
              break;
            case 'Y':
              ops->at(dim - (index++) - 1) = PauliOp::Y;
              break;
            case 'Z':
              ops->at(dim - (index++) - 1) = PauliOp::Z;
              break;
            case 'I':
              ops->at(dim - (index++) - 1) = PauliOp::I;
              break;
            
            default:
              assert(false);
              break;
            }
          }
          setExpectationVariables();
        } 

        PauliOperator& operator=(const PauliOperator& other) {
          dim = other.dim;
          ops = std::make_shared<std::vector<PauliOp> >(
            other.ops.get()->begin(), other.ops.get()->end()); 
          coeff = other.coeff;
          non_trivial = other.non_trivial;
          xmask = other.xmask;
          zmask = other.zmask;
          phase = other.phase;
          sign = other.sign;
          x_indices = other.x_indices;

          // assert(coeff. > 0.0);
          return *this;
        }
        PauliOperator(const PauliOperator& other) {
          dim = other.dim;
          ops = std::make_shared<std::vector<PauliOp> >(
            other.ops.get()->begin(), other.ops.get()->end()); 
          coeff = other.coeff;
          non_trivial = other.non_trivial;
          xmask = other.xmask;
          zmask = other.zmask;
          phase = other.phase;
          sign = other.sign;
          x_indices = other.x_indices;
          // assert(coeff.real() > 0.0 || coeff.imag() > 0.0);
        }
        PauliOperator conj() const {
          return PauliOperator(*ops, std::complex<ValType>(coeff.real(), -coeff.imag()));
        }
        PauliOperator operator*(const PauliOperator& other) const {
          std::complex<ValType> new_coeff = coeff * other.coeff;
          IdxType dim1 = ops->size();
          IdxType dim2 = other.ops->size();
          if (dim1 != dim2) {
            fprintf(stderr, "Pauli strings of different size: %s %s\n", pauliToString().c_str(), other.pauliToString().c_str());
            assert(false);
          }
          std::vector<PauliOp> opstring (dim1);

          IdxType minval = std::min(dim1, dim2);
          for (IdxType i = 0; i < minval; i++) {
            // Ensure self * other ordering
            PauliOp op1 = ops.get()->at(i);
            PauliOp op2 = other.ops.get()->at(i);
            opstring[i] = products[op1][op2];

            std::complex<ValType> new_contribution;
            if (commutationRelations[op1][op2]) {
              new_contribution = std::complex<ValType> (signRelations[op1][op2], 0.0);
            } else {
              new_contribution = std::complex<ValType> (0.0, signRelations[op1][op2]);
            }
            new_coeff *= new_contribution;
          }
          return PauliOperator(opstring, new_coeff);
          
        }
        PauliOperator& operator*=(ValType scalar){
          coeff *= scalar;
          return *this;
        }
        PauliOperator operator*(ValType scalar) const {
          PauliOperator newop (*this);
          newop.coeff *= scalar;
          return newop;
        }
        const std::vector<IdxType>& get_xindices() const {
          return x_indices;
        }
        const IdxType get_xmask() const {
          return xmask;
        }
        const IdxType get_zmask() const {
          return zmask;
        }
        // Dump the Pauli operator to string
        std::string pauliToString() const {
          std::stringstream ss;
          std::vector<PauliOp> reversed(dim);
          if (coeff != 1.0) {
            ss << "(" << coeff.real() << " ";
            if (coeff.imag() >= 0) {
              ss << "+ " << coeff.imag() << "i)";
            } else {
              ss << "- " << fabs(coeff.imag()) << "i)";
            }
            // ss
          }
          // ss << "*";
          // Reverse to ensure correct order (little endian)
          std::reverse_copy(ops->begin(), ops->end(), reversed.begin());
          for (int pauli_index: reversed) {
            ss << PAULI_OP_NAMES[pauli_index];
          }
          return ss.str();
        }
        bool parity(IdxType other, ValType& sign) const {
          bool parity_xor = 0;
          sign = 1;
          for (IdxType i = 0; i < dim; i++) {
            PauliOp op = ops->at(i);
            if (op == PauliOp::I) {
              continue;
            }
            // Y Pauli diagonalized via an S gate w/ a Hadamard
            bool bit = (other & (1 << i)) > 0;
            if (op == PauliOp::Y && bit) {
              // parity_xor ^= 1;
              sign *= -1;
            }
            // assert (op == PauliOp::Z);
            parity_xor ^= bit;
          }
          return parity_xor;
        }
        bool parity(const PauliOperator& other) const {
          IdxType mindim = std::min(other.dim, dim);
          bool parity_xor = 0;
          for (IdxType i = 0; i < mindim; i++) {
            parity_xor = parity_xor ^ commutationRelations[ops.get()->at(i)][other.ops.get()->at(i)];
          }
          return parity_xor;
        }
        bool isNonTrivial() const { return non_trivial;}
        bool QWC(PauliOperator& other) {
          IdxType mindim = std::min(other.dim, dim);
          bool qwc = 1;
          for (IdxType i = 0; i < mindim; i++) {
            qwc = qwc && commutationRelations[ops.get()->at(i)][other.ops.get()->at(i)];
          }
          return qwc;
        }
        bool GC(PauliOperator& other) {
          return parity(other);
        }
        std::complex<ValType> getCoeff() const {return coeff;}
        void setCoeff (std::complex<ValType> new_coeff) {coeff = new_coeff;}
        const std::shared_ptr<std::vector<PauliOp> > getOps () const {return ops;};
        const PauliOp& at (IdxType _ix) const {return ops->at(_ix);}
        bool operator==(const PauliOperator& other) const {
          if (other.dim != dim) {
            // NOTE: This disregards redundant identities
            return false;
          }
          for (size_t i = 0; i < dim; i++) {
            if (ops->at(i) != other.ops->at(i)) {
              return false;
            }
          }
          return true;
        }
        ValType sample_expectation(std::unordered_map<IdxType, IdxType> counts) const {
          ValType expect = 0.0;
          IdxType shots = 0;
          for (auto& i: counts) {
            IdxType observation = i.first;
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
        std::size_t prodval = 1lu;
        for (PauliOp op: *op.getOps()) {
          prodval *= static_cast<size_t>(op);
        }
        return std::hash<int>{} (prodval);
      }
    };
    std::ostream& operator<<(std::ostream& out, const PauliOperator& op);

    using PauliMap = std::unordered_map<PauliOperator, std::complex<ValType>, PauliHash>;
  };// namespace vqe
};// namespace nwqsim
#endif