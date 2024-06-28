#ifndef MEASUREMENT
#define MEASUREMENT
// #include "hamiltonian.hh"
#include <vector>
#include <string>
#include "environment.hpp"
#include "gate.hpp"
#include "observable/pauli_operator.hpp"
#include "transform/transform.hpp"
#include "circuit.hpp"
#include <unordered_set>
/**
 * Generic Measurement Class (assuming single operator) 
*/

namespace NWQSim {
  namespace VQE {
    class Measurement: public Circuit {
      protected:
        PauliOperator op;
        PauliOperator diagonal_op;
        double sign;
      public:
        Measurement(const PauliOperator& _op): op(_op), Circuit(_op.numQubits()) {
          size_t index = 0;
          sign = 1.0;
          IdxType dim = _op.get_dim();
          IdxType xmask = _op.get_xmask();
          IdxType zmask = _op.get_zmask();

          for (IdxType index = 0; index < dim; index++) {
            IdxType xbit = (xmask >> (1 << index)) >> index;
            IdxType zbit = (zmask >> (1 << index)) >> index;
            switch (xbit + 2 * zbit)
            {
            case 1: // X
              /* code */
              H(index);
              break;
            case 3: // Y
              RZ(-PI/2, index);
              H(index);
            
            default:
              break;
            }
          }
          
          diagonal_op = PauliOperator(0, op.get_xmask() | op.get_zmask(), op.get_dim());
          std::cout << diagonal_op << std::endl;
        };
        ValType operatorExpectation(std::unordered_map<IdxType, IdxType> counts) {
          ValType v = diagonal_op.sample_expectation(counts);
          return sign * v;
        }
    };
  };// namespace vqe
};// namespace nwqsim
#endif