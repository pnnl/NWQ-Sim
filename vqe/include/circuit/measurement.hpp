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
          std::vector<PauliOp> opstring;
          opstring.reserve(op.numQubits());
          for (PauliOp opval: *op.getOps()) {
            if (opval == PauliOp::I) {
              opstring.push_back(PauliOp::I);
              index++;
              continue;
            }
            if (opval == PauliOp::X) {
              H(index);
            }
            if (opval == PauliOp::Y) {
              S(index);
              H(index);
              sign *= -1;
            }
            opstring.push_back(PauliOp::Z);
            index++;
          }
          diagonal_op = PauliOperator(opstring);
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