#include "observable/pauli_operator.hpp"
#include "utils.hpp"
#include <string>
#include <iostream>

namespace NWQSim{
namespace VQE {
  std::ostream& operator<<(std::ostream& out, const PauliOperator& op) {
        out << op.pauliToString();
        return out;
      };
  std::string to_binary_string(NWQSim::IdxType val, NWQSim::IdxType n_qubits) {
    std::stringstream ss;
    for (IdxType i = n_qubits-1; i >= 0; i--) {
      ss << ((val & (1 << i)) ? "1" : "0");
    }

    return ss.str();
  }
};// namespace NWQSim

};// namespace VQE
