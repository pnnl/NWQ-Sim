#ifndef ANSATZ
#define ANSATZ
// #include "hamiltonian.hh"
#include <vector>
#include <string>
#include "environment.hpp"
#include "gate.hpp"
#include "observable/fermionic_operator.hpp"
#include "observable/pauli_operator.hpp"
#include "transform/transform.hpp"
#include "circuit.hpp"
#include <unordered_set>
/**
 * Generic Ansatz Class (to be inherited by specific types, e.g. UCCSD) 
*/

namespace NWQSim {
  namespace VQE {
    class Ansatz: public Circuit {
      protected:
        std::shared_ptr<std::vector<ValType> > theta;
        std::shared_ptr<std::vector<IdxType> > parameterized_gates; // Indices of parameterized gates
        std::shared_ptr<std::vector<IdxType> > gate_parameter_pointers; // Gate parameter indices
        std::shared_ptr<std::vector<ValType> > gate_coefficients; // Gate coefficients

      public:
        Ansatz(IdxType n_qubits): Circuit(n_qubits) {
          theta = std::make_shared<std::vector<ValType> >();
          gate_coefficients = std::make_shared<std::vector<ValType> >();
          parameterized_gates = std::make_shared<std::vector<IdxType> >();
          gate_parameter_pointers = std::make_shared<std::vector<IdxType> >();
          };
        ~Ansatz() {};
        void assignGateParam(IdxType gate_index, ValType param) {
          Gate& g = gates->at(gate_index);
          g.theta = param;
          /* The `theta` field is used by all of the single qubit gates currently supported,
              RI, RX, RY, RZ, P
          */ 
          // 
        }
        void setParams(const std::vector<ValType> params) {
          assert (params.size() == theta->size());
          std::copy(params.begin(), params.end(), theta->begin());
          IdxType index = 0;
          for (IdxType gate_id: *parameterized_gates) {
            ValType paramval =  theta->at(gate_parameter_pointers->at(index));
            paramval *= gate_coefficients->at(index);
            assignGateParam(gate_id, paramval);
            index++;
          }
        }
        std::string toQASM3() {
          std::stringstream outstream;
          outstream << "OPENQASM 3.0;\ninclude \"stdgates.inc\";\n";
          outstream << "qubit[" << num_qubits() << "] q;\n";
          for (size_t i = 0; i < theta->size(); i++) {
            outstream << "input float[64] theta_" << i << ";\n";
          }
          size_t param_ix = 0;
          size_t gate_ix = 0;
          for (Gate& g: *gates) {
            std::string lower_name = OP_NAMES[g.op_name];
            if (lower_name == "expect") {
              continue;
            }
            std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(),
              [](unsigned char c){ return std::tolower(c); });

            outstream << lower_name;
            // std::cout << gate_ix 
            if (param_ix < parameterized_gates->size() and gate_ix == parameterized_gates->at(param_ix)) {
              outstream << "(" << gate_coefficients->at(param_ix) << "*theta_" << gate_parameter_pointers->at(param_ix) << ")";
              param_ix++;
            } else {
              switch(g.op_name) {
                case OP::RX:
                case OP::RY:
                case OP::RZ:
                case OP::RI:
                case OP::P:
                  outstream << "(" << g.theta << ")";
                  break;
                case OP::CX: {
                  outstream << " q[" << g.ctrl << "],";
                  break;
                }
                default:
                break;
              }
            }
            outstream << " q[" << g.qubit << "];\n";
            gate_ix++;
          }

          return outstream.str();
        }
        virtual std::vector<std::string> getFermionicOperatorStrings() const {
          std::vector<std::string> result;
          throw std::runtime_error("Fermionic operators not specified for this ansatz\n");
          return result;
        };

        Ansatz compose(const Circuit& other, std::vector<IdxType>& qubit_mapping) {
          Ansatz composition = Ansatz(*this);
          for (IdxType i = 0; i < num_qubits(); i++) {
            if (i != qubit_mapping[i]) {
              composition.SWAP(i, qubit_mapping[i]);
            }
          }
          // std::concate
          composition.gates->insert(composition.gates->end(), other.gates->begin(), other.gates->end());
          // for (size_t i = 0; i < other->gates.size(); i++) {

          // }
          // std::copy(other.gates->begin(), other.gates->end(), composition.gates->begin() + gates->size());
          return composition;
        }
        // Accessors
        const IdxType numParams() const { return theta->size(); };
        const std::shared_ptr<std::vector<IdxType> > getParamGateIndices() const {return parameterized_gates;}
        const std::shared_ptr<std::vector<IdxType> > getParamGatePointers() const {return gate_parameter_pointers;}
        std::vector<ValType> getGateParams() const {
          std::vector<ValType> result(parameterized_gates->size());
          // Get the parameter values for each gate
          std::generate(result.begin(), result.end(), 
          [&, n = 0] () mutable {
            return gates->at(parameterized_gates->at(n++)).theta;
          } );
          return result;}
        // const/Non-const access
        const std::shared_ptr<std::vector<ValType> > getParams() const {return theta;}
        std::vector<ValType>* getParams() {return theta.get();}

        void OneParamGate(enum OP _op_name,
              IdxType _qubit,
              IdxType _ctrl = -1,
              IdxType _n_qubits = 1,
              IdxType _param_idx = -1,
              ValType _coeff = 1.0,
              ValType _param = 0,
              IdxType _repetition = 0) {
          IdxType index = gates->size();
          IdxType param_index = _param_idx;
          ValType param_value = _param;
          if (param_index == -1) {
            param_index = theta->size();
            theta->push_back(param_value);
          } else {
            param_value = theta->at(param_index);
          }
          gate_coefficients->push_back(_coeff);
          gate_parameter_pointers->push_back(param_index);
          parameterized_gates->push_back(index);
          switch(_op_name) {
            case OP::RI: {
              RI(param_value, _qubit);
              break;
            }
            case OP::RX: {
              RX(param_value, _qubit);
              break;
            }
            case OP::RY: {

              RY(param_value, _qubit);
              break;
            }
            case OP::RZ: {
              RZ(param_value, _qubit);
              break;
            }
            case OP::P: {
              P(param_value, _qubit);
              break;
            }
            default:
              fprintf(stderr, "%s is not a supported 1-Param Gate\n", OP_NAMES[(int)_op_name]);
              assert(false);
          }
        };   
        void ExponentialGate(const PauliOperator& pauli_op,
                             enum OP _op_name, 
                             ValType _coeff = 1.0, 
                             ValType _param = 0.0, 
                             IdxType param_index = -1) {
          std::vector<IdxType> non_trivial;
          IdxType index = 0;
          for (PauliOp op: *pauli_op.getOps()) {
            if (op != PauliOp::I) {
              non_trivial.push_back(index);
              if (op == PauliOp::X) {
                H(index);
              }
              if (op == PauliOp::Y) {
                RZ(-PI/2, index);
                H(index);
              }
            }
            index++;
          }
          for (IdxType i = non_trivial.size() - 1; i >= 1; i--) {
            CX(non_trivial[i], non_trivial[i-1]);
          }
          OneParamGate(_op_name, non_trivial.front(), -1, 1, param_index, _coeff, _param, 0);
          for (IdxType i = 0; i < non_trivial.size()-1; i++) {
            CX(non_trivial[i+1], non_trivial[i]);
          }
          for (IdxType idx: non_trivial) {
            PauliOp op = pauli_op.at(idx);
            if (op != PauliOp::I) {
              if (op == PauliOp::X) {
                H(idx);
              }
              if (op == PauliOp::Y) {
                H(idx);
                RZ(PI / 2, idx);
              }
            }
          }
        }
    };
    class UCCSD: public Ansatz {
      protected:
        const MolecularEnvironment& env;
        IdxType n_singles;
        IdxType n_doubles;
        IdxType trotter_n;
        std::vector<std::vector<FermionOperator> > fermion_operators;
        void getFermionOps();
        void buildAnsatz(std::vector<std::vector<PauliOperator> > pauli_oplist);
      public:
        UCCSD(const MolecularEnvironment& _env, Transformer transform, IdxType _trotter_n = 1): 
                                  env(_env),
                                  trotter_n(_trotter_n),
                                  Ansatz(2 * _env.n_spatial) {
          n_singles = 2 * env.n_occ * env.n_virt;
          n_doubles = 2 * env.n_occ * env.n_occ * env.n_virt * env.n_virt;
          fermion_operators.reserve(n_singles + n_doubles);
          getFermionOps();
          std::vector<std::vector<PauliOperator> > pauli_ops;
          pauli_ops.reserve(4 * n_singles + 16 * n_doubles);
          transform(env, fermion_operators, pauli_ops, true);  
          buildAnsatz(pauli_ops);
        };
        virtual std::vector<std::string> getFermionicOperatorStrings() const override {
          std::vector<std::string> result;
          result.reserve(fermion_operators.size());
          for (auto& oplist : fermion_operators) {
            std::ostringstream opstream;
            bool first = true;
            for (auto& op: oplist) {
              if (!first) {
                opstream << " ";
              } else {
                first = false;
              }
              opstream << op.toString(env.n_occ, env.n_virt);
            }
            result.push_back(opstream.str());
          }
          return result;
        };
        const MolecularEnvironment& getEnv() const {return env;};
    };
  };// namespace vqe
};// namespace nwqsim
#endif