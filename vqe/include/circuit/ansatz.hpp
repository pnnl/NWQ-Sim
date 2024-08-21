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
#include "utils.hpp"
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
        std::shared_ptr<std::vector<std::vector<std::pair<IdxType, ValType> > > > gate_parameter_pointers; // Gate parameter indices
        std::shared_ptr<std::vector<ValType> > gate_coefficients; // Gate coefficients
        std::unordered_map<std::string, IdxType> excitation_index_map;
      public:
        Ansatz(IdxType n_qubits): Circuit(n_qubits) {
          theta = std::make_shared<std::vector<ValType> >();
          gate_coefficients = std::make_shared<std::vector<ValType> >();
          parameterized_gates = std::make_shared<std::vector<IdxType> >();
          gate_parameter_pointers = std::make_shared<std::vector<std::vector<std::pair<IdxType, ValType> > > >();
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
        virtual 
        void setParams(const std::vector<ValType>& params) {
          assert (params.size() == theta->size());
          std::copy(params.begin(), params.end(), theta->begin());
          IdxType index = 0;
          for (IdxType gate_id: *parameterized_gates) {
            ValType paramval = 0;
            for (auto pair: gate_parameter_pointers->at(index)) {
              paramval += theta->at(pair.first) * pair.second;
            }
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
            if (g.op_name == OP::EXPECT) {
              continue;
            }
            std::transform(lower_name.begin(), lower_name.end(), lower_name.begin(),
              [](unsigned char c){ return std::tolower(c); });

            outstream << lower_name;
            if (false && param_ix < parameterized_gates->size() && gate_ix == parameterized_gates->at(param_ix)) {
              outstream << "(" << gate_coefficients->at(param_ix) << "*(";
              std::vector<std::pair<IdxType, ValType> > expr_vec = gate_parameter_pointers->at(param_ix);
              for (IdxType i = 0; i < expr_vec.size(); i++) {
                outstream << expr_vec.at(i).second << "*theta_" << expr_vec.at(i).first;
                if (i < expr_vec.size() - 1) {
                  outstream << " * ";
                } else {
                  outstream << ")";
                }
              }
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
                default:
                if (g.n_qubits == 2) {
                  outstream << " q[" << g.ctrl << "],";
                }
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
        };
        virtual std::vector<std::pair<std::string, ValType> > getFermionicOperatorParameters() const {
          std::vector<std::pair<std::string, ValType> > result;
          throw std::runtime_error("Fermionic operators not specified for this ansatz\n");
        }
        void compose(const Circuit& other, std::vector<IdxType>& qubit_mapping) {
          // Ansatz composition = Ansatz(*this);
          for (IdxType i = 0; i < num_qubits(); i++) {
            if (i != qubit_mapping[i]) {
              SWAP(i, qubit_mapping[i]);
            }
          }
          // std::concate
          gates->insert(gates->end(), other.gates->begin(), other.gates->end());
          // for (size_t i = 0; i < other->gates.size(); i++) {

          // }
          // std::copy(other.gates->begin(), other.gates->end(), composition.gates->begin() + gates->size());
          // return composition;
        }
        // Accessors
        virtual IdxType numParams() const { return theta->size(); };
        const std::shared_ptr<std::vector<IdxType> > getParamGateIndices() const {return parameterized_gates;}
        const std::shared_ptr<std::vector<std::vector<std::pair<IdxType, ValType> > > > getParamGatePointers() const {return gate_parameter_pointers;}
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
        std::vector<ValType>& getParamRef() {return *theta.get();}

        void OneParamGate(enum OP _op_name,
              IdxType _qubit,
              IdxType _ctrl = -1,
              IdxType _n_qubits = 1,
              std::vector<std::pair<IdxType, ValType> > _param_ids = {{-1, 0.0}},
              ValType _coeff = 1.0,
              IdxType _repetition = 0) {
          IdxType index = gates->size();
          auto param_index = _param_ids;
          ValType param_value = 0.0;
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
                             std::vector<std::pair<IdxType, ValType> > _param,
                             ValType _coeff = 1.0) {
          std::vector<IdxType> non_trivial;
          IdxType dim = pauli_op.get_dim();
          IdxType xmask = pauli_op.get_xmask();
          IdxType zmask = pauli_op.get_zmask();

          for (IdxType index = 0; index < dim; index++) {
            IdxType xbit = (xmask >> index) & 1;
            IdxType zbit = (zmask >> index) & 1;
            switch (xbit + 2 * zbit)
            {
            case 1: // X
              /* code */
              non_trivial.push_back(index);
              H(index);
              break;
            case 2: // Z
              non_trivial.push_back(index);
              break;

            case 3: // Y
              non_trivial.push_back(index);
              RZ(-PI/2, index);
              H(index);
              break;
            
            default:
              break;
            }
          }
          for (IdxType i = non_trivial.size() - 1; i >= 1; i--) {
            CX(non_trivial[i], non_trivial[i-1]);
          }
          OneParamGate(_op_name, non_trivial.front(), -1, 1, _param, _coeff, 0);
          for (IdxType i = 0; i < non_trivial.size()-1; i++) {
            CX(non_trivial[i+1], non_trivial[i]);
          }
          for (IdxType index: non_trivial) {
            IdxType xbit = (xmask >> index) & 1;
            IdxType zbit = (zmask >> index) & 1;
            switch (xbit + 2 * zbit)
            {
            case 1: // X
              /* code */
              H(index);
              break;
            case 2: // Z
              break;

            case 3: // Y
              H(index);
              RZ(PI/2, index);
              break;
      
            default:
              break;
            }
          }
        }
        virtual void buildAnsatz() {};
    };
    class UCCSD: public Ansatz {
      protected:
        const MolecularEnvironment& env;
        IdxType n_singles;
        IdxType n_doubles;
        IdxType trotter_n;
        IdxType unique_params;
        Transformer qubit_transform;
        // bool enforce_symmetries;
        /** 
         * Enforce symmetries for each term. Each fermionic term will have one symmetry entry. If no symmetries are enforced, 
         * symmetries[i] = {{i, 1.0}}; Otherwise, symmetries[i] = {{j, 1.0}, {k, -1.0}} denotes that theta_i must be equal to theta_j - theta_k
         */
        std::vector<std::vector<std::pair<IdxType, ValType> > > symmetries;
        std::vector<IdxType> fermion_ops_to_params; // map from fermion operators to parameters (used in update)
        std::vector<std::vector<FermionOperator> > fermion_operators;
        virtual void getFermionOps();
      public:
        UCCSD(const MolecularEnvironment& _env, Transformer _qubit_transform, IdxType _trotter_n = 1): 
                                  env(_env),
                                  trotter_n(_trotter_n),
                                  qubit_transform(_qubit_transform),
                                  Ansatz(2 * _env.n_spatial) {
          n_singles = 2 * env.n_occ * env.n_virt;
          IdxType c2virtual = choose2(env.n_virt);
          IdxType c2occupied = choose2(env.n_occ);
          n_doubles = 2 * (env.n_occ) * c2virtual + 2 * (env.n_virt) * c2occupied + env.n_occ * env.n_virt +\
              c2occupied * c2virtual * 4;
          fermion_operators.reserve(n_singles + n_doubles);
          symmetries = std::vector<std::vector<std::pair<IdxType, ValType> > >((n_singles + n_doubles));
          fermion_ops_to_params.resize(n_doubles + n_singles);
          std::fill(fermion_ops_to_params.begin(), fermion_ops_to_params.end(), -1);
          unique_params = 0;
          
        };
        virtual void buildAnsatz() override;
        virtual std::vector<std::string> getFermionicOperatorStrings() const override {
          std::vector<std::string> result;
          result.reserve(fermion_operators.size());
          for (auto& oplist : fermion_operators) {
            std::string opstring = "";
            bool first = true;
            for (auto& op: oplist) {
              if (!first) {
                opstring = " " + opstring;
              } else {
                first = false;
              }
              opstring = op.toString(env.n_occ, env.n_virt) + opstring;
            }
            result.push_back(opstring);
          }
          return result;
        };
        virtual std::vector<std::pair<std::string, ValType> > getFermionicOperatorParameters() const override {
          std::vector<std::pair<std::string, ValType> > result;
          result.reserve(fermion_operators.size());
          for (size_t i = 0; i < fermion_operators.size(); i++) {
            const auto &oplist = fermion_operators.at(i);
            const std::vector<std::pair<IdxType, ValType> > &param_expr = symmetries[i];
            ValType param = 0.0;
            for (auto& i: param_expr) {
              param += i.second * theta->at(fermion_ops_to_params[i.first]);
            }
            std::string opstring = "";
            bool first = true;
            bool is_distinct = true;
            bool is_mixed = true;
            std::vector<IdxType> indices_seen(env.n_spatial, 0);
            for (auto& op: oplist) {
              if (!first) {
                opstring = " " + opstring;
              } else {
                first = false;
                is_mixed = op.getSpin() != oplist[1].getSpin();
              }
              if (indices_seen[op.getOrbitalIndex(env)]) {
                is_distinct = false;
              }
              indices_seen[op.getOrbitalIndex(env)] += 1;
              opstring = op.toString(env.n_occ, env.n_virt) + opstring;
            }
            result.push_back(std::make_pair(opstring, param));
            if (is_distinct && oplist.size() == 4 && is_mixed) {
              first = true;
              opstring = "";
              for (auto& op: oplist) {
                if (!first) {
                  opstring = " " + opstring;
                } else {
                  first = false;
                }
                opstring = op.spinReversed().toString(env.n_occ, env.n_virt) + opstring;
              }
              result.push_back(std::make_pair(opstring, param));
            }
          }
          return result;
        };
        
        const MolecularEnvironment& getEnv() const {return env;};
        virtual IdxType numParams() const override { return unique_params; };
    };
    
  };// namespace vqe
};// namespace nwqsim
#endif