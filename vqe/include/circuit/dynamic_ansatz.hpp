#include "circuit/ansatz.hpp"
#ifndef DYNAMIC_ANSATZ_HH
#define DYNAMIC_ANSATZ_HH
namespace NWQSim { 
  namespace VQE {
    using FermionOplist = std::vector<FermionOperator>;
    using SymmetryGroup = std::vector<FermionOplist>; // for now, just assume we're preserving spin-conjugated symmetries

    class DynamicAnsatz: public Ansatz {
      protected:
        const MolecularEnvironment& env;
        std::vector<SymmetryGroup> operator_pool; // list of Fermion operators, grouped according to spin symmetry
        std::vector<std::pair<IdxType, ValType> > symmetry_terms; // list of Fermion operators, grouped according to spin symmetry
        std::vector< std::vector<std::vector<FermionOperator> > > fermi_op_pool;
        std::vector< std::vector<std::vector<FermionOperator> > > fermi_op_selections;
        std::vector<std::vector<PauliOperator> > pauli_op_pool;
      
      public:
        DynamicAnsatz(const MolecularEnvironment& _env): env(_env), Ansatz(_env.n_spatial * 2) {}
        virtual void buildAnsatz() override;
        void make_op_pool(Transformer tf) {
          generate_excitations(fermi_op_pool, env);
          pauli_op_pool.resize(fermi_op_pool.size());
          for (size_t i = 0; i < fermi_op_pool.size(); i++) {

            std::vector<std::vector<PauliOperator> > temp_pool;
            tf(env, fermi_op_pool[i], temp_pool, true);
            for (auto oplist: temp_pool) {
              for (auto op: oplist) {
                pauli_op_pool[i].push_back(op);
              }
            }
          } 
        }
        const std::vector<std::vector<PauliOperator> >& get_pauli_op_pool() {
          return pauli_op_pool;
        }
        virtual void add_operator(IdxType operator_index, ValType param) {
          theta->push_back(param);
          std::vector<std::pair<IdxType, ValType > > idxvals = {{theta->size() - 1, 1.0}};
          // std::cout << "Adding operator ";
          // size_t index = 0;
          // size_t num_opstrings = fermi_op_pool[operator_index].size();
          // for (auto opstring: fermi_op_pool[operator_index]) {
          //   size_t num_ops = opstring.size();
          //   size_t o_index = 0;
          //   for (auto single_qubit_op: opstring) {
          //     std::cout << single_qubit_op.toString(env.n_occ, env.n_virt);
          //     if ((o_index++) < num_ops - 1)
          //       std::cout << " ";
          //   }
          //   if ((index ++) < num_opstrings - 1) {
          //     std::cout << ", ";
          //   }
          // }
          // std::cout << std::endl;
          for (PauliOperator pauli: pauli_op_pool[operator_index]) {
            if (pauli.isNonTrivial() && abs(pauli.getCoeff().real()) > 1e-10) {
              ExponentialGate(pauli, OP::RZ, idxvals, 2 * pauli.getCoeff().real());
            }
          }
          fermi_op_selections.push_back(fermi_op_pool[operator_index]);
        }; 
      virtual std::vector<std::string> getFermionicOperatorStrings() const override {
          std::vector<std::string> result;
          result.reserve(fermi_op_selections.size());
          for (size_t i = 0; i < fermi_op_selections.size(); i++) {
            const auto &opgroup = fermi_op_selections.at(i);
            std::string opstring = "";
            bool first = true;
            for (auto oplist: opgroup) {
              for (auto& op: oplist) {
                if (!first) {
                  opstring = " " + opstring;
                } else {
                  first = false;
                }
                opstring = op.toString(env.n_occ, env.n_virt) + opstring;
              }
              if (opgroup.size() > 1)
                opstring = opstring + ", ";
            }
            result.push_back(opstring);
          }
          return result;
        };
        virtual std::vector<std::pair<std::string, ValType> > getFermionicOperatorParameters() const override {
          std::vector<std::pair<std::string, ValType> > result;
          result.reserve(fermi_op_selections.size());
          for (size_t i = 0; i < fermi_op_selections.size(); i++) {
            const auto &opgroup = fermi_op_selections.at(i);
            ValType param = theta->at(i);
            std::string opstring = "";
            bool first = true;
            for (auto oplist: opgroup) {
              if (!first) {
                opstring = opstring + " ";
              } else {
                first = false;
              }
              std::string term = "";
              size_t index = 0;
              for (auto& op: oplist) {
                term = op.toString(env.n_occ, env.n_virt) + " " + term;
              }
              opstring += "(" + term + ")";
              if ((index++) < oplist.size() - 1)
                opstring = opstring + ",";
            }
            result.push_back(std::make_pair(opstring, param));
          }
          return result;
        };
        
    };
  };
};

#endif