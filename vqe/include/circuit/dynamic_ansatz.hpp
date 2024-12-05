#include "circuit/ansatz.hpp"
#ifndef DYNAMIC_ANSATZ_HH
#define DYNAMIC_ANSATZ_HH
namespace NWQSim { 
  namespace VQE {
    using FermionOplist = std::vector<FermionOperator>;
    using SymmetryGroup = std::vector<FermionOplist>; // for now, just assume we're preserving spin-conjugated symmetries
    enum class PoolType {
      Fermionic,
      SingletGSD,
      Pauli,
      MinimalPauli,
      FermionicOrigin
    };
    class DynamicAnsatz: public Ansatz {
      protected:
        const MolecularEnvironment& env;
        std::vector<SymmetryGroup> operator_pool; // list of Fermion operators, grouped according to spin symmetry
        std::vector<std::pair<IdxType, ValType> > symmetry_terms; // list of Fermion operators, grouped according to spin symmetry
        std::vector< std::vector<std::vector<FermionOperator> > > fermi_op_pool;
        std::vector< std::vector<std::vector<FermionOperator> > > fermi_op_selections;
        std::vector< PauliOperator > pauli_op_selections;
        std::vector<std::vector<PauliOperator> > pauli_op_pool;
        PoolType pool_type;

      public:
        DynamicAnsatz(const MolecularEnvironment& _env, PoolType _pool_type = PoolType::Fermionic): env(_env), Ansatz(_env.n_spatial * 2), pool_type(_pool_type) {
          ansatz_name = "Dynamic Ansatz"; // MZ: placeholder, will specify the type of ansatz in the make_op_pool()
        }
        virtual void buildAnsatz() override;
        void make_op_pool(Transformer tf, IdxType seed, IdxType pool_size) {
          if ((pool_type == PoolType::Fermionic) || (pool_type == PoolType::SingletGSD) || (pool_type == PoolType::FermionicOrigin)) {
            if (pool_type == PoolType::Fermionic) {
              generate_fermionic_excitations(fermi_op_pool, env);
              ansatz_name = "UCCSD Minimal (Dynamic)";
            } else if (pool_type == PoolType::SingletGSD) {
              generate_singletGSD_excitations(fermi_op_pool, env);
              ansatz_name = "Singlet GSD (Dynamic)";
            } else if (pool_type == PoolType::FermionicOrigin) {
              generate_fermionic_excitations_origin(fermi_op_pool, env);
              ansatz_name = "UCCSD Original (Dynamic)";
            }
            std::cout << "DEBUG: Pool Size " << fermi_op_pool.size() << std::endl;
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
          // } else if (pool_type == PoolType::SingletGSD) {
          //   generate_singletGSD_excitations(fermi_op_pool, env);
          //   pauli_op_pool.resize(fermi_op_pool.size());
          //   for (size_t i = 0; i < fermi_op_pool.size(); i++) {
          //     std::vector<std::vector<PauliOperator> > temp_pool;
          //     tf(env, fermi_op_pool[i], temp_pool, true);
          //     for (auto oplist: temp_pool) {
          //       for (auto op: oplist) {
          //         pauli_op_pool[i].push_back(op);
          //       }
          //     }
          //   }
          } else if (pool_type == PoolType::Pauli) {
            generate_pauli_excitations(pauli_op_pool, env, pool_size, seed);
            ansatz_name = "UCCSD Qubit (Dynamic)";
          } else if (pool_type == PoolType::MinimalPauli) {
            generate_minimal_pauli_excitations(pauli_op_pool, env);
            ansatz_name = "UCCSD Qubit (Minimal Pauli, Dynamic)";
          } else {
            throw std::runtime_error("Invalid pool type");
          }
          
        }
        const std::vector<std::vector<PauliOperator> >& get_pauli_op_pool() {
          return pauli_op_pool;
        }
        virtual std::string get_operator_string(IdxType operator_index) {
          std::string opstring = "";
          if (pool_type == PoolType::Fermionic) {
                size_t index = 0;
                auto symmetry_group = fermi_op_pool[operator_index];
                bool first = true;
                for (auto oplist: symmetry_group) {
                  if (!first) {
                    opstring = opstring + " ";
                  } else {
                    first = false;
                  }
                  std::string term = "";
                  for (auto& op: oplist) {
                    term = op.toString(env.n_occ, env.n_virt) + " " + term;
                  }
                  opstring += "(" + term + ")";
                  if ((index++) < oplist.size() - 1)
                    opstring = opstring + ",";
                }
          } else {
            // since we know that each operator sum in the pool is a single Pauli string
            PauliOperator pauli = pauli_op_pool[operator_index][0];
            opstring = pauli.pauliToString(false);
          }
          return opstring;
        }
        virtual void add_operator(IdxType operator_index, ValType param) {
          theta->push_back(param);
          std::vector<std::pair<IdxType, ValType > > idxvals = {{theta->size() - 1, 1.0}};
          if (pool_type == PoolType::Fermionic) {
            // For each Pauli term in the mapped pool, evolve using the same parameter
            for (PauliOperator pauli: pauli_op_pool[operator_index]) {
              if (pauli.isNonTrivial() && abs(pauli.getCoeff().real()) > 1e-10) {
                ExponentialGate(pauli, OP::RZ, idxvals, 2 * pauli.getCoeff().real());
              }
            }
            // slap that bad boy on the operator selection list
            fermi_op_selections.push_back(fermi_op_pool[operator_index]);
          } else {
            // since we know that each operator sum in the pool is a single Pauli string
            PauliOperator pauli = pauli_op_pool[operator_index][0];
            // First add the Pauli evolution to the ansatz
            ExponentialGate(pauli, OP::RZ, idxvals, 2 * pauli.getCoeff().real());
            // Append to the Pauli operator pool
            pauli_op_selections.push_back(pauli);
          }

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
          if (pool_type == PoolType::Fermionic) {
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
          } else {
            result.reserve(pauli_op_selections.size());
            for (size_t i = 0; i < pauli_op_selections.size(); i++) {
              const auto &opgroup = pauli_op_selections.at(i);
              ValType param = theta->at(i);
              std::string opstring = pauli_op_selections[i].pauliToString(false);
              result.push_back(std::make_pair(opstring, param));
            }
          }
          return result;
        };
        
    };
  };
};

#endif