#ifndef HAMIL
#define HAMIL
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <memory>
#include "observable/fermionic_operator.hpp"
#include "observable/pauli_operator.hpp"
#include "transform/transform.hpp"
#include "utils.hpp" 
#include "state.hpp" 



namespace NWQSim {
  namespace VQE {
    typedef std::unordered_map<PauliOperator, ValType, PauliHash> ExpectationMap;
    class Hamiltonian {
      protected:
        MolecularEnvironment env;
        std::vector<std::vector<FermionOperator> > fermi_operators;
        std::vector<std::vector<PauliOperator> > pauli_operators;
        IdxType n_ops;
      
      public:
        Hamiltonian(std::string input_path, IdxType n_particles, bool xacc_scheme,
                    Transformer transform = getJordanWignerTransform);
        Hamiltonian(const std::vector<std::pair<std::string, std::complex<double>>>& input_ops, 
                    IdxType n_particles, bool xacc_scheme,
                    Transformer transform);
        Hamiltonian(MolecularEnvironment _env,
                    const std::vector<std::vector<FermionOperator> > _fermion_operators, 
                    Transformer transform = getJordanWignerTransform): env(_env), fermi_operators(_fermion_operators) {
          n_ops = 0;
          std::vector<PauliOperator> single_oplist;
          transform(env, fermi_operators, single_oplist, false);
          std::list<std::vector<IdxType>> pauli_cliques;
          sorted_insertion(single_oplist, pauli_cliques, false);
          pauli_operators.resize(pauli_cliques.size());
          IdxType index = 0;
          for (auto& clique: pauli_cliques) {
            pauli_operators[index].reserve(clique.size());
            for (IdxType pauli_idx: clique) {
              pauli_operators[index].push_back(single_oplist[pauli_idx]);
            }
            index++;
          }
          for (auto& i : pauli_operators) {
            n_ops += i.size();
          }
        }
        Hamiltonian(MolecularEnvironment _env,
                    const std::vector<std::vector<PauliOperator> > _pauli_operators, 
                    Transformer transform = getJordanWignerTransform): env(_env), pauli_operators(_pauli_operators) {
          n_ops = 0;
          for (auto& i : pauli_operators) {
            n_ops += i.size();
          }
        }
        void
        get_pauli_coeffs(std::vector<double> & result) const {
          // result.resize(pauli_operators.size;
          result.reserve(pauli_operators.size());
          for (auto& pauli_list: pauli_operators) {
            for (auto& pauli: pauli_list) {
              result.push_back(pauli.getCoeff().real());
            }
          }
        };
        ValType expectation(const ExpectationMap& expectations) const {
          ValType expect = env.constant;
          for (auto& pauli_list: pauli_operators) {
            for (auto& pauli: pauli_list) {
              if (expectations.find(pauli) == expectations.end()) {
                throw std::runtime_error("Pauli expectation not provided");
              }
              expect += pauli.getCoeff().real() * expectations.at(pauli);
            }
          }
          return expect;
        }
        IdxType num_ops () const {
          IdxType result = 0;
          for (auto& i : pauli_operators) {
            result += i.size();
          }
          return result;
        }
        const MolecularEnvironment& getEnv () const {return env;}
        const std::vector<std::vector<PauliOperator> >& getPauliOperators() const {return pauli_operators;};
        const std::vector<std::vector<FermionOperator> >& getFermionicOperators() const {return fermi_operators;};
    };
    void make_commutators(std::shared_ptr<Hamiltonian> hamil,
                          const std::vector<std::vector<PauliOperator> >& pauli_op_pool,
                          std::vector<std::vector<PauliOperator> >& commutator_list);
  };// namespace vqe
};// namespace nwqsim
#endif