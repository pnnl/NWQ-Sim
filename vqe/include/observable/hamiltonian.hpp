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



namespace NWQSim {
  namespace VQE {
    typedef std::unordered_map<PauliOperator, ValType, PauliHash> ExpectationMap;
    class Hamiltonian {
      protected:
        MolecularEnvironment env;
        std::vector<std::vector<FermionOperator> > fermi_operators;
        std::vector<std::vector<PauliOperator> > pauli_operators;
        IdxType n_ops;
        Transformer qubit_transform;
      
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
          transform(env, fermi_operators, pauli_operators, false);
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
        Hamiltonian(const std::vector<std::vector<PauliOperator> > _pauli_operators, 
                    IdxType n_particles,
                    bool use_xacc,
                    Transformer transform = getJordanWignerTransform):  pauli_operators(_pauli_operators) {
          IdxType num_qubits = pauli_operators.front().front().get_dim();
          IdxType n_spatial = num_qubits / 2;
          env = MolecularEnvironment(n_spatial, n_particles, use_xacc);
          n_ops = 0;
          for (auto& i : pauli_operators) {
            n_ops += i.size();
          }
        }
        Hamiltonian(MolecularEnvironment _env,
                    const std::vector<std::vector<PauliOperator> > _pauli_operators): env(_env), pauli_operators(_pauli_operators) {
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
        Transformer getTransformer () const {return qubit_transform;}
        const std::vector<std::vector<PauliOperator> >& getPauliOperators() const {return pauli_operators;};
        const std::vector<std::vector<FermionOperator> >& getFermionicOperators() const {return fermi_operators;};
    };
  };// namespace vqe
};// namespace nwqsim
#endif