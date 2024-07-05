#ifndef VQE_ADAPT_HH
#define VQE_ADAPT_HH
#include "vqe_state.hpp"
#include "circuit/dynamic_ansatz.hpp"
namespace NWQSim {
  namespace VQE {
    class AdaptVQE {
      protected:
        std::shared_ptr<DynamicAnsatz> a;
        std::shared_ptr<VQEState> state;
        std::vector< std::vector<FermionOperator> > fermi_op_pool;
        std::vector<std::vector<PauliOperator> > pauli_op_pool;
        std::vector<ValType> gradient_magnitudes;
        std::vector<std::shared_ptr<Hamiltonian> > hamil;
        std::vector<std::shared_ptr<Hamiltonian> > commutators;
    public:
      AdaptVQE(std::shared_ptr<VQEState> backend, std::shared_ptr<Hamiltonian> _hamil): state(backend), hamil(_hamil) {}
      void make_commutators() {
        std::shared_ptr<Hamiltonian> hamil = state->get_hamiltonian();
        const auto& pauli_strings = hamil->getPauliOperators()[0];
        commutators.resize(fermi_op_pool.size());
        for (size_t i = 0; i < pauli_op_pool.size(); i++) {
          std::unordered_map<PauliOperator,  std::complex<double>, PauliHash> pmap;
          std::vector<PauliOperator> oplist = pauli_op_pool[i];
          for (auto p1: pauli_strings) {
            for (auto p2: oplist) {
              auto p12 = p1 * p2;
              auto p21 = (p2 * p1);
              p21 *= -1.0;
              if (pmap.find(p12) == pmap.end()) {
                pmap[p12] = p12.getCoeff();
              } else {
                pmap[p12] += p12.getCoeff();
              }
              if (pmap.find(p21) == pmap.end()) {
                pmap[p21] = p21.getCoeff();
              } else {
                pmap[p21] += p21.getCoeff();
              }
            }
          }
          std::vector<std::vector<PauliOperator> > comm_ops(1);
          comm_ops[0].reserve(pmap.size());
          for (auto pair: pmap) {
            if (abs(pair.second.real()) > 1e-10 || abs(pair.second.imag()) > 1e-10) {
              comm_ops.push_back(pair.first);
            }
          }
          commutators[i] = std::make_shared<Hamiltonian>(hamil->getEnv(), comm_ops);
        }
      }
      void optimize(std::vector<double>& parameters, ValType& ene, IdxType maxiter, ValType reltol = 1e-4) {
        ene = hamil->get_constant();
        IdxType iter = 0;
        ValType prev_ene = 1 + ene;
        make_commutators();
        while(abs((prev_ene - ene) / prev_ene) >= reltol && iter < maxiter) {
          prev_ene = ene;
          a->getFermionicOperatorParameters();
          IdxType max_ind = 0;
          double max_ene = -MAXFLOAT;
          for (size_t i = 0; i < pauli_op_pool.size(); i++) {
            a->add_operator(pauli_op_pool[i]);
            state->set_hamil(commutators[i]);
            state->set_ansatz(a);
            state->initialize();
            double ene = state->energy();
            if (ene > max_ene) {
              max_ind = i;
              max_ene = ene;
            }
            a->remove_last_operator();
          }
          a->add_operator(pauli_op_pool[i]);
          state->set_ansatz(a);
          state->set_hamil(hamil);
          parameters.push_back(0.0);
          state->optimize(parameters, ene);
          iter++;
        }
        
      }

    };
  };
};


#endif