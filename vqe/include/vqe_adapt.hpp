#ifndef VQE_ADAPT_HH
#define VQE_ADAPT_HH
#include "state.hpp"
#include "vqe_state.hpp"
#include "environment.hpp"
#include "utils.hpp"
#include "circuit/dynamic_ansatz.hpp"
#include <cstddef>
namespace NWQSim {
  namespace VQE {
    class AdaptVQE {
      protected:
        std::shared_ptr<DynamicAnsatz> ansatz;
        std::shared_ptr<VQEState> state;
        std::vector<ValType> gradient_magnitudes;
        std::shared_ptr<Hamiltonian> hamil;
        std::shared_ptr<Ansatz> gradient_measurement;
        const std::complex<ValType> imag = {0, 1.0};
        std::vector<std::vector<std::vector<double> > > commutator_coeffs;
        std::vector<std::vector<std::vector<IdxType> > > commutator_zmasks;
        std::vector<ObservableList*> gradient_observables;
        std::vector<IdxType> observable_sizes;
    public:
      AdaptVQE(std::shared_ptr<DynamicAnsatz> _ans, std::shared_ptr<VQEState> backend, std::shared_ptr<Hamiltonian> _hamil): 
            ansatz(_ans), 
            state(backend), 
            hamil(_hamil) {
        gradient_measurement = std::make_shared<Ansatz>(_ans->num_qubits());
      }
      ~AdaptVQE() {
        size_t index = 0;
        for (auto i: gradient_observables)
          state->delete_observables(i, observable_sizes[index++]);
      }
      void commutator(std::vector<PauliOperator>& oplist1, 
                      std::vector<PauliOperator>& oplist2, 
                      PauliMap& summation) {
        for (auto p1: oplist1) {
              for (auto p2: oplist2) {
                auto p12 = (p1 * p2);
                auto p21 = (p2 * p1);
                p21 *= -1.0 * imag;
                p12 *= imag;
                if (summation.find(p12) == summation.end()) {
                  summation[p12] = p12.getCoeff();
                } else {
                  summation[p12] += p12.getCoeff();
                }
                if (summation.find(p21) == summation.end()) {
                  summation[p21] = p21.getCoeff();
                } else {
                  summation[p21] += p21.getCoeff();
                }
              }
            }
      }
      void make_commutators() {
        std::shared_ptr<Hamiltonian> hamil = state->get_hamiltonian();
        const auto& pauli_strings = hamil->getPauliOperators();

        const auto& pauli_op_pool = ansatz->get_pauli_op_pool();

        commutator_coeffs.resize(pauli_op_pool.size());
        // gradient_observables.reserve(pauli_op_pool.size());
        commutator_zmasks.resize(pauli_op_pool.size());
        gradient_magnitudes.resize(pauli_op_pool.size());
        gradient_observables.resize(pauli_op_pool.size());
        observable_sizes.resize(pauli_op_pool.size());
        for (size_t i = 0; i < pauli_op_pool.size(); i++) {
          std::unordered_map<PauliOperator,  std::complex<double>, PauliHash> pmap;
          std::vector<PauliOperator> oplist = pauli_op_pool[i];
          for (auto hamil_oplist: pauli_strings) {
            commutator(hamil_oplist, oplist, pmap);
          }
          std::vector<PauliOperator> comm_ops;
          comm_ops.reserve(pmap.size());
          for (auto pair: pmap) {
            if (abs(pair.second.real()) > 1e-10 || abs(pair.second.imag()) > 1e-10) {
              PauliOperator op(pair.first);
              op.setCoeff(pair.second);
              comm_ops.push_back(op);
            }
          }
          // Create commuting groups using the (nonoverlapping) Sorted Insertion heuristic (see Crawford et. al 2021)
          std::list<std::vector<IdxType>> cliques;
          sorted_insertion(comm_ops, cliques, false);
          auto cliqueiter = cliques.begin();
          commutator_coeffs[i].resize(cliques.size());
          commutator_zmasks[i].resize(cliques.size());
          state->allocate_observables(gradient_observables[i], cliques.size());
          std::vector<IdxType> qubit_mapping (ansatz->num_qubits());
          std::iota(qubit_mapping.begin(), qubit_mapping.end(), 0);
          observable_sizes[i] = cliques.size();
          // For each clique, construct a measurement circuit and append
          for (size_t j = 0; j < cliques.size(); j++) {
            std::vector<IdxType>& clique = *cliqueiter;
            std::vector<PauliOperator> commuting_group (clique.size());
            std::transform(clique.begin(), clique.end(),
              commuting_group.begin(), [&] (IdxType ind) {return comm_ops[ind];});
            // NOTE: IN PROCESS OF API UPDATE!!!! PPOD!!!!!
 
            PauliOperator common = make_common_op(commuting_group, 
                                                  commutator_zmasks[i][j], 
                                                  commutator_coeffs[i][j]);
            
            Measurement circ1 (common, false);
            gradient_measurement->compose(circ1, qubit_mapping);            
            state->set_exp_gate(gradient_measurement, &gradient_observables[i][j], commutator_zmasks[i][j], commutator_coeffs[i][j]);
            Measurement circ2 (common, true);
            gradient_measurement->compose(circ2, qubit_mapping);      
            cliqueiter++;  
          }
          // commutators[i] = std::make_shared<Hamiltonian>(hamil->getEnv(), comm_ops_grouped);
        }
      }
      void optimize(std::vector<double>& parameters, ValType& ene, IdxType maxiter, ValType reltol = 1e-5, ValType reltol_fval = 1e-7) {
        ene = hamil->getEnv().constant;
        state->initialize();
        ValType constant = ene;
        IdxType iter = 0;
        ValType prev_ene = 1 + ene;
        const auto& pauli_op_pool = ansatz->get_pauli_op_pool();
        maxiter = 50;
        std::uniform_real_distribution<ValType> dist(0, 2 * PI);
        std::mt19937_64 rng (2423);
        while(iter < maxiter) {
          prev_ene = ene;
          IdxType max_ind = 0; 
          double max_ene = -MAXFLOAT;
          // std::fill(gradient_magnitudes.begin(), gradient_magnitudes.end(), 0);
          state->call_simulator(gradient_measurement);
          std::fill(gradient_magnitudes.begin(), gradient_magnitudes.end(), 0);
          state->get_exp_values(gradient_observables, observable_sizes, gradient_magnitudes);
          double grad_norm = std::accumulate(gradient_magnitudes.begin(), gradient_magnitudes.end(), 0.0, [] (ValType a, ValType b) {
            return a + b * b;
          });
          max_ind = std::max_element(gradient_magnitudes.begin(),
                                     gradient_magnitudes.end(),
                                     [] (ValType a, ValType b) {return abs(a) < abs(b);}) - gradient_magnitudes.begin();
          for (auto i: gradient_magnitudes) {
            std::cout << i << " ";
          }
          std::cout << std::endl;
          // std::cout << gradient_magnitudes << std::endl;
          if (std::sqrt(grad_norm) < reltol) {
            break;
          }
          ValType paramval = 0.0;//dist(rng);
          ansatz->add_operator(max_ind, paramval);
          // state->swap_hamil(hamil);
          parameters.push_back(paramval);
          state->optimize(parameters, ene);
          ValType denom = abs(prev_ene) > 0 ? prev_ene : 1.0;
          if (abs((ene - prev_ene) / denom) < reltol_fval) {
            break;
          }
          iter++;
        }
        
      }

    };
  };
};


#endif