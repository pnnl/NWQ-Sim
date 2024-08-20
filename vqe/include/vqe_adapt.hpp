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
    const std::complex<ValType> imag = {0, 1.0}; // constant `i` value
    /**
    * @brief  Implementation of ADAPT-VQE optimizer
    * @note   Uses commutator-based gradient calculation with either Fermionic or Pauli operator pools
    */
    class AdaptVQE {
      protected:
        std::shared_ptr<DynamicAnsatz> ansatz; // Iterative Ansatz being constructed
        std::shared_ptr<VQEState> state; // VQE state for energy calculation and ansatz optimization
        std::vector<ValType> gradient_magnitudes; // vector of gradient magnitudes
        std::shared_ptr<Hamiltonian> hamil; // Hamiltonian observable
        std::shared_ptr<Ansatz> gradient_measurement; // Measurement circuit for gradient estimation
        std::vector<std::vector<std::vector<double> > > commutator_coeffs; // Coefficients for commutator operators
        std::vector<std::vector<std::vector<IdxType> > > commutator_zmasks; // Zmasks for commutator operators
        std::vector<ObservableList*> gradient_observables; // Vector of structure pointers for measurement circuit
        std::vector<IdxType> observable_sizes; // Stores the number of commuting cliques for each commutator  
    public:

      //Ctor
      AdaptVQE(std::shared_ptr<DynamicAnsatz> _ans, std::shared_ptr<VQEState> backend, std::shared_ptr<Hamiltonian> _hamil): 
            ansatz(_ans), 
            state(backend), 
            hamil(_hamil) {
        gradient_measurement = std::make_shared<Ansatz>(_ans->num_qubits());
      }
      //Dtor
      ~AdaptVQE() {
        size_t index = 0;
        for (auto i: gradient_observables)
          state->delete_observables(i, observable_sizes[index++]);
      }

     /**
      * @brief  Calculate the commutator of two sums over Pauli strings
      * @note   Uses a std::unordered_map to keep track of coefficients to avoid duplicate/zero terms. Computes [oplist1, oplist2]
      * @param  oplist1: Sum over Pauli terms, first operator in commutator
      * @param  oplist2: Sum over Pauli terms, second operator in commutator
      * @param  summation: std::unordered_map to track coefficient sums
      * @retval None
      */
      void commutator(std::vector<PauliOperator>& oplist1, 
                      std::vector<PauliOperator>& oplist2, 
                      PauliMap& summation) {
        for (auto p1: oplist1) {
              for (auto p2: oplist2) {
                auto p12 = (p1 * p2);
                auto p21 = (p2 * p1);
                // multiply by an imaginary factor 
                p21 *= -1.0 * imag; // 
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
        IdxType poolsize = pauli_op_pool.size();
        commutator_coeffs.resize(poolsize);
        // gradient_observables.reserve(pauli_op_pool.size());
        commutator_zmasks.resize(poolsize);
        gradient_magnitudes.resize(poolsize);
        gradient_observables.resize(poolsize);
        observable_sizes.resize(poolsize);
        size_t num_pauli_terms_total = 0;
        for (size_t i = 0; i <poolsize; i++) {
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
          num_pauli_terms_total += comm_ops.size();
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
            PauliOperator common = make_common_op(commuting_group, 
                                                  commutator_zmasks[i][j], 
                                                  commutator_coeffs[i][j]);
            
            Measurement circ1 (common, false); // QWC measurement circuit $U_M$
            gradient_measurement->compose(circ1, qubit_mapping);         // add to gradient measurement
            // add a gate to compute the expectation values   
            state->set_exp_gate(gradient_measurement, gradient_observables[i] + j, commutator_zmasks[i][j], commutator_coeffs[i][j]);
            Measurement circ2 (common, true); // inverse of the measurement circuit $U_M^\dagger$
            gradient_measurement->compose(circ2, qubit_mapping);  // add the inverse
            cliqueiter++;  
          }
          // commutators[i] = std::make_shared<Hamiltonian>(hamil->getEnv(), comm_ops_grouped);
        }
        if (state->get_process_rank() == 0)
          std::cout << "Generated " << pauli_op_pool.size() << " commutators with " << num_pauli_terms_total << " (possibly degenerate) Individual Pauli Strings" << std::endl;
      }
      void optimize(std::vector<double>& parameters, ValType& ene, IdxType maxiter, ValType abstol = 1e-5, ValType fvaltol = 1e-7) {
        ene = hamil->getEnv().constant;
        state->initialize();
        ValType constant = ene;
        IdxType iter = 0;
        ValType prev_ene = 1 + ene;
        const auto& pauli_op_pool = ansatz->get_pauli_op_pool();
        while(iter < maxiter) {
          prev_ene = ene;
          IdxType max_ind = 0; 
          double max_ene = -MAXFLOAT;
          // std::fill(gradient_magnitudes.begin(), gradient_magnitudes.end(), 0);
          state->call_simulator(gradient_measurement);
          std::fill(gradient_magnitudes.begin(), gradient_magnitudes.end(), 0);
          state->get_exp_values(gradient_observables, observable_sizes, gradient_magnitudes);
          double grad_norm = std::sqrt(std::accumulate(gradient_magnitudes.begin(), gradient_magnitudes.end(), 0.0, [] (ValType a, ValType b) {
            return a + b * b;
          }));
          if (grad_norm < abstol) {
            break;
          }
          max_ind = std::max_element(gradient_magnitudes.begin(),
                                     gradient_magnitudes.end(),
                                     [] (ValType a, ValType b) {return abs(a) < abs(b);}) - gradient_magnitudes.begin();

          ValType paramval = 0.0;//dist(rng);
          ansatz->add_operator(max_ind, paramval);
          // state->swap_hamil(hamil);
          parameters.push_back(paramval);
          state->optimize(parameters, ene);
          if (state->get_process_rank() == 0) {
            std::cout << "ADAPT Iteration " << iter << ", Fval = " << ene << std::endl;
            std::cout << "\tSelected Operator: " << ansatz->get_operator_string(max_ind) << ", Current gradient norm = " << grad_norm << std::endl;
          }
          if (abs((ene - prev_ene)) < fvaltol) {
            break;
          }
          iter++;
        }
        
      }

    };
  };
};


#endif