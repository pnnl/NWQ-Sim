#include "observable/hamiltonian.hpp"
#include "circuit/ansatz.hpp"
#include "circuit/measurement.hpp"
#include "environment.hpp"
#include "observable/fermionic_operator.hpp"
#include "observable/pauli_operator.hpp"
#include "state.hpp"
#include "utils.hpp"
#include "vqe_state.hpp"
#include <fstream>
#include <sstream>
#include <regex>
#include <queue>
#include <list>
#include <vector>
namespace NWQSim{
  namespace VQE {
    struct FermiArgs {
      IdxType qubit;
      std::complex<ValType> coeff;
      FermionOpType typeval;
      IdxType term_index;
    };
    const std::regex pattern("\\(\\s*([\\d\\.e\\+-]+),\\s*([\\d\\.]+)\\)([\\d^\\s]+)");

    void read_fermion_operators(std::string input_path,
                                std::vector<std::vector<FermionOperator>>& fermi_operators,
                                MolecularEnvironment& env,
                                size_t n_particles,
                                bool xacc_scheme) {
       std::ifstream input_file(input_path);
      if (!input_file.is_open()) {
        throw std::runtime_error("Could not open file");
      }
      std::string line;
      std::istringstream linestream;

      std::string throwaway;
      std::complex<ValType> coeff;
      std::smatch match;
      IdxType max_index = 0;
        
      std::queue<FermiArgs> arglist;
      IdxType term_counter = 0;
      ValType energy_constant = 0;
      while (std::getline(input_file, line)) {
        if (line.length() == 0 || line[0] == 'c') {
          continue;
        }
        
        if (std::regex_search(line, match, pattern)) {
          if (match.size() != 4) {
            continue;
          }
          coeff = {std::stod(match.str(1)), std::stod(match.str(2))};
          if (match.str(3) == "") {
            energy_constant = coeff.real();
          }
          std::string fermi_ops = match.str(3);
          linestream = std::istringstream(fermi_ops);
          std::vector<std::string> term_ops;
          while(std::getline(linestream, throwaway, ' ')) {
            FermionOpType optype;
            IdxType qubit_index;
            if (throwaway == "" || throwaway == "+")  {
              continue;
            }
            if (throwaway.back() == '^') {
              optype = Creation;
              qubit_index = std::stoll(throwaway.substr(0, throwaway.length() - 1));
            } else {
              optype = Annihilation;
              qubit_index = std::stoll(throwaway);
            }
            max_index = std::max(max_index, qubit_index);

            term_ops.push_back(throwaway);
            FermiArgs args;
            args.coeff = coeff;
            args.qubit = qubit_index;
            args.term_index = term_counter;
            args.typeval = optype;
            arglist.push(args);
          }
          term_counter ++;
        } else {

        std::cout << "No match " <<  line << std::endl;
        }
      }
      if ((max_index + 1) % 2) {
        throw std::runtime_error("Odd number of spin orbitals");  
      }
        
      fermi_operators.resize(term_counter);
      IdxType n_spatial_orbitals = (max_index + 1) / 2;
      IdxType n_occ = n_particles / 2;
      IdxType n_virt = n_spatial_orbitals - n_occ;
      env = MolecularEnvironment(n_spatial_orbitals, 
                                 n_particles,
                                 xacc_scheme,
                                 energy_constant);
      IdxType coeff_index = 0;
      IdxType otypeval, spinval, orbital_index;
      Spin spin;
      OrbitalType orbitaltype;
      
      while(!arglist.empty()) {
        auto args = arglist.front();
        arglist.pop();
        
        getFermiInfoFromQubit(args.qubit, orbital_index, spinval, otypeval, n_occ, n_virt, xacc_scheme);
        spin = spinval ? Down : Up;
        orbitaltype = otypeval ? Virtual : Occupied;
        fermi_operators[args.term_index].push_back(
          FermionOperator(orbital_index, 
                          orbitaltype, 
                          spin, 
                          args.typeval, 
                          xacc_scheme,
                          args.coeff));
      };
    }


    Hamiltonian::Hamiltonian(std::string input_path, IdxType n_particles, bool xacc_scheme,
                    Transformer transform) {
      read_fermion_operators(input_path, fermi_operators, env, n_particles, xacc_scheme);
      transform(env, fermi_operators, pauli_operators, false);
      n_ops = 0;
      for (auto& i : pauli_operators) {
        n_ops += i.size();
      }
// #ifndef NDEBUG
//       std::cout << "Hamiltonian Pauli Strings:" << std::endl;
//       for (auto& plist: pauli_operators) {
//         for (auto& p: plist) {
//           std::cout << p << "\n";
//         }
//       }
// #endif
    };
    Hamiltonian::Hamiltonian(const std::vector<std::pair<std::string, std::complex<double>>>& input_ops, 
                             IdxType n_particles, 
                             bool xacc_scheme,
                             Transformer transform) {


      std::string throwaway;
      std::complex<ValType> coeff;
      IdxType max_index = 0;
      std::queue<FermiArgs> arglist;
      std::istringstream linestream;
      IdxType term_counter = 0;
      ValType energy_constant = 0;
      for (auto& pair: input_ops) {
        

          coeff = pair.second;
          if (pair.first == "") {
            energy_constant = coeff.real();
          }
          std::string fermi_ops = pair.first;
          linestream = std::istringstream(fermi_ops);
          std::vector<std::string> term_ops;
          while(std::getline(linestream, throwaway, ' ')) {
            FermionOpType optype;
            IdxType qubit_index;
            if (throwaway == "" || throwaway == "+")  {
              continue;
            }
            if (throwaway.back() == '^') {
              optype = Creation;
              qubit_index = std::stoll(throwaway.substr(0, throwaway.length() - 1));
            } else {
              optype = Annihilation;
              qubit_index = std::stoll(throwaway);
            }
            max_index = std::max(max_index, qubit_index);

            term_ops.push_back(throwaway);
            FermiArgs args;
            args.coeff = coeff;
            args.qubit = qubit_index;
            args.term_index = term_counter;
            args.typeval = optype;
            arglist.push(args);
          }
          term_counter ++;
      }
      if ((max_index + 1) % 2) {
        throw std::runtime_error("Odd number of spin orbitals");  
      }
      fermi_operators.resize(term_counter);
      IdxType n_spatial_orbitals = (max_index + 1) / 2;
      IdxType n_occ = n_particles / 2;
      IdxType n_virt = n_spatial_orbitals - n_occ;
      env = MolecularEnvironment(n_spatial_orbitals, 
                                 n_particles,
                                 xacc_scheme,
                                 energy_constant);
      IdxType coeff_index = 0;
      IdxType otypeval, spinval, orbital_index;
      Spin spin;
      OrbitalType orbitaltype;
      
      while(!arglist.empty()) {
        auto args = arglist.front();
        arglist.pop();
        
        getFermiInfoFromQubit(args.qubit, orbital_index, spinval, otypeval, n_occ, n_virt, xacc_scheme);
        spin = spinval ? Down : Up;
        orbitaltype = otypeval ? Virtual : Occupied;
        fermi_operators[args.term_index].push_back(
          FermionOperator(orbital_index, 
                          orbitaltype, 
                          spin, 
                          args.typeval, 
                          xacc_scheme,
                          args.coeff));
      };
      transform(env, fermi_operators, pauli_operators, false);
      n_ops = 0;
      for (auto& i : pauli_operators) {
        n_ops += i.size();
      }
// #ifndef NDEBUG
//       std::cout << "Hamiltonian Pauli Strings:" << std::endl;
//       for (auto& plist: pauli_operators) {
//         for (auto& p: plist) {
//           std::cout << p << "\n";
//         }
//       }
// #endif
    };
  void commutator(std::vector<PauliOperator>& oplist1, 
                        std::vector<PauliOperator>& oplist2, 
                        PauliMap& summation) {
          const std::complex<double> imag = {0, 1.0};
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
    void make_commutators(std::shared_ptr<Hamiltonian> hamil,
                          const std::vector<std::vector<PauliOperator> >& pauli_op_pool,
                          std::vector<std::vector<PauliOperator> >& commutator_list) {
        const auto& pauli_strings = hamil->getPauliOperators();
        commutator_list.reserve(pauli_op_pool.size());
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
          commutator_list.push_back(comm_ops);
        }
      }
    PauliOperator make_common_op(const std::vector<PauliOperator>& pauli_list, 
                               std::vector<IdxType>& zmasks,
                               std::vector<ValType>& coeffs) {
      zmasks.reserve(pauli_list.size());
      coeffs.reserve(pauli_list.size());
      IdxType composite_xmask = 0;
      IdxType composite_zmask = 0;
      IdxType dim = -1;
      for (const PauliOperator& pauli: pauli_list) {
        dim = std::max(dim, pauli.get_dim());
        composite_xmask |= pauli.get_xmask();
        coeffs.push_back(pauli.getCoeff().real());
        coeffs.back() *= (pauli.count_y() % 2) ? -1.0 : 1.0;
        zmasks.push_back(pauli.get_zmask() | pauli.get_xmask());
        composite_zmask |= pauli.get_zmask();
      }
      return PauliOperator(composite_xmask, composite_zmask, dim);
      
    }
    void make_group_measurement_circuits(const std::vector<std::vector<PauliOperator> >& commutator_list,
                                           std::shared_ptr<Ansatz> measurement_circuit,
                                           std::vector<std::vector<IdxType> >& commutator_zmasks,
                                           std::vector<std::vector<ValType> >& commutator_coeffs,
                                           std::vector<std::vector<ObservableList*> >& gradient_observables,
                                           VQEState* state) {
        IdxType num_qubits = commutator_list.front().front().get_dim();
        gradient_observables.resize(commutator_list.size());
        commutator_coeffs.resize(commutator_list.size());
        commutator_zmasks.resize(commutator_list.size());
        // Create commuting groups using the (nonoverlapping) Sorted Insertion heuristic (see Crawford et. al 2021)
        for (size_t i = 0; i < commutator_list.size(); i ++) {
          std::list<std::vector<IdxType>> cliques;
          auto comm_ops = commutator_list[i];
          sorted_insertion(comm_ops, cliques, false);
          auto cliqueiter = cliques.begin();
          gradient_observables[i].resize(cliques.size());
          std::vector<IdxType> qubit_mapping (num_qubits);
          std::iota(qubit_mapping.begin(), qubit_mapping.end(), 0);
          // For each clique, construct a measurement circuit and append
          for (size_t j = 0; j < cliques.size(); j++) {
            std::vector<IdxType>& clique = *cliqueiter;
            std::vector<PauliOperator> commuting_group (clique.size());
            std::transform(clique.begin(), clique.end(),
              commuting_group.begin(), [&] (IdxType ind) {return comm_ops.at(ind);});
            // NOTE: IN PROCESS OF API UPDATE!!!! PPOD!!!!!
            PauliOperator common = make_common_op(commuting_group, 
                                                  commutator_zmasks[i], 
                                                  commutator_coeffs[i]);
            
            Measurement circ1 (common, false);
            measurement_circuit->compose(circ1, qubit_mapping);            
            state->set_exp_gate(measurement_circuit, gradient_observables[i][j], commutator_zmasks[i], commutator_coeffs[i]);
            Measurement circ2 (common, true);
            measurement_circuit->compose(circ2, qubit_mapping);      
            cliqueiter++;  
          } 
        }
         
          // commutators[i] = std::make_shared<Hamiltonian>(hamil->getEnv(), comm_ops_grouped);
      }
  }; // namespace VQE
}; // namespace NWQSim