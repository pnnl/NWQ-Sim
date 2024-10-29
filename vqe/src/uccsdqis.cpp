#include "circuit/ansatz.hpp"
#include "circuit/dynamic_ansatz.hpp"
#include "observable/fermionic_operator.hpp"
#include "utils.hpp"


#include <set>
#include <tuple>
// #include <map>




namespace NWQSim {
  namespace VQE {

    class UCCSDQis: public Ansatz {
      protected:
        const MolecularEnvironment& env;
        IdxType n_singles;
        IdxType n_doubles;
        IdxType trotter_n;
        IdxType unique_params;
        IdxType symm_level;
        Transformer qubit_transform;
        /** 
         * Enforce symmetries for each term. Each fermionic term will have one symmetry entry. If no symmetries are enforced, 
         * symmetries[i] = {{i, 1.0}}; Otherwise, symmetries[i] = {{j, 1.0}, {k, -1.0}} denotes that theta_i must be equal to theta_j - theta_k
         */
        std::vector<std::vector<std::pair<IdxType, ValType> > > symmetries;
        std::vector<IdxType> fermion_ops_to_params; // map from fermion operators to parameters (used in update)
        std::vector<std::vector<FermionOperator> > fermion_operators;
        std::set<std::tuple< IdxType, IdxType, IdxType, IdxType >> existing_tuples; // MZ: for recording symmetry
        // std::map<std::tuple<IdxType, IdxType, IdxType, IdxType>, std::vector<std::pair<IdxType, double>>> existing_terms; // MZ: for recording symmetry

      public:
        UCCSDQis(const MolecularEnvironment& _env, Transformer _qubit_transform, IdxType _trotter_n = 1, IdxType _symm_level = 3): 
                                  env(_env),
                                  trotter_n(_trotter_n),
                                  symm_level(_symm_level),
                                  qubit_transform(_qubit_transform),
                                  Ansatz(2 * _env.n_spatial) {
          // ----- MZ added ------
          // Single exictations: 1 alpha (beta) electron to 1 virtual alpha (beta) orbital
          // So # alpha-alpha single is n_occ * n_virt, same thing for beta, so total number is 2 * n_occ * n_virt
          //
          // Double exictations have 3 cases: alpha-alpha, beta-beta, alpha-beta
          // => # alpha-alpha double is nChooseK(n_virt, 2)*nChooseK(n_occ, 2), same thing for beta-beta
          // => # alpha-beta double is like (# alpha-alpha single)*(# beta-beta single), i.e.,  (n_occ * n_virt)^2
          // So total number of double excitations is 2*(nChooseK(n_virt, 2)*nChooseK(n_occ, 2)) + (n_occ * n_virt)^2
          // ---------------------
          n_singles = 2 * env.n_occ * env.n_virt;
          IdxType c2virtual = choose2(env.n_virt);
          IdxType c2occupied = choose2(env.n_occ);
          n_doubles = 2*c2virtual*c2occupied + (env.n_occ) * (env.n_virt) * (env.n_occ) * (env.n_virt); // MZ: this is exact, no need for 10*

          fermion_operators.reserve(n_singles + n_doubles);
          symmetries = std::vector<std::vector<std::pair<IdxType, ValType> > >((n_singles + n_doubles));
          fermion_ops_to_params.resize(n_doubles + n_singles);
          std::fill(fermion_ops_to_params.begin(), fermion_ops_to_params.end(), -1);
          unique_params = 0;
          existing_tuples = {};
          // existing_terms = {};
        };
 
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
          }
          return result;
        };
        
        const MolecularEnvironment& getEnv() const {return env;};
        virtual IdxType numParams() const override { return unique_params; };
        virtual IdxType numOps() const override { return fermion_operators.size(); };

    void add_double_excitation(FermionOperator i, FermionOperator j, FermionOperator r, FermionOperator s,  const std::vector<std::pair<IdxType, double>>& symm_expr, bool param) {

        // use this index as the unique parameter to create the symmetry
        symmetries[fermion_operators.size()] = symm_expr;
        fermion_operators.push_back({i, j, r, s});
        // record the string::parameter mapping
        if (param) {
          fermion_ops_to_params[fermion_operators.size()-1] = unique_params++;
          excitation_index_map[to_fermionic_string(fermion_operators.back(), env)] = unique_params-1;
        }
    }
    void add_double_excitation(FermionOperator i, FermionOperator j, FermionOperator r, FermionOperator s) {

        // use this index as the unique parameter to create the symmetry
        symmetries[fermion_operators.size()] = {{fermion_operators.size(), 1.0}};
        fermion_operators.push_back({i, j, r, s});
        // record the string::parameter mapping
        fermion_ops_to_params[fermion_operators.size()-1] = unique_params++;
        excitation_index_map[to_fermionic_string(fermion_operators.back(), env)] = unique_params-1;
    }

   /**
    * @brief  Generate Fermionic operators for UCCSD
    * @note   Symmetry-linked operators (e.g. by anticommutation, spin reversal) share parameters
    * @retval None
    */
    void  getFermionOps() {
      // Alpha Single Excitations
      for (IdxType p = 0; p < env.n_occ; p++) {
        FermionOperator occupied_annihilation_up (p, Occupied, Up, Annihilation, env.xacc_scheme);
        for (IdxType q = 0; q < env.n_virt; q++) {
          // creation operator
          FermionOperator virtual_creation_up (q, Virtual, Up, Creation, env.xacc_scheme);
          // use this index as the unique parameter to create the symmetry
          symmetries[fermion_operators.size()] = {{fermion_operators.size(), 1.0}};
          fermion_ops_to_params[fermion_operators.size()] = unique_params++;
          fermion_operators.push_back({occupied_annihilation_up, virtual_creation_up});
          // record the string::parameter mapping
          excitation_index_map[to_fermionic_string(fermion_operators.back(), env)] = unique_params-1;
        }
      }
      // Beta Single Excitations
      // In a closed shell system, the beta single excitations has the same parameters as the alpha single excitations
      for (IdxType p = 0; p < env.n_occ; p++) {
        FermionOperator occupied_annihilation_down (p, Occupied, Down, Annihilation, env.xacc_scheme);
        for (IdxType q = 0; q < env.n_virt; q++) {
          FermionOperator virtual_creation_down (q, Virtual, Down, Creation, env.xacc_scheme);
          if (symm_level >= 1) {
            // Add a pointer to the corresponding Up spin term to share parameters
            symmetries[fermion_operators.size()] = {{fermion_operators.size() - env.n_occ * env.n_virt, 1.0}};
            fermion_operators.push_back({occupied_annihilation_down, virtual_creation_down});
            excitation_index_map[to_fermionic_string(fermion_operators.back(), env)] = fermion_operators.size() - env.n_occ * env.n_virt - 1;
          } else {
            symmetries[fermion_operators.size()] = {{fermion_operators.size(), 1.0}};
            fermion_ops_to_params[fermion_operators.size()] = unique_params++;
            fermion_operators.push_back({occupied_annihilation_down, virtual_creation_down});
            excitation_index_map[to_fermionic_string(fermion_operators.back(), env)] = unique_params-1;
            // record the string::parameter mapping
          }
        }
      }
    /*===========Double Excitations===========*/
    // MZ: re-write the original code, my python head cannot get it
    // alpha-alpha and beta-beta

    for (IdxType i = 0; i < env.n_occ; i++) {
        FermionOperator i_occ_ann_up (i, Occupied, Up, Annihilation, env.xacc_scheme);
        FermionOperator i_occ_ann_dw (i, Occupied, Down, Annihilation, env.xacc_scheme);
        for (IdxType r = 0; r < env.n_virt; r++) {
            FermionOperator r_virt_cre_up (r, Virtual, Up, Creation, env.xacc_scheme);
            FermionOperator r_virt_cre_dw (r, Virtual, Down, Creation, env.xacc_scheme);
            for (IdxType j = i+1; j < env.n_occ; j++) {
                FermionOperator j_occ_ann_up (j, Occupied, Up, Annihilation, env.xacc_scheme); 
                FermionOperator j_occ_ann_dw (j, Occupied, Down, Annihilation, env.xacc_scheme);
                for (IdxType s = r+1; s < env.n_virt; s++) {
                    FermionOperator s_virt_cre_dw (s, Virtual, Down, Creation, env.xacc_scheme);
                    FermionOperator s_virt_cre_up (s, Virtual, Up, Creation, env.xacc_scheme);
                    if (symm_level >= 1) {
                        IdxType term1 = fermion_operators.size();
                        add_double_excitation(i_occ_ann_up, j_occ_ann_up, r_virt_cre_up, s_virt_cre_up, {{term1, 1.0}}, true); // MZ: alpha-alpha
                        add_double_excitation(i_occ_ann_dw, j_occ_ann_dw, r_virt_cre_dw, s_virt_cre_dw, {{term1, 1.0}}, false); // MZ: beta-beta
                    } else {
                        add_double_excitation(i_occ_ann_up, j_occ_ann_up, r_virt_cre_up, s_virt_cre_up); // MZ: alpha-alpha
                        add_double_excitation(i_occ_ann_dw, j_occ_ann_dw, r_virt_cre_dw, s_virt_cre_dw); // MZ: beta-beta
                    }
                }
            }
        }
    }

    // // mixed
    // For some reason, accuracy is even worse at Symm 3
    // for (IdxType i = 0; i < env.n_occ; i++) {
    //     FermionOperator i_occ_ann_up (i, Occupied, Up, Annihilation, env.xacc_scheme);
    //     FermionOperator i_occ_ann_dw (i, Occupied, Down, Annihilation, env.xacc_scheme);
    //     for (IdxType r = 0; r < env.n_virt; r++) {
    //         FermionOperator r_virt_cre_up (r, Virtual, Up, Creation, env.xacc_scheme);
    //         FermionOperator r_virt_cre_dw (r, Virtual, Down, Creation, env.xacc_scheme);
    //         for (IdxType j = 0; j < env.n_occ; j++) {
    //             FermionOperator j_occ_ann_dw (j, Occupied, Down, Annihilation, env.xacc_scheme);
    //             FermionOperator j_occ_ann_up (j, Occupied, Up, Annihilation, env.xacc_scheme);
    //             for (IdxType s = 0; s < env.n_virt; s++) {
    //                 FermionOperator s_virt_cre_dw (s, Virtual, Down, Creation, env.xacc_scheme);
    //                 FermionOperator s_virt_cre_up (s, Virtual, Up, Creation, env.xacc_scheme);
    //                 if ( (symm_level < 2) || (i == j && r == s) ) {
    //                   add_double_excitation(i_occ_ann_up, j_occ_ann_dw, r_virt_cre_up, s_virt_cre_dw); // MZ: alpha-beta, not need for beta-alpha
    //                 } else {
    //                   std::tuple<IdxType, IdxType, IdxType, IdxType> new_tuple = {i,j,r,s};
    //                   if (existing_tuples.find(new_tuple) != existing_tuples.end()) {
    //                     // The tuple exist in the set
    //                     add_double_excitation(j_occ_ann_up, i_occ_ann_dw, s_virt_cre_up, r_virt_cre_dw, existing_terms.at(new_tuple), false);
    //                     existing_tuples.erase(new_tuple);
    //                     existing_terms.erase(new_tuple);
    //                   } else {
    //                     // The tuple does not exist in the set
    //                     IdxType term2 = fermion_operators.size();
    //                     std::vector<std::pair<IdxType, double>> sym_term= {{term2, 1.0}};
    //                     add_double_excitation(i_occ_ann_up, j_occ_ann_dw, r_virt_cre_up, s_virt_cre_dw, sym_term, true); // MZ: alpha-beta, not need for beta-alpha
    //                     existing_tuples.insert({j, i, s, r});
    //                     existing_terms.insert({ {j, i, s, r}, sym_term });
    //                   }
    //                 }

    //             }
    //         }
    //     }
    // }


    // It seems like the trottertor error is larger (extra 1e-3 hartree) at Symm 3
    for (IdxType i = 0; i < env.n_occ; i++) {
        FermionOperator i_occ_ann_up (i, Occupied, Up, Annihilation, env.xacc_scheme);
        FermionOperator i_occ_ann_dw (i, Occupied, Down, Annihilation, env.xacc_scheme);
        for (IdxType r = 0; r < env.n_virt; r++) {
            FermionOperator r_virt_cre_up (r, Virtual, Up, Creation, env.xacc_scheme);
            FermionOperator r_virt_cre_dw (r, Virtual, Down, Creation, env.xacc_scheme);
            for (IdxType j = 0; j < env.n_occ; j++) {
                FermionOperator j_occ_ann_dw (j, Occupied, Down, Annihilation, env.xacc_scheme);
                FermionOperator j_occ_ann_up (j, Occupied, Up, Annihilation, env.xacc_scheme);
                for (IdxType s = 0; s < env.n_virt; s++) {
                    FermionOperator s_virt_cre_dw (s, Virtual, Down, Creation, env.xacc_scheme);
                    FermionOperator s_virt_cre_up (s, Virtual, Up, Creation, env.xacc_scheme);
                    if ( (symm_level < 2) || (i == j && r == s) ) {
                      add_double_excitation(i_occ_ann_up, j_occ_ann_dw, r_virt_cre_up, s_virt_cre_dw); // MZ: alpha-beta, not need for beta-alpha
                    } else {
                      std::tuple<IdxType, IdxType, IdxType, IdxType> new_tuple = {i,j,r,s};
                      if (existing_tuples.find(new_tuple) != existing_tuples.end()) {
                        // The tuple exist in the set, so we skip the term and erase the tuple
                        existing_tuples.erase(new_tuple);
                      } else {
                        // The tuple does not exist in the set, so we add the term
                        IdxType term2 = fermion_operators.size();
                        add_double_excitation(i_occ_ann_up, j_occ_ann_dw, r_virt_cre_up, s_virt_cre_dw, {{term2, 1.0}}, true); // MZ: alpha-beta, not need for beta-alpha
                        add_double_excitation(j_occ_ann_up, i_occ_ann_dw, s_virt_cre_up, r_virt_cre_dw, {{term2, 1.0}}, false);
                        existing_tuples.insert({j, i, s, r});
                      }
                    }

                }
            }
        }
    }
      
    #ifndef NDEBUG
      // Print out the generated Fermionic operators
      for (size_t i = 0; i < fermion_operators.size(); i++) {
        std::string fermi_op = "";
        for (auto i: fermion_operators[i]) {
          fermi_op += i.toString(env.n_occ, env.n_virt) + " ";
        }
        std::cout << i << ": " << fermi_op << " || [";
        for (auto sym: symmetries[i]) {
          std::cout << "{" <<sym.first << ", " << fermion_ops_to_params[sym.first] << ", " << sym.second << "}" << ", ";
        }
        std::cout << "]  " << fermion_ops_to_params[i] << std::endl;
        // if (symmetries[i].size())
        // std::cout  << std::endl;
      }
      std::cout << fermion_operators.size() << " " << unique_params << std::endl;
    #endif

    } // getFermionOps


    void  buildAnsatz()  override {
      getFermionOps();
      // assert((n_doubles + n_singles) == fermion_operators.size());
      std::cout << n_singles << " " << n_doubles << std::endl;
      std::cout << "Generated " << fermion_operators.size() << " operators." << std::endl;
      theta->resize(unique_params * trotter_n);
      // exit(0);
      std::vector<std::vector<PauliOperator> > pauli_oplist;
      if (env.xacc_scheme) {
        for (IdxType i = 0; i < env.n_occ; i++) {
          X(i);
          X(i+env.n_spatial);
        }
      } else {
        for (IdxType i = 0; i < 2 * env.n_occ; i++) {
          X(i);
        }
      }
      pauli_oplist.reserve(4 * n_singles + 16 * n_doubles);
      qubit_transform(env, fermion_operators, pauli_oplist, true);  
      IdxType index = 0; // parameter index, shares parameters for Pauli evolution gates corresponding to the same Fermionic operator within the same Trotter step
      for (auto& fermionic_group: pauli_oplist) {
        bool wasused = 0;
        for (auto& pauli: fermionic_group) {
          assert (pauli.getCoeff().imag() == 0.0);
          double coeff = pauli.getCoeff().real(); 
          
          std::vector<std::pair<IdxType, ValType> > idxvals(symmetries[index].size());

          std::transform(symmetries[index].begin(), symmetries[index].end(),idxvals.begin(), 
            [&] (std::pair<IdxType, ValType> val) {
              return std::make_pair(fermion_ops_to_params[val.first], val.second);
            } );

          if (pauli.isNonTrivial() && abs(coeff) > 1e-10) { 
            ExponentialGate(pauli, OP::RZ, idxvals, 2 * coeff);
            wasused = 1;
          }
        }
        if(!wasused) {
          printf("%lld parameter not used for operator\n", index);
        }
        index++;
      }
      for (IdxType i = 0; i < trotter_n - 1; i++) {
        for (auto& fermionic_group: pauli_oplist) {
          for (auto& pauli: fermionic_group) {
            double coeff = pauli.getCoeff().real();
            if (pauli.isNonTrivial() && abs(coeff) > 1e-10)  {  
            std::vector<std::pair<IdxType, ValType> > idxvals = symmetries[index];
            std::transform(symmetries[index].begin(), symmetries[index].end(),idxvals.begin(), 
            [&, i] (std::pair<IdxType, ValType> val) {
              return std::make_pair(fermion_ops_to_params[val.first] + (i + 1) * unique_params, val.second);
            } );
            ExponentialGate(pauli, OP::RZ, idxvals, 2 * coeff);
            }
          }
          index++;
        }
      }
    }


    }; // class UCCSDQis
  
  
  };// namespace vqe
};// namespace nwqsim