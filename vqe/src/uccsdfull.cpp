#include "circuit/ansatz.hpp"
#include "circuit/dynamic_ansatz.hpp"

namespace NWQSim {
  namespace VQE {
    class UCCSDFull: public Ansatz {
      protected:
        const MolecularEnvironment& env;
        IdxType n_singles;
        IdxType n_doubles;
        IdxType trotter_n;
        IdxType unique_params;
        bool symm_enforce;
        Transformer qubit_transform;
        /** 
         * Enforce symmetries for each term. Each fermionic term will have one symmetry entry. If no symmetries are enforced, 
         * symmetries[i] = {{i, 1.0}}; Otherwise, symmetries[i] = {{j, 1.0}, {k, -1.0}} denotes that theta_i must be equal to theta_j - theta_k
         */
        std::vector<std::vector<std::pair<IdxType, ValType> > > symmetries;
        std::vector<IdxType> fermion_ops_to_params; // map from fermion operators to parameters (used in update)
        std::vector<std::vector<FermionOperator> > fermion_operators;
      public:
        UCCSDFull(const MolecularEnvironment& _env, Transformer _qubit_transform, IdxType _trotter_n = 1, bool _symm_enforce = true): 
                                  env(_env),
                                  trotter_n(_trotter_n),
                                  symm_enforce(_symm_enforce),
                                  qubit_transform(_qubit_transform),
                                  Ansatz(2 * _env.n_spatial) {
          n_singles = 2 * env.n_occ * env.n_virt;
          IdxType c2virtual = choose2(env.n_virt);
          IdxType c2occupied = choose2(env.n_occ);
          n_doubles = 10 * (env.n_occ) * (env.n_virt) * (env.n_occ) * (env.n_virt);
          fermion_operators.reserve(n_singles + n_doubles);
          symmetries = std::vector<std::vector<std::pair<IdxType, ValType> > >((n_singles + n_doubles));
          fermion_ops_to_params.resize(n_doubles + n_singles);
          std::fill(fermion_ops_to_params.begin(), fermion_ops_to_params.end(), -1);
          unique_params = 0;
          
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
            // if (is_distinct && oplist.size() == 4 && is_mixed) {
            //   first = true;
            //   opstring = "";
            //   for (auto& op: oplist) {
            //     if (!first) {
            //       opstring = " " + opstring;
            //     } else {
            //       first = false;
            //     }
            //     opstring = op.spinReversed().toString(env.n_occ, env.n_virt) + opstring;
            //   }
            //   result.push_back(std::make_pair(opstring, param));
            // }
          }
          return result;
        };
        
        const MolecularEnvironment& getEnv() const {return env;};
        virtual IdxType numParams() const override { return unique_params; };
    

    /**
      * @brief  Generate Fermionic operators for UCCSD
      * @note   Symmetry-linked operators (e.g. by anticommutation, spin reversal) share parameters
      * @retval None
      */
      void getFermionOps()
      {
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
        for (IdxType p = 0; p < env.n_occ; p++) {
          FermionOperator occupied_annihilation_down (p, Occupied, Down, Annihilation, env.xacc_scheme);
          for (IdxType q = 0; q < env.n_virt; q++) {
            FermionOperator virtual_creation_down (q, Virtual, Down, Creation, env.xacc_scheme);
            if (symm_enforce) {
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
        for (IdxType i = 0; i < env.n_occ; i++) {
          FermionOperator occupied_annihilation_up_1 (i, Occupied, Up, Annihilation, env.xacc_scheme);
          FermionOperator occupied_annihilation_down_1 (i, Occupied, Down, Annihilation, env.xacc_scheme);
          for (IdxType j = 0; j < env.n_occ; j++) {
            FermionOperator occupied_annihilation_up_2 (j, Occupied, Up, Annihilation, env.xacc_scheme);
            FermionOperator occupied_annihilation_down_2 (j, Occupied, Down, Annihilation, env.xacc_scheme);
            for (IdxType r = 0; r < env.n_virt; r++) {
              FermionOperator occupied_virtual_up_1 (r, Virtual, Up, Creation, env.xacc_scheme);
              FermionOperator occupied_virtual_down_1 (r, Virtual, Down, Creation, env.xacc_scheme);
              
              for (IdxType s = 0; s < env.n_virt; s++) {
                FermionOperator occupied_virtual_up_2 (s, Virtual, Up, Creation, env.xacc_scheme);
                FermionOperator occupied_virtual_down_2 (s, Virtual, Down, Creation, env.xacc_scheme);
                if (r != s && i != j) {
                  symmetries[fermion_operators.size()] = {{fermion_operators.size(), 1.0}};
                  fermion_ops_to_params[fermion_operators.size()] = unique_params++;
                  fermion_operators.push_back({occupied_annihilation_up_1,
                                              occupied_annihilation_up_2,
                                              occupied_virtual_up_1,
                                              occupied_virtual_up_2});
                  symmetries[fermion_operators.size()] = {{fermion_operators.size(), 1.0}};
                  fermion_ops_to_params[fermion_operators.size()] = unique_params++;
                  excitation_index_map[to_fermionic_string(fermion_operators.back(), env)] = unique_params-1;
                  fermion_operators.push_back({occupied_annihilation_down_1,
                                              occupied_annihilation_down_2,
                                              occupied_virtual_down_1,
                                              occupied_virtual_down_2});
                }

                symmetries[fermion_operators.size()] = {{fermion_operators.size(), 1.0}};
                fermion_ops_to_params[fermion_operators.size()] = unique_params++;
                excitation_index_map[to_fermionic_string(fermion_operators.back(), env)] = unique_params-1;
                fermion_operators.push_back({occupied_annihilation_up_1,
                                            occupied_annihilation_down_2,
                                            occupied_virtual_down_1,
                                            occupied_virtual_up_2});

                symmetries[fermion_operators.size()] = {{fermion_operators.size(), 1.0}};
                fermion_ops_to_params[fermion_operators.size()] = unique_params++;
                excitation_index_map[to_fermionic_string(fermion_operators.back(), env)] = unique_params-1;
                fermion_operators.push_back({occupied_annihilation_down_1,
                                            occupied_annihilation_up_2,
                                            occupied_virtual_up_1,
                                            occupied_virtual_down_2});
                

              }
            }
          }
          }
        /*===========Double Excitations===========*/
        // Mixed-spin symmetries, all distinct spatial orbitals
        // for (IdxType i = 0; i < env.n_occ; i++) {
        //   FermionOperator occ_down_1 (i, Occupied, Down, Annihilation, env.xacc_scheme);
        //   FermionOperator occ_up_1 (i, Occupied, Up, Annihilation, env.xacc_scheme);
        //   // Strictly greater than i
        //   for (IdxType j = i+1; j < env.n_occ; j++) {
        //     if (i == j) {
        //       continue;
        //     }
        //     FermionOperator occ_down_2 (j, Occupied, Down, Annihilation, env.xacc_scheme);
        //     FermionOperator occ_up_2 (j, Occupied, Up, Annihilation, env.xacc_scheme);
        //     for (IdxType r = 0; r < env.n_virt; r++) {
        //       FermionOperator virt_down_3 (r, Virtual, Down, Creation, env.xacc_scheme);
        //       FermionOperator virt_up_3 (r, Virtual, Up, Creation, env.xacc_scheme);
        //       // Strictly greater than r 
        //       for (IdxType s = r+1; s < env.n_virt; s++) {
        //         if (r == s) {
        //           continue;
        //         }
        //         FermionOperator virt_down_4 (s, Virtual, Down, Creation, env.xacc_scheme);
        //         FermionOperator virt_up_4 (s, Virtual, Up, Creation, env.xacc_scheme);
        //         IdxType alpha_term = fermion_operators.size();
        //         // All alpha excitation: alpha_4*alpha_3*alpha_2*alpha_1
        //         fermion_operators.push_back({
        //               occ_down_1,
        //               occ_down_2,
        //               virt_down_3,
        //               virt_down_4});
        //         IdxType beta_term = fermion_operators.size();
        //         // All beta excitation: beta_4*beta_3*beta_2*beta_1
        //         fermion_operators.push_back({
        //               occ_up_1, 
        //               occ_up_2,
        //               virt_up_3,
        //               virt_up_4});
        //         IdxType mixed_term1 = fermion_operators.size();
        //         // Mixed term 1: beta_4*alpha_3*alpha_2*beta_1
        //         fermion_operators.push_back({
        //               occ_up_1,
        //               occ_down_2,
        //               virt_down_3,
        //               virt_up_4});
        //         IdxType mixed_term2 = fermion_operators.size();
        //         // Mixed term 2: alpha_4*beta_3*beta_2*alpha_1
        //         fermion_operators.push_back({
        //               occ_down_1,
        //               occ_up_2,
        //               virt_up_3,
        //               virt_down_4});
        //         IdxType mixed_term3 = fermion_operators.size();
        //         // Mixed term 3: beta_4*alpha_3*alpha_1*beta_2
        //         fermion_operators.push_back({
        //               occ_up_2,
        //               occ_down_1,
        //               virt_down_3,
        //               virt_up_4});
        //         IdxType mixed_term4 = fermion_operators.size();
        //         // Mixed term 4: alpha_4*beta_3*beta_1*alpha_2
        //         fermion_operators.push_back({
        //               occ_down_2,
        //               occ_up_1,
        //               virt_up_3,
        //               virt_down_4});
        //         IdxType mixed_term5 = fermion_operators.size();
        //         // Mixed term 3: beta_4*alpha_3*alpha_1*beta_2
        //         fermion_operators.push_back({
        //               occ_up_2,
        //               occ_down_1,
        //               virt_down_4,
        //               virt_up_3});
        //         IdxType mixed_term6 = fermion_operators.size();
        //         // Mixed term 4: alpha_4*beta_3*beta_1*alpha_2
        //         fermion_operators.push_back({
        //               occ_down_2,
        //               occ_up_1,
        //               virt_up_4,
        //               virt_down_3});
        //         IdxType mixed_term7 = fermion_operators.size();
        //         // Mixed term 3: beta_4*alpha_3*alpha_1*beta_2
        //         fermion_operators.push_back({
        //               occ_up_1,
        //               occ_down_2,
        //               virt_down_4,
        //               virt_up_3});
        //         IdxType mixed_term8 = fermion_operators.size();
        //         // Mixed term 4: alpha_4*beta_3*beta_1*alpha_2
        //         fermion_operators.push_back({
        //               occ_down_1,
        //               occ_up_2,
        //               virt_up_4,
        //               virt_down_3});
              
        //         // Add the parameter pointers, all 4 values determined by two parameters
        //         if (symm_enforce) {
        //           symmetries[alpha_term] = {{mixed_term2, 1.0}, {mixed_term4, -1.0}};
        //           symmetries[beta_term] = {{mixed_term2, 1.0}, {mixed_term4, -1.0}};
        //           symmetries[mixed_term1] = {{mixed_term2, 1.0}};
        //           symmetries[mixed_term2] = {{mixed_term2, 1.0}};
        //           symmetries[mixed_term3] = {{mixed_term4, 1.0}};
        //           symmetries[mixed_term4] = {{mixed_term4, 1.0}};

        //           fermion_ops_to_params[mixed_term2] = unique_params++;
        //           fermion_ops_to_params[mixed_term4] = unique_params++;
        //           excitation_index_map[to_fermionic_string(fermion_operators[mixed_term2], env)] = unique_params - 2;
        //           excitation_index_map[to_fermionic_string(fermion_operators[mixed_term4], env)] = unique_params - 1;
        //         } else {

        //           symmetries[alpha_term] = {{alpha_term, 1.0}};
        //           symmetries[beta_term] = {{beta_term, 1.0}};
        //           symmetries[mixed_term1] = {{mixed_term1, 1.0}};
        //           symmetries[mixed_term2] = {{mixed_term2, 1.0}};
        //           symmetries[mixed_term3] = {{mixed_term3, 1.0}};
        //           symmetries[mixed_term4] = {{mixed_term4, 1.0}};
        //           symmetries[mixed_term5] = {{mixed_term5, 1.0}};
        //           symmetries[mixed_term6] = {{mixed_term6, 1.0}};
        //           symmetries[mixed_term7] = {{mixed_term7, 1.0}};
        //           symmetries[mixed_term8] = {{mixed_term8, 1.0}};

        //           fermion_ops_to_params[alpha_term] = unique_params++;
        //           fermion_ops_to_params[beta_term] = unique_params++;
        //           fermion_ops_to_params[mixed_term1] = unique_params++;
        //           fermion_ops_to_params[mixed_term2] = unique_params++;
        //           fermion_ops_to_params[mixed_term3] = unique_params++;
        //           fermion_ops_to_params[mixed_term4] = unique_params++;
        //           fermion_ops_to_params[mixed_term5] = unique_params++;
        //           fermion_ops_to_params[mixed_term6] = unique_params++;
        //           fermion_ops_to_params[mixed_term7] = unique_params++;
        //           fermion_ops_to_params[mixed_term8] = unique_params++;
        //           excitation_index_map[to_fermionic_string(fermion_operators[alpha_term], env)] = unique_params - 10;
        //           excitation_index_map[to_fermionic_string(fermion_operators[beta_term], env)] = unique_params - 9;
        //           excitation_index_map[to_fermionic_string(fermion_operators[mixed_term1], env)] = unique_params - 8;
        //           excitation_index_map[to_fermionic_string(fermion_operators[mixed_term2], env)] = unique_params - 7;
        //           excitation_index_map[to_fermionic_string(fermion_operators[mixed_term3], env)] = unique_params - 6;
        //           excitation_index_map[to_fermionic_string(fermion_operators[mixed_term4], env)] = unique_params - 5;
        //           excitation_index_map[to_fermionic_string(fermion_operators[mixed_term5], env)] = unique_params - 4;
        //           excitation_index_map[to_fermionic_string(fermion_operators[mixed_term6], env)] = unique_params - 3;
        //           excitation_index_map[to_fermionic_string(fermion_operators[mixed_term7], env)] = unique_params - 2;
        //           excitation_index_map[to_fermionic_string(fermion_operators[mixed_term8], env)] = unique_params - 1;
        //         }
        //       }
        //     }
        //   }
        // }
        // Degenerate occupied excitations (WITHOUT the degenerate occupied/degenerate virtual)
        // for (IdxType i = 0; i < env.n_occ; i++) {
        //   FermionOperator occ_down_1 (i, Occupied, Down, Annihilation, env.xacc_scheme);
        //   FermionOperator occ_up_1 (i, Occupied, Up, Annihilation, env.xacc_scheme);
        //     for (IdxType r = 1; r < env.n_virt; r++) {
        //     FermionOperator virt_down_2 (r, Virtual, Down, Creation, env.xacc_scheme);
        //     FermionOperator virt_up_2 (r, Virtual, Up, Creation, env.xacc_scheme);
        //       // disallow r == s
        //       for (IdxType s = 0; s < env.n_virt; s++) {
        //         if (s == r) {
        //           continue;
        //         }
        //       FermionOperator virt_down_3 (s, Virtual, Down, Creation, env.xacc_scheme);
        //       FermionOperator virt_up_3 (s, Virtual, Up, Creation, env.xacc_scheme);
        //         IdxType term = fermion_operators.size();
        //         fermion_operators.push_back({
        //               occ_up_1,
        //               occ_down_1,
        //               virt_down_2,
        //               virt_up_3});
        //         fermion_ops_to_params[term] = unique_params++;
        //         symmetries[term] = {{term, 1.0}};
        //         term++;
        //         fermion_operators.push_back({
        //               occ_down_1,
        //               occ_up_1,
        //               virt_up_2,
        //               virt_down_3});
        //       if (symm_enforce) {
        //         symmetries[term] = {{term - 1, 1.0}};
        //         excitation_index_map[to_fermionic_string(fermion_operators[term], env)] = unique_params - 1;
        //       } else {
        //         fermion_ops_to_params[term] = unique_params++;
        //         symmetries[term] = {{term, 1.0}};
        //         excitation_index_map[to_fermionic_string(fermion_operators[term-1], env)] = unique_params - 2;
        //         excitation_index_map[to_fermionic_string(fermion_operators[term], env)] = unique_params - 1;
        //       }
              //   term++;
              //   fermion_operators.push_back({
              //         occ_down_1,
              //         occ_up_1,
              //         virt_up_3,
              //         virt_down_2});
              // fermion_ops_to_params[term] = unique_params++;
              // symmetries[term] = {{term, 1.0}};
              //   term++;
              //   fermion_operators.push_back({
              //         occ_down_1,
              //         occ_up_1,
              //         virt_up_2,
              //         virt_down_3});
              // fermion_ops_to_params[term] = unique_params++;
              // symmetries[term] = {{term, 1.0}};
              // excitation_index_map[to_fermionic_string(fermion_operators[term-3], env)] = unique_params - 4;
              // excitation_index_map[to_fermionic_string(fermion_operators[term-2], env)] = unique_params - 3;
        //     }
        //   }
        // }

        // Degenerate virtual excitations (WITH the degenerate occupied/degenerate virtual)
        // for (IdxType i = 0; i < env.n_virt; i++) {
        //   FermionOperator virt_down_3 (i, Virtual, Down, Creation, env.xacc_scheme);
        //   FermionOperator virt_up_3 (i, Virtual, Up, Creation, env.xacc_scheme);
        //     for (IdxType r = 0; r < env.n_occ; r++) {
        //   FermionOperator occ_down_1 (r, Occupied, Down, Annihilation, env.xacc_scheme);
        //   FermionOperator occ_up_1 (r, Occupied, Up, Annihilation, env.xacc_scheme);
        //       // allow r == s
        //       for (IdxType s = 0; s < env.n_occ; s++) {
        //       FermionOperator occ_down_2 (s, Occupied, Down, Annihilation, env.xacc_scheme);
        //       FermionOperator occ_up_2 (s, Occupied, Up, Annihilation, env.xacc_scheme);
                
        //         IdxType term = fermion_operators.size();
        //         fermion_operators.push_back({
        //               occ_up_1,
        //               occ_down_2,
        //               virt_down_3,
        //               virt_up_3});
        //         fermion_ops_to_params[term] = unique_params++;
        //         symmetries[term] = {{term, 1.0}};
        //         excitation_index_map[to_fermionic_string(fermion_operators[term], env)] = unique_params - 1;
        //         if (r < s) {
        //           term++;
        //          fermion_operators.push_back({
        //               occ_down_1,
        //               occ_up_2,
        //               virt_up_3,
        //               virt_down_3});
        //         if (symm_enforce) {
        //           symmetries[term] = {{term - 1, 1.0}};
        //           excitation_index_map[to_fermionic_string(fermion_operators[term], env)] = unique_params - 1;
        //         } else {
        //           fermion_ops_to_params[term] = unique_params++;
        //           symmetries[term] = {{term, 1.0}};
        //           excitation_index_map[to_fermionic_string(fermion_operators[term-1], env)] = unique_params - 2;
        //           excitation_index_map[to_fermionic_string(fermion_operators[term], env)] = unique_params - 1;
        //         }
        //           // fermion_operators.push_back({
        //           //       occ_down_2,
        //           //       occ_up_1,
        //           //       virt_up_3,
        //           //       virt_down_3});
        //           // fermion_ops_to_params[term] = unique_params++;
        //           // symmetries[term] = {{term, 1.0}};
        //           // excitation_index_map[to_fermionic_string(fermion_operators[term], env)] = unique_params - 1;
        //           // term++;
        //           // fermion_operators.push_back({
        //           //       occ_down_1,
        //           //       occ_up_2,
        //           //       virt_up_3,
        //           //       virt_down_3});
        //           // fermion_ops_to_params[term] = unique_params++;
        //           // symmetries[term] = {{term, 1.0}};
        //           // excitation_index_map[to_fermionic_string(fermion_operators[term], env)] = unique_params - 1;
        //         }
                
        //     }
        //   }
        // }
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
      };
      void buildAnsatz() override {
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

  };
}; //namespace VQE
}; //namespace NWQSim