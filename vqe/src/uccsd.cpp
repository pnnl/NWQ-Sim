#include "circuit/ansatz.hpp"
#include "circuit/dynamic_ansatz.hpp"
#include "observable/fermionic_operator.hpp"
#include "utils.hpp"
namespace NWQSim {
  namespace VQE {
    void UCCSD::add_double_excitation(FermionOperator i, FermionOperator j, FermionOperator r, FermionOperator s,  const std::vector<std::pair<IdxType, double>>& symm_expr, bool param) {

        // use this index as the unique parameter to create the symmetry
        symmetries[fermion_operators.size()] = symm_expr;
        fermion_operators.push_back({i, j, r, s});
        // record the string::parameter mapping
        if (param) {
          fermion_ops_to_params[fermion_operators.size()-1] = unique_params++;
          excitation_index_map[to_fermionic_string(fermion_operators.back(), env)] = unique_params-1;
        }
    }
    void UCCSD::add_double_excitation(FermionOperator i, FermionOperator j, FermionOperator r, FermionOperator s) {

        // use this index as the unique parameter to create the symmetry
        symmetries[fermion_operators.size()] = {{fermion_operators.size(), 1.0}};
        fermion_operators.push_back({i, j, r, s});
        // record the string::parameter mapping
        fermion_ops_to_params[fermion_operators.size()-1] = unique_params++;
        excitation_index_map[to_fermionic_string(fermion_operators.back(), env)] = unique_params-1;

    }
   /**
    * @brief  Generate 4-term excitation
    * @note   operator is a_s^\dagger a_r^\dagger a_j a_i
    * @param  r: index of the first 
    * @retval None
    */
    void UCCSD::generate_mixed_excitation(IdxType i, IdxType j, IdxType r, IdxType s) {
      // Make the necessary single-orbital operators
      FermionOperator j_alpha (j, Occupied, Up, Annihilation, env.xacc_scheme);
      FermionOperator j_beta (j, Occupied, Down, Annihilation, env.xacc_scheme);
      FermionOperator s_alpha (s, Virtual, Up, Creation, env.xacc_scheme);
      FermionOperator s_beta (s, Virtual, Down, Creation, env.xacc_scheme);
      FermionOperator i_alpha (i, Occupied, Up, Annihilation, env.xacc_scheme);
      FermionOperator i_beta (i, Occupied, Down, Annihilation, env.xacc_scheme);
      FermionOperator r_alpha (r, Virtual, Up, Creation, env.xacc_scheme);
      FermionOperator r_beta (r, Virtual, Down, Creation, env.xacc_scheme);
      if ((i == j && r != s) || (i != j && r == s)) {
        if (symm_level >= 3) {
          IdxType term1 = fermion_operators.size();
          // only one parameter needed to enforce the symmetry 
          add_double_excitation(i_alpha, j_beta, r_beta, s_alpha, {{term1, 1.0}}, true);
          add_double_excitation(j_alpha, i_beta, s_beta, r_alpha, {{term1, 1.0}}, false);
          add_double_excitation(i_beta, j_alpha, s_alpha, r_beta, {{term1, 1.0}}, false);
          add_double_excitation(j_beta, i_alpha, r_alpha, s_beta, {{term1, 1.0}}, false);
        } else {
          add_double_excitation(i_alpha, j_beta, r_beta, s_alpha);
          add_double_excitation(i_alpha, j_beta, s_beta, r_alpha);
          add_double_excitation(i_beta, j_alpha, s_alpha, r_beta);
          add_double_excitation(j_beta, i_alpha, r_alpha, s_beta);
        }
        return;
      } else if (i == j && r == s) {
        // double excitation
        if (symm_level >= 2) {
          // only one parameter needed to enforce the symmetry 
          IdxType term1 = fermion_operators.size();
          add_double_excitation(i_alpha, j_beta, r_beta, s_alpha, {{term1, 1.0}}, true);
          add_double_excitation(i_beta, j_alpha, r_alpha, s_beta, {{term1, 1.0}}, false);
        } else {
          add_double_excitation(i_alpha, j_beta, r_beta, s_alpha);
          add_double_excitation(i_beta, j_alpha, r_alpha, s_beta);
        }
        return;

      }
      // use the mixed excitation terms as the free variables (alpha_r beta_s beta_j alpha_i - alpha_s beta_r beta_i alpha_j)
      IdxType mixed_term1 = fermion_operators.size();
      IdxType mixed_term2 = fermion_operators.size() + 1;
      if (symm_level >= 3) {
        // s _r _j i
        add_double_excitation(i_alpha, j_beta, r_beta, s_alpha, {{mixed_term1, 1.0}}, true);
        // s _r _i j
        add_double_excitation(j_alpha, i_beta, r_beta, s_alpha, {{mixed_term2, 1.0}}, true);
        // now add the pure alpha/beta terms
        add_double_excitation(j_alpha, i_alpha, r_alpha, s_alpha, {{mixed_term1, -1.0}, {mixed_term2, 1.0}}, false);
        add_double_excitation(j_alpha, i_alpha, s_alpha, r_alpha, {{mixed_term1, 1.0},  {mixed_term2, -1.0}}, false);
        add_double_excitation(i_alpha, j_alpha, r_alpha, s_alpha, {{mixed_term1, 1.0},  {mixed_term2, -1.0}}, false);
        add_double_excitation(i_alpha, j_alpha, s_alpha, r_alpha, {{mixed_term1, -1.0},  {mixed_term2, 1.0}}, false);
        add_double_excitation(j_beta, i_beta, r_beta, s_beta, {{mixed_term1, -1.0}, {mixed_term2, 1.0}}, false);
        add_double_excitation(j_beta, i_beta, s_beta, r_beta, {{mixed_term1, 1.0},  {mixed_term2, -1.0}}, false);
        add_double_excitation(i_beta, j_beta, r_beta, s_beta, {{mixed_term1, 1.0},  {mixed_term2, -1.0}}, false);
        add_double_excitation(i_beta, j_beta, s_beta, r_beta, {{mixed_term1, -1.0},  {mixed_term2, 1.0}}, false);

        // now for the mixed shenanigans
        // _s r j _i
        add_double_excitation(i_beta, j_alpha, r_alpha, s_beta, {{mixed_term1, 1.0}}, false);
        // _s r i _j
        add_double_excitation(j_beta, i_alpha, r_alpha, s_beta, {{mixed_term2, 1.0}}, false);
        // r _s _j i
        add_double_excitation(i_alpha, j_beta, s_beta, r_alpha, {{mixed_term2, 1.0}}, false);
        // _r s j _i = s _r _i j
        add_double_excitation(i_alpha, j_beta, s_beta, r_alpha, {{mixed_term2, 1.0}}, false);
        // r _s _i j = s _r _j i
        add_double_excitation(i_alpha, j_beta, s_beta, r_alpha, {{mixed_term1, 1.0}}, false);
        // _r s i _j = s _r _j i
        add_double_excitation(i_beta, j_alpha, s_alpha, r_beta, {{mixed_term1, 1.0}}, false);
      } else {
        // s _r _j i
        add_double_excitation(i_alpha, j_beta, r_beta, s_alpha);
        // s _r _i j
        add_double_excitation(j_alpha, i_beta, r_beta, s_alpha);
        // now add the pure alpha/beta terms
        add_double_excitation(j_alpha, i_alpha, r_alpha, s_alpha);
        add_double_excitation(j_alpha, i_alpha, s_alpha, r_alpha);
        add_double_excitation(i_alpha, j_alpha, r_alpha, s_alpha);
        add_double_excitation(i_alpha, j_alpha, s_alpha, r_alpha);
        add_double_excitation(j_beta, i_beta, r_beta, s_beta);
        add_double_excitation(j_beta, i_beta, s_beta, r_beta);
        add_double_excitation(i_beta, j_beta, r_beta, s_beta);
        add_double_excitation(i_beta, j_beta, s_beta, r_beta);

        // now for the mixed shenanigans
        // _s r j _i
        add_double_excitation(i_beta, j_alpha, r_alpha, s_beta);
        // _s r i _j
        add_double_excitation(j_beta, i_alpha, r_alpha, s_beta);
        // r _s _j i
        add_double_excitation(i_alpha, j_beta, s_beta, r_alpha);
        // _r s j _i
        add_double_excitation(i_alpha, j_beta, s_beta, r_alpha);
        // r _s _i j
        add_double_excitation(i_alpha, j_beta, s_beta, r_alpha);
        // _r s i _j
        add_double_excitation(i_beta, j_alpha, s_alpha, r_beta);
      }
      


    }
   /**
    * @brief  Generate Fermionic operators for UCCSD
    * @note   Symmetry-linked operators (e.g. by anticommutation, spin reversal) share parameters
    * @retval None
    */
    void UCCSD::getFermionOps()
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
      // Mixed-spin symmetries, all distinct spatial orbitals
      for (IdxType i = 0; i < env.n_occ; i++) {
        FermionOperator occ_down_1 (i, Occupied, Down, Annihilation, env.xacc_scheme);
        FermionOperator occ_up_1 (i, Occupied, Up, Annihilation, env.xacc_scheme);
        // Strictly greater than i
        for (IdxType j = i+1; j < env.n_occ; j++) {
          if (i == j) {
            continue;
          }
          FermionOperator occ_down_2 (j, Occupied, Down, Annihilation, env.xacc_scheme);
          FermionOperator occ_up_2 (j, Occupied, Up, Annihilation, env.xacc_scheme);
          for (IdxType r = 0; r < env.n_virt; r++) {
            FermionOperator virt_down_3 (r, Virtual, Down, Creation, env.xacc_scheme);
            FermionOperator virt_up_3 (r, Virtual, Up, Creation, env.xacc_scheme);
            // Strictly greater than r 
            for (IdxType s = r+1; s < env.n_virt; s++) {
              if (r == s)
                continue;
              
              generate_mixed_excitation(i, j, r, s);
            }
          }
        }
      }
      // Degenerate occupied excitations (WITHOUT the degenerate occupied/degenerate virtual)
      for (IdxType i = 0; i < env.n_occ; i++) {
        FermionOperator occ_down_1 (i, Occupied, Down, Annihilation, env.xacc_scheme);
        FermionOperator occ_up_1 (i, Occupied, Up, Annihilation, env.xacc_scheme);
          for (IdxType r = 1; r < env.n_virt; r++) {
          FermionOperator virt_down_2 (r, Virtual, Down, Creation, env.xacc_scheme);
          FermionOperator virt_up_2 (r, Virtual, Up, Creation, env.xacc_scheme);
            // disallow r == s
            for (IdxType s = 0; s < r; s++) {
              if (s == r) {
                continue;
              }
             generate_mixed_excitation(i, i, r, s);
          }
        }
      }

      // Degenerate virtual excitations (WITH the degenerate occupied/degenerate virtual)
      for (IdxType i = 0; i < env.n_virt; i++) {
        FermionOperator virt_down_3 (i, Virtual, Down, Creation, env.xacc_scheme);
        FermionOperator virt_up_3 (i, Virtual, Up, Creation, env.xacc_scheme);
          for (IdxType r = 0; r < env.n_occ; r++) {
        FermionOperator occ_down_1 (r, Occupied, Down, Annihilation, env.xacc_scheme);
        FermionOperator occ_up_1 (r, Occupied, Up, Annihilation, env.xacc_scheme);
            // allow r == s
            for (IdxType s = 0; s < r+1; s++) {
                generate_mixed_excitation(r, s, i, i);
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
    };
    void UCCSD::buildAnsatz() {
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
    void DynamicAnsatz::buildAnsatz() {
      // exit(0);
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

      
    }
}; //namespace VQE
}; //namespace NWQSim