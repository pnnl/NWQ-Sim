#include "circuit/ansatz.hpp"
#include "circuit/dynamic_ansatz.hpp"
namespace NWQSim {
  namespace VQE {
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
          excitation_index_map[to_fermionic_string(fermion_operators.back(), env)] = unique_params;
        }
      }
      // Beta Single Excitations
      for (IdxType p = 0; p < env.n_occ; p++) {
        FermionOperator occupied_annihilation_down (p, Occupied, Down, Annihilation, env.xacc_scheme);
        for (IdxType q = 0; q < env.n_virt; q++) {
          FermionOperator virtual_creation_down (q, Virtual, Down, Creation, env.xacc_scheme);
          // Add a pointer to the corresponding Up spin term to share parameters
          symmetries[fermion_operators.size()] = {{fermion_operators.size() - env.n_occ * env.n_virt, 1.0}};
          fermion_operators.push_back({occupied_annihilation_down, virtual_creation_down});
          excitation_index_map[to_fermionic_string(fermion_operators.back(), env)] = fermion_operators.size() - env.n_occ * env.n_virt - 1;
        }
      }
      /*===========Double Excitations===========*/
      // Mixed-spin symmetries, all distinct spatial orbitals
      for (IdxType i = 0; i < env.n_occ; i++) {
        FermionOperator occ_down_1 (i, Occupied, Down, Annihilation, env.xacc_scheme);
        FermionOperator occ_up_1 (i, Occupied, Up, Annihilation, env.xacc_scheme);
        // Strictly greater than i
        for (IdxType j = i+1; j < env.n_occ; j++) {
          FermionOperator occ_down_2 (j, Occupied, Down, Annihilation, env.xacc_scheme);
          FermionOperator occ_up_2 (j, Occupied, Up, Annihilation, env.xacc_scheme);
          for (IdxType r = 0; r < env.n_virt; r++) {
          FermionOperator virt_down_3 (r, Virtual, Down, Creation, env.xacc_scheme);
          FermionOperator virt_up_3 (r, Virtual, Up, Creation, env.xacc_scheme);
            // Strictly greater than r 
            for (IdxType s = r+1; s < env.n_virt; s++) {
              FermionOperator virt_down_4 (s, Virtual, Down, Creation, env.xacc_scheme);
              FermionOperator virt_up_4 (s, Virtual, Up, Creation, env.xacc_scheme);
              IdxType alpha_term = fermion_operators.size();
              // All alpha excitation: alpha_4*alpha_3*alpha_2*alpha_1
              fermion_operators.push_back({
                    occ_down_1,
                    occ_down_2,
                    virt_down_3,
                    virt_down_4});
              IdxType beta_term = fermion_operators.size();
              // All beta excitation: beta_4*beta_3*beta_2*beta_1
              fermion_operators.push_back({
                    occ_up_1, 
                    occ_up_2,
                    virt_up_3,
                    virt_up_4});
              IdxType mixed_term1 = fermion_operators.size();
              // Mixed term 1: beta_4*alpha_3*alpha_2*beta_1
              fermion_operators.push_back({
                    occ_up_1,
                    occ_down_2,
                    virt_down_3,
                    virt_up_4});
              IdxType mixed_term2 = fermion_operators.size();
              // Mixed term 2: beta_4*alpha_3*beta_2*alpha_1
              fermion_operators.push_back({
                    occ_down_1,
                    occ_up_2,
                    virt_down_3,
                    virt_up_4});
              // Add the parameter pointers, all 4 values determined by two parameters
              symmetries[alpha_term] = {{mixed_term1, 1.0}, {mixed_term2, 1.0}};
              symmetries[beta_term] = {{mixed_term1, 1.0}, {mixed_term2, 1.0}};
              symmetries[mixed_term1] = {{mixed_term1, 1.0}};
              symmetries[mixed_term2] = {{mixed_term2, 1.0}};

              fermion_ops_to_params[mixed_term1] = unique_params++;
              fermion_ops_to_params[mixed_term2] = unique_params++;
              excitation_index_map[to_fermionic_string(fermion_operators[mixed_term1], env)] = unique_params - 2;
              excitation_index_map[to_fermionic_string(fermion_operators[mixed_term2], env)] = unique_params - 1;
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
            FermionOperator virt_down_3 (s, Virtual, Down, Creation, env.xacc_scheme);
            FermionOperator virt_up_3 (s, Virtual, Up, Creation, env.xacc_scheme);
              IdxType term = fermion_operators.size();
              fermion_operators.push_back({
                    occ_up_1,
                    occ_down_1,
                    virt_down_2,
                    virt_up_3});
              fermion_ops_to_params[term] = unique_params++;
              symmetries[term] = {{term, 1.0}};
              term++;
              fermion_operators.push_back({
                    occ_up_1,
                    occ_down_1,
                    virt_down_3,
                    virt_up_2});
            fermion_ops_to_params[term] = unique_params++;
            symmetries[term] = {{term, 1.0}};
            excitation_index_map[to_fermionic_string(fermion_operators[term-1], env)] = unique_params - 2;
            excitation_index_map[to_fermionic_string(fermion_operators[term], env)] = unique_params - 1;
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
            for (IdxType s = 0; s < r + 1; s++) {
            FermionOperator occ_down_2 (s, Occupied, Down, Annihilation, env.xacc_scheme);
            FermionOperator occ_up_2 (s, Occupied, Up, Annihilation, env.xacc_scheme);
              
              IdxType term = fermion_operators.size();
              fermion_operators.push_back({
                    occ_up_1,
                    occ_down_2,
                    virt_down_3,
                    virt_up_3});
              fermion_ops_to_params[term] = unique_params++;
              symmetries[term] = {{term, 1.0}};
              excitation_index_map[to_fermionic_string(fermion_operators[term], env)] = unique_params - 1;

              if (r > s) {
                term++;
                fermion_operators.push_back({
                      occ_up_2,
                      occ_down_1,
                      virt_down_3,
                      virt_up_3});
                fermion_ops_to_params[term] = unique_params++;
                symmetries[term] = {{term, 1.0}};
                excitation_index_map[to_fermionic_string(fermion_operators[term], env)] = unique_params - 1;
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
    };
    void UCCSD::buildAnsatz() {
      getFermionOps();
      assert((n_doubles + n_singles) == fermion_operators.size());
      std::cout << "Generated " << n_doubles + n_singles << " operators." << std::endl;
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