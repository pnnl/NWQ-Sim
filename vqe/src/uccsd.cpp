#include "circuit/ansatz.hpp"
#include "environment.hpp"
#include "observable/fermionic_operator.hpp"
#include "utils.hpp"
#include <unordered_map>
namespace NWQSim {
  namespace VQE {
    void UCCSD::getFermionOps()
    {
      // Single excitation
      for (IdxType p = 0; p < env.n_occ; p++) {
        FermionOperator occupied_annihilation_up (p, Occupied, Up, Annihilation, env.xacc_scheme);
        for (IdxType q = 0; q < env.n_virt; q++) {
          FermionOperator virtual_creation_up (q, Virtual, Up, Creation, env.xacc_scheme);
          symmetries[fermion_operators.size()] = {{fermion_operators.size(), 1.0}};
          fermion_ops_to_params[fermion_operators.size()] = unique_params++;
          fermion_operators.push_back({occupied_annihilation_up, virtual_creation_up});
// #ifndef NDEBUG
                // std::cout << "SE: " << virtual_creation_down.qubitIndex(env.n_occ, env.n_virt) << " " << 
                //              occupied_annihilation_down.qubitIndex(env.n_occ, env.n_virt) << std::endl;
                // std::cout << "SE: " <<  virtual_creation_up.qubitIndex(env.n_occ, env.n_virt) << " " << 
                //              occupied_annihilation_up.qubitIndex(env.n_occ, env.n_virt) << std::endl;
// #endif
        }
      }
      for (IdxType p = 0; p < env.n_occ; p++) {
        FermionOperator occupied_annihilation_down (p, Occupied, Down, Annihilation, env.xacc_scheme);
        for (IdxType q = 0; q < env.n_virt; q++) {
          FermionOperator virtual_creation_down (q, Virtual, Down, Creation, env.xacc_scheme);
          symmetries[fermion_operators.size()] = {{fermion_operators.size() - env.n_occ * env.n_virt, 1.0}};
          fermion_operators.push_back({occupied_annihilation_down, virtual_creation_down});
// #ifndef NDEBUG
                // std::cout << "SE: " << virtual_creation_down.qubitIndex(env.n_occ, env.n_virt) << " " << 
                //              occupied_annihilation_down.qubitIndex(env.n_occ, env.n_virt) << std::endl;
                // std::cout << "SE: " <<  virtual_creation_up.qubitIndex(env.n_occ, env.n_virt) << " " << 
                //              occupied_annihilation_up.qubitIndex(env.n_occ, env.n_virt) << std::endl;
// #endif
        }
      }
      // Double excitation
    for (IdxType i = 0; i < env.n_occ; i++) {
      FermionOperator occ_down_1 (i, Occupied, Down, Annihilation, env.xacc_scheme);
      FermionOperator occ_up_1 (i, Occupied, Up, Annihilation, env.xacc_scheme);
      for (IdxType j = i+1; j < env.n_occ; j++) {
        FermionOperator occ_down_2 (j, Occupied, Down, Annihilation, env.xacc_scheme);
        FermionOperator occ_up_2 (j, Occupied, Up, Annihilation, env.xacc_scheme);
        for (IdxType r = 0; r < env.n_virt; r++) {
        FermionOperator virt_down_1 (r, Virtual, Down, Creation, env.xacc_scheme);
        FermionOperator virt_up_1 (r, Virtual, Up, Creation, env.xacc_scheme);
          for (IdxType s = r+1; s < env.n_virt; s++) {
            FermionOperator virt_down_2 (s, Virtual, Down, Creation, env.xacc_scheme);
            FermionOperator virt_up_2 (s, Virtual, Up, Creation, env.xacc_scheme);
            IdxType alpha_term = fermion_operators.size();
            fermion_operators.push_back({
                  occ_down_1,
                  occ_down_2,
                  virt_down_2,
                  virt_down_1});
            IdxType beta_term = fermion_operators.size();
            fermion_operators.push_back({
                  occ_up_1, 
                  occ_up_2,
                  virt_up_2,
                  virt_up_1});
            IdxType mixed_term1 = fermion_operators.size();
            fermion_operators.push_back({
                  occ_up_1,
                  occ_down_2,
                  virt_down_2,
                  virt_up_1});
            IdxType mixed_term2 = fermion_operators.size();
            fermion_operators.push_back({
                  occ_down_1,
                  occ_up_2,
                  virt_down_2,
                  virt_up_1});
            symmetries[alpha_term] = {{mixed_term1, 1.0}, {mixed_term2, -1.0}};
            symmetries[beta_term] = {{mixed_term1, 1.0}, {mixed_term2, -1.0}};
            symmetries[mixed_term1] = {{mixed_term1, 1.0}};
            symmetries[mixed_term2] = {{mixed_term2, 1.0}};

          fermion_ops_to_params[mixed_term1] = unique_params++;
          fermion_ops_to_params[mixed_term2] = unique_params++;
          }
        }
      }
    }
    for (IdxType i = 0; i < env.n_occ; i++) {
      FermionOperator occ_down_1 (i, Occupied, Down, Annihilation, env.xacc_scheme);
      FermionOperator occ_up_1 (i, Occupied, Up, Annihilation, env.xacc_scheme);
      for (IdxType j = 0; j < i + 1; j++) {
        FermionOperator occ_down_2 (j, Occupied, Down, Annihilation, env.xacc_scheme);
        FermionOperator occ_up_2 (j, Occupied, Up, Annihilation, env.xacc_scheme);
        for (IdxType r = 0; r < env.n_virt; r++) {
        FermionOperator virt_down_1 (r, Virtual, Down, Creation, env.xacc_scheme);
        FermionOperator virt_up_1 (r, Virtual, Up, Creation, env.xacc_scheme);
          for (IdxType s = 0; s < r + 1; s++) {
          FermionOperator virt_down_2 (s, Virtual, Down, Creation, env.xacc_scheme);
          FermionOperator virt_up_2 (s, Virtual, Up, Creation, env.xacc_scheme);
            
            IdxType term = fermion_operators.size();
            fermion_operators.push_back({
                  occ_up_1,
                  occ_down_2,
                  virt_down_2,
                  virt_up_1});
          fermion_ops_to_params[term] = unique_params++;
          symmetries[term] = {{term, 1.0}};
          if (i != j || r != s) {
            symmetries.at(fermion_operators.size()) = {{fermion_operators.size(), 1.0}};
            fermion_ops_to_params.at(fermion_operators.size()) = unique_params++;
            fermion_operators.push_back({
                  occ_down_1,
                  occ_up_2,
                  virt_up_2,
                  virt_down_1});
          }
            
          }
        }
      }
    }
#ifndef NDEBUG
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
    struct
    FermiHash {
      size_t operator()(const std::vector<FermionOperator>& op) const {
        size_t value = 1;
        for (auto i: op) {
          value *= std::hash<IdxType>{}(i.getType() + 1 + i.getOrbitalType() * 2 + i.getType()*4 + i.getSpin()*8);
        }
        return value;
      }
    };
    struct
    FermiEq {
      bool operator()(const std::vector<FermionOperator>& op, const std::vector<FermionOperator>& op2) const {
        if (op.size() != op2.size()) {
          return false;
        }
        auto it1 = op.begin();
        auto it2 = op2.begin();
        for (size_t i = 0; i < op.size(); i++) {
          if (*it1 != *it2) {
            return false;
          }
        }
        return true;
      }
    };
    using FermiMap = std::unordered_map<std::vector<FermionOperator>, double, FermiHash, FermiEq>;
    void UCCSD::loadParameters(std::string amplitude_path) {
      std::vector<std::vector<FermionOperator>> provided_excitations;
      MolecularEnvironment env2 = env;
      read_fermion_operators(amplitude_path, provided_excitations, env2, env.n_part, env.xacc_scheme);
      FermiMap mp;
      for (auto i: provided_excitations) {
        mp[i] = i.front().getCoeff().real();
      }
      IdxType index = 0;
      for (auto i: fermion_operators) {
        if (mp.find(i) != mp.end()) {
          theta->assign(index++, mp[i]);
        }
      }

    }
    void UCCSD::buildAnsatz(std::vector<std::vector<PauliOperator> > pauli_oplist) {
      unique_params *= trotter_n;
      for (IdxType i = 0; i < env.n_occ; i++) {
        X(i);
        X(i+env.n_spatial);
      }
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
      IdxType params_per_rep = unique_params / trotter_n;
      for (IdxType i = 0; i < trotter_n - 1; i++) {
        for (auto& fermionic_group: pauli_oplist) {
          for (auto& pauli: fermionic_group) {
            double coeff = pauli.getCoeff().real();
            if (pauli.isNonTrivial() && abs(coeff) > 1e-10)  {  
            std::vector<std::pair<IdxType, ValType> > idxvals = symmetries[index];
            std::transform(symmetries[index].begin(), symmetries[index].end(),idxvals.begin(), 
            [params_per_rep, i] (std::pair<IdxType, ValType> val) {
              return std::make_pair(val.first + (i + 1) * params_per_rep, val.second);
            } );
            ExponentialGate(pauli, OP::RZ, idxvals, 2 * coeff);
            }
          }
          index++;
        }
      }
    }
}; //namespace VQE
}; //namespace NWQSim