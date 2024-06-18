#include "circuit/ansatz.hpp"
namespace NWQSim {
  namespace VQE {
    void UCCSD::getFermionOps()
    {
      // Single excitation
      for (IdxType p = 0; p < env.n_occ; p++) {
        FermionOperator occupied_annihilation_down (p, Occupied, Down, Annihilation);
        FermionOperator occupied_annihilation_up (p, Occupied, Up, Annihilation);
        for (IdxType q = 0; q < env.n_virt; q++) {
          FermionOperator virtual_creation_down (q, Virtual, Down, Creation);
          FermionOperator virtual_creation_up (q, Virtual, Up, Creation);
          fermion_operators.push_back({occupied_annihilation_down, virtual_creation_down});
          fermion_operators.push_back({occupied_annihilation_up, virtual_creation_up});
// #ifndef NDEBUG
                // std::cout << "SE: " << virtual_creation_down.qubitIndex(env.n_occ, env.n_virt) << " " << 
                //              occupied_annihilation_down.qubitIndex(env.n_occ, env.n_virt) << std::endl;
                // std::cout << "SE: " <<  virtual_creation_up.qubitIndex(env.n_occ, env.n_virt) << " " << 
                //              occupied_annihilation_up.qubitIndex(env.n_occ, env.n_virt) << std::endl;
// #endif
        }
      }

      // Double excitation
      std::vector<FermionOperator> occ_1, occ_2, virt_1, virt_2;
      for (IdxType p = 0; p < env.n_virt; p++) {
        FermionOperator virt_down_1 (p, Virtual, Down, Creation);
        FermionOperator virt_up_1 (p, Virtual, Up, Creation);
        for (IdxType r = 0; r < env.n_occ; r++) {
          FermionOperator occ_down_1 (r, Occupied, Down, Annihilation);
          FermionOperator occ_up_1 (r, Occupied, Up, Annihilation);
          for (IdxType q = 0; q < env.n_virt; q++) {
            FermionOperator virt_down_2 (q, Virtual, Down, Creation);
            FermionOperator virt_up_2 (q, Virtual, Up, Creation);
              for (IdxType s = 0; s < env.n_occ; s++) {
              FermionOperator occ_down_2 (s, Occupied, Down, Annihilation);
              FermionOperator occ_up_2 (s, Occupied, Up, Annihilation);
              if (r > s && p > q) {
                // Down spin all together
                // std::cout << virt_down_1.qubitIndex(env.n_occ, env.n_virt) << " " << 
                //              virt_down_2.qubitIndex(env.n_occ, env.n_virt) << " " << 
                //              occ_down_1.qubitIndex(env.n_occ, env.n_virt) << " " << 
                //              occ_down_2.qubitIndex(env.n_occ, env.n_virt) << std::endl;
                fermion_operators.push_back({
                  occ_down_1,
                  occ_down_2,
                  virt_down_1,
                  virt_down_2,
                                             });
// #ifndef NDEBUG
                // std::cout <<  "DE: " << virt_up_1.qubitIndex(env.n_occ, env.n_virt) << " " << 
                //              virt_up_2.qubitIndex(env.n_occ, env.n_virt) << " " << 
                //              occ_up_1.qubitIndex(env.n_occ, env.n_virt) << " " << 
                //              occ_up_2.qubitIndex(env.n_occ, env.n_virt) << std::endl;
// #endif
                // Up spin all together
                fermion_operators.push_back({
                                             occ_up_1,
                                             occ_up_2,
                                             virt_up_1,
                                             virt_up_2,});
              }

// #ifndef NDEBUG
              // std::cout <<  "DE: " << virt_down_1.qubitIndex(env.n_occ, env.n_virt) << " " << 
              //               virt_up_2.qubitIndex(env.n_occ, env.n_virt) << " " << 
              //               occ_up_1.qubitIndex(env.n_occ, env.n_virt) << " " << 
              //               occ_down_2.qubitIndex(env.n_occ, env.n_virt) << std::endl;
// #endif
              fermion_operators.push_back({
                                           occ_up_1,
                                           occ_down_2,
                                           virt_down_1,
                                           virt_up_2,});
            }
              
          }
          }
        }
    };
    void UCCSD::buildAnsatz(std::vector<std::vector<PauliOperator> > pauli_oplist) {
      theta->resize(pauli_oplist.size() * trotter_n);
      for (IdxType i = 0; i < env.n_occ ; i++) {
        X(i);
        X(i + env.n_spatial);
      }
      IdxType index = 0; // parameter index, shares parameters for Pauli evolution gates corresponding to the same Fermionic operator within the same Trotter step
      for (auto& fermionic_group: pauli_oplist) {
        for (auto& pauli: fermionic_group) {
          assert (pauli.getCoeff().imag() == 0.0);
          double coeff = pauli.getCoeff().real();
          if (pauli.isNonTrivial() && abs(coeff) > 1e-10) { 
            ExponentialGate(pauli, OP::RZ, 2 * coeff, 0.0, index);
          }
        }
        index++;
      }
      for (IdxType i = 0; i < trotter_n - 1; i++) {
        for (auto& fermionic_group: pauli_oplist) {
          for (auto& pauli: fermionic_group) {
            double coeff = pauli.getCoeff().real();
            if (pauli.isNonTrivial() && abs(coeff) > 1e-10)  {  
              ExponentialGate(pauli, OP::RZ, 2 * coeff, 0.0, index);
            }
          }
          index++;
        }
      }
    }
}; //namespace VQE
}; //namespace NWQSim