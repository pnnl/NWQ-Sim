#include "observable/pauli_operator.hpp"
#include "observable/fermionic_operator.hpp"
#include "transform/transform.hpp"
#include "utils.hpp"
#include <random>

namespace NWQSim{
namespace VQE {
/**
 * @brief  Generate a set of single/double Fermionic excitations
 * @note   Same set of operators evolved in a UCCSD Ansatz
 * @param  fermion_operators: Output vector to store Fermionic operator products
 * @param  env: Struct containing molecular information
 * @retval None
 */
void generate_fermionic_excitations(std::vector<std::vector<std::vector<FermionOperator> > >& fermion_operators,
                                    const MolecularEnvironment& env) {
    // Pre-calculate expected sizes to avoid reallocations
    const size_t n_singles = env.n_occ * env.n_virt * 2;  // single excitations (alpha and beta)
    const size_t n_same_spin = env.n_occ * (env.n_occ - 1) * env.n_virt * (env.n_virt - 1) / 4;  // same-spin doubles
    const size_t n_mixed_spin = env.n_occ * env.n_occ * env.n_virt * env.n_virt / 2;  // mixed-spin doubles
    const size_t expected_size = n_singles + n_same_spin + n_mixed_spin;
    
    fermion_operators.clear();
    fermion_operators.reserve(expected_size);

    std::vector<std::vector<FermionOperator>> single_excitation(2);
    single_excitation[0].reserve(2);
    single_excitation[1].reserve(2);
    
    // Generate single excitations
    for (IdxType p = 0; p < env.n_occ; p++) {
        for (IdxType q = 0; q < env.n_virt; q++) {
            const FermionOperator occ_ann_up(p, Occupied, Up, Annihilation, env.xacc_scheme);
            const FermionOperator occ_ann_down(p, Occupied, Down, Annihilation, env.xacc_scheme);
            const FermionOperator virt_cre_up(q, Virtual, Up, Creation, env.xacc_scheme);
            const FermionOperator virt_cre_down(q, Virtual, Down, Creation, env.xacc_scheme);
            fermion_operators.push_back({
                {occ_ann_up, virt_cre_up},
                {occ_ann_down, virt_cre_down}
            });
        }
    }

    std::vector<std::vector<FermionOperator>> double_excitation(2);
    double_excitation[0].reserve(4);
    double_excitation[1].reserve(4);

    // Generate same-spin double excitations
    for (IdxType i = 0; i < env.n_occ; i++) {
        for (IdxType r = 0; r < env.n_virt; r++) {
            const FermionOperator i_occ_ann_up(i, Occupied, Up, Annihilation, env.xacc_scheme);
            const FermionOperator i_occ_ann_dw(i, Occupied, Down, Annihilation, env.xacc_scheme);
            const FermionOperator r_virt_cre_up(r, Virtual, Up, Creation, env.xacc_scheme);
            const FermionOperator r_virt_cre_dw(r, Virtual, Down, Creation, env.xacc_scheme);
            for (IdxType j = i + 1; j < env.n_occ; j++) {
                const FermionOperator j_occ_ann_up(j, Occupied, Up, Annihilation, env.xacc_scheme);
                const FermionOperator j_occ_ann_dw(j, Occupied, Down, Annihilation, env.xacc_scheme);
                for (IdxType s = r + 1; s < env.n_virt; s++) {
                    const FermionOperator s_virt_cre_up(s, Virtual, Up, Creation, env.xacc_scheme);
                    const FermionOperator s_virt_cre_dw(s, Virtual, Down, Creation, env.xacc_scheme);
                    fermion_operators.push_back({
                        {i_occ_ann_up, j_occ_ann_up, r_virt_cre_up, s_virt_cre_up},
                        {i_occ_ann_dw, j_occ_ann_dw, r_virt_cre_dw, s_virt_cre_dw}
                    });
                }
            }
        }
    }

    // Generate mixed-spin double excitations
    for (IdxType i = 0; i < env.n_occ; i++) {
        for (IdxType r = 0; r < env.n_virt; r++) {
            const FermionOperator i_occ_ann_up(i, Occupied, Up, Annihilation, env.xacc_scheme);
            const FermionOperator r_virt_cre_up(r, Virtual, Up, Creation, env.xacc_scheme);
            for (IdxType j = 0; j < env.n_occ; j++) {
                if (j < i) continue;
                const FermionOperator j_occ_ann_dw(j, Occupied, Down, Annihilation, env.xacc_scheme);
                for (IdxType s = 0; s < env.n_virt; s++) {
                    if (j == i && s < r) continue;
                    const FermionOperator s_virt_cre_dw(s, Virtual, Down, Creation, env.xacc_scheme);
                    if (i == j && r == s) {
                        fermion_operators.push_back({
                            {i_occ_ann_up, j_occ_ann_dw, r_virt_cre_up, s_virt_cre_dw}
                        });
                    } else {
                        fermion_operators.push_back({
                            {i_occ_ann_up, j_occ_ann_dw, r_virt_cre_up, s_virt_cre_dw},
                            {FermionOperator(j, Occupied, Up, Annihilation, env.xacc_scheme),
                             FermionOperator(i, Occupied, Down, Annihilation, env.xacc_scheme),
                             FermionOperator(s, Virtual, Up, Creation, env.xacc_scheme),
                             FermionOperator(r, Virtual, Down, Creation, env.xacc_scheme)}
                        });
                    }
                }
            }
        }
    }
}

/**
 * @brief  Generate a set of single/double singlet and triplet excitations
 * @note   Same set of operators evolved in a singlet GSD Ansatz for ADAPT-VQE in Qubit-ADAPT-VQE paper
 * @param  fermion_operators: Output vector to store Fermionic operator products
 * @param  env: Struct containing molecular information
 * @retval None
 */
void generate_singlet_gsd_excitations(std::vector<std::vector<std::vector<FermionOperator> > >& fermion_operators,
                                    const MolecularEnvironment& env) {
    const int n_singles = (env.n_spatial * (env.n_spatial - 1));
    const int n_doubles = counting_doubles(env.n_spatial);
    const size_t total_size = n_singles + n_doubles;
    
    fermion_operators.clear();
    fermion_operators.reserve(total_size);

    std::vector<std::vector<FermionOperator>> single_excitation(2);
    single_excitation[0].reserve(2);
    single_excitation[1].reserve(2);
    
    //===========Single Excitations===========
    for (IdxType p = 0; p < env.n_spatial; p++) {
        for (IdxType q = p+1; q < env.n_spatial; q++) {
            const IdxType pi = spind_to_ind(p, env.n_occ);
            const IdxType qi = spind_to_ind(q, env.n_occ);
            const auto pov = occ_or_vir(p, env.n_occ);
            const auto qov = occ_or_vir(q, env.n_occ);
            
            const FermionOperator virtual_creation_up(pi, pov, Up, Creation, env.xacc_scheme);
            const FermionOperator virtual_creation_down(pi, pov, Down, Creation, env.xacc_scheme);
            const FermionOperator occupied_annihilation_up(qi, qov, Up, Annihilation, env.xacc_scheme);
            const FermionOperator occupied_annihilation_down(qi, qov, Down, Annihilation, env.xacc_scheme);
            
            fermion_operators.push_back({
                {occupied_annihilation_up, virtual_creation_up},
                {occupied_annihilation_down, virtual_creation_down}
            });
        }
    }

    //===========Double Excitations===========
    std::vector<std::vector<FermionOperator>> double_excitation(6);  // Max size needed for triplet case
    for (auto& vec : double_excitation) {
        vec.reserve(4);
    }

    for (IdxType r = 0; r < env.n_spatial; r++) {
        for (IdxType s = r; s < env.n_spatial; s++) {
            const IdxType ri = spind_to_ind(r, env.n_occ);
            const IdxType si = spind_to_ind(s, env.n_occ);
            const auto rov = occ_or_vir(r, env.n_occ);
            const auto sov = occ_or_vir(s, env.n_occ);
            const FermionOperator ra(ri, rov, Up, Annihilation, env.xacc_scheme);
            const FermionOperator rb(ri, rov, Down, Annihilation, env.xacc_scheme);
            const FermionOperator sa(si, sov, Up, Annihilation, env.xacc_scheme);
            const FermionOperator sb(si, sov, Down, Annihilation, env.xacc_scheme);
            for (IdxType p = 0; p < env.n_spatial; p++) {
                const IdxType pi = spind_to_ind(p, env.n_occ);
                const auto pov = occ_or_vir(p, env.n_occ);
                const FermionOperator pa(pi, pov, Up, Creation, env.xacc_scheme);
                const FermionOperator pb(pi, pov, Down, Creation, env.xacc_scheme);
                for (IdxType q = p; q < env.n_spatial; q++) {
                    // Skip invalid combinations early
                    if ((p == r) && (q == s)) continue;
                    const IdxType qi = spind_to_ind(q, env.n_occ);
                    const auto qov = occ_or_vir(q, env.n_occ);
                    const FermionOperator qa(qi, qov, Up, Creation, env.xacc_scheme);
                    const FermionOperator qb(qi, qov, Down, Creation, env.xacc_scheme);
                    if (p != q) {
                        if (r == s) {
                            // Only singlet with two terms
                            fermion_operators.push_back({
                                {rb, ra, qb, pa},
                                {-1.0 * rb, ra, qa, pb}
                            });
                        } else if (((r != s) && (q != r)) || ((p == r) && (q != s) && (r != s)) || 
                                 ((p == s) && (q != r) && (r != s)) || ((q == s) && (p != r)) || 
                                 ((q == r) && (p != s))) {
                            // Triplet case
                            fermion_operators.push_back({
                                {2.0 * sa, ra, qa, pa},
                                {sb, ra, qb, pa},
                                {sa, rb, qb, pa},
                                {qa, pb, sb, ra},
                                {qa, pb, sa, rb},
                                {2.0 * sb, rb, qb, pb}
                            });
                        }
                    } else {  // p == q
                        if ((q != r) && (r != s)) {
                            // Group 3 - only singlet
                            fermion_operators.push_back({
                                {sb, ra, pb, pa},
                                {sa, rb, pb, pa}
                            });
                        } else if ((q == r) && (r != s)) {
                            // Group 4 - only singlet
                            fermion_operators.push_back({
                                {sb, ra, pa, pb},
                                {sa, rb, pb, pa}
                            });
                        } else if ((q != r) && (r == s)) {
                            // Group 5 - only singlet
                            fermion_operators.push_back({
                                {rb, ra, pb, pa}
                            });
                        }
                    }
                }
            }
        }
    }
};




// MZ: The problems in the following code is the error on the symmetry.
void generate_fermionic_excitations_origin(std::vector<std::vector<std::vector<FermionOperator> > >& fermion_operators,
                                    const MolecularEnvironment& env) {
    // Single excitation
    for (IdxType p = 0; p < env.n_occ; p++) {
      FermionOperator occupied_annihilation_up (p, Occupied, Up, Annihilation, env.xacc_scheme);
      FermionOperator occupied_annihilation_down (p, Occupied, Down, Annihilation, env.xacc_scheme);
      for (IdxType q = 0; q < env.n_virt; q++) {
        FermionOperator virtual_creation_up (q, Virtual, Up, Creation, env.xacc_scheme);
        FermionOperator virtual_creation_down (q, Virtual, Down, Creation, env.xacc_scheme);
        fermion_operators.push_back({{occupied_annihilation_up, virtual_creation_up},
                                      {occupied_annihilation_down, virtual_creation_down}});
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
        FermionOperator virt_down_3 (r, Virtual, Down, Creation, env.xacc_scheme);
        FermionOperator virt_up_3 (r, Virtual, Up, Creation, env.xacc_scheme);
          for (IdxType s = r+1; s < env.n_virt; s++) {
            FermionOperator virt_down_4 (s, Virtual, Down, Creation, env.xacc_scheme);
            FermionOperator virt_up_4 (s, Virtual, Up, Creation, env.xacc_scheme);
            IdxType alpha_term = fermion_operators.size();
            fermion_operators.push_back(
               {{occ_down_1,
                occ_down_2,
                virt_down_3,
                virt_down_4},
               {occ_up_1, 
                occ_up_2,
                virt_up_3,
                virt_up_4},
                {occ_up_1,
                 occ_down_2, 
                 virt_down_3,
                 virt_up_4},
                {occ_down_1,
                 occ_up_2,
                 virt_down_3,
                 virt_up_4}}
                );
          }
        }
      }
    }
    for (IdxType i = 0; i < env.n_occ; i++) {
      FermionOperator occ_down_1 (i, Occupied, Down, Annihilation, env.xacc_scheme);
      FermionOperator occ_up_1 (i, Occupied, Up, Annihilation, env.xacc_scheme);
        for (IdxType r = 1; r < env.n_virt; r++) {
        FermionOperator virt_down_2 (r, Virtual, Down, Creation, env.xacc_scheme);
        FermionOperator virt_up_2 (r, Virtual, Up, Creation, env.xacc_scheme);
          for (IdxType s = 0; s < r; s++) {
          FermionOperator virt_down_3 (s, Virtual, Down, Creation, env.xacc_scheme);
          FermionOperator virt_up_3 (s, Virtual, Up, Creation, env.xacc_scheme);
            IdxType term = fermion_operators.size();
            fermion_operators.push_back(
                 {{occ_up_1,
                  occ_down_1,
                  virt_down_2,
                  virt_up_3},
                  {occ_up_1,
                  occ_down_1,
                  virt_down_3,
                  virt_up_2}});
        }
      }
    }

    for (IdxType i = 0; i < env.n_virt; i++) {
      FermionOperator virt_down_3 (i, Virtual, Down, Creation, env.xacc_scheme);
      FermionOperator virt_up_3 (i, Virtual, Up, Creation, env.xacc_scheme);
        for (IdxType r = 0; r < env.n_occ; r++) {
      FermionOperator occ_down_1 (r, Occupied, Down, Annihilation, env.xacc_scheme);
      FermionOperator occ_up_1 (r, Occupied, Up, Annihilation, env.xacc_scheme);
          for (IdxType s = 0; s < r + 1; s++) {
          FermionOperator occ_down_2 (s, Occupied, Down, Annihilation, env.xacc_scheme);
          FermionOperator occ_up_2 (s, Occupied, Up, Annihilation, env.xacc_scheme);
            if (r > s) {
              fermion_operators.push_back({
                    {occ_up_1,
                    occ_down_2,
                    virt_down_3,
                    virt_up_3},
                    {occ_up_2,
                    occ_down_1,
                    virt_down_3,
                    virt_up_3}});
            } else {
              fermion_operators.push_back({
                    {occ_up_1,
                    occ_down_2,
                    virt_down_3,
                    virt_up_3}});
            }
        }
        }
      }
};




/**
 * @brief  Generate Pauli Operator Pool
 * @note   Same operators as UCCSD, just with single Pauli Strings. TODO: Remove redundant Paulis
 * @param  pauli_operators: Output vector of observables
 * @param  env: Molecular environment structure
 * @retval None
 */
void generate_pauli_excitations(std::vector<std::vector<PauliOperator> >& pauli_operators,
                                    const MolecularEnvironment& env,
                                    IdxType subsample,
                                    IdxType seed) {
  IdxType n_singles = env.n_occ * env.n_virt;
  IdxType n_doubles = env.n_occ * (env.n_occ) * env.n_virt * (env.n_virt) +\
              choose2(env.n_occ) * choose2(env.n_virt) * 2; 
  std::vector<std::vector<std::vector<FermionOperator> > > fermion_operators;
  fermion_operators.reserve(n_singles);
  generate_fermionic_excitations(fermion_operators, env);

  std::vector<PauliOperator> temporary_storage;
  if (subsample > 0) {
    temporary_storage.reserve(4 * n_singles + 16 * n_doubles);
  } else {
    pauli_operators.reserve(4 * n_singles + 16 * n_doubles);
  }
    
  IdxType index = 0;
  // NOTE: Not the most efficient way of doing this, kind of a workaround due to the data structures used elsewhere
  for (size_t i = 0; i < fermion_operators.size(); i++) {
  // Perform the JordanWigner mapping for each Fermionic operator product
    std::vector<std::vector<PauliOperator> > mapper_temp;
    getJordanWignerTransform(env, fermion_operators[i], mapper_temp, true);
    // Iterate ove the JW mapped Paulis and slap them on the operator pool
    for (auto pauli_list: mapper_temp) {
      for (auto pauli: pauli_list) {
        // If we're not subsampling or if this is a selected index
        if (subsample < 0) {
          pauli_operators.push_back({pauli});\
        } else {
          temporary_storage.push_back(pauli);
        }
        index++; // keep count of indices
      }
    }
  }
  std::vector<IdxType> selection;
  // If we're subsampling, generate the list of operator indices to choose a priori
  if (subsample > 0) {
    std::vector<IdxType> indices(temporary_storage.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::mt19937_64 random_engine (seed);
    std::shuffle(indices.begin(), indices.end(), random_engine);
    // take the first N indices, where N is either the number of samples or the number of operators (whichever is smaller)
    selection = std::vector<IdxType>(
      indices.begin(),
      indices.begin() + std::min(subsample, (IdxType)indices.size()));
    pauli_operators.reserve(selection.size());
    for (IdxType i : selection) {
      pauli_operators.push_back({temporary_storage[i]});
    }
  }
  std::vector<IdxType>::iterator iter = selection.begin();
};

 /**
  * @brief  Construct the minimal operator pool using the strategy desc. in Appdx. C of Tang et al. 2021
  * @note   
  * @param  _pauli_operators: Output vector of Pauli operators (each element is a singleton)
  * @param  _env: Molecular environment
  * @retval None
  */
  void generate_minimal_pauli_excitations(std::vector<std::vector<PauliOperator > >& _pauli_operators,
                                  const MolecularEnvironment& _env) {
    IdxType n_qubits = _env.n_spatial * 2;
    _pauli_operators.reserve(2 * n_qubits - 2);
    for (size_t i = 0; i < n_qubits - 1; i++) {
      // First construct the operator of the form II...Z_{i+1}Y_i...III
      IdxType xmask_1 = 1 << i;
      IdxType zmask_1 = (1 << i) + (1 << (i + 1));
      PauliOperator pauli_1(xmask_1, zmask_1, n_qubits);
      _pauli_operators.push_back({pauli_1});
      // Now construct the operator of the form II...Y_{i+1}...III
      IdxType xmask_2 = (1 << (i + 1));
      IdxType zmask_2 = (1 << (i + 1));
      PauliOperator pauli_2(xmask_2, zmask_2, n_qubits);
      _pauli_operators.push_back({pauli_2});
    }
  };
};// namespace VQE
};// namespace NWQSim