#include "observable/pauli_operator.hpp"
#include "observable/fermionic_operator.hpp"
#include "transform/transform.hpp"
#include "utils.hpp"
#include <algorithm>
#include <random>
#include <stdexcept>
#include <string>
#include <iostream>
#include <sstream>
#include <queue>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <regex>

// // Muqing: for modifying generate_fermionic_excitations()
// #include <set>
// #include <tuple>

// // Sean: Modify QWC or GC commutation grouping in sorted_insertion()

namespace NWQSim{
namespace VQE {
  std::ostream& operator<<(std::ostream& out, const PauliOperator& op) {
        out << op.pauliToString();
        return out;
      };
  std::string to_binary_string(NWQSim::IdxType val, NWQSim::IdxType n_qubits) {
    std::stringstream ss;
    for (IdxType i = n_qubits-1; i >= 0; i--) {
      ss << ((val & (1ll << i)) ? "1" : "0");
    }

    return ss.str();

  }
  struct BFSFrame {
  IdxType node;
  IdxType hops;
};
void all_pairs_shortest_paths(const std::vector<std::list<IdxType> >& adj_list, std::vector<IdxType>& distances) {
  const IdxType white = 0; // not yet visited
  const IdxType black = 1; // visited
  const IdxType N = adj_list.size();
  distances.resize(N * N);
  std::fill(distances.begin(), distances.end(), -1);
  // assume graph is connected
  for (IdxType node = 0; node < N; node++) {
    std::vector<IdxType> colors(N, white);
    colors[node] = black;
    std::queue<BFSFrame> bfs_queue;
    bfs_queue.push({node, 0});
    IdxType queue_node, hops;
    while(!bfs_queue.empty()) {
      BFSFrame frame = bfs_queue.front();
      bfs_queue.pop();
      queue_node = frame.node;
      hops = frame.hops;

      distances[node * N + queue_node] = hops;
      distances[queue_node * N + node] = hops;

      for (auto neighbor: adj_list[queue_node]) {
        if (colors[neighbor] == white) {
          colors[neighbor] = black;
          bfs_queue.push({neighbor, hops+1});
        }
      }
    }
}
}
void readgraph(std::string fpath, std::vector<std::list<IdxType> >& adj_list) {
  std::ifstream infile;
  infile.open(fpath);
  if (!infile.is_open()) {
    throw std::runtime_error("File " + fpath + " not found\n");
  }
  IdxType u, v;
  IdxType nodes, edges;
  infile >> nodes >> edges;
  adj_list.resize(nodes);
  std::istringstream linestream;
  std::string line;
  while (std::getline(infile, line)) {
    if (line.length() < 2) {
      continue;
    }
    linestream = std::istringstream(line);
    linestream >> u;
    linestream >> v;
    adj_list[u].push_back(v);
    adj_list[v].push_back(u);
  }
}

/*QWC or GC selected in this function!*/
void sorted_insertion(const std::vector<PauliOperator>& paulilist, std::list<std::vector<IdxType> >& cliques, bool overlap) {
std::vector<IdxType> sorted_nodes;
    sorted_nodes.reserve(paulilist.size());
    std::vector<IdxType> argorder(paulilist.size());
    std::iota(argorder.begin(), argorder.end(), 0);
    // sort the nodelist in descending coeff order
    std::sort(argorder.begin(), argorder.end(), [&paulilist] (IdxType i, IdxType j) {
        return abs(paulilist[i].getCoeff().real()) > abs(paulilist[j].getCoeff().real());});
    for (auto i: argorder) {
        sorted_nodes.push_back(i);
    }
    std::vector<IdxType> added_already;
    std::vector<bool> added_flag (sorted_nodes.size(), 0);
    if (overlap)
        added_already.reserve(sorted_nodes.size());
    while (sorted_nodes.size() > 0) {
        std::list<IdxType> clique;
        std::vector<IdxType> not_added;
        not_added.reserve(sorted_nodes.size());
        for (IdxType node: sorted_nodes) {
            if (clique.empty()) {
                clique.push_back(node);
                continue;
            }
            
            bool commutes = true;
            for (auto other: clique) {
                if (!paulilist[node].QWC(paulilist[other])) {
                    commutes = false;
                    break;
                }
            }
            if (commutes) {
                clique.push_back(node);
            } else {
                not_added.push_back(node);
            }
        }
        if (overlap) {
            for (IdxType other_node: added_already) {
                bool commutes = true;
                for (auto resident: clique) {
                    if (!paulilist[other_node].QWC(paulilist[resident])) {
                        commutes = false;
                        break;
                    }
                } 
                if (commutes) {
                    clique.push_back(other_node);
                } 
            }
            for (IdxType i: clique) {
                if (added_flag[i] == 0) {
                    added_already.push_back(i);
                    added_flag[i] = 1;
                }
            }
        }
        if (!clique.empty()) {
          cliques.push_back(std::vector<IdxType>(clique.begin(), clique.end()));
        }
        sorted_nodes = not_added;
    }
  std::cout << "Cliques: " << cliques.size() << std::endl;
}
/**
 * @brief Make an XACC-formatted excitation string for a product of fermionic operators 
 * @note   
 * @param  product: Vector of individual annihilation or creation operators
 * @retval 
 */
std::string to_fermionic_string(const std::vector<FermionOperator>& product, const MolecularEnvironment& env){
    std::string opstring = "";
    bool first = true;
    IdxType prodterms = product.size();
    IdxType index = 0;
    for (auto& op: product) {
      if (!first) {
        opstring = " " + opstring;
      } 
      first = false;
      opstring = op.toString(env.n_occ, env.n_virt) + opstring;
    }
    return opstring;
};


const std::regex pattern("^\\s*(\\d+)\\s+(\\d+)\\s+(?:(\\d+)\\s+(\\d+)){0,1}\\s+([\\d\\.\\+e-]+)", std::regex_constants::multiline);
/**
 * @brief Read in a set of amplitudes
 * @note   
 * @param  product: Vector of individual annihilation or creation operators
 * @retval
 * */ 
void read_amplitudes(std::string fpath, std::vector<ValType>& params, const std::unordered_map<std::string, IdxType>& idx_map){
    std::ifstream infile;
    infile.open(fpath);
    if (!infile.is_open()) {
      std::cout << "Could not open amplitude file " << fpath << std::endl;
    }
    std::vector<bool> read_in(params.size(), 0);
    std::string line;
    std::smatch match;
    while(std::getline(infile, line)) {
        if (line.length() == 0) {
          continue;
        }

        if (std::regex_search(line, match, pattern)) {
          ValType amplitude = std::stod(match.str(5));
          std::string key;
          // if single excitation
          
          if (match.str(3) == "") {
            key = match.str(1) + "^ " + match.str(2);
          } else {
            key = match.str(1) + "^ " + match.str(2) + "^ " + match.str(3) + " " + match.str(4);
          }
          if (idx_map.find(key) == idx_map.end()) {
            continue;
          }
          params.at(idx_map.at(key)) = amplitude;
          read_in.at(idx_map.at(key)) = 1;
        }
    }
    IdxType index = 0;
    for (auto i: read_in) {
      if (!i) {
        // throw std::runtime_error("Parameter " + std::to_string(index) + " not provided\n");
        params.at(index) = 0;
      }
      index++;
    }
};

/* -------------------------------- MZ: MOVED TO ansatz_pool.cpp --------------------------------*/
// MZ: for better naming and the rest functions has nothing to do these pool generation functions
// /**
//  * @brief  Generate a set of single/double Fermionic excitations
//  * @note   Same set of operators evolved in a UCCSD Ansatz
//  * @param  fermion_operators: Output vector to store Fermionic operator products
//  * @param  env: Struct containing molecular information
//  * @retval None
//  */
// // void generate_singlet_gsd_excitations(std::vector<std::vector<std::vector<FermionOperator> > >& fermion_operators,
// //                                     const MolecularEnvironment& env) {
// void generate_fermionic_excitations(std::vector<std::vector<std::vector<FermionOperator> > >& fermion_operators,
//                                     const MolecularEnvironment& env) {
//     // int single_counter = 0;
//     // int double_counter = 0;
//     std::set<std::tuple< IdxType, IdxType, IdxType, IdxType >> existing_tuples; // MZ: for recording symmetry
//     // MZ: could use a better indexing iteration scheme to avoid the use of sets for efficiency
//     // but usually this is not a bottleneck
//     /*===========Single Excitations===========*/
//     for (IdxType p = 0; p < env.n_occ; p++) {
//       FermionOperator occupied_annihilation_up (p, Occupied, Up, Annihilation, env.xacc_scheme);
//       FermionOperator occupied_annihilation_down (p, Occupied, Down, Annihilation, env.xacc_scheme);
//       for (IdxType q = 0; q < env.n_virt; q++) {
//         // creation operator
//         FermionOperator virtual_creation_up (q, Virtual, Up, Creation, env.xacc_scheme);
//         FermionOperator virtual_creation_down (q, Virtual, Down, Creation, env.xacc_scheme);
//         fermion_operators.push_back({{occupied_annihilation_up, virtual_creation_up},
//                                       {occupied_annihilation_down, virtual_creation_down}});
//         // single_counter += 1;
//       }
//     }
//     /*===========Double Excitations===========*/
//     // alpha-alpha and beta-beta
//     for (IdxType i = 0; i < env.n_occ; i++) {
//         FermionOperator i_occ_ann_up (i, Occupied, Up, Annihilation, env.xacc_scheme);
//         FermionOperator i_occ_ann_dw (i, Occupied, Down, Annihilation, env.xacc_scheme);
//         for (IdxType r = 0; r < env.n_virt; r++) {
//             FermionOperator r_virt_cre_up (r, Virtual, Up, Creation, env.xacc_scheme);
//             FermionOperator r_virt_cre_dw (r, Virtual, Down, Creation, env.xacc_scheme);
//             for (IdxType j = i+1; j < env.n_occ; j++) {
//                 FermionOperator j_occ_ann_up (j, Occupied, Up, Annihilation, env.xacc_scheme); 
//                 FermionOperator j_occ_ann_dw (j, Occupied, Down, Annihilation, env.xacc_scheme);
//                 for (IdxType s = r+1; s < env.n_virt; s++) {
//                     FermionOperator s_virt_cre_dw (s, Virtual, Down, Creation, env.xacc_scheme);
//                     FermionOperator s_virt_cre_up (s, Virtual, Up, Creation, env.xacc_scheme);
//                     fermion_operators.push_back({ 
//                       {i_occ_ann_up, j_occ_ann_up, r_virt_cre_up, s_virt_cre_up},
//                       {i_occ_ann_dw, j_occ_ann_dw, r_virt_cre_dw, s_virt_cre_dw}
//                     });
//                     // double_counter += 1;
//                 }
//             }
//         }
//     }
//     // alpha-beta
//     for (IdxType i = 0; i < env.n_occ; i++) {
//         FermionOperator i_occ_ann_up (i, Occupied, Up, Annihilation, env.xacc_scheme);
//         FermionOperator i_occ_ann_dw (i, Occupied, Down, Annihilation, env.xacc_scheme);
//         for (IdxType r = 0; r < env.n_virt; r++) {
//             FermionOperator r_virt_cre_up (r, Virtual, Up, Creation, env.xacc_scheme);
//             FermionOperator r_virt_cre_dw (r, Virtual, Down, Creation, env.xacc_scheme);
//             for (IdxType j = 0; j < env.n_occ; j++) {
//                 FermionOperator j_occ_ann_dw (j, Occupied, Down, Annihilation, env.xacc_scheme);
//                 FermionOperator j_occ_ann_up (j, Occupied, Up, Annihilation, env.xacc_scheme);
//                 for (IdxType s = 0; s < env.n_virt; s++) {
//                     FermionOperator s_virt_cre_dw (s, Virtual, Down, Creation, env.xacc_scheme);
//                     FermionOperator s_virt_cre_up (s, Virtual, Up, Creation, env.xacc_scheme);
//                     if (i == j && r == s) {
//                       fermion_operators.push_back({{i_occ_ann_up, j_occ_ann_dw, r_virt_cre_up, s_virt_cre_dw}});
//                       // double_counter += 1;
//                     } else {
//                       std::tuple<IdxType, IdxType, IdxType, IdxType> new_tuple = {i,j,r,s};
//                       if (existing_tuples.find(new_tuple) != existing_tuples.end()) {
//                         // The tuple exist in the set, so we skip the term and erase the tuple
//                         existing_tuples.erase(new_tuple);
//                       } else {
//                         // The tuple does not exist in the set, so we add the term
//                         fermion_operators.push_back({
//                           {i_occ_ann_up, j_occ_ann_dw, r_virt_cre_up, s_virt_cre_dw},
//                           {j_occ_ann_up, i_occ_ann_dw, s_virt_cre_up, r_virt_cre_dw}
//                         });
//                         existing_tuples.insert({j, i, s, r});
//                         // double_counter += 1;
//                       }
//                     }
//                 }
//             }
//         }
//     }
// };



// /**
//  * @brief  Generate a set of singlet single/double Fermionic excitations
//  * @note   Same set of operators evolved in a singlet GSD Ansatz
//  * @param  fermion_operators: Output vector to store Fermionic operator products
//  * @param  env: Struct containing molecular information
//  * @retval None
//  */
// // void generate_fermionic_excitations(std::vector<std::vector<std::vector<FermionOperator> > >& fermion_operators,
// //                                     const MolecularEnvironment& env) {
// void generate_singlet_gsd_excitations(std::vector<std::vector<std::vector<FermionOperator> > >& fermion_operators,
//                                     const MolecularEnvironment& env) {
//     std::cout << ">>>> Select SingletGSD Ansatz <<<<\n" << std::endl;
//     int single_counter = 0;
//     int double_counter = 0;
//     int doublt_term1_counter = 0;
//     int doublt_term2_counter = 0;
//     int doublt_term4_counter = 0;
//     int doublt_term6_counter = 0;
//     std::set<std::tuple< IdxType, IdxType, IdxType, IdxType >> existing_tuples; // MZ: for recording symmetry
//     // MZ: could use a better indexing iteration scheme to avoid the use of sets for efficiency
//     // but usually this is not a bottleneck
//     //===========Single Excitations===========
//      for (IdxType p = 0; p < env.n_spatial; p++) {
//       IdxType pi = spind_to_ind(p, env.n_occ);
//       auto pov = occ_or_vir(p, env.n_occ);
//       FermionOperator virtual_creation_up (pi, pov, Up, Creation, env.xacc_scheme);
//       FermionOperator virtual_creation_down (pi, pov, Down, Creation, env.xacc_scheme);
//       for (IdxType q = p+1; q < env.n_spatial; q++) {
//         IdxType qi = spind_to_ind(q, env.n_occ);
//         auto qov = occ_or_vir(q, env.n_occ);
//         // creation operator
//         FermionOperator occupied_annihilation_up (qi, qov, Up, Annihilation, env.xacc_scheme);
//         FermionOperator occupied_annihilation_down (qi, qov, Down, Annihilation, env.xacc_scheme);
//         fermion_operators.push_back({{occupied_annihilation_up, virtual_creation_up},
//                                     {occupied_annihilation_down, virtual_creation_down}});
//         single_counter += 1;
//       }
//     }
//     //===========Double Excitations===========
//     // Singlets and Triplets
//     int rs = -1;
//     for (IdxType r = 0; r < env.n_spatial; r++) {
//       IdxType ri = spind_to_ind(r, env.n_occ);
//       auto rov = occ_or_vir(r, env.n_occ);
//       FermionOperator ra (ri, rov, Up, Annihilation, env.xacc_scheme);
//       FermionOperator rb (ri, rov, Down, Annihilation, env.xacc_scheme);
//       for (IdxType s = r; s < env.n_spatial; s++) {
//         IdxType si = spind_to_ind(s, env.n_occ);
//         auto sov = occ_or_vir(s, env.n_occ);
//         FermionOperator sa (si, sov, Up, Annihilation, env.xacc_scheme);
//         FermionOperator sb (si, sov, Down, Annihilation, env.xacc_scheme);
//         rs += 1;
//         int pq = -1;
//         for (IdxType p = 0; p < env.n_spatial; p++) {
//           IdxType pi = spind_to_ind(p, env.n_occ);
//           auto pov = occ_or_vir(p, env.n_occ);
//           FermionOperator pa (pi, pov, Up, Creation, env.xacc_scheme);
//           FermionOperator pb (pi, pov, Down, Creation, env.xacc_scheme);
//           for (IdxType q = p; q < env.n_spatial; q++) {
//             IdxType qi = spind_to_ind(q, env.n_occ);
//             auto qov = occ_or_vir(q, env.n_occ);
//             FermionOperator qa (qi, qov, Up, Creation, env.xacc_scheme);
//             FermionOperator qb (qi, qov, Down, Creation, env.xacc_scheme);
//             pq += 1;
//             if (rs > pq) {
//               continue; 
//             }
//             if ( (p == r) && (q == s) ) {
//               continue;
//             }
//             if (p!=q) {
//               if (r == s) { // only singlet with two terms, not sure Group 3 or Group 4 cases
//                 fermion_operators.push_back({
//                   {0.5*rb, ra,qb, pa},
//                   {-0.5*rb, ra,qa, pb}
//                 });
//                 double_counter += 1;
//                 doublt_term2_counter += 1;
//                 continue;
//               }
//               if ( ((r!=s)&&(q!=r)) || ((p==r)&&(q!=s)&(r!=s)) || ((p==s)&&(q!=r)&(r!=s)) || ((q==s)&&(p!=r)) || ((q==r)&&(p!=s))) {
//                 // Trtiplet
//                 fermion_operators.push_back({
//                   {2.0/sqrt(24.0)*sa,ra,qa,pa},
//                   {1.0/sqrt(24.0)*sb,ra,qb,pa},
//                   {1.0/sqrt(24.0)*sa,rb,qb,pa,},
//                   {1.0/sqrt(24.0)*qa,pb,sb,ra,},
//                   {1.0/sqrt(24.0)*qa,pb,sa,rb,},
//                   {2.0/sqrt(24.0)*sb,rb,qb,pb}
//                 });
//                 double_counter += 1;
//                 doublt_term6_counter += 1;
//                 continue;
//               }
//             }
//             // Group 3 to 5
//             if (p == q) {
//               // Group 3
//               if ((q != r)&&(r!=s)) {
//                 // only singlet
//                 fermion_operators.push_back({
//                   {0.5*sb,ra,pb,pa},
//                   {0.5*sa,rb,pb,pa}
//                 });
//                 double_counter += 1;
//                 doublt_term2_counter += 1;
//                 continue;
//               }
//               // Group 4
//               if ((q == r)&&(r!=s)) {
//                 // only singlet
//                 fermion_operators.push_back({
//                   {0.5*sb,ra,pa,pb},
//                   {0.5*sa,rb,pb,pa}
//                 });
//                 double_counter += 1;
//                 doublt_term2_counter += 1;
//                 continue;
//               }
//               // Group 5
//               if ((q != r)&&(r==s)) {
//                 // only singlet
//                 fermion_operators.push_back({
//                   {1.0/sqrt(2.0)*rb,ra,pb,pa}
//                 });
//                 double_counter += 1;
//                 doublt_term1_counter += 1;
//                 continue;
//               }
//             } // p == q

//           } // q
//         } // p
//       } // s
//     } // r

//     int rs_t = -1;
//     for (IdxType r = 0; r < env.n_spatial; r++) {
//       IdxType ri = spind_to_ind(r, env.n_occ);
//       auto rov = occ_or_vir(r, env.n_occ);
//       FermionOperator ra (ri, rov, Up, Annihilation, env.xacc_scheme);
//       FermionOperator rb (ri, rov, Down, Annihilation, env.xacc_scheme);
//       for (IdxType s = r; s < env.n_spatial; s++) {
//         IdxType si = spind_to_ind(s, env.n_occ);
//         auto sov = occ_or_vir(s, env.n_occ);
//         FermionOperator sa (si, sov, Up, Annihilation, env.xacc_scheme);
//         FermionOperator sb (si, sov, Down, Annihilation, env.xacc_scheme);
//         rs_t += 1;
//         int pq_t = -1;
//         for (IdxType p = 0; p < env.n_spatial; p++) {
//           IdxType pi = spind_to_ind(p, env.n_occ);
//           auto pov = occ_or_vir(p, env.n_occ);
//           FermionOperator pa (pi, pov, Up, Creation, env.xacc_scheme);
//           FermionOperator pb (pi, pov, Down, Creation, env.xacc_scheme);
//           for (IdxType q = p; q < env.n_spatial; q++) {
//             IdxType qi = spind_to_ind(q, env.n_occ);
//             auto qov = occ_or_vir(q, env.n_occ);
//             FermionOperator qa (qi, qov, Up, Creation, env.xacc_scheme);
//             FermionOperator qb (qi, qov, Down, Creation, env.xacc_scheme);
//             pq_t += 1;
//             if (rs_t > pq_t) {
//               continue; 
//             }
//             if ( (p == r) && (q == s) ) {
//               continue;
//             }
//             if ((p == q) || (r==s)) {
//               continue;
//             } // singlet cases

//             if ( ((r!=s)&&(q!=r)) || ((p==r)&&(q!=s)&(r!=s)) || ((p==s)&&(q!=r)&(r!=s)) || ((q==s)&&(p!=r)) || ((q==r)&&(p!=s))) {
//                 fermion_operators.push_back({
//                   {0.5/sqrt(2.0)*sb,ra,qb,pa},
//                   {-0.5/sqrt(2.0)*sa,rb,qb,pa},
//                   {-0.5/sqrt(2.0)*sb,ra,qa,pb},
//                   {0.5/sqrt(2.0)*sa,rb,qa,pb}
//                 });
//               double_counter += 1;
//               doublt_term4_counter += 1;
//               continue;
//             }

//           } // q
//         } // p
//       } // s
//     } // r
// };

// // MZ: The problems in the following code is the error on the symmetry.
// void generate_fermionic_excitations_origin(std::vector<std::vector<std::vector<FermionOperator> > >& fermion_operators,
//                                     const MolecularEnvironment& env) {
//     // Single excitation
//     for (IdxType p = 0; p < env.n_occ; p++) {
//       FermionOperator occupied_annihilation_up (p, Occupied, Up, Annihilation, env.xacc_scheme);
//       FermionOperator occupied_annihilation_down (p, Occupied, Down, Annihilation, env.xacc_scheme);
//       for (IdxType q = 0; q < env.n_virt; q++) {
//         FermionOperator virtual_creation_up (q, Virtual, Up, Creation, env.xacc_scheme);
//         FermionOperator virtual_creation_down (q, Virtual, Down, Creation, env.xacc_scheme);
//         fermion_operators.push_back({{occupied_annihilation_up, virtual_creation_up},
//                                       {occupied_annihilation_down, virtual_creation_down}});
//       }
//     }
//     // Double excitation
//     for (IdxType i = 0; i < env.n_occ; i++) {
//       FermionOperator occ_down_1 (i, Occupied, Down, Annihilation, env.xacc_scheme);
//       FermionOperator occ_up_1 (i, Occupied, Up, Annihilation, env.xacc_scheme);
//       for (IdxType j = i+1; j < env.n_occ; j++) {
//         FermionOperator occ_down_2 (j, Occupied, Down, Annihilation, env.xacc_scheme);
//         FermionOperator occ_up_2 (j, Occupied, Up, Annihilation, env.xacc_scheme);
//         for (IdxType r = 0; r < env.n_virt; r++) {
//         FermionOperator virt_down_3 (r, Virtual, Down, Creation, env.xacc_scheme);
//         FermionOperator virt_up_3 (r, Virtual, Up, Creation, env.xacc_scheme);
//           for (IdxType s = r+1; s < env.n_virt; s++) {
//             FermionOperator virt_down_4 (s, Virtual, Down, Creation, env.xacc_scheme);
//             FermionOperator virt_up_4 (s, Virtual, Up, Creation, env.xacc_scheme);
//             IdxType alpha_term = fermion_operators.size();
//             fermion_operators.push_back(
//                {{occ_down_1,
//                 occ_down_2,
//                 virt_down_3,
//                 virt_down_4},
//                {occ_up_1, 
//                 occ_up_2,
//                 virt_up_3,
//                 virt_up_4},
//                 {occ_up_1,
//                  occ_down_2, 
//                  virt_down_3,
//                  virt_up_4},
//                 {occ_down_1,
//                  occ_up_2,
//                  virt_down_3,
//                  virt_up_4}}
//                 );
//           }
//         }
//       }
//     }
//     for (IdxType i = 0; i < env.n_occ; i++) {
//       FermionOperator occ_down_1 (i, Occupied, Down, Annihilation, env.xacc_scheme);
//       FermionOperator occ_up_1 (i, Occupied, Up, Annihilation, env.xacc_scheme);
//         for (IdxType r = 1; r < env.n_virt; r++) {
//         FermionOperator virt_down_2 (r, Virtual, Down, Creation, env.xacc_scheme);
//         FermionOperator virt_up_2 (r, Virtual, Up, Creation, env.xacc_scheme);
//           for (IdxType s = 0; s < r; s++) {
//           FermionOperator virt_down_3 (s, Virtual, Down, Creation, env.xacc_scheme);
//           FermionOperator virt_up_3 (s, Virtual, Up, Creation, env.xacc_scheme);
//             IdxType term = fermion_operators.size();
//             fermion_operators.push_back(
//                  {{occ_up_1,
//                   occ_down_1,
//                   virt_down_2,
//                   virt_up_3},
//                   {occ_up_1,
//                   occ_down_1,
//                   virt_down_3,
//                   virt_up_2}});
//         }
//       }
//     }

//     for (IdxType i = 0; i < env.n_virt; i++) {
//       FermionOperator virt_down_3 (i, Virtual, Down, Creation, env.xacc_scheme);
//       FermionOperator virt_up_3 (i, Virtual, Up, Creation, env.xacc_scheme);
//         for (IdxType r = 0; r < env.n_occ; r++) {
//       FermionOperator occ_down_1 (r, Occupied, Down, Annihilation, env.xacc_scheme);
//       FermionOperator occ_up_1 (r, Occupied, Up, Annihilation, env.xacc_scheme);
//           for (IdxType s = 0; s < r + 1; s++) {
//           FermionOperator occ_down_2 (s, Occupied, Down, Annihilation, env.xacc_scheme);
//           FermionOperator occ_up_2 (s, Occupied, Up, Annihilation, env.xacc_scheme);
//             if (r > s) {
//               fermion_operators.push_back({
//                     {occ_up_1,
//                     occ_down_2,
//                     virt_down_3,
//                     virt_up_3},
//                     {occ_up_2,
//                     occ_down_1,
//                     virt_down_3,
//                     virt_up_3}});
//             } else {
//               fermion_operators.push_back({
//                     {occ_up_1,
//                     occ_down_2,
//                     virt_down_3,
//                     virt_up_3}});
//             }
//         }
//         }
//       }
// };

// /**
//  * @brief  Generate Pauli Operator Pool
//  * @note   Same operators as UCCSD, just with single Pauli Strings. TODO: Remove redundant Paulis
//  * @param  pauli_operators: Output vector of observables
//  * @param  env: Molecular environment structure
//  * @retval None
//  */
// void generate_pauli_excitations(std::vector<std::vector<PauliOperator> >& pauli_operators,
//                                     const MolecularEnvironment& env,
//                                     IdxType subsample,
//                                     IdxType seed) {
//   IdxType n_singles = env.n_occ * env.n_virt;
//   IdxType n_doubles = env.n_occ * (env.n_occ) * env.n_virt * (env.n_virt) +\
//               choose2(env.n_occ) * choose2(env.n_virt) * 2; 
//   std::vector<std::vector<std::vector<FermionOperator> > > fermion_operators;
//   fermion_operators.reserve(n_singles);
//   generate_fermionic_excitations(fermion_operators, env);

//   std::vector<PauliOperator> temporary_storage;
//   if (subsample > 0) {
//     temporary_storage.reserve(4 * n_singles + 16 * n_doubles);
//   } else {
//     pauli_operators.reserve(4 * n_singles + 16 * n_doubles);
//   }
    
//   IdxType index = 0;
//   // NOTE: Not the most efficient way of doing this, kind of a workaround due to the data structures used elsewhere
//   for (size_t i = 0; i < fermion_operators.size(); i++) {
//   // Perform the JordanWigner mapping for each Fermionic operator product
//     std::vector<std::vector<PauliOperator> > mapper_temp;
//     getJordanWignerTransform(env, fermion_operators[i], mapper_temp, true);
//     // Iterate ove the JW mapped Paulis and slap them on the operator pool
//     for (auto pauli_list: mapper_temp) {
//       for (auto pauli: pauli_list) {
//         // If we're not subsampling or if this is a selected index
//         if (subsample < 0) {
//           pauli_operators.push_back({pauli});\
//         } else {
//           temporary_storage.push_back(pauli);
//         }
//         index++; // keep count of indices
//       }
//     }
//   }
//   std::vector<IdxType> selection;
//   // If we're subsampling, generate the list of operator indices to choose a priori
//   if (subsample > 0) {
//     std::vector<IdxType> indices(temporary_storage.size());
//     std::iota(indices.begin(), indices.end(), 0);
//     std::mt19937_64 random_engine (seed);
//     std::shuffle(indices.begin(), indices.end(), random_engine);
//     // take the first N indices, where N is either the number of samples or the number of operators (whichever is smaller)
//     selection = std::vector<IdxType>(
//       indices.begin(),
//       indices.begin() + std::min(subsample, (IdxType)indices.size()));
//     pauli_operators.reserve(selection.size());
//     for (IdxType i : selection) {
//       pauli_operators.push_back({temporary_storage[i]});
//     }
//   }
//   std::vector<IdxType>::iterator iter = selection.begin();

  
// };

//  /**
//   * @brief  Construct the minimal operator pool using the strategy desc. in Appdx. C of Tang et al. 2021
//   * @note   
//   * @param  _pauli_operators: Output vector of Pauli operators (each element is a singleton)
//   * @param  _env: Molecular environment
//   * @retval None
//   */
//   void generate_minimal_pauli_excitations(std::vector<std::vector<PauliOperator > >& _pauli_operators,
//                                   const MolecularEnvironment& _env) {
//     IdxType n_qubits = _env.n_spatial * 2;
//     _pauli_operators.reserve(2 * n_qubits - 2);
//     for (size_t i = 0; i < n_qubits - 1; i++) {
//       // First construct the operator of the form II...Z_{i+1}Y_i...III
//       IdxType xmask_1 = 1 << i;
//       IdxType zmask_1 = (1 << i) + (1 << (i + 1));
//       PauliOperator pauli_1(xmask_1, zmask_1, n_qubits);
//       _pauli_operators.push_back({pauli_1});
//       // Now construct the operator of the form II...Y_{i+1}...III
//       IdxType xmask_2 = (1 << (i + 1));
//       IdxType zmask_2 = (1 << (i + 1));
//       PauliOperator pauli_2(xmask_2, zmask_2, n_qubits);
//       _pauli_operators.push_back({pauli_2});
//     }
//   };
/* --------------------------------------------------------------------------------------------------------*/

 /**
  * @brief  Make a common measurement operator for a QWC group
  * @note   `pauli_list` must consist of QWC operators
  * @param  pauli_list: list of PauliOperators
  * @param  zmasks: output for the observable zmasks to measure each PauliOperator after diagonalization
  * @param  coeffs: coefficients for each Pauli operator (sign corrected) following diagonalization
  * @retval PauliOperator
  */
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
};// namespace VQE
};// namespace NWQSim

std::ostream& operator<<(std::ostream& out, const std::vector<double>& target) {
  out << "[";
  size_t len = target.size();
  if (len > 0) {
      for (size_t i = 0; i < len - 1; i++) {
          out << target[i] << ", ";
      }
      out << target[len-1];
  }
  out << "]";
  return out;
}