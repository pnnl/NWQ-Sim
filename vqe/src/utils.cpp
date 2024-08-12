#include "observable/pauli_operator.hpp"
#include "observable/fermionic_operator.hpp"
#include "utils.hpp"
#include <string>
#include <iostream>
#include <sstream>
#include <queue>
#include <fstream>
#include <vector>

namespace NWQSim{
namespace VQE {
  std::ostream& operator<<(std::ostream& out, const PauliOperator& op) {
        out << op.pauliToString();
        return out;
      };
  std::string to_binary_string(NWQSim::IdxType val, NWQSim::IdxType n_qubits) {
    std::stringstream ss;
    for (IdxType i = n_qubits-1; i >= 0; i--) {
      ss << ((val & (1 << i)) ? "1" : "0");
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
}
void generate_fermionic_excitations(std::vector<std::vector<std::vector<FermionOperator> > >& fermion_operators,
                                    const MolecularEnvironment& env) {
  // Single excitation
      for (IdxType p = 0; p < env.n_occ; p++) {
        FermionOperator occupied_annihilation_up (p, Occupied, Up, Annihilation, env.xacc_scheme);
        FermionOperator occupied_annihilation_down (p, Occupied, Down, Annihilation, env.xacc_scheme);
        for (IdxType q = 0; q < env.n_virt; q++) {
          FermionOperator virtual_creation_up (q, Virtual, Up, Creation, env.xacc_scheme);
          FermionOperator virtual_creation_down (q, Virtual, Down, Creation, env.xacc_scheme);
          fermion_operators.push_back({{occupied_annihilation_up, virtual_creation_up}, {occupied_annihilation_down, virtual_creation_down}});
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
            fermion_operators.push_back({{
                  occ_down_1,
                  occ_down_2,
                  virt_down_2,
                  virt_down_1},
                  {occ_up_1,
                   occ_up_2,
                   virt_up_2,
                   virt_up_1},
                   {
                  occ_up_1,
                  occ_down_2,
                  virt_down_2,
                  virt_up_1},
                  {
                  occ_down_1 * -1.0,
                  occ_up_2,
                  virt_down_2,
                  virt_up_1}
                  });
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
            if (i != j && r != s) {
              fermion_operators.push_back({
                    {occ_up_1,
                    occ_down_2,
                    virt_down_2,
                    virt_up_1},
                    {occ_down_1,
                    occ_up_2,
                    virt_up_2,
                    virt_down_1}});
            } else {
              fermion_operators.push_back({
                    {occ_down_1,
                    occ_up_2,
                    virt_up_2,
                    virt_down_1}});

            }
          }
            
        }
      }
    }
};

void generate_pauli_excitations(std::vector<std::vector<PauliOperator> >& pauli_operators,
                                    const MolecularEnvironment& env) {
  IdxType n_singles = env.n_occ * env.n_virt;
  IdxType n_doubles = env.n_occ * env.n_virt;
  std::vector<std::vector<std::vector<FermionOperator> > > fermion_operators;
  fermion_operators.reserve(n_singles);
  // Single excitation
      for (IdxType p = 0; p < env.n_occ; p++) {
        FermionOperator occupied_annihilation_up (p, Occupied, Up, Annihilation, env.xacc_scheme);
        FermionOperator occupied_annihilation_down (p, Occupied, Down, Annihilation, env.xacc_scheme);
        for (IdxType q = 0; q < env.n_virt; q++) {
          FermionOperator virtual_creation_up (q, Virtual, Up, Creation, env.xacc_scheme);
          FermionOperator virtual_creation_down (q, Virtual, Down, Creation, env.xacc_scheme);
          fermion_operators.push_back({{occupied_annihilation_up, virtual_creation_up}, {occupied_annihilation_down, virtual_creation_down}});
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
            fermion_operators.push_back({{
                  occ_down_1,
                  occ_down_2,
                  virt_down_2,
                  virt_down_1},
                  {occ_up_1,
                   occ_up_2,
                   virt_up_2,
                   virt_up_1},
                   {
                  occ_up_1,
                  occ_down_2,
                  virt_down_2,
                  virt_up_1},
                  {
                  occ_down_1 * -1.0,
                  occ_up_2,
                  virt_down_2,
                  virt_up_1}
                  });
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
            if (i != j && r != s) {
              fermion_operators.push_back({
                    {occ_up_1,
                    occ_down_2,
                    virt_down_2,
                    virt_up_1},
                    {occ_down_1,
                    occ_up_2,
                    virt_up_2,
                    virt_down_1}});
            } else {
              fermion_operators.push_back({
                    {occ_down_1,
                    occ_up_2,
                    virt_up_2,
                    virt_down_1}});

            }
          }
            
        }
      }
    }
};
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