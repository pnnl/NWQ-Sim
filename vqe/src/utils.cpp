#include "observable/pauli_operator.hpp"
#include "utils.hpp"
#include <string>
#include <iostream>
#include <sstream>
#include <queue>

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