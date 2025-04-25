#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <memory>
#include <cmath>
#include <sstream>
#include "../include/backendManager.hpp"
#include "../include/state.hpp"
#include "../include/circuit.hpp"
#include "../include/nwq_util.hpp"

using Coord = std::pair<int, int>;
using Triangle = std::vector<Coord>;

// Hex lattice directions (for 3 edges around a vertex)
const std::vector<Coord> hex_neighbors = {
    {1, 0}, {0, 1}, {-1, 1},
    {-1, 0}, {0, -1}, {1, -1}
};

// Triangles made from 3 connected neighbors
std::vector<Triangle> generate_triangles(const std::map<Coord, int>& qubits) {
    std::vector<Triangle> triangles;
    for (const auto& [coord, idx] : qubits) {
        int x = coord.first, y = coord.second;
        for (int i = 0; i < 6; i++) {
            Coord a = {x + hex_neighbors[i].first, y + hex_neighbors[i].second};
            Coord b = {x + hex_neighbors[(i + 1) % 6].first, y + hex_neighbors[(i + 1) % 6].second};

            if (qubits.count(a) && qubits.count(b)) {
                // Avoid duplicates by sorting the triangle
                std::vector<Coord> tri = {coord, a, b};
                std::sort(tri.begin(), tri.end());
                if (std::find(triangles.begin(), triangles.end(), tri) == triangles.end()) {
                    triangles.push_back(tri);
                }
            }
        }
    }
    return triangles;
}

void add_stabilizers(std::shared_ptr<NWQSim::Circuit> circuit,
                     const std::vector<Triangle>& triangles,
                     const std::map<Coord, int>& qubit_coords,
                     int& ancilla_start) {
    for (const auto& tri : triangles) {
        std::vector<int> q;
        for (const auto& c : tri) q.push_back(qubit_coords.at(c));

        int x_anc = ancilla_start++;
        circuit->H(x_anc);
        for (int qi : q) circuit->CX(x_anc, qi);
        circuit->H(x_anc);
        circuit->M(x_anc);

        int z_anc = ancilla_start++;
        for (int qi : q) circuit->CX(qi, z_anc);
        circuit->M(z_anc);
    }
}

int main() {
    for (int d = 3; d <= 9; d += 2) {
        std::map<Coord, int> qubit_coords;
        int radius = d;
        int q_index = 0;

        // Create hexagonal lattice within given radius
        for (int x = -radius; x <= radius; x++) {
            for (int y = -radius; y <= radius; y++) {
                if (std::abs(x + y) > radius) continue;
                qubit_coords[{x, y}] = q_index++;
            }
        }

        auto triangles = generate_triangles(qubit_coords);
        int n_data = qubit_coords.size();
        int n_anc = 2 * triangles.size();
        int n_qubits = n_data + n_anc;

        auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
        int ancilla_start = n_data;

        add_stabilizers(circuit, triangles, qubit_coords, ancilla_start);

        double timer = 0;
        auto state = BackendManager::create_state("cpu", n_qubits, "stab");
        state->sim(circuit, timer);

        std::ostringstream filename;
        filename << "/people/garn195/NWQ-Sim/stabilizer/hex_color_bench/" << d << ".txt";
        std::ofstream outfile(filename.str());
        outfile << "cpu\n" << timer / 1000.0 << "\n" << d << "\n"
                << triangles.size() << "\n" << n_qubits << "\n";
    }

    return 0;
}
