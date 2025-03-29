#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>
#include "../include/backendManager.hpp"
#include "../include/state.hpp"
#include "../include/circuit.hpp"
#include "../include/nwq_util.hpp"

// Apply LDPC stabilizers and measure syndromes
void apply_ldpc_stabilizers(std::shared_ptr<NWQSim::Circuit> circuit, int distance) {
    int num_qubits = 2 * distance;
    int num_ancillas = distance;
    
    // Apply stabilizers
    for (int i = 0; i < distance; i++) {
        circuit->CX(i, (i + 2) % num_qubits);
        circuit->H((i + 1) % num_qubits);
        circuit->CX((i + 1) % num_qubits, (i + 3) % num_qubits);
        circuit->H((i + 1) % num_qubits);
    }
    
    // Measure stabilizers using ancilla qubits
    for (int i = 0; i < distance; i++) {
        int ancilla = num_qubits + i;
        circuit->CX(i, ancilla);
        circuit->CX((i + 1) % num_qubits, ancilla);
        circuit->M(ancilla);
    }
}

int main() {
    for (int d = 3; d < 52; d += 2) {
        int distance = d;
        int n_qubits = 2 * distance + distance;  // Data + ancilla qubits
        int shots = 10;
        int rounds = 1;
        auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
        
        // Apply LDPC stabilizers and syndrome measurements
        for (int i = 0; i < rounds; i++) {
            apply_ldpc_stabilizers(circuit, distance);
        }

        std::string backend = "nvgpu";
        std::string sim_method = "stab";
        double timer = 0;

        std::cout << "Creating state" << std::endl;
        auto state = BackendManager::create_state(backend, n_qubits, sim_method);
        
        std::cout << "Starting sim" << std::endl;
        state->sim(circuit, timer);
        
        std::ostringstream filename;
        filename << "/people/garn195/NWQ-Sim/stabilizer/ldpc_code_test/" << backend << "_" << sim_method << "_" << distance << ".txt";
        std::ofstream outfile(filename.str());
        if (!outfile) {
            std::cerr << "Error opening file: " << filename.str() << std::endl;
        }
        
        outfile << "stab" << std::endl;
        outfile << timer / 1000.0 << std::endl;
        outfile << distance << std::endl;
        outfile << rounds << std::endl;
        outfile << n_qubits << std::endl;
    }
    return 0;
}
