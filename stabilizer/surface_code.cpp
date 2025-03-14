#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../include/backendManager.hpp"
#include "../include/state.hpp"
#include "../include/circuit.hpp"
#include "../include/nwq_util.hpp"

// Convert 2D (i, j) grid index to 1D qubit index
inline int qubit_index(int i, int j, int grid_size) {
    return i * grid_size + j;
}

// X stabilizer measurements (alternating plaquettes)
void measure_x_stabilizers(std::shared_ptr<NWQSim::Circuit> circuit, int distance) {
    int grid_size = distance + 1;
    
    for (int i = 0; i < distance; i++) {
        for (int j = 0; j < distance; j++) {
            if ((i + j) % 2 == 1) { // X-stabilizer condition
                int ancilla = qubit_index(i, j, grid_size);
                
                circuit->H(ancilla);  // Prepare ancilla in |+>

                // Apply CNOTs with data neighbors
                if (i+1 < grid_size) circuit->CX(ancilla, qubit_index(i+1, j, grid_size));
                if (j+1 < grid_size) circuit->CX(ancilla, qubit_index(i, j+1, grid_size));
                if (i-1 >= 0) circuit->CX(ancilla, qubit_index(i-1, j, grid_size));
                if (j-1 >= 0) circuit->CX(ancilla, qubit_index(i, j-1, grid_size));

                circuit->H(ancilla);  // Convert back to computational basis
                circuit->M(ancilla);  // Measure stabilizer
            }
        }
    }
}

// Z stabilizer measurements (alternating plaquettes)
void measure_z_stabilizers(std::shared_ptr<NWQSim::Circuit> circuit, int distance) {
    int grid_size = distance + 1;
    
    for (int i = 0; i < distance; i++) {
        for (int j = 0; j < distance; j++) {
            if ((i + j) % 2 == 0) { // Z-stabilizer condition
                int ancilla = qubit_index(i, j, grid_size);

                // Apply CNOTs with data neighbors
                if (i+1 < grid_size) circuit->CX(qubit_index(i+1, j, grid_size), ancilla);
                if (j+1 < grid_size) circuit->CX(qubit_index(i, j+1, grid_size), ancilla);
                if (i-1 >= 0) circuit->CX(qubit_index(i-1, j, grid_size), ancilla);
                if (j-1 >= 0) circuit->CX(qubit_index(i, j-1, grid_size), ancilla);

                circuit->M(ancilla);  // Measure stabilizer
            }
        }
    }
}


int main()
{



    // Main Circuit
    measure_x_stabilizers();
    measure_z_stabilizers();

    int n_qubits;
    int shots = 10;

    int rounds = 1;
    NWQSim::IdxType S_count = 0;
    NWQSim::IdxType H_count = rounds * n_qubits;
    NWQSim::IdxType CX_count = 0;
    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);

    std::string backend = "nvgpu";
    std::string sim_method = "stab";
    double timer = 0;
    
    /*Create T and Measurement Tableaus with only stabilizers. T starts empty, M starts as identity.*/
    std::cout << "Creating state" << std::endl;
    auto state = BackendManager::create_state(backend, n_qubits, sim_method);

    
    // std::vector<std::shared_ptr<Circuit>> circuit2D = {circuit, circuit};

    std::cout << "Starting sim" << std::endl;

    // state->sim(circuit, timer);
    state->sim(circuit, timer);
    state->print_res_state();
    return 0;
}