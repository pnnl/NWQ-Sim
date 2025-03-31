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
void measure_x_stabilizers(std::shared_ptr<NWQSim::Circuit> circuit, std::shared_ptr<NWQSim::Circuit> circuit2D, std::vector<int>& gate_chunks, int distance) {
    int grid_size = distance + 1;
    std::vector<int> ancillas;

    gate_chunks.push_back(0);    
    // Apply all Hadamards first
    for (int i = 0; i < distance; i++) {
        for (int j = 0; j < distance; j++) {
            if ((i + j) % 2 == 1) { // X-stabilizer condition
                int ancilla = qubit_index(i, j, grid_size);
                circuit2D->H(ancilla);  // Prepare ancilla in |+>
                ancillas.push_back(ancilla);
                gate_chunks.back()++;
            }
        }
    }
    
    gate_chunks.push_back(0);    
    // Apply CNOTs in separate loops
    for (int ancilla : ancillas) {
        int i = ancilla / grid_size;
        int j = ancilla % grid_size;
        if (i+1 < grid_size)
        {
            circuit2D->CX(ancilla, qubit_index(i+1, j, grid_size));
            gate_chunks.back()++;
        }
    }
    
    gate_chunks.push_back(0);    
    for (int ancilla : ancillas) {
        int i = ancilla / grid_size;
        int j = ancilla % grid_size;
        if (j+1 < grid_size){
            circuit2D->CX(ancilla, qubit_index(i, j+1, grid_size));
            gate_chunks.back()++;
        }
    }

    gate_chunks.push_back(0);    
    for (int ancilla : ancillas) {
        int i = ancilla / grid_size;
        int j = ancilla % grid_size;
        if (i-1 >= 0){
            circuit2D->CX(ancilla, qubit_index(i-1, j, grid_size));
            gate_chunks.back()++;
        }
    }

    gate_chunks.push_back(0);
    for (int ancilla : ancillas) {
        int i = ancilla / grid_size;
        int j = ancilla % grid_size;
        if (j-1 >= 0){
            circuit2D->CX(ancilla, qubit_index(i, j-1, grid_size));
            gate_chunks.back()++;
        }
    }
    
    gate_chunks.push_back(0);    
    // Apply all Hadamards again
    for (int ancilla : ancillas) {
        circuit2D->H(ancilla);
        gate_chunks.back()++;
    }

    
    // // Measure all ancillas
    // for (int ancilla : ancillas) {
    //     circuit->M(ancilla);
    // }
}

// Z stabilizer measurements (alternating plaquettes)
void measure_z_stabilizers(std::shared_ptr<NWQSim::Circuit> circuit, std::shared_ptr<NWQSim::Circuit> circuit2D, std::vector<int>& gate_chunks, int distance) {
    int grid_size = distance + 1;
    std::vector<int> ancillas;

    gate_chunks.push_back(0);
    // Collect all ancillas
    for (int i = 0; i < distance; i++) {
        for (int j = 0; j < distance; j++) {
            if ((i + j) % 2 == 0) { // Z-stabilizer condition
                int ancilla = qubit_index(i, j, grid_size);
                ancillas.push_back(ancilla);
            }
        }
    }

    gate_chunks.push_back(0);
    // Apply CNOTs in separate loops
    for (int ancilla : ancillas) {
        int i = ancilla / grid_size;
        int j = ancilla % grid_size;
        if (i+1 < grid_size){
            circuit2D->CX(qubit_index(i+1, j, grid_size), ancilla);
            gate_chunks.back()++;
        }
    }

    gate_chunks.push_back(0);
    for (int ancilla : ancillas) {
        int i = ancilla / grid_size;
        int j = ancilla % grid_size;
        if (j+1 < grid_size){
            circuit2D->CX(qubit_index(i, j+1, grid_size), ancilla);
            gate_chunks.back()++;
        }
    }

    gate_chunks.push_back(0);
    for (int ancilla : ancillas) {
        int i = ancilla / grid_size;
        int j = ancilla % grid_size;
        if (i-1 >= 0){
            circuit2D->CX(qubit_index(i-1, j, grid_size), ancilla);
            gate_chunks.back()++;
        }
    }

    gate_chunks.push_back(0);
    for (int ancilla : ancillas) {
        int i = ancilla / grid_size;
        int j = ancilla % grid_size;
        if (j-1 >= 0){
            circuit2D->CX(qubit_index(i, j-1, grid_size), ancilla);
            gate_chunks.back()++;
        }
    }
    
    // // Measure all ancillas
    // for (int ancilla : ancillas) {
    //     circuit->M(ancilla);
    // }
}


int main()
{
    for(int d = 3; d < 100; d+=2)
    {
    int distance = d;
    int n_qubits = 2 * pow(distance, 2) + 1;
    int shots = 10;
    int rounds = 1;
    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
    auto circuit2D = std::make_shared<NWQSim::Circuit>(n_qubits);

    std::string backend = "nvgpu";
    std::string sim_method = "stab";
    double timer = 0;
    
    /*Create T and Measurement Tableaus with only stabilizers. T starts empty, M starts as identity.*/
    std::cout << "Creating state" << std::endl;
    auto state = BackendManager::create_state(backend, n_qubits, sim_method);
    std::vector<int> gate_chunks;
    
    //Add surface code routines to the circuit
    for(int i = 0; i < rounds; i++)
    {
        measure_x_stabilizers(circuit, circuit2D, gate_chunks, distance);
        state->sim2D(circuit2D, gate_chunks, timer);
        // state->sim(circuit, timer);

        // circuit->clear();
        circuit2D->clear();

        gate_chunks.clear();


        measure_z_stabilizers(circuit, circuit2D, gate_chunks, distance);
        state->sim2D(circuit2D, gate_chunks, timer);
        // state->sim(circuit, timer);
        // circuit->clear();
        circuit2D->clear();

        gate_chunks.clear();



        // Test for desync
        // for(int n = 0; n < n_qubits; n+=2)
        // {
        //     circuit->H(n);
        //     circuit->S(n);
        //     circuit->S(n);
        //     circuit->H(n);
        //     circuit->M(n);
        // }
        // for(int n = 1; n < n_qubits; n+=2)
        // {
        //     circuit->H(n);
        //     circuit->M(n);
        // }
    }

    // state->simBitwise(circuit, timer);
    // state->print_res_state();
    // NWQSim::IdxType *results = state->measure_all(shots);
    // for(int i = 0; i < shots; i++)
    //     std::cout << "Result " << i << ": " << results[i] << std::endl;


    std::ostringstream filename;
    filename << "/people/garn195/NWQ-Sim/stabilizer/initialize_code/" << backend << "2D_" << sim_method << "_" << distance << ".txt";
    std::ofstream outfile(filename.str());
    if (!outfile) {
        std::cerr << "Error opening file: " << filename.str() << std::endl;
    }

    outfile << "stab2D" << std::endl;
    outfile << timer/1000.0 << std::endl;
    outfile << distance << std::endl;
    outfile << rounds << std::endl;
    outfile << n_qubits << std::endl;
    
    }

    return 0;
}