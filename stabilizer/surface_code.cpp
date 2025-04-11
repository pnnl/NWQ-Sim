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
    int grid_size = (2 * distance) - 1;
    
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            if ((i%2 == 1) && (j%2 == 0))
            {
                int ancilla = qubit_index(i, j, grid_size);
                
                circuit->H(ancilla); 

                if (i+1 < grid_size) circuit->CX(ancilla, qubit_index(i+1, j, grid_size));
                if (i-1 >= 0) circuit->CX(ancilla, qubit_index(i-1, j, grid_size));
                if (j+1 < grid_size) circuit->CX(ancilla, qubit_index(i, j+1, grid_size));
                if (j-1 >= 0) circuit->CX(ancilla, qubit_index(i, j-1, grid_size));

                circuit->H(ancilla); 
                circuit->M(ancilla); 
            }
        }
    }
}

void measure_z_stabilizers(std::shared_ptr<NWQSim::Circuit> circuit, int distance) {
    int grid_size = (2 * distance) - 1;
    
    for (int i = 0; i < grid_size; i++) {
        for (int j = 0; j < grid_size; j++) {
            if ((i%2 == 0) && (j%2 == 1))
            {
                int ancilla = qubit_index(i, j, grid_size);

                if (j+1 < grid_size) circuit->CX(qubit_index(i, j+1, grid_size), ancilla);
                if (j-1 >= 0) circuit->CX(qubit_index(i, j-1, grid_size), ancilla);
                if (i+1 < grid_size) circuit->CX(qubit_index(i+1, j, grid_size), ancilla);
                if (i-1 >= 0) circuit->CX(qubit_index(i-1, j, grid_size), ancilla);

                circuit->M(ancilla); 
            }
        }
    }
}


int main()
{
    for(int d = 81; d < 153; d+=2)
    {
    int distance = d;
    int n_qubits = pow((2 * distance) - 1, 2);
    int shots = 10;
    int rounds = 1;
    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
    
    //Add surface code routines to the circuit
    for(int i = 0; i < rounds; i++)
    {
        measure_x_stabilizers(circuit, distance);
        measure_z_stabilizers(circuit, distance);
    }

    
    std::string backend = "cpu";
    std::string sim_method = "stab";
    double timer = 0;
    
    /*Create T and Measurement Tableaus with only stabilizers. T starts empty, M starts as identity.*/
    std::cout << "Creating state" << std::endl;
    auto state = BackendManager::create_state(backend, n_qubits, sim_method);

    
    // std::vector<std::shared_ptr<Circuit>> circuit2D = {circuit, circuit};

    std::cout << "Starting sim" << std::endl;

    state->sim(circuit, timer);
    // state->simBitwise(circuit, timer);
    // state->print_res_state();
    // NWQSim::IdxType *results = state->measure_all(shots);
    // for(int i = 0; i < shots; i++)
    //     std::cout << "Result " << i << ": " << results[i] << std::endl;


    std::ostringstream filename;
    filename << "/people/garn195/NWQ-Sim/stabilizer/fowler_surface_code/" << backend << "_" << sim_method << "_" << distance << ".txt";
    std::ofstream outfile(filename.str());
    if (!outfile) {
        std::cerr << "Error opening file: " << filename.str() << std::endl;
    }

    outfile << "cpu" << std::endl;
    outfile << timer/1000.0 << std::endl;
    outfile << distance << std::endl;
    outfile << rounds << std::endl;
    outfile << n_qubits << std::endl;
    
    }

    return 0;
}