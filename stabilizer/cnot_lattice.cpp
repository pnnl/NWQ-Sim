#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>
#include <filesystem>
#include <fstream>

#include "../include/backendManager.hpp"
#include "../include/state.hpp"
#include "../include/circuit.hpp"
#include "../include/nwq_util.hpp"

// #include "src/T_separation.hpp"
#include "src/qasm_extraction.hpp"

int main() {
    std::string folder_path = "/Users/garn195/Project Repositories/NWQ-Sim/stabilizer/stim_to_qasm_files";
    int n_qubits = 4;

    std::string backend = "cpu";
    std::string sim_method = "stab";

    for(const auto& entry : std::filesystem::directory_iterator(folder_path)) 
    {
        if(entry.is_regular_file() && entry.path().extension() == ".qasm") 
        {
            std::string input = entry.path().string();
            std::cout << "Processing: " << input << std::endl;

            auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
            if(appendQASMToCircuit(circuit, input, n_qubits))
            {
                // std::vector<NWQSim::Gate> gates = circuit->get_gates();
                // std::cout << circuit->to_string() << std::endl;

                double timer = 0;
                auto state = BackendManager::create_state(backend, n_qubits, sim_method);
                state->sim(circuit, timer);

                std::cout << "Sim time: " << timer/1000.0 << "s" << std::endl;

                std::string name = "";
                std::ostringstream filename;
                filename << "/Users/garn195/Project Repositories/NWQ-Sim/stabilizer/surface_operation_bench/" << backend << "_" << sim_method << "_" << n_qubits << ".txt";
                std::ofstream outfile(filename.str());
                if (!outfile) 
                {
                    std::cerr << "Error opening file: " << filename.str() << std::endl;
                }
                
                outfile << backend << std::endl;
                outfile << timer/1000.0 << std::endl;
                outfile << n_qubits << std::endl;
                outfile << circuit->num_gates() << std::endl;
            }
        }
    }
    return 0;
}