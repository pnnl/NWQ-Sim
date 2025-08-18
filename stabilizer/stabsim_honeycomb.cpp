#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "../include/backendManager.hpp"
#include "../include/state.hpp"
#include "../include/circuit.hpp"
#include "../include/nwq_util.hpp"

#include "src/qasm_extraction.hpp"

int main() {
    std::string folder_path = "/Users/garn195/Project Repositories/NWQ-Sim/stabilizer/surface_operation_bench/";
    std::string output_path = "/Users/garn195/Project Repositories/NWQ-Sim/stabilizer/surface_operation_bench/";
    
    std::string backend = "nvgpu";
    std::string sim_method = "stab";

    // Iterate through all honeycomb QASM files
    for(const auto& entry : std::filesystem::directory_iterator(folder_path)) 
    {
        if(entry.is_regular_file() && entry.path().extension() == ".qasm") 
        {
            std::string filename = entry.path().filename().string();
            
            // Check if this is a honeycomb file (format: honeycomb_d{distance}.qasm)
            if(filename.find("honeycomb_d") == 0) 
            {
                std::string input = entry.path().string();
                std::cout << "Processing: " << input << std::endl;

                // Extract distance from filename
                // Format: honeycomb_d{distance}.qasm
                std::string base_name = filename.substr(0, filename.find_last_of("."));
                size_t d_pos = base_name.find("_d") + 2;
                std::string distance_str = base_name.substr(d_pos);
                int distance = std::stoi(distance_str);

                // Count qubits from QASM file
                std::ifstream qasm_file(input);
                std::string line;
                int n_qubits = 0;
                
                while(std::getline(qasm_file, line)) {
                    if(line.find("qreg") != std::string::npos) {
                        size_t bracket_pos = line.find("[");
                        size_t close_bracket = line.find("]");
                        if(bracket_pos != std::string::npos && close_bracket != std::string::npos) {
                            std::string qubit_count = line.substr(bracket_pos + 1, close_bracket - bracket_pos - 1);
                            n_qubits = std::stoi(qubit_count);
                            break;
                        }
                    }
                }
                qasm_file.close();

                auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
                if(appendQASMToCircuit(circuit, input, n_qubits))
                {
                    double timer = 0;
                    auto state = BackendManager::create_state(backend, n_qubits, sim_method);
                    state->sim(circuit, timer);

                    std::cout << "Sim time: " << timer/1000.0 << "s" << std::endl;

                    // Write benchmark results
                    std::ostringstream output_filename;
                    output_filename << output_path << backend << "_" << sim_method 
                                  << "_honeycomb_d" << distance << ".txt";
                    
                    std::ofstream outfile(output_filename.str());
                    if (!outfile) 
                    {
                        std::cerr << "Error opening file: " << output_filename.str() << std::endl;
                        continue;
                    }
                    
                    outfile << backend << std::endl;
                    outfile << timer/1000.0 << std::endl;
                    outfile << n_qubits << std::endl;
                    outfile << circuit->num_gates() << std::endl;
                    outfile << distance << std::endl;
                    
                    outfile.close();
                }
                else
                {
                    std::cerr << "Failed to load QASM file: " << input << std::endl;
                }
            }
        }
    }
    return 0;
}