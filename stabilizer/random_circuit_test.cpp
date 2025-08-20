#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <regex>

#include "../include/backendManager.hpp"
#include "../include/state.hpp"
#include "../include/circuit.hpp"
#include "../include/nwq_util.hpp"

#include "/people/garn195/NWQ-Sim/stabilizer/src/qasm_extraction.hpp"

int main() {
    std::string folder_path = "/people/garn195/NWQ-Sim/stabilizer/random_bench_dense";
    std::string output_path = "/people/garn195/NWQ-Sim/stabilizer/random_bench_dense/";
    
    std::string backend = "cpu";
    std::string sim_method = "stab";

    // Iterate through all random QASM files
    for(const auto& entry : std::filesystem::directory_iterator(folder_path)) 
    {
        if(entry.is_regular_file() && entry.path().extension() == ".qasm") 
        {
            std::string filename = entry.path().filename().string();
            
            // Check if this is a random circuit file (format: random_q{qubits}_r{rounds}.qasm)
            if(filename.find("random_q") == 0) 
            {
                std::string input = entry.path().string();
                std::cout << "Processing: " << input << std::endl;

                // Extract qubits and rounds from filename
                std::regex re("random_q(\\d+)_r(\\d+).qasm");
                std::smatch match;
                int n_qubits = 0;
                int rounds = 0;
                if (std::regex_search(filename, match, re) && match.size() > 2) {
                    n_qubits = std::stoi(match.str(1));
                    rounds = std::stoi(match.str(2));
                } else {
                    std::cerr << "Could not parse qubits and rounds from filename: " << filename << std::endl;
                    continue;
                }

                auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
                if(appendQASMToCircuit(circuit, input, n_qubits))
                {
                    double timer = 0;
                    auto state = BackendManager::create_state(backend, n_qubits, sim_method);
                    state->allocate_measurement_buffers(1000000);
                    state->sim(circuit, timer);

                    std::cout << "Sim time: " << timer/1000.0 << "s" << std::endl;

                    // Write benchmark results
                    std::ostringstream output_filename;
                    output_filename << output_path << backend << "_" << sim_method 
                                  << "_random_q" << n_qubits << "_r" << rounds << ".txt";
                    
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
                    outfile << rounds << std::endl;
                    
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
