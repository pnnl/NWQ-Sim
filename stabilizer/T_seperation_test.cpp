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

#include "src/T_separation.hpp"
#include "src/qasm_extraction.hpp"

namespace fs = std::filesystem;

int main()
{
    std::string inputFolder = "/people/garn195/NWQ-Sim/stabilizer/qft_sk";//Input folder path
    std::string outputFolder = "/people/garn195/NWQ-Sim/stabilizer/QFT_transpilation_data";//Output folder path

    // Create the output folder if it doesn't exist
    fs::create_directories(outputFolder);

    // Loop through all files in the input folder
    for (const auto& entry : fs::directory_iterator(inputFolder)) 
    {
        if (entry.is_regular_file()) 
        {
            fs::path inputPath = entry.path();
            std::string filename = inputPath.filename().string();

            fs::path outputPath = fs::path(outputFolder) / (inputPath.stem().string() + ".txt");

            // Open the input and output files
            std::ifstream infile(inputPath);
            std::ofstream outfile(outputPath);

            if (!infile.is_open() || !outfile.is_open()) 
            {
                std::cerr << "Could not open " << inputPath << " or " << outputPath << std::endl;
                continue;
            } 
            /* Small test */
            int n_qubits = 4;
            auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);

            appendQASMToCircuit(circuit, inputPath, n_qubits);
            NWQSim::IdxType total_gates = (circuit->get_gates()).size();


            //Measurement circuit will be filled in the passthrough function
            auto M_circ = std::make_shared<NWQSim::Circuit>(n_qubits);


            std::string backend = "CPU";
            std::string sim_method = "stab";
            double timer = 0;
            
            /*Create T and Measurement Tableaus with only stabilizers. T starts empty, M starts as identity.*/
            auto T_tab = BackendManager::create_state(backend, n_qubits, sim_method);
            T_tab->delete_all_rows();
            T_tab->remove_destabilizers();

            auto M_tab = BackendManager::create_state(backend, n_qubits, "stab");
            std::chrono::duration<long long, std::ratio<1, 1000000>> proc_time;
            T_process(circuit, M_circ, T_tab, proc_time);
            /*Run T passthrough after tableaus have been prepared and measurement operations absorbed*/
            M_tab->sim(M_circ, timer);
            M_tab->remove_destabilizers();
            // M_tab->print_res_state();
            T_passthrough(T_tab, M_tab, outfile, proc_time, 10);


            //Put the M circuit back in forward time after the new Clifford gates have been appended
            std::cout << "---- T tableau -----" << std::endl;

            // T_tab->print_res_state();
            std::cout << "---- M tableau -----" << std::endl;

            // M_tab->print_res_state();
            // std::cout << "M tableau measurement results: " << (M_tab->measure_all(10))[0] << std::endl;
            /*Measurement Tableau*/

            outfile << "Total gates: " << total_gates << std::endl;
            outfile.close(); // Close the file


        }
    }
    return 0;
}