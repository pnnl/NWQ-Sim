#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../include/backendManager.hpp"
#include "../include/state.hpp"
#include "../include/circuit.hpp"
#include "../include/nwq_util.hpp"

int extractQubitIndex(const std::string& qubitStr) 
{
    std::regex qubitRegex("q\\w*\\[(\\d+)\\]");
    std::smatch match;
    if (std::regex_search(qubitStr, match, qubitRegex) && match.size() > 1) 
    {
        return std::stoi(match.str(1));
    }
    else
    {
        std::cerr << "Gate called but no qubit index!" << std::endl;
        return -1;
    }
    
}
int extractNumQubit(const std::string& qubitStr) 
{
    size_t start = qubitStr.find('[');
    size_t end = qubitStr.find(']');
    
    if (start != std::string::npos && end != std::string::npos && start + 1 < end) {
        return std::stoi(qubitStr.substr(start + 1, end - start - 1));
    }
    
    return -1; // Return -1 if no valid integer is found
}

void appendQASMToCircuit(std::shared_ptr<NWQSim::Circuit>& circuit, const std::string& filename, int& n_qubits) 
{
    std::ifstream file(filename);
    if (!file.is_open()) 
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    int tCount = 0;
    while (std::getline(file, line)) 
    {
        std::string gate;
        int qubit1, qubit2;

        if (line.empty() || line.find("include") != std::string::npos || line.find("gate") != std::string::npos || line.find("barrier") != std::string::npos || line.find("OPENQASM") != std::string::npos || line.find("creg") != std::string::npos || line.find("#include") != std::string::npos)
            continue;

        std::istringstream lineStream(line);
        lineStream >> gate;

        if(line.find("qreg") != std::string::npos)
        {
            n_qubits = extractNumQubit(line);
            circuit->set_num_qubits(n_qubits);
            // std::cout << n_qubits << std::endl;
        }
        else if(gate == "tdg")
        {
            tCount++;
            std::string qubitStr;
            lineStream >> qubitStr;
            qubit1 = extractQubitIndex(qubitStr);
            if (qubit1 != -1) 
            {
                circuit->TDG(qubit1);
            }
            // std::cout << gate << qubit1 << ", " << qubit2 << std::endl;
        }
        else if(gate == "h")
        {
            std::string qubitStr;
            lineStream >> qubitStr;
            qubit1 = extractQubitIndex(qubitStr);
            if (qubit1 != -1) 
            {
                circuit->H(qubit1);
            }
            // std::cout << gate << qubit1 << ", " << qubit2 << std::endl;
        }
        else if(gate == "s")
        {
            std::string qubitStr;
            lineStream >> qubitStr;
            qubit1 = extractQubitIndex(qubitStr);
            if (qubit1 != -1) 
            {
                circuit->S(qubit1);
            }
            // std::cout << gate << qubit1 << ", " << qubit2 << std::endl;
        }
        else if(gate == "sdg")
        {
            std::string qubitStr;
            lineStream >> qubitStr;
            qubit1 = extractQubitIndex(qubitStr);
            if (qubit1 != -1) 
            {
                circuit->SDG(qubit1);
            }
            // std::cout << gate << qubit1 << ", " << qubit2 << std::endl;
        }
        else if(gate == "t")
        {
            tCount++;

            std::string qubitStr;
            lineStream >> qubitStr;
            qubit1 = extractQubitIndex(qubitStr);
            if (qubit1 != -1) 
            {
                circuit->T(qubit1);
            }
            // std::cout << gate << qubit1 << ", " << qubit2 << std::endl;
        }
        else if(gate == "m")
        {
            std::string qubitStr;
            lineStream >> qubitStr;
            qubit1 = extractQubitIndex(qubitStr);
            if (qubit1 != -1) 
            {
                circuit->M(qubit1);
            }
            // std::cout << gate << qubit1 << ", " << qubit2 << std::endl;
        }
        else if(gate == "reset")
        {
            std::string qubitStr;
            lineStream >> qubitStr;
            qubit1 = extractQubitIndex(qubitStr);
            if (qubit1 != -1) 
            {
                circuit->RESET(qubit1);
            }
            // std::cout << gate << qubit1 << ", " << qubit2 << std::endl;
        }
        else if(gate == "cx")
        {
            std::string qubitStr1, qubitStr2;
            //Read the full line after cx and split it by the comma
            std::getline(lineStream, qubitStr1, ',');
            std::getline(lineStream, qubitStr2);

            qubitStr1.erase(0, qubitStr1.find_first_not_of(" \t"));
            qubitStr1.erase(qubitStr1.find_last_not_of(" \t") + 1);
            qubitStr2.erase(0, qubitStr2.find_first_not_of(" \t"));
            qubitStr2.erase(qubitStr2.find_last_not_of(" \t") + 1);

            qubit1 = extractQubitIndex(qubitStr1);
            qubit2 = extractQubitIndex(qubitStr2);

            if (qubit1 != -1 && qubit2 != -1) 
            {
                circuit->CX(qubit1, qubit2);
            }
            // std::cout << gate << qubit1 << ", " << qubit2 << std::endl;
        }
        else if(gate == "cxyz")
        {
            std::string qubitStr1, qubitStr2;
            //Read the full line after cxyz and split it by the comma
            std::getline(lineStream, qubitStr1, ',');
            std::getline(lineStream, qubitStr2);

            qubitStr1.erase(0, qubitStr1.find_first_not_of(" \t"));
            qubitStr1.erase(qubitStr1.find_last_not_of(" \t") + 1);
            qubitStr2.erase(0, qubitStr2.find_first_not_of(" \t"));
            qubitStr2.erase(qubitStr2.find_last_not_of(" \t") + 1);

            qubit1 = extractQubitIndex(qubitStr1);
            qubit2 = extractQubitIndex(qubitStr2);

            if (qubit1 != -1 && qubit2 != -1) 
            {
                circuit->CX(qubit1, qubit2);
                circuit->CY(qubit1, qubit2);
                circuit->CZ(qubit1, qubit2);
            }
            // std::cout << gate << qubit1 << ", " << qubit2 << std::endl;
        }
        else
        {
            std::cout << gate << " does not match a gate."; 
        }
    }
    std::cout << "------\n\n\n TCount in qasm: " << tCount << " \n\n\n-----" << std::endl;
}

// Create a circuit with 2 qubits
int main(){
    std::vector<int> qubit_test = {4, 8, 16, 32, 64, 128, 256, 512, 1028, 2048, 4096, 8192, 16384};
    for(int i = 0; i < qubit_test.size(); i++)
    {
        std::cout << "Starting program" << std::endl;
        int n_qubits = qubit_test[i];
        int shots = 10;

        // NWQSim::IdxType S_count = 100000;
        // NWQSim::IdxType H_count = 100000;
        // NWQSim::IdxType CX_count = 100000;

        auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);

        // std::string inFile = "/Users/garn195/Project Repositories/NWQ-Sim/stabilizer/T_transpilation_test/adder_n10.qasm";
        // if(inFile != "")
        //     appendQASMToCircuit(circuit, inFile, n_qubits);

        // circuit->H(0);
        // circuit->S(0);
        // circuit->S(0);
        // circuit->H(0);
        // circuit->H(1);
        // circuit->S(1);
        // circuit->S(1);
        // circuit->H(1);
        // circuit->CX(1,0);
        // circuit->H(0);
        // circuit->H(3);
        // circuit->CX(2,3);
        // circuit->CX(2,3);
        // circuit->CX(1,2);
        // circuit->CX(3,0);
        // circuit->CX(0,1);
        // circuit->CX(0,1);
        // circuit->CX(2,3);
        // circuit->CX(2,3);
        // circuit->S(3);
        // circuit->CX(3,0);
        // circuit->H(3);


        // std::srand(std::time(nullptr));  // Seed random number generator

        // for(int i = 0; i < 100000; i++) 
        // {
        //     for(int j = 0; j < n_qubits; j++)
        //     {
        //         circuit->H((std::rand() % (n_qubits-1)));
        //         circuit->CX((std::rand() % (n_qubits-1)),(std::rand() % (n_qubits)));
        //         circuit->S((std::rand() % (n_qubits-1)));
        //     }
        // }
        for(int i = 0; i < 100000; i++) 
        {
            for(int j = 0; j < n_qubits; j++)
            {
                circuit->H((std::rand() % (n_qubits-1)));
                circuit->CX((std::rand() % (n_qubits-1)),(std::rand() % (n_qubits)));
                circuit->S((std::rand() % (n_qubits-1)));
            }
        }

        std::string backend = "nvgpu";
        std::string sim_method = "stab";
        double timer = 0;
        
        /*Create T and Measurement Tableaus with only stabilizers. T starts empty, M starts as identity.*/
        std::cout << "Creating state" << std::endl;
        auto state = BackendManager::create_state(backend, n_qubits, sim_method);


        std::cout << "Starting sim" << std::endl;
        state->sim(circuit, timer);
        // state->print_res_state();
        // NWQSim::IdxType* results = state->measure_all(shots);

        // for(int i = 0; i < shots; i++)
        //     std::cout << "Result " << i << ": " << results[i] << std::endl;

        std::cout << "Sim time: " << timer/1000.0 << "s" << std::endl;

        NWQSim::IdxType gate_count = S_count + H_count + CX_count;

        std::string name = "";
        std::ostringstream filename;
        filename << "/people/garn195/NWQ-Sim/stabilizer/sim_bench/" << backend << "_" << sim_method << "_" << n_qubits << ".txt";
        std::ofstream outfile(filename.str());
        if (!outfile) {
            std::cerr << "Error opening file: " << filename.str() << std::endl;
        }

        outfile << sim_method << std::endl;
        outfile << timer/1000.0 << std::endl;
        outfile << n_qubits << std::endl;
        outfile << S_count << std::endl;
        outfile << H_count << std::endl;
        outfile << CX_count << std::endl;

        outfile.close(); // Close the file
    }

    return 0;
}