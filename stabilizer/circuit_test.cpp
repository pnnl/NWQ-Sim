#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../include/backendManager.hpp"
#include "../include/state.hpp"
#include "../include/circuit.hpp"
#include "../include/nwq_util.hpp"

namespace NWQSim{


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
        
        return -1;
    }

    void appendQASMToCircuit(std::shared_ptr<Circuit>& circuit, const std::string& filename, int& n_qubits) 
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
}
int main()
{

    std::vector<int> qubit_test;
    for(int i = 4; i < 200; i+=4)
        qubit_test.push_back(pow(i,2));
    for(int i = 0; i < qubit_test.size(); i++)
    {
        std::cout << "Starting program" << std::endl;
        int n_qubits = qubit_test[i];
        int shots = 10;

        int layers = std::log2(n_qubits);
        auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
        auto m_circuit = std::make_shared<NWQSim::Circuit>(n_qubits);

        std::cout << "Building circuit" << std::endl;

        std::string backend = "cpu";
        std::string sim_method = "stab";
        double timer = 0;
        
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, 1);
        std::uniform_real_distribution<> prob(0.0, 1.0);

        /*Create T and Measurement Tableaus with only stabilizers. T starts empty, M starts as identity.*/
        std::cout << "Creating state" << std::endl;
        auto state = BackendManager::create_state(backend, n_qubits, sim_method);

        for(int k = 0; k < layers; k++)
        {
            int half_qubit = qubit_test[i]/2;
            for(int j = 0; j < half_qubit; j++) 
            {
                if (dis(gen))
                    circuit->H(j);
                else
                    circuit->S(j);

                circuit->CX(j, j+half_qubit);
                
                if (prob(gen) < 0.2)
                    circuit->M(j+half_qubit);
            }    
        }


        state->sim(circuit, timer);

        std::cout << "Sim time: " << timer/1000.0 << "s" << std::endl;

        std::string name = "";
        std::ostringstream filename;
        filename << "/Users/garn195/Project Repositories/NWQ-Sim/stabilizer/surface_code_laptop/" << backend << "_" << sim_method << "_" << n_qubits << ".txt";
        std::ofstream outfile(filename.str());
        if (!outfile) {
            std::cerr << "Error opening file: " << filename.str() << std::endl;
        }

        outfile << backend << std::endl;
        outfile << timer/1000.0 << std::endl;
        outfile << n_qubits << std::endl;

        outfile.close();
    }

    return 0;
}
