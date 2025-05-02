#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>
#include <regex>

#include "../../include/backendManager.hpp"
#include "../../include/state.hpp"
#include "../../include/circuit.hpp"
#include "../../include/nwq_util.hpp"

int extractQubitIndex(const std::string& qubitStr) 
{
    std::regex qubitRegex(R"((\w+)\[(\d+)\])");  // match format like qregless[0]
    std::smatch match;
    if (std::regex_search(qubitStr, match, qubitRegex) && match.size() >= 3) 
    {
        try {
            return std::stoi(match.str(2));  // Get the number inside the brackets
        } catch (const std::exception& e) {
            std::cerr << "Error parsing qubit index from: " << qubitStr << " - " << e.what() << std::endl;
            return -1;
        }
    }
    else
    {
        std::cerr << "Gate called but no qubit index! String: " << qubitStr << std::endl;
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

std::vector<std::string> split_by_semicolon(const std::string& str) {
    std::vector<std::string> tokens;
    std::istringstream ss(str);
    std::string token;
    while (std::getline(ss, token, ';')) {
        size_t start = token.find_first_not_of(" \t");
        size_t end = token.find_last_not_of(" \t");
        if (start != std::string::npos)
            tokens.push_back(token.substr(start, end - start + 1));
    }
    return tokens;
}

bool appendQASMToCircuit(std::shared_ptr<NWQSim::Circuit>& circuit, const std::string& filename, int& n_qubits) 
{
    std::ifstream file(filename);
    if (!file.is_open()) 
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return false;
    }

    std::string line;
    int tCount = 0;

    while (std::getline(file, line)) 
    {
        if (line.empty() ||
            line.find("include") != std::string::npos ||
            line.find("gate") != std::string::npos ||
            line.find("barrier") != std::string::npos ||
            line.find("OPENQASM") != std::string::npos ||
            line.find("creg") != std::string::npos ||
            line.find("#include") != std::string::npos)
            continue;

        auto instructions = split_by_semicolon(line);

        for (const std::string& instr : instructions)
        {
            std::istringstream lineStream(instr);
            std::string gate;
            int qubit1 = -1, qubit2 = -1;

            lineStream >> gate;

            if(instr.find("//") != std::string::npos)
            {
                continue;
            }
            else if(instr.find("qreg") != std::string::npos)
            {
                n_qubits = extractNumQubit(instr);
                if(n_qubits > 30000)
                    return false;
                circuit->set_num_qubits(n_qubits);
                std::cout << n_qubits << std::endl;
                
            }
            else if(gate == "tdg")
            {
                tCount++;
                std::string qubitStr;
                lineStream >> qubitStr;
                qubit1 = extractQubitIndex(qubitStr);
                if (qubit1 != -1) circuit->TDG(qubit1);
            }
            else if(gate == "t")
            {
                tCount++;
                std::string qubitStr;
                lineStream >> qubitStr;
                qubit1 = extractQubitIndex(qubitStr);
                if (qubit1 != -1) circuit->T(qubit1);
            }
            else if(gate == "sdg")
            {
                std::string qubitStr;
                lineStream >> qubitStr;
                qubit1 = extractQubitIndex(qubitStr);
                if (qubit1 != -1) circuit->SDG(qubit1);
            }
            else if(gate == "s")
            {
                std::string qubitStr;
                lineStream >> qubitStr;
                qubit1 = extractQubitIndex(qubitStr);
                if (qubit1 != -1) circuit->S(qubit1);
            }
            else if(gate == "h")
            {
                std::string qubitStr;
                lineStream >> qubitStr;
                qubit1 = extractQubitIndex(qubitStr);
                if (qubit1 != -1) circuit->H(qubit1);
            }
            else if(gate == "m" || gate == "measure")
            {
                std::string qubitStr;
                lineStream >> qubitStr;
                qubit1 = extractQubitIndex(qubitStr);
                if (qubit1 != -1) circuit->M(qubit1);
            }
            else if(gate == "reset")
            {
                std::string qubitStr;
                lineStream >> qubitStr;
                qubit1 = extractQubitIndex(qubitStr);
                if (qubit1 != -1) circuit->RESET(qubit1);
            }
            else if(gate == "cx" || gate == "cxyz")
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
                    if (gate == "cxyz")
                    {
                        circuit->CY(qubit1, qubit2);
                        circuit->CZ(qubit1, qubit2);
                    }
                }
            }
            // else if (!gate.empty())
            // {
            //     std::cout << gate << " does not match a gate.\n";
            // }
        }
    }  
    return true;

    // std::cout << "------\n\n\n TCount in qasm: " << tCount << " \n\n\n-----" << std::endl;
}