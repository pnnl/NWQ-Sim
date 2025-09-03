#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>
#include <regex>

#include "../../backendManager.hpp"
#include "../../state.hpp"
#include "../../circuit.hpp"
#include "../../nwq_util.hpp"

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


std::vector<double> extract_params(const std::string& params_str) {
    std::vector<double> params;
    std::string content = params_str.substr(params_str.find('(') + 1, params_str.find(')') - params_str.find('(') - 1);
    std::stringstream ss(content);
    std::string item;
    while (std::getline(ss, item, ',')) {
        params.push_back(std::stod(item));
    }
    return params;
}

std::string extract_base_gate(const std::string& gate_str) {
    size_t paren_pos = gate_str.find('(');
    if (paren_pos != std::string::npos) {
        return gate_str.substr(0, paren_pos);
    }
    return gate_str;
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
            std::string gate_full;
            int qubit1 = -1, qubit2 = -1;

            lineStream >> gate_full;
            std::string gate = extract_base_gate(gate_full);

            if(instr.find("//") != std::string::npos)
            {
                continue;
            }
            else if(instr.find("qreg") != std::string::npos)
            {
                n_qubits = extractNumQubit(instr);
                circuit->set_num_qubits(n_qubits);
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
            else if(gate == "damp")
            {
                auto params = extract_params(gate_full);
                std::string qubitStr;
                lineStream >> qubitStr;
                qubit1 = extractQubitIndex(qubitStr);
                if (qubit1 != -1) circuit->DAMP(qubit1, params[0], params[1]);
            }
            else if (gate == "dep1")
            {
                auto params = extract_params(gate_full);
                std::string qubitStr;
                lineStream >> qubitStr;
                qubit1 = extractQubitIndex(qubitStr);
                if (qubit1 != -1) circuit->DEP1(qubit1, params[0]);
            }
            else if (gate == "dep2")
            {
                auto params = extract_params(gate_full);
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
                    circuit->DEP2(qubit1, qubit2, params[0]);
                }
            }
            else if (gate == "chan1")  // NEW: CHAN1 – single qubit, 3 parameters
            {
                std::vector<double> params = extract_params(gate_full);
                std::string qubitStr;
                lineStream >> qubitStr;
                qubit1 = extractQubitIndex(qubitStr);
                // std::cerr << "[qasm_extraction] chan1 q=" << qubit1
                //           << " params=" << params.size() << std::endl; // stderr
                if (qubit1 != -1)
                {
                    // std::cout<< "CHAN1 reached in extract!: "  << params[1] << std::endl;
                    circuit->CHAN1(qubit1, params); // expects vector<double>{px,py,pz}
                }
            }
            else if (gate == "chan2")  // NEW: CHAN2 – two qubits, up to 15 parameters
            {
                std::vector<double> params = extract_params(gate_full);
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
                    // Forward the full vector<double> of probabilities (length can be 15)
                    circuit->CHAN2(qubit1, qubit2, params);
                }
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
        }
    }  
    return true;
}
