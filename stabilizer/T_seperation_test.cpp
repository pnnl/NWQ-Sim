#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../include/backendManager.hpp"
#include "../include/state.hpp"
#include "../include/circuit.hpp"
#include "../include/nwq_util.hpp"

#include "src/T_separation.hpp"

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


int main(){

    /* Small test */
    int n_qubits = 4;
    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);

    std::string inFile = "/Users/garn195/Project Repositories/NWQ-Sim/stabilizer/T_transpilation_test/qft_n18_iter2.qasm";
    std::string outFile = "/Users/garn195/Project Repositories/NWQ-Sim/stabilizer/stab_T_bench/qft_n18_iter2.txt";
    appendQASMToCircuit(circuit, inFile, n_qubits);


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
    M_tab->print_res_state();
    T_passthrough(T_tab, M_tab, outFile, proc_time, 10);


    //Put the M circuit back in forward time after the new Clifford gates have been appended
    std::cout << "---- M tableau -----" << std::endl;

    M_tab->print_res_state();
    std::cout << "M tableau measurement results: " << (M_tab->measure_all(10))[0] << std::endl;
    /*Measurement Tableau*/



    return 0;
}
 