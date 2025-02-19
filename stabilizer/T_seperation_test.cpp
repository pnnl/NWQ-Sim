#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../include/backendManager.hpp"
#include "../include/state.hpp"
#include "../include/circuit.hpp"
#include "../include/nwq_util.hpp"

#include "src/T_seperation.hpp"

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

void appendQASMToCircuit(std::shared_ptr<NWQSim::Circuit>& circuit, const std::string& filename) 
{
    std::ifstream file(filename);
    if (!file.is_open()) 
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) 
    {
        std::string gate;
        int qubit1, qubit2;

        if (line.empty() || line.find("include") != std::string::npos || line.find("gate") != std::string::npos || line.find("barrier") != std::string::npos || line.find("OPENQASM") != std::string::npos || line.find("qreg") != std::string::npos || line.find("creg") != std::string::npos || line.find("#include") != std::string::npos)
            continue;

        std::istringstream lineStream(line);
        lineStream >> gate;
        if(gate == "tdg")
        {
            std::string qubitStr;
            lineStream >> qubitStr;
            qubit1 = extractQubitIndex(qubitStr);
            if (qubit1 != -1) 
            {
                circuit->TDG(qubit1);
            }
            std::cout << gate << qubit1 << ", " << qubit2 << std::endl;
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
            std::cout << gate << qubit1 << ", " << qubit2 << std::endl;
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
            std::cout << gate << qubit1 << ", " << qubit2 << std::endl;
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
            std::cout << gate << qubit1 << ", " << qubit2 << std::endl;
        }
        else if(gate == "t")
        {
            std::string qubitStr;
            lineStream >> qubitStr;
            qubit1 = extractQubitIndex(qubitStr);
            if (qubit1 != -1) 
            {
                circuit->T(qubit1);
            }
            std::cout << gate << qubit1 << ", " << qubit2 << std::endl;
        }
        else if(gate == "measure")
        {
            std::string qubitStr;
            lineStream >> qubitStr;
            qubit1 = extractQubitIndex(qubitStr);
            if (qubit1 != -1) 
            {
                circuit->M(qubit1);
            }
            std::cout << gate << qubit1 << ", " << qubit2 << std::endl;
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
            std::cout << gate << qubit1 << ", " << qubit2 << std::endl;
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
            std::cout << gate << qubit1 << ", " << qubit2 << std::endl;
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
            std::cout << gate << qubit1 << ", " << qubit2 << std::endl;
        }
        else
        {
            std::cout << gate << " does not match a gate."; 
        }
    }
}


int main(){

    /* Small test */
    int n_qubits = 18;
    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
    // circuit->T(1);
    // circuit->T(1);
    // circuit->T(1);

    // circuit->H(3);
    // circuit->CX(2,0);
    // circuit->CX(0,1);
    // circuit->H(3);
    // circuit->CX(0,2);
    // circuit->CX(0,2);
    // circuit->T(1);
    // circuit->CX(3,2);
    // circuit->CX(2,0);
    // circuit->H(3);
    // circuit->T(2);
    // circuit->CX(0,1);
    // circuit->CX(1,0);
    // circuit->CX(3,2);
    // circuit->CX(1,0);

    /*adder_n10*/
    // int n_qubits = 10;
    // auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
    // circuit->CX(0, 4);
    // circuit->CX(1, 5);
    // circuit->CX(2, 6);
    // circuit->CX(3, 7);
    // circuit->CX(0, 8);
    // circuit->H(0);
    // circuit->CX(4, 0);
    // circuit->S(0);
    // circuit->S(0);
    // circuit->S(0);
    // circuit->T(0);
    // circuit->CX(8, 0);
    // circuit->T(0);
    // circuit->CX(4, 0);
    // circuit->S(0);
    // circuit->S(0);
    // circuit->S(0);
    // circuit->T(0);
    // circuit->T(4);
    // circuit->CX(8, 0);
    // circuit->T(0);
    // circuit->H(0);
    // circuit->CX(1, 0);
    // circuit->H(1);
    // circuit->CX(5, 1);
    // circuit->S(1);
    // circuit->S(1);
    // circuit->S(1);
    // circuit->T(1);
    // circuit->CX(0, 1);
    // circuit->T(1);
    // circuit->CX(5, 1);
    // circuit->S(1);
    // circuit->S(1);
    // circuit->S(1);
    // circuit->T(1);
    // circuit->CX(0, 1);
    // circuit->T(1);
    // circuit->H(1);
    // circuit->CX(2, 1);
    // circuit->H(2);
    // circuit->T(5);
    // circuit->CX(0, 5);
    // circuit->T(0);
    // circuit->S(5);
    // circuit->S(5);
    // circuit->S(5);
    // circuit->T(5);
    // circuit->CX(0, 5);
    // circuit->H(5);
    // circuit->T(5);
    // circuit->T(5);
    // circuit->T(5);
    // circuit->T(5);
    // circuit->H(5);
    // circuit->CX(0, 5);
    // circuit->CX(6, 2);
    // circuit->S(2);
    // circuit->S(2);
    // circuit->S(2);
    // circuit->T(2);
    // circuit->CX(1, 2);
    // circuit->T(2);
    // circuit->CX(6, 2);
    // circuit->S(2);
    // circuit->S(2);
    // circuit->S(2);
    // circuit->T(2);
    // circuit->CX(1, 2);
    // circuit->T(2);
    // circuit->H(2);
    // circuit->CX(3, 2);
    // circuit->H(3);
    // circuit->T(6);
    // circuit->CX(1, 6);
    // circuit->T(1);
    // circuit->S(6);
    // circuit->S(6);
    // circuit->S(6);
    // circuit->T(6);
    // circuit->CX(1, 6);
    // circuit->H(6);
    // circuit->T(6);
    // circuit->T(6);
    // circuit->T(6);
    // circuit->T(6);
    // circuit->H(6);
    // circuit->CX(1, 6);
    // circuit->CX(7, 3);
    // circuit->S(3);
    // circuit->S(3);
    // circuit->S(3);
    // circuit->T(3);
    // circuit->CX(2, 3);
    // circuit->T(3);
    // circuit->CX(7, 3);
    // circuit->S(3);
    // circuit->S(3);
    // circuit->S(3);
    // circuit->T(3);
    // circuit->CX(2, 3);
    // circuit->T(3);
    // circuit->H(3);
    // circuit->T(7);
    // circuit->CX(2, 7);
    // circuit->T(2);
    // circuit->S(7);
    // circuit->S(7);
    // circuit->S(7);
    // circuit->T(7);
    // circuit->CX(2, 7);
    // circuit->H(7);
    // circuit->T(7);
    // circuit->T(7);
    // circuit->T(7);
    // circuit->T(7);
    // circuit->H(7);
    // circuit->CX(2, 7);
    // circuit->CX(8, 4);
    // circuit->S(4);
    // circuit->S(4);
    // circuit->S(4);
    // circuit->T(4);
    // circuit->T(8);
    // circuit->CX(8, 4);
    // circuit->H(4);
    // circuit->T(4);
    // circuit->T(4);
    // circuit->T(4);
    // circuit->T(4);
    // circuit->H(4);
    // circuit->CX(8, 4);
    // circuit->CX(3, 9);
    // circuit->H(3);
    // circuit->CX(7, 3);
    // circuit->S(3);
    // circuit->S(3);
    // circuit->S(3);
    // circuit->T(3);
    // circuit->CX(2, 3);
    // circuit->T(3);
    // circuit->CX(7, 3);
    // circuit->S(3);
    // circuit->S(3);
    // circuit->S(3);
    // circuit->T(3);
    // circuit->CX(2, 3);
    // circuit->T(3);
    // circuit->H(3);
    // circuit->T(7);
    // circuit->CX(2, 7);
    // circuit->T(2);
    // circuit->S(7);
    // circuit->S(7);
    // circuit->S(7);
    // circuit->T(7);
    // circuit->CX(2, 7);
    // circuit->CX(3, 2);
    // circuit->H(2);
    // circuit->CX(6, 2);
    // circuit->S(2);
    // circuit->S(2);
    // circuit->S(2);
    // circuit->T(2);
    // circuit->CX(1, 2);
    // circuit->T(2);
    // circuit->CX(6, 2);
    // circuit->S(2);
    // circuit->S(2);
    // circuit->S(2);
    // circuit->T(2);
    // circuit->CX(1, 2);
    // circuit->T(2);
    // circuit->H(2);
    // circuit->T(6);
    // circuit->CX(1, 6);
    // circuit->T(1);
    // circuit->S(6);
    // circuit->S(6);
    // circuit->S(6);
    // circuit->T(6);
    // circuit->CX(1, 6);
    // circuit->CX(2, 1);
    // circuit->H(1);
    // circuit->CX(5, 1);
    // circuit->S(1);
    // circuit->S(1);
    // circuit->S(1);
    // circuit->T(1);
    // circuit->CX(0, 1);
    // circuit->T(1);
    // circuit->CX(5, 1);
    // circuit->S(1);
    // circuit->S(1);
    // circuit->S(1);
    // circuit->T(1);
    // circuit->CX(0, 1);
    // circuit->T(1);
    // circuit->H(1);
    // circuit->T(5);
    // circuit->CX(0, 5);
    // circuit->T(0);
    // circuit->S(5);
    // circuit->S(5);
    // circuit->S(5);
    // circuit->T(5);
    // circuit->CX(0, 5);
    // circuit->CX(1, 0);
    // circuit->H(0);
    // circuit->CX(4, 0);
    // circuit->S(0);
    // circuit->S(0);
    // circuit->S(0);
    // circuit->T(0);
    // circuit->H(5);
    // circuit->T(5);
    // circuit->T(5);
    // circuit->T(5);
    // circuit->T(5);
    // circuit->H(5);
    // circuit->CX(1, 5);
    // circuit->H(6);
    // circuit->T(6);
    // circuit->T(6);
    // circuit->T(6);
    // circuit->T(6);
    // circuit->H(6);
    // circuit->CX(2, 6);
    // circuit->H(7);
    // circuit->T(7);
    // circuit->T(7);
    // circuit->T(7);
    // circuit->T(7);
    // circuit->H(7);
    // circuit->CX(3, 7);
    // circuit->CX(8, 0);
    // circuit->T(0);
    // circuit->CX(4, 0);
    // circuit->S(0);
    // circuit->S(0);
    // circuit->S(0);
    // circuit->T(0);
    // circuit->T(4);
    // circuit->CX(8, 0);
    // circuit->T(0);
    // circuit->H(0);
    // circuit->CX(8, 4);
    // circuit->S(4);
    // circuit->S(4);
    // circuit->S(4);
    // circuit->T(4);
    // circuit->T(8);
    // circuit->CX(8, 4);
    // circuit->CX(0, 8);
    // circuit->H(4);
    // circuit->T(4);
    // circuit->T(4);
    // circuit->T(4);
    // circuit->T(4);
    // circuit->H(4);
    // circuit->CX(0, 4);


    appendQASMToCircuit(circuit, "/people/garn195/NWQ-Sim/stabilizer/qft_test/qft_n18_iter1.qasm");

    //Measurement circuit will be filled in the passthrough function
    auto M_circ = std::make_shared<NWQSim::Circuit>(n_qubits);

    std::string backend = "CPU";
    std::string sim_method = "stab";
    double timer = 0;
    
    /*Create T and Measurement Tableaus with only stabilizers. T starts empty, M starts as identity.*/
    auto T_tab = BackendManager::create_state(backend, n_qubits, sim_method);
    T_tab->delete_all_rows();
    T_tab->remove_destabilizers();


    /*Run T passthrough after tableaus have been prepared*/
    T_passthrough(circuit, T_tab, M_circ, 10);

    /*T Tableau*/
    std::cout << "---- T tableau -----" << std::endl;
    T_tab->print_res_state();
    /*T Tableau*/

    /*Measurement Tableau*/
    circuit_reverse(M_circ);

    //After all M_tab gates and additional Clifford gates from optimization
    //Simulate the M circuit evolution after construction
    auto M_tab = BackendManager::create_state(backend, n_qubits, "stab");
    M_tab->sim(M_circ, timer);
    std::cout << "---- M tableau -----" << std::endl;

    M_tab->print_res_state();
    std::cout << "M tableau measurement results: " << (M_tab->measure_all(10))[0] << std::endl;
    /*Measurement Tableau*/



    return 0;
}
 