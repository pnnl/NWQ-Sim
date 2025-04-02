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
#include "src/qasm_extraction.hpp"



int main(){

    /* Small test */
    int n_qubits = 4;
    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);

    std::string inFile = "/Users/garn195/Project Repositories/NWQ-Sim/stabilizer/T_transpilation_test/multiplier_n45.qasm";
    std::string outFile = "/Users/garn195/Project Repositories/NWQ-Sim/stabilizer/stab_T_bench/multiplier_n45.txt";
    std::ofstream outfile(outFile); // Open a file for writing
    appendQASMToCircuit(circuit, inFile, n_qubits);
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
    M_tab->print_res_state();
    T_passthrough(T_tab, M_tab, outfile, proc_time, 10);


    //Put the M circuit back in forward time after the new Clifford gates have been appended
    std::cout << "---- M tableau -----" << std::endl;

    M_tab->print_res_state();
    std::cout << "M tableau measurement results: " << (M_tab->measure_all(10))[0] << std::endl;
    /*Measurement Tableau*/

    outfile << total_gates << std::endl;
    outfile.close(); // Close the file

    return 0;
}