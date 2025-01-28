#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../backendManager.hpp"
#include "../state.hpp"
#include "../circuit.hpp"
#include "../nwq_util.hpp"

#include "T_seperation.hpp"

int main(){

    /* Small test */
    // int n_qubits = 4;
    // auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
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
    int n_qubits = 10;
    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
    circuit->CX(0, 4);
    circuit->CX(1, 5);
    circuit->CX(2, 6);
    circuit->CX(3, 7);
    circuit->CX(0, 8);
    circuit->H(0);
    circuit->CX(4, 0);
    circuit->S(0);
    circuit->S(0);
    circuit->S(0);
    circuit->T(0);
    circuit->CX(8, 0);
    circuit->T(0);
    circuit->CX(4, 0);
    circuit->S(0);
    circuit->S(0);
    circuit->S(0);
    circuit->T(0);
    circuit->T(4);
    circuit->CX(8, 0);
    circuit->T(0);
    circuit->H(0);
    circuit->CX(1, 0);
    circuit->H(1);
    circuit->CX(5, 1);
    circuit->S(1);
    circuit->S(1);
    circuit->S(1);
    circuit->T(1);
    circuit->CX(0, 1);
    circuit->T(1);
    circuit->CX(5, 1);
    circuit->S(1);
    circuit->S(1);
    circuit->S(1);
    circuit->T(1);
    circuit->CX(0, 1);
    circuit->T(1);
    circuit->H(1);
    circuit->CX(2, 1);
    circuit->H(2);
    circuit->T(5);
    circuit->CX(0, 5);
    circuit->T(0);
    circuit->S(5);
    circuit->S(5);
    circuit->S(5);
    circuit->T(5);
    circuit->CX(0, 5);
    circuit->H(5);
    circuit->T(5);
    circuit->T(5);
    circuit->T(5);
    circuit->T(5);
    circuit->H(5);
    circuit->CX(0, 5);
    circuit->CX(6, 2);
    circuit->S(2);
    circuit->S(2);
    circuit->S(2);
    circuit->T(2);
    circuit->CX(1, 2);
    circuit->T(2);
    circuit->CX(6, 2);
    circuit->S(2);
    circuit->S(2);
    circuit->S(2);
    circuit->T(2);
    circuit->CX(1, 2);
    circuit->T(2);
    circuit->H(2);
    circuit->CX(3, 2);
    circuit->H(3);
    circuit->T(6);
    circuit->CX(1, 6);
    circuit->T(1);
    circuit->S(6);
    circuit->S(6);
    circuit->S(6);
    circuit->T(6);
    circuit->CX(1, 6);
    circuit->H(6);
    circuit->T(6);
    circuit->T(6);
    circuit->T(6);
    circuit->T(6);
    circuit->H(6);
    circuit->CX(1, 6);
    circuit->CX(7, 3);
    circuit->S(3);
    circuit->S(3);
    circuit->S(3);
    circuit->T(3);
    circuit->CX(2, 3);
    circuit->T(3);
    circuit->CX(7, 3);
    circuit->S(3);
    circuit->S(3);
    circuit->S(3);
    circuit->T(3);
    circuit->CX(2, 3);
    circuit->T(3);
    circuit->H(3);
    circuit->T(7);
    circuit->CX(2, 7);
    circuit->T(2);
    circuit->S(7);
    circuit->S(7);
    circuit->S(7);
    circuit->T(7);
    circuit->CX(2, 7);
    circuit->H(7);
    circuit->T(7);
    circuit->T(7);
    circuit->T(7);
    circuit->T(7);
    circuit->H(7);
    circuit->CX(2, 7);
    circuit->CX(8, 4);
    circuit->S(4);
    circuit->S(4);
    circuit->S(4);
    circuit->T(4);
    circuit->T(8);
    circuit->CX(8, 4);
    circuit->H(4);
    circuit->T(4);
    circuit->T(4);
    circuit->T(4);
    circuit->T(4);
    circuit->H(4);
    circuit->CX(8, 4);
    circuit->CX(3, 9);
    circuit->H(3);
    circuit->CX(7, 3);
    circuit->S(3);
    circuit->S(3);
    circuit->S(3);
    circuit->T(3);
    circuit->CX(2, 3);
    circuit->T(3);
    circuit->CX(7, 3);
    circuit->S(3);
    circuit->S(3);
    circuit->S(3);
    circuit->T(3);
    circuit->CX(2, 3);
    circuit->T(3);
    circuit->H(3);
    circuit->T(7);
    circuit->CX(2, 7);
    circuit->T(2);
    circuit->S(7);
    circuit->S(7);
    circuit->S(7);
    circuit->T(7);
    circuit->CX(2, 7);
    circuit->CX(3, 2);
    circuit->H(2);
    circuit->CX(6, 2);
    circuit->S(2);
    circuit->S(2);
    circuit->S(2);
    circuit->T(2);
    circuit->CX(1, 2);
    circuit->T(2);
    circuit->CX(6, 2);
    circuit->S(2);
    circuit->S(2);
    circuit->S(2);
    circuit->T(2);
    circuit->CX(1, 2);
    circuit->T(2);
    circuit->H(2);
    circuit->T(6);
    circuit->CX(1, 6);
    circuit->T(1);
    circuit->S(6);
    circuit->S(6);
    circuit->S(6);
    circuit->T(6);
    circuit->CX(1, 6);
    circuit->CX(2, 1);
    circuit->H(1);
    circuit->CX(5, 1);
    circuit->S(1);
    circuit->S(1);
    circuit->S(1);
    circuit->T(1);
    circuit->CX(0, 1);
    circuit->T(1);
    circuit->CX(5, 1);
    circuit->S(1);
    circuit->S(1);
    circuit->S(1);
    circuit->T(1);
    circuit->CX(0, 1);
    circuit->T(1);
    circuit->H(1);
    circuit->T(5);
    circuit->CX(0, 5);
    circuit->T(0);
    circuit->S(5);
    circuit->S(5);
    circuit->S(5);
    circuit->T(5);
    circuit->CX(0, 5);
    circuit->CX(1, 0);
    circuit->H(0);
    circuit->CX(4, 0);
    circuit->S(0);
    circuit->S(0);
    circuit->S(0);
    circuit->T(0);
    circuit->H(5);
    circuit->T(5);
    circuit->T(5);
    circuit->T(5);
    circuit->T(5);
    circuit->H(5);
    circuit->CX(1, 5);
    circuit->H(6);
    circuit->T(6);
    circuit->T(6);
    circuit->T(6);
    circuit->T(6);
    circuit->H(6);
    circuit->CX(2, 6);
    circuit->H(7);
    circuit->T(7);
    circuit->T(7);
    circuit->T(7);
    circuit->T(7);
    circuit->H(7);
    circuit->CX(3, 7);
    circuit->CX(8, 0);
    circuit->T(0);
    circuit->CX(4, 0);
    circuit->S(0);
    circuit->S(0);
    circuit->S(0);
    circuit->T(0);
    circuit->T(4);
    circuit->CX(8, 0);
    circuit->T(0);
    circuit->H(0);
    circuit->CX(8, 4);
    circuit->S(4);
    circuit->S(4);
    circuit->S(4);
    circuit->T(4);
    circuit->T(8);
    circuit->CX(8, 4);
    circuit->CX(0, 8);
    circuit->H(4);
    circuit->T(4);
    circuit->T(4);
    circuit->T(4);
    circuit->T(4);
    circuit->H(4);
    circuit->CX(0, 4);



    int shots = 10;

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
    auto start = std::chrono::high_resolution_clock::now();
    T_passthrough(circuit, T_tab, M_circ, 10);
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time elapsed: " << duration.count() << " milliseconds" << std::endl;


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
 