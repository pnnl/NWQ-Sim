#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../backendManager.hpp"
#include "../state.hpp"
#include "../circuit.hpp"
#include "../nwq_util.hpp"

// Create a circuit with 2 qubits
int main(){
    int n_qubits = 4;
    int shots = 10;

    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);

    // circuit->H(3);
    // circuit->CX(2,0);
    // circuit->CX(0,1);
    // circuit->H(3);
    // circuit->CX(0,2);
    // circuit->CX(0,2);
    // //circuit->T(1);
    // circuit->CX(3,2);
    // circuit->CX(2,0);
    // circuit->H(3);
    // //circuit->T(2);
    // circuit->CX(0,1);
    // circuit->CX(1,0);
    // circuit->CX(3,2);
    // circuit->CX(1,0);

    circuit -> H(0);
    circuit -> CX(0,1);

    std::string backend = "CPU";
    std::string sim_method = "stab";
    double timer = 0;
    
    /*Create T and Measurement Tableaus with only stabilizers. T starts empty, M starts as identity.*/
    auto state = BackendManager::create_state(backend, n_qubits, sim_method);

    state->sim(circuit, timer);
    state->print_res_state();
    std::cout << "M tableau measurement results: " << (state->measure_all(10))[0] << std::endl;
}