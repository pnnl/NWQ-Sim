#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../../../include/backendManager.hpp"
#include "../../../include/state.hpp"
#include "../../../include/circuit.hpp"
#include "../../../include/nwq_util.hpp"

int main()
{

    std::cout << "Starting program" << std::endl;
    int n_qubits = 3;
    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
    std::cout << "Building circuit" << std::endl;

    circuit->H(0);
    circuit->S(1);
    circuit->CX(0, 1);
    circuit->CX(1, 2);

    circuit->H(2);
    circuit->SDG(0);
    circuit->CX(2, 0);
    circuit->H(1);

    circuit->S(2);
    circuit->CX(0,2);
    circuit->SDG(1);
    circuit->H(0);

    circuit->CX(1, 0);
    circuit->CX(2, 1);
    circuit->S(0);
    circuit->H(2);

    std::string backend = "cpu";
    std::string sim_method = "stab";
    double timer = 0;
    
    /*Create T and Measurement Tableaus with only stabilizers. T starts empty, M starts as identity.*/
    std::cout << "Creating state" << std::endl;
    auto state = BackendManager::create_state(backend, n_qubits, sim_method);

    std::cout << "Starting sim bitwise" << std::endl;

    state->sim(circuit, timer);
    state->print_res_state();

    return 0;
}