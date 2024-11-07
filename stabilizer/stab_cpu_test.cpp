#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../include/backendManager.hpp"
#include "../include/state.hpp"
#include "../include/circuit.hpp"
#include "../include/nwq_util.hpp"

int main(){

    // Create a circuit with 2 qubits
    int n_qubits = 3;
    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);

    // Add some gates to the circuit
    //circuit->RZ(0.125, 0);

    std::string backend = "CPU";
    std::string sim_method = "stab";
    double timer;
    auto state = BackendManager::create_state(backend, n_qubits, sim_method);
    std::cout << "TEST" << std::endl;

    state->sim(circuit, timer);
    int shots = 5;

    NWQSim::IdxType* result = state->measure_all(shots);
    for(int i = 0; i < (1 << n_qubits); i++)
    {
        std::cout << "Result " << i << ": " << result[i] << std::endl;
    }


	return 0;
}	