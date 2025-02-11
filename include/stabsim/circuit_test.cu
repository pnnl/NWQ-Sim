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

    circuit -> H(0);
    circuit -> S(0);
    circuit -> S(0);
    circuit -> H(0);

    std::string backend = "NVGPU";
    std::string sim_method = "stab";
    double timer = 0;
    
    /*Create T and Measurement Tableaus with only stabilizers. T starts empty, M starts as identity.*/
    auto state = BackendManager::create_state(backend, n_qubits, sim_method);

    state->sim(circuit, timer);
    state->print_res_state();
    NWQSim::IdxType* results = state->measure_all(shots);

    for(int i = 0; i < shots; i++)
        std::cout << "Result " << i << ": " << results[i] << std::endl;
}