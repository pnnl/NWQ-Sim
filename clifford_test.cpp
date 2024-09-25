#include <string>
#include <complex>
#include <iostream>
#include <vector>

#include "NWQ-Sim/include/backendManager.hpp"
#include "NWQ-Sim/include/state.hpp"
#include "NWQ-Sim/include/circuit.hpp"
#include "NWQ-Sim/include/nwq_util.hpp"

int main(){
    // Create a circuit with 2 qubits
	int n_qubits = 2;
	auto circuit = std::make_shared<Circuit>(n_qubits);

	// Add some gates to the circuit
	circuit->H(0);
	circuit->CX(0, 1);
	//circuit->RZ(0.125, 0);
	std::string backend = "CPU";
	std::string sim_method = "sv";
	auto state = BackendManager::create_state(backend, n_qubits, sim_method);

	// Add measurement operation to circuit then simulate
	int shots = 1024;
	//circuit->MA(shots);
	state->sim(circuit);
	//long long int *result = state->get_results();


	return 0;
}