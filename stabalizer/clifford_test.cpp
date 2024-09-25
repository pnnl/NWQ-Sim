#include <string>
#include <complex>
#include <iostream>
#include <vector>

#include "../include/backendManager.hpp"
#include "../include/state.hpp"
#include "../include/circuit.hpp"
#include "../include/nwq_util.hpp"

int main(){
    // Create a circuit with 2 qubits
	int n_qubits = 2;
	auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);

	// Add some gates to the circuit
	circuit->H(0);
	circuit->CX(0, 1);
	//circuit->RZ(0.125, 0);
	std::string backend = "CPU";
	std::string sim_method = "sv";
	auto state = BackendManager::create_state(backend, n_qubits, sim_method);

	// Add measurement operation to circuit then simulate
	int shots = 1;
	circuit->MA(shots);
	state->sim(circuit);
	long long int *results = state->get_results();

	for(int i = 0; i < shots; i++)
	{
		std::cout << results[0] << std::endl;
	}

	

	return 0;
}