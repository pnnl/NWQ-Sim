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
	circuit->M(0);

	std::string backend = "CPU";
	std::string sim_method = "sv";
	
	// Add measurement operation to circuit then simulate
	int shots = 10;
	auto state1 = BackendManager::create_state(backend, n_qubits, sim_method);
	state1->sim(circuit);
	long long int *results1 = state1->get_results();

	for(int i = 0; i < shots; i++)
	{
		std::cout << results1[i] << std::endl;
	}

	auto state2 = BackendManager::create_state(backend, n_qubits, sim_method);
	state2->clifford_sim(circuit);
	long long int *results2 = state2->get_results();

	for(int i = 0; i < shots; i++)
	{
		std::cout << results2[i] << std::endl;
	}

	return 0;
}