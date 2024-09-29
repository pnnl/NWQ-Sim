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
	int n_qubits = 5;
	int circuit_length = 1;
	auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);

	// Add some gates to the circuit
	for(int i = 0; i < circuit_length; i++)
	{
		circuit ->H(0);
		circuit ->CX(0, 1);
		circuit ->CX(0, 2);
		circuit ->CX(0, 3);
		circuit ->CX(0, 4);

		// circuit ->H(0);
		circuit ->H(1);
		circuit ->H(2);
		circuit ->H(3);
		circuit ->H(4);

		// circuit ->S(0);
		// circuit ->S(1);
		// circuit ->S(2);
		// circuit ->S(3);
		// circuit ->S(4);

		// circuit ->CX(1, 2);
		// circuit ->CX(1, 3);
		// circuit ->CX(1, 4);
	}
	circuit->M(0);
	circuit->M(1);
	circuit->M(2);
	circuit->M(3);
	circuit->M(4);
	


	std::string backend = "CPU";
	std::string sim_method = "sv";
	
	// Add measurement operation to circuit then simulate

	auto state1 = BackendManager::create_state(backend, n_qubits, sim_method);
	state1->sim(circuit);
	long long int *results1 = state1->get_results();

	std::cout << "NWQSim: " << results1[0] << "\n" << std::endl;

	auto state2 = BackendManager::create_state(backend, n_qubits, sim_method);
	state2->clifford_sim(circuit);
	long long int *results2 = state2->get_results();

	std::cout << "Stabalizer sim: " << results2[0] << std::endl;

	std::cout << "Done" << std::endl;

	return 0;
}