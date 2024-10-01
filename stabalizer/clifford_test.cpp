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
	std::ofstream outfile("qubits_time.txt");
	//std::ofstream outfile("qubits_memory.txt");

	for(int rounds = 0; rounds < 1; rounds++){
		for(int n_qubits = 10; n_qubits < 12; n_qubits++)
		{
			auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
			for(int j = 0; j < rounds; j++)
			{
				// Add some gates to the circuit
				for(int i = 0; i < n_qubits; i++)
				{
					circuit -> H(i);
				}
				for(int i = 0; i < n_qubits-1; i++)
				{
					circuit -> CX(i,i+1);
				}
				for(int i = 0; i < n_qubits; i++)
				{
					circuit -> S(i);
				}
			}

			for(int i = 0; i < n_qubits; i++)
			{
				circuit -> M(i);
			}

			std::string backend = "CPU";
			std::string sim_method = "sv";
			
			// Add measurement operation to circuit then simulate

			auto state1 = BackendManager::create_state(backend, n_qubits, sim_method);
			double state1_time;
			state1->sim(circuit, state1_time);
			long long int *results1 = state1->get_results();

			std::cout << "NWQSim: " << results1[0] << "\n" << std::endl;

			auto state2 = BackendManager::create_state(backend, n_qubits, sim_method);
			double state2_time;
			state2->clifford_sim(circuit, state2_time);
			long long int *results2 = state2->get_results();

			std::cout << "Stabalizer sim: " << results2[0] << std::endl;

			outfile << std::setw(5) << rounds << " ";
			outfile << std::setw(5) << n_qubits << " ";
			outfile << std::setw(5) << state1_time << " ";
			outfile << std::setw(5) << state2_time << "\n";

			std::cout << "Done" << std::endl;
		}
	}

	outfile.close();

	return 0;
}