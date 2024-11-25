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
	
	int number_of_qubits = 4;
	int num_rounds = 1;
	for(int rounds = num_rounds; rounds <= num_rounds; rounds++){
		for(int n_qubits = number_of_qubits; n_qubits <= number_of_qubits; n_qubits++)
		{
			auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
			circuit -> H(0);
			circuit -> CX(0,1);

			// circuit -> CX(0,1);
			// circuit -> CX(1,2);
			// for(int j = 0; j < rounds; j++)
			// {
			// 	Add some gates to the circuit
			// 	for(int i = 0; i < 1; i++)
			// 	{
			// 		circuit -> H(i);
			// 	}
			// 	for(int i = 0; i < n_qubits-1; i++)
			// 	{
			// 		circuit -> CX(i,i+1);
			// 	}
			// 	for(int i = 0; i < n_qubits; i++)
			// 	{
			// 		circuit -> S(i);
			// 	}
			// }

			// for(int i = n_qubits-1; i <= n_qubits; i++)
			// {
			// 	circuit -> M(i);
			// }
			circuit -> M(1);
			// circuit -> M(1);
			// circuit -> M(2);
			// circuit -> M(3);

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

			std::cout << "stabilizer sim: " << results2[0] << std::endl;

			outfile << std::setw(5) << rounds << " ";
			outfile << std::setw(5) << n_qubits << " ";
			outfile << std::setw(5) << state1_time << " ";
			outfile << std::setw(5) << state2_time << "\n";

			std::cout << "Done\n\n\n" << std::endl;
		}//num qubits
	}//num rounds

	outfile.close();

	auto circuit = std::make_shared<NWQSim::Circuit>(number_of_qubits);
	circuit -> H(0);
	circuit -> CX(0,1);
	circuit -> M(1);

	NWQSim::Tableau test = NWQSim::Tableau(circuit, number_of_qubits);
	test.simulate();

    std::cout << "\n------ initial stabilizers: " << std::endl;
	std::vector<std::string> paulistrings = test.get_stabilizers();
	for(int i = 0; i < paulistrings.size(); i++)
	{
		std::cout << paulistrings[i] << std::endl;
	}
	std::cout << "------" << std::endl; 

	test.add_stabilizer("IXYZ");
	//test.simulate();

	std::cout << "\n------ destabilizers after addition: " << std::endl;
	paulistrings = test.get_destabilizers();
	for(int i = 0; i < paulistrings.size(); i++)
	{
		std::cout << paulistrings[i] << std::endl;
	}
	std::cout << "------" << std::endl; 

	std::cout << "\n------ stabilizers after addition: " << std::endl;
	paulistrings = test.get_stabilizers();
	for(int i = 0; i < paulistrings.size(); i++)
	{
		std::cout << paulistrings[i] << std::endl;
	}
	std::cout << "------" << std::endl; 

	return 0;
}