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

    int number_of_qubits = 4;

	auto circuit = std::make_shared<NWQSim::Circuit>(number_of_qubits);
	

	NWQSim::Tableau test = NWQSim::Tableau(number_of_qubits);


	test.replace_stabilizer("XXII", 0);
	test.replace_stabilizer("ZZII", 1);
	test.replace_stabilizer("IIXX", 2);
	test.replace_stabilizer("IIZZ", 3);

	std::cout << "\n------ stabilizers: " << std::endl;
	std::vector<std::string> paulistrings = test.get_stabilizers();
	for(int i = 0; i < paulistrings.size(); i++)
	{
		std::cout << paulistrings[i] << std::endl;
	}
	std::cout << "------" << std::endl;

	test.print_table(test);

	// std::vector<std::vector<uint>> testGraph = test.convert_to_graph();

	// test.print_table(testGraph);

	ComplexMatrix DM = test.get_DM();

	std::cout << DM << std::endl;

	return 0;
}	