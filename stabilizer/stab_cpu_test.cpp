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
    int n_qubits = 2;
    int shots = 1;

    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);

    // for(int i = 0; i < 500000; i++)
    // {
    //     circuit->CX(0,1);
    //     circuit->CX(1,2);
    //     circuit->CX(2,3);
    // }
    circuit->X(0);
    circuit->X(1);
    // circuit->M(0);
    // circuit->M(1);

    std::string backend = "CPU";
    std::string sim_method = "stab";
    double timer = 0;
    auto state = BackendManager::create_state(backend, n_qubits, sim_method);
    
    state->sim(circuit, timer);

    // std::cout << "Result: " << *(state->get_results()) << std::endl;

    std::cout << "Result: " << *(state->measure_all(shots)) << std::endl;



    // NWQSim::IdxType* result = state->measure_all(shots);

    // for(int i = 0; i < shots; i++)
    // {
    //     std::cout << "Result " << i << ": " << result[i] << std::endl;
    // }


    // state->replace_stabilizer("XXII", 0);
	// state->replace_stabilizer("ZZII", 1);
	// state->replace_stabilizer("IIXX", 2);
	// state->replace_stabilizer("IIZZ", 3);
    // std::cout << "Actual result: " << *(state->get_results()) << std::endl;

    //std::vector<std::vector<double>> DM = state->get_density_matrix();
    // for(int i = 0; i < DM.size(); i++)
    // {
    //     for(int j = 0; j < DM.size(); j++)
    //     {
    //         std::cout << DM[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // std::vector<std::string> stabs = state->get_stabilizers();
    // for(int i = 0; i < stabs.size(); i++)
    // {
    //     std::cout << stabs[i] << std::endl;
    // }

    //std::vector<std::vector<int>> GM = state->get_graph_matrix();
    // for(int i = 0; i < GM.size(); i++)
    // {
    //     for(int j = 0; j < GM.size(); j++)
    //     {
    //         std::cout << GM[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // sim_method = "sv";
    // auto state2 = BackendManager::create_state(backend, n_qubits, sim_method);
    // double timer2;

    // state2->sim(circuit, timer2);
    
    // long long int* result2 = state2->measure_all(shots);

    // for(int i = 0; i < shots; i++)
    // {
    //     std::cout << "Result " << i << ": " << result2[i] << std::endl;
    // }

	return 0;
}	