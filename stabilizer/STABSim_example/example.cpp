#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../../include/backendManager.hpp"
#include "../../include/state.hpp"
#include "../../include/circuit.hpp"
#include "../../include/nwq_util.hpp"

#include "../src/qasm_extraction.hpp"

int main()
{
    /*Number of qubits in the simulation. Qasm files will automatically change this, 
    however, a vector can make dynamically generating circuits of different sizes easier*/
    std::vector<int> circuit_size = {128}; 

    for(int i = 0; i < circuit_size.size(); i++)
    {
        std::cout << "Starting program" << std::endl;

        //Define the number of test qubits and final measurement shots
        int n_qubits = circuit_size[i];
        int shots = 10;
        //If you want a repeated circuit in one run, define repition rounds
        int repetition_rounds = 1;
        //Create a circuit object
        auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);


        /*Append a test circuit from a QASM file to the NWQSim circuit*/
        // std::string inFile = "/Users/garn195/Project Repositories/NWQ-Sim/stabilizer/T_transpilation_test/adder_n10.qasm";
        // if(inFile != "")
        //     appendQASMToCircuit(circuit, inFile, n_qubits);


        /*Or generate a circuit*/

        //Random direct additions to the circuit
        // std::srand(std::time(nullptr));
        // for(int j = 0; j < 100000; j++) 
        // {
        //     circuit->H((std::rand() % (n_qubits-1)));
        //     circuit->CX((std::rand() % (n_qubits-1)),(std::rand() % (n_qubits)));
        //     circuit->S((std::rand() % (n_qubits-1)));
        // }

        /*Or make a layered circuit*/

        //Gate additions for a layered circuit
        // std::vector<NWQSim::Gate> gate_layer;
        // for(int k = 0; k < n_qubits; k++)
        // {
        //     NWQSim::Gate G(NWQSim::OP::H, (std::rand() % (n_qubits-1)));
        //     gate_layer.push_back(G);
        // }
        // std::vector<NWQSim::Gate> full_circuit;

        // for(int j = 0; j < repetition_rounds; j++)
        // {
        //     full_circuit.insert(full_circuit.end(),gate_layer.begin(),gate_layer.end());
        // }
        // circuit->set_gates(full_circuit);


        /*Tiny circuit*/
        circuit -> H(63);
        circuit -> S(63);
        circuit -> S(63);
        circuit -> H(63);


        /*For STAB_GPU-Sim2D, independent gates should be user defined by gate chunks.
        Gate chunks is an array of the sequential number of gates that can be run 
        at a time without overwriting unfinished memory.*/
        std::vector<int> gate_chunks (repetition_rounds, n_qubits);


        std::cout << "Circuit built" << std::endl;


        //"nvgpu" for gpu sim, "cpu" for cpu sim
        std::string backend = "cpu";
        //"stab", "sv", "dm", etc.
        std::string sim_method = "stab";
        
        std::cout << "Creating state" << std::endl;
        double timer = 0;
        auto state = BackendManager::create_state(backend, n_qubits, sim_method);

        std::cout << "Starting sim" << std::endl;

        /*Once a 'state' obejct has been created with the appropriate backend, qubit size, and sim method, its time to simulate*/
        /*For most sim_methods including 1D stabilizer parallelization, use sim(circuit, timer); 
        For parallelized stabilizers and gates on STABSim, use sim2D(circuit, gate_chunks,timer);*/
        // state->sim(circuit, timer);
        state->sim(circuit, timer);

        std::cout << "Simulation complete" << std::endl;

        /*View the tableau of a stabilizer simulation by calling state->print_res_state();*/
        // state->print_res_state();

        /*View a snapshot of the results over some number of shots (without changing the state) by calling state->measure_all(shots);*/
        NWQSim::IdxType* results = state->measure_all(shots);

        /*View each shot result*/
        for(int i = 0; i < shots; i++)
        {
            std::cout << "Result " << i << ": " << results[i] << std::endl;
        }

        /*View the simulation time*/
        std::cout << "Sim time: " << timer/1000.0 << "s" << std::endl;

    }
    return 0;
}
