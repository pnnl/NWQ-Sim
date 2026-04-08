#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>
#include <random>
#include <algorithm> // For std::sort

#include "../include/backendManager.hpp"
#include "../include/state.hpp"
#include "../include/circuit.hpp"
#include "../include/nwq_util.hpp"
#include "../include/stabsim/stab_cpu.hpp"
#ifdef CUDA_ENABLED
#include "../include/stabsim/stab_cuda.cuh"
#endif

int main() 
{
    int n_qubits = 3; //15 // A few hundred qubits 
    int rounds = 1; //3 // Number of rounds to simulate
    
    double timer_cpu = 0;
    double timer_cuda = 0;
    int num_measurements = 0;

    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);

    circuit->H(0);
    circuit->S(0);
    circuit->RZ(5*PI/4, 0);
    circuit->H(0);
    circuit->M(0);


    auto cpu_state = BackendManager::create_state("cpu", n_qubits, "stab");

    // Simulate on both backends
    std::cout << "Simulating on CPU..." << std::endl;
    double m_results = 0;
    
    std::uniform_int_distribution<int> dist(1, 999999999);
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 engine(seed);
    for(int i = 0; i < 1000; i++)
    {
        // cpu_state->print_res_state();

        cpu_state->set_seed(dist(engine)+i);
        cpu_state->sim(circuit, timer_cpu);
        m_results += cpu_state->get_measurement_results()[0];
        // std::cout << "first: " << cpu_state->get_measurement_results()[0] << std::endl;
        // std::cout << "second: "  << cpu_state->get_measurement_results()[1] << std::endl;
        // std::cout << "third: "  << cpu_state->get_measurement_results()[2] << std::endl;


        // cpu_state->print_res_state();

        cpu_state->reset_state();
    }
    

    std::cout << "CPU sim time: " << timer_cpu / 1000.0 << "s" << std::endl;
    std::cout << m_results << std::endl;

    return 0;
}