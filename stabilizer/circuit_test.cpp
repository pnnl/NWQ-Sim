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
#include "../include/stabsim/stab_cuda.cuh"

int main() {
    int n_qubits = 7; // A few hundred qubits
    int rounds = 1; // Number of rounds to simulate
    
    double timer_cpu = 0;
    double timer_cuda = 0;
    int num_measurements = 0;

    // Create a deterministic circuit
    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
    for(int i = 0; i < rounds; ++i) 
    {
        for (int i = 0; i < n_qubits; i+=2) 
        {
            circuit->M(i%5);
            num_measurements++;

            circuit->H(i);
            circuit->S(i);
            circuit->S(i);
            circuit->H(i);
            circuit->CX(2, i);
            circuit->CX(0, 1);
            
            circuit->RESET(i%3);
            circuit->H(i%7);

        }
        // circuit->H(5);
        // circuit->M(5);
        // Add some measurement gates to test that part of the logic
        for (int i = 0; i < n_qubits; ++i) 
        {
            circuit->CX(1, 2);
            circuit->S(i);
            circuit->M(i);
            circuit->CX(1, 0);
            circuit->RESET(i%10);
            num_measurements++;
        }
    }
    


    // Create CPU and GPU states
    auto cpu_state = BackendManager::create_state("cpu", n_qubits, "stab");
    auto cuda_state = BackendManager::create_state("nvgpu", n_qubits, "stab");

    // Allocate measurement buffers for CUDA
    cuda_state->allocate_measurement_buffers(num_measurements);

    // Simulate on both backends
    std::cout << "Simulating on CPU..." << std::endl;
    cpu_state->sim(circuit, timer_cpu);
    std::cout << "CPU sim time: " << timer_cpu / 1000.0 << "s" << std::endl;

    std::cout << "Simulating on GPU..." << std::endl;
    cuda_state->sim(circuit, timer_cuda);
    std::cout << "CUDA sim time: " << timer_cuda / 1000.0 << "s" << std::endl;

    // Compare the results
    std::cout << "\nComparing final states..." << std::endl;

    auto cpu_stabilizers = cpu_state->get_stabilizers();
    auto cuda_stabilizers = cuda_state->get_stabilizers();
    auto cpu_measurements = cpu_state->get_measurement_results();
    auto cuda_measurements = cuda_state->get_measurement_results();

    // cpu_state->print_res_state();
    // cuda_state->print_res_state();
    

    // Sort both vectors to ensure order doesn't matter
    std::sort(cpu_stabilizers.begin(), cpu_stabilizers.end());
    std::sort(cuda_stabilizers.begin(), cuda_stabilizers.end());

    bool match = true;
    if (cpu_stabilizers.size() != cuda_stabilizers.size()) {
        std::cerr << "Mismatch in number of stabilizers: "
                    << "CPU=" << cpu_stabilizers.size() << ", "
                    << "CUDA=" << cuda_stabilizers.size() << std::endl;
        match = false;
    } else {
        for (size_t i = 0; i < cpu_stabilizers.size(); ++i) {
            if (cpu_stabilizers[i] != cuda_stabilizers[i]) {
                std::cerr << "Mismatch in stabilizer at index " << i << ":\n"
                            << "  CPU: " << cpu_stabilizers[i] << "\n"
                            << "  CUDA:" << cuda_stabilizers[i] << std::endl;
                match = false;
            }
        }
    }

    if (cpu_measurements.size() != cuda_measurements.size()) {
        std::cerr << "Mismatch in number of measurements: "
                    << "CPU=" << cpu_measurements.size() << ", "
                    << "CUDA=" << cuda_measurements.size() << std::endl;
        match = false;
    } else {
        for (size_t i = 0; i < cpu_measurements.size(); ++i) {
            if (cpu_measurements[i] != cuda_measurements[i]) {
                std::cerr << "Mismatch in measurement at index " << i << ": "
                            << "CPU=" << cpu_measurements[i] << ", "
                            << "CUDA=" << cuda_measurements[i] << std::endl;
                match = false;
            }
        }
    }

    // for (int i = 0; i < n_qubits; ++i) {
    //     std::cout << "Match at " << i << ": "
    //                 << "CPU=" << cpu_measurements[i] << ", "
    //                 << "CUDA=" << cuda_measurements[i] << std::endl;
    // }

    if (match) {
        std::cout << "Validation PASSED: Stabilizers and measurements match." << std::endl;
    } else {
        std::cout << "Validation FAILED: Stabilizers and/or measurements do NOT match." << std::endl;
    }
    

    // Free measurement buffers for CUDA
    cuda_state->free_measurement_buffers();

    return 0;
}