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
    int n_qubits = 1; //15 // A few hundred qubits 
    int rounds = 1; //3 // Number of rounds to simulate
    
    double timer_cpu = 0;
    double timer_cuda = 0;
    int num_measurements = 0;

    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);

    circuit->H(0);
    circuit->RZ(PI/4, 0);
    circuit->S(0);
    circuit->S(0); 
    circuit->H(0);
    circuit->M(0);

    // Random circuit generator (reproducible)
    std::mt19937 rng(12345); // change seed to vary the circuit
    std::uniform_int_distribution<int> gate_dist(0, 4);          // 0:H,1:S,2:CX,3:RESET,4:M
    std::uniform_int_distribution<int> qubit_dist(0, n_qubits-1);

    // Number of random operations per round (tune as desired)
    // int ops_per_round = n_qubits;

    // for (int round = 0; round < rounds; ++round) 
    // {
    //     for (int g = 0; g < ops_per_round; ++g) 
    //     {
    //         int gate = gate_dist(rng);
    //         if (gate == 2 && n_qubits >= 2) 
    //         {
    //             // CX: pick distinct ctrl/target
    //             int ctrl = qubit_dist(rng);
    //             int target = qubit_dist(rng);
    //             while (target == ctrl) target = qubit_dist(rng);
    //             circuit->CX(ctrl, target);
    //         } 
    //         else 
    //         {
    //             int q = qubit_dist(rng);
    //             switch (gate) 
    //             {
    //                 case 0: 
    //                 {
    //                     circuit->M(q);
    //                     num_measurements++; 
    //                     break;
    //                 }
    //                 case 1: circuit->S(q); break;
    //                 // case 2: circuit->RESET(q); break;
    //                 case 3: 
    //                 {
    //                     circuit->S(q); 
    //                     break;
    //                 }
    //                 default: 
    //                 {
    //                     circuit->H(q);
    //                     circuit->S(q);
    //                     circuit->S(q);
    //                     circuit->H(q);
    //                     break; // fallback
    //                 }
    //             }
    //         }
    //     }
    // }




    // Create CPU and GPU states
    auto cpu_state = BackendManager::create_state("cpu", n_qubits, "stab");

    // auto cpu_state2 = BackendManager::create_state("cpu", n_qubits, "stab");


    // auto cuda_state = BackendManager::create_state("nvgpu", n_qubits, "stab");

    // Allocate measurement buffers for CUDA
    // cuda_state->allocate_measurement_buffers(num_measurements);

    // Simulate on both backends
    std::cout << "Simulating on CPU..." << std::endl;
    double m_results = 0;
    
    std::uniform_int_distribution<int> dist(1, 999999999);

    for(int i = 0; i < 1; i++)
    {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::mt19937 engine(seed);
        cpu_state->set_seed(dist(engine));
        cpu_state->sim(circuit, timer_cpu);
        m_results += cpu_state->get_measurement_results()[0];
        std::cout << cpu_state->get_measurement_results()[0] << std::endl;
        // cpu_state->print_res_state();

        cpu_state->reset_state();
    }
    
    // cpu_state2->sim(circuit, timer_cpu);

    std::cout << "CPU sim time: " << timer_cpu / 1000.0 << "s" << std::endl;
    // auto cpu_measurements = cpu_state->get_measurement_results();
    // cpu_state->print_res_state();
    std::cout << m_results << std::endl;

    // std::cout << "Simulating on GPU..." << std::endl;
    // cuda_state->sim(circuit, timer_cuda);
    // std::cout << "CUDA sim time: " << timer_cuda / 1000.0 << "s" << std::endl;

    // // Compare the results
    // std::cout << "\nComparing final states..." << std::endl;

    // auto cpu_stabilizers = cpu_state->get_stabilizers();
    // auto cpu_stabilizers = cpu_state2->get_stabilizers();

    // auto cuda_stabilizers = cuda_state->get_stabilizers();
    // auto cpu_measurements = cpu_state->get_measurement_results();
    // auto cuda_measurements = cuda_state->get_measurement_results();

    // std::cout << circuit->to_string() << std::endl;
    // cpu_state2->print_res_state();

    // cuda_state->print_res_state();
    



    // Sort both vectors to ensure order doesn't matter
    // std::sort(cpu_stabilizers.begin(), cpu_stabilizers.end());
    // std::sort(cuda_stabilizers.begin(), cuda_stabilizers.end());

    // bool match = true;
    // std::string stab_test(n_qubits, 'I');
    // for (size_t i = 0; i < cuda_stabilizers.size(); ++i) 
    // {
    //     if (cuda_stabilizers[i] ==  stab_test) {
    //         std::cerr << "Indentity Stabilizer in CUDA doesn't make sense" << std::endl;
    //         match = false;
    //         break;
    //     }
    // }


    // if (cpu_stabilizers.size() != cuda_stabilizers.size()) {
    //     std::cerr << "Mismatch in number of stabilizers: "
    //                 << "CPU=" << cpu_stabilizers.size() << ", "
    //                 << "CUDA=" << cuda_stabilizers.size() << std::endl;
    //     match = false;
    // } 
    // else
    // {
    //     for (size_t i = 0; i < cpu_stabilizers.size(); ++i) {
    //         if (cpu_stabilizers[i] != cuda_stabilizers[i]) {
    //             std::cerr << "Mismatch in stabilizer at index " << i << ":\n"
    //                         << "  CPU: " << cpu_stabilizers[i] << "\n"
    //                         << "  CUDA:" << cuda_stabilizers[i] << std::endl;
    //             match = false;
    //         }
    //     }
    // }

    // if (cpu_measurements.size() != cuda_measurements.size()) {
    //     std::cerr << "Mismatch in number of measurements: "
    //                 << "CPU=" << cpu_measurements.size() << ", "
    //                 << "CUDA=" << cuda_measurements.size() << std::endl;
    //     match = false;
    // } 
    // else 
    // {
    //     for (size_t i = 0; i < cpu_measurements.size(); ++i) {
    //         if (c != cuda_measurements[i]) {
    //             std::cerr << "Mismatch in measurement at index " << i << ": "
    //                         << "CPU=" << cpu_measurements[i] << ", "
    //                         << "CUDA=" << cuda_measurements[i] << std::endl;
    //             match = false;
    //         }
    //     }
    // }

    // for (int i = 0; i < cpu_measurements.size(); ++i) {
    //     std::cout << "M" << i << ": "
    //                 << "CPU=" << cpu_measurements[i] << ", "
    //                 << "CUDA=" << cuda_measurements[i] << std::endl;
    // }

    // if (match) {
    //     std::cout << "Validation PASSED: Stabilizers and measurement numbers match." << std::endl;
    // } else {
    //     std::cout << "Validation FAILED: Stabilizers and/or measurement numbers do NOT match." << std::endl;
    // }
    

    // Free measurement buffers for CUDA
    // cuda_state->free_measurement_buffers();

    return 0;
}