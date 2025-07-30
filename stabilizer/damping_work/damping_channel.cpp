#define EIGEN


#include <string>
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>
#include <ranges>

#include "../../include/backendManager.hpp"
#include "../../include/state.hpp"
#include "../../include/circuit.hpp"
#include "../../include/nwq_util.hpp"
#include "../src/qasm_extraction.hpp"

void print_dm(double num, ComplexMatrix rho, double total_time) 
{
    std::cout << "T2=" << num << std::endl;
    for (int i = 0; i < rho.rows(); ++i) {
        for (int j = 0; j < rho.cols(); ++j) {
            std::cout << rho(i, j).real() << " " << rho(i, j).imag() << " ";
        }
    }
    std::cout << std::endl;
    std::cout << "Time=" << total_time << std::endl;
}

int main(int argc, char* argv[]) {
    // Expecting 2 arguments: tau and T1
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <qubits> <tau> <T1> <T2> <shots> <qasm filepath>\n";
        return 1;
    }

    // Convert command-line arguments to doubles
    int n_qubits = std::atof(argv[1]);
    double tau = std::atof(argv[2]);
    double T1 = std::atof(argv[3]);
    double T2 = std::atof(argv[4]);
    double iters = std::atof(argv[5]);


    double total_time = 0;
    double lambda = 1/T2 - 1/(2*T1);
    double p_amp = 1 - exp(-tau/T1);
    double p_phase = .5*(1 - exp(-tau*lambda));
    

    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);

    std::string backend = "cpu";
    std::string sim_method = "stab";

    //Get the base circuit from a qasm file
    if(argv[6])
    {
        std::string qasmfile = argv[6];
        appendQASMToCircuit(circuit, qasmfile, n_qubits);
        for(int n = 0; n < n_qubits; n++)
            circuit->DAMP(n, p_phase, p_amp);
    }
    //Use a manually built circuit
    else
    {
        circuit->H(0);
        circuit->DAMP(0, p_phase, p_amp);
    }

    
    std::vector<ComplexMatrix> densityMatrices;
    ComplexMatrix avgDM = ComplexMatrix::Zero(1<<1, 1<<1);
    auto stab_state = BackendManager::create_state(backend, n_qubits, sim_method);
    for(int num_iters = 0; num_iters < iters; num_iters++)
    {
        double timer = 0;
        stab_state->sim(circuit, timer);
        total_time += timer;  //in seconds
        stab_state->reset_state();
        // avgDM += (stab_state->get_density_matrix());
        // NWQSim::IdxType* results = stab_state->measure_all(shots);
        // std::cout << "Results: ";
        // for(int i = 0; i < shots; i++)
        //     std::cout << results[i] << " ";
    }

    // avgDM /= static_cast<Complex>(iters);

    // std::cout << "Average density matrix:\n" << avgDM << "\n";

    // std::cout << "Trace: " << avgDM.trace() << "\n";

    print_dm(T2, avgDM, total_time/1000); 

    // sim_method = "dm";
    // auto dmCircuit = std::make_shared<NWQSim::Circuit>(n_qubits);
    // dmCircuit->H(0);
    // dmCircuit->CX(0,1);
    // std::vector<NWQSim::IdxType> qubits = {0, 1};
    // // for(int each : std::views::iota(0, n_qubits))
    // // {
    // //     qubits.push_back(each);
    // // }
    // const std::vector<NWQSim::IdxType> qubit_list = qubits;
    // dmCircuit->MOD_NOISE("SET", "T1", 100, qubit_list);
    // dmCircuit->MOD_NOISE("SET", "T2", 100, qubit_list);
    // auto DM_state = BackendManager::create_state(backend, n_qubits, sim_method);
    // std::cout << "Starting DM Sim" << std::endl;
    // DM_state->sim(circuit, timer);
    // DM_state->print_res_state();
    return 0;
}
