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
void print_time(double total_time) 
{
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
    double t = std::atof(argv[2]);
    double T1 = std::atof(argv[3]);
    double T2 = std::atof(argv[4]);
    double iters = std::atof(argv[5]);

    double total_time;
    double lambda = 1/T2 - 1/(2*T1);
    double p_amp = 1 - exp(-t/T1);
    double p_phase = 1 - exp(-t*lambda);
    

    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);

    std::string backend = "cpu";
    std::string sim_method = "stab";
    bool qasm_exists = (argv[6] != nullptr) && (argv[6][0] != '\0');

    if(qasm_exists)
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
        for(int q = 0; q < n_qubits-1; q++)
        {
            circuit->CX(q, q+1);
        }
        for(int q = 0; q < n_qubits; q++)
        {
            circuit->DAMP(q, p_phase, p_amp);
        }
    }

    std::vector<ComplexMatrix> densityMatrices;
    ComplexMatrix avgDM;
    if(n_qubits < 8 && !qasm_exists)
        avgDM = ComplexMatrix::Zero(1<<n_qubits, 1<<n_qubits);
    auto stab_state = BackendManager::create_state(backend, n_qubits, sim_method);
    double timer;
    for(int num_iters = 0; num_iters < iters; num_iters++)
    {
        timer = 0;
        stab_state->sim(circuit, timer);
        total_time += timer;
        if(n_qubits < 8 && !qasm_exists)
            avgDM += stab_state->get_density_matrix();
        stab_state->reset_state();
    }

    //If manually built, return the DM
    if(!qasm_exists)
    {
        avgDM /= static_cast<Complex>(iters);
        print_dm(T2, avgDM, total_time);
    }
    else
        print_time(total_time);

    return 0;
}