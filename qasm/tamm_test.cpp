#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <csignal>
#include "../include/backendManager.hpp"
#include "../include/state.hpp"
#include "../include/circuit.hpp"

using namespace NWQSim;
using ValType = double;

void segfault_handler(int) {
    std::cerr << "[FATAL] Segmentation fault\n";
    std::exit(1);
}

int main(int argc, char** argv) {
    std::signal(SIGSEGV, segfault_handler);
    tamm::initialize(argc, argv);

    int n_qubits = 10;
    int bd_min   = 1;
    int bd_max   = 20;
    int SHOTS    = 8192;

    if (argc > 1) n_qubits = std::stoi(argv[1]);
    if (argc > 2) bd_min   = std::stoi(argv[2]);
    if (argc > 3) bd_max   = std::stoi(argv[3]);
    if (argc > 4) SHOTS    = std::stoi(argv[4]);
    if (bd_min > bd_max) {
        std::cerr << "[ERROR] bd_min (" << bd_min
                  << ") > bd_max (" << bd_max << ")\n";
        return 1;
    }

    std::cerr << "[INFO] n_qubits=" << n_qubits
              << "  bd_min="   << bd_min
              << "  bd_max="   << bd_max
              << "  shots="    << SHOTS << "\n";

    auto circuit = std::make_shared<Circuit>(n_qubits);
    circuit->H(0);
    for (int i = 0; i < n_qubits - 1; ++i) {
        circuit->CX(i, i + 1);
    }

    long long idx_zero = 0LL;
    long long idx_one  = (n_qubits < 63)
                       ? ((1LL << n_qubits) - 1)
                       : -1LL;

    ValType inv_sqrt2 = 1.0/std::sqrt(2.0);

    std::ofstream csv("bond_vs_fidelity.csv");
    if (!csv) {
        std::cerr << "[ERROR] Cannot open bond_vs_fidelity.csv\n";
        return 1;
    }
    csv << "bond_dimension,fidelity\n";
    csv << std::fixed << std::setprecision(8);

    for (int bd = bd_min; bd <= bd_max; ++bd) {
        std::cerr << "[INFO] TN sampling bond_dim=" << bd << "\n";

        auto tn_state = BackendManager::create_state("CPU",
                                                     n_qubits,
                                                     "TN",
                                                     bd);
        tn_state->sim(circuit);
        long long* samples = tn_state->measure_all(SHOTS);

        long long cnt0 = 0, cnt1 = 0;
        for (int i = 0; i < SHOTS; ++i) {
            if (samples[i] == idx_zero)     ++cnt0;
            else if (samples[i] == idx_one) ++cnt1;
        }

        ValType sum   = ValType(cnt0 + cnt1);
        ValType beta0 = std::sqrt(cnt0 / sum);
        ValType beta1 = std::sqrt(cnt1 / sum);

        ValType overlap = inv_sqrt2*(beta0 + beta1);
        ValType fidelity = overlap*overlap;

        csv << bd << "," << fidelity << "\n";
        std::cerr << "[RESULT] bd=" << bd
                  << "  fidelity=" << fidelity << "\n";
    }

    tamm::finalize();
    csv.close();
    return 0;
}

